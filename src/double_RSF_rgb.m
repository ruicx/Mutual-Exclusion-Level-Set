function [u1, u2] = double_RSF_rgb(varargin)

if length(varargin) == 3
    u1        = varargin{1};
    u2        = varargin{2};
    Img       = varargin{3};
    param1.lambda1   = 1.5 * 1e-1;
    param1.lambda2   = 1.3 * 1e-1;
    param1.nu        = 0.008*255*2555; % coefficient of the length term
    param1.timestep  = .1;             % time step
    param1.mu        = 1;              % coefficient of the level set (distance) regularization term P(\phi)
    param1.epsilon   = 1;              % the papramater in the definition of smoothed Dirac function
    param1.sigma     = 5.0;            % Gaussian kernel size
    param1.alpha     = 5.4;            % coefficient of prior term
    param1.eta       = 60;             % coefficient of mutual exclusion term
    param1.draw_step = 0;              % draw step, default 0 is not draw
    param2           = param1;
elseif length(varargin) == 6
    u1 = varargin{1};
    u2 = varargin{2};
    Img = varargin{3};
    prior1 = varargin{4};
    prior2 = varargin{5};
    param1 = varargin{6};
    param2 = param1;
elseif length(varargin) == 7
    u1 = varargin{1};
    u2 = varargin{2};
    Img = varargin{3};
    prior1 = varargin{4};
    prior2 = varargin{5};
    param1 = varargin{6};
    if iscell(varargin{7})
        % double_RSF(u1, u2, Img, prior1, prior2, param, callbacks)
        param2 = param1;
        callbacks = varargin{7};
    else
        % double_RSF(u1, u2, Img, prior1, prior2, param1, param2)
        param2 = varargin{7};
    end
elseif length(varargin) == 8
    % double_RSF(u1, u2, Img, prior1, prior2, param1, param2, callbacks)
    u1 = varargin{1};
    u2 = varargin{2};
    Img = varargin{3};
    prior1 = varargin{4};
    prior2 = varargin{5};
    param1 = varargin{6};
    param2 = varargin{7};
    callbacks = varargin{8};
else
    error('Number of input cannot be resove')
end

% warning for inequality param
numIter = param1.numIter;
if param1.numIter ~= param2.numIter
    warning('param1.numIter ~= param2.numIter using param1.numIter')
end

% default not draw
if isfield(param1, 'draw_step')
    draw_step = param1.draw_step;
    if param1.draw_step ~= param2.draw_step
        warning('param1.draw_step ~= param2.draw_step using param1.draw_step')
    end
else
    draw_step = 0;
end

sigma = param1.sigma;
param1.Ksigma = fspecial('gaussian', round(2 * sigma) * 2 + 1, sigma); % the Gaussian kernel
param1.KI = split_conv3(Img, param1.Ksigma);
param1.KONE = split_conv3(ones(size(Img)), param1.Ksigma);
sigma = param2.sigma;
param2.Ksigma = fspecial('gaussian', round(2 * sigma) * 2 + 1, sigma); % the Gaussian kernel
param2.KI = split_conv3(Img, param2.Ksigma);
param2.KONE = split_conv3(ones(size(Img)), param2.Ksigma);


% u = u0;
for aa = 1:numIter

    % RSF moethd
    if exist('callbacks', 'var')
        [u1, ds1] = mutual_rsf(u1, u2, prior1, Img, param1);
        [u2, ds2] = mutual_rsf(u2, u1, prior2, Img, param2);
    else
        u1 = mutual_rsf(u1, u2, prior1, Img, param1);
        u2 = mutual_rsf(u2, u1, prior2, Img, param2);
    end

    % plot contour
    if draw_step > 0 && mod(aa, draw_step) == 0
        if ~exist('fhandle1', 'var')
            fhandle1 = figure;
        end
        figure(fhandle1);
        % imagesc(Img, [min(Img(:)), max(Img(:))]);
        imshow(uint8(Img));
        colormap(gray); hold on; axis off, axis equal
        contour(u1, [0 0], 'r');
        contour(u2, [0 0], 'b');
        title([num2str(aa), ' iterations']);
        hold off;
    end

    % callbacks
    if exist('callbacks', 'var')
        for cb_cell = callbacks
            cb_fun = cb_cell{:};
            cb_fun(u1, ds1, u2, ds2, aa);
        end
    end

end % ! END for

end % ! END RSF


function varargout = mutual_rsf(u, phi, prior, Img, param)
%mutual_rsf - single step of RSF method with region mutual exclusion term
%
% Syntax: varargout = mutual_rsf(u, phi, prior, Img, lambda1, lambda2, nu, timestep, mu, epsilon, alpha, eta, KI, KONE, Ksigma)
%
% single step of RSF method with region mutual exclusion term

    lambda1   = param.lambda1; 
    lambda2   = param.lambda2; 
    nu        = param.nu;      % coefficient of the length term
    timestep  = param.timestep;% time step
    mu        = param.mu;      % coefficient of the level set (distance) regularization term P(\phi)
    epsilon   = param.epsilon; % the papramater in the definition of smoothed Dirac function
    alpha     = param.alpha;   % coefficient of prior term
    eta       = param.eta;     % coefficient of mutual exclusion term
    KI        = param.KI;
    KONE      = param.KONE;
    Ksigma    = param.Ksigma;  % Gaussian kernel size

    u = NeumannBoundCond(u);
    K = curvature_central(u);

    DrcU = (epsilon/pi) ./ (epsilon^2.+u.^2); % eq.(9)

    [f1, f2] = localBinaryFit(Img, u, KI, KONE, Ksigma, epsilon);


    % compute lambda1*e1-lambda2*e2
    s1 = lambda1 .* f1.^2 - lambda2 .* f2.^2;      % compute lambda1*e1-lambda2*e2 in the 1st term in eq. (15) in IEEE TIP 08
    s2 = lambda1 .* f1- lambda2.*f2;
    dataForce = (lambda1 -lambda2) * KONE .* Img .* Img + split_conv3(s1, Ksigma)-2.*Img.*split_conv3(s2,Ksigma);
                                                   % eq.(15)
    A = -DrcU .* sum(dataForce, 3);                % 1st term in eq. (15)
    P = mu * (4 * del2(u) - K);                    % 3rd term in eq. (15), where 4*del2(u) computes the laplacian (d^2u/dx^2 + d^2u/dy^2)
    L = nu .* DrcU .* K;                           % 2nd term in eq. (15)
    AT = -alpha * (u - prior);                     % prior constract term
    MET = -DrcU .* eta .* Heaviside(phi, epsilon); % mutual exclusion term
    % fprintf('%.4f %.4f\n', sum(MET(:)), sum(AT(:)));
    u = u + timestep * (L + P + A + AT + MET);     % eq.(15)

    if nargout == 1
        varargout{1} = u;
    else
        ds.L = L;
        ds.P = P;
        ds.A = A;
        ds.AT = AT;
        ds.MET = MET;
        ds.DrcU = DrcU;
        ds.Heaviside_phi = Heaviside(phi, epsilon);
        varargout = {u, ds};
    end

end %! END mutual_rsf


function [f1, f2] = localBinaryFit(Img, u, KI, KONE, Ksigma, epsilon)
    % compute f1 and f2
    Hu = 0.5*(1+(2/pi)*atan(u./epsilon)); % eq.(8)

    I = Img.*Hu;
    c1 = split_conv3(Hu,Ksigma);
    c2 = split_conv3(I,Ksigma);          % the numerator of eq.(14) for i = 1
    f1 = c2./(c1);                        % compute f1 according to eq.(14) for i = 1
    f2 = (KI-c2)./(KONE-c1);              % compute f2 according to the formula in Section IV-A,
                                          % which is an equivalent expression of eq.(14) for i = 2.
end % ! END localBinaryFit


function g = NeumannBoundCond(f)
    % Neumann boundary condition
    [nrow,ncol] = size(f);
    g = f;
    g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);
    g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);
    g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);
end % ! END NeumannBoundCond


function k = curvature_central(u)
    % compute curvature
    [ux,uy] = gradient(u);
    normDu = sqrt(ux.^2+uy.^2+1e-10);   % the norm of the gradient plus a small possitive number
                                        % to avoid division by zero in the following computation.
    Nx = ux./normDu;
    Ny = uy./normDu;
    [nxx,~] = gradient(Nx);
    [~,nyy] = gradient(Ny);
    k = nxx+nyy;                        % compute divergence
end % ! END curvature_central

function rsl = split_conv3(left, right)
    rsl = ones(size(left));
    for aa = 1:size(left, 3)
        rsl(:,:,aa) = conv2(left(:,:,aa), right, 'same');
    end
end % ! split_conv3 END