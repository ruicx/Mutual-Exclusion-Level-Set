function u_cell = multi_RSF(varargin)

if length(varargin) == 4
    % double_RSF(Img, init_cell, prior_cell, param_cell)
    Img        = varargin{1};
    u_cell     = varargin{2};
    prior_cell = varargin{3};
    param_cell = varargin{4};
elseif length(varargin) == 5
    % double_RSF(Img, init_cell, prior_cell, param_cell, callbacks)
    Img        = varargin{1};
    u_cell     = varargin{2};
    prior_cell = varargin{3};
    param_cell = varargin{4};
    callbacks  = varargin{5};
else
    error('Number of input cannot be resove')
end

if length(u_cell) ~= length(prior_cell)
    error('Number of initial contour and prior are inequal')
end

% padding param
if length(param_cell) == 1
    param_cell = repmat(param_cell, length(u_cell), 1);
elseif length(u_cell) ~= length(param_cell)
    error('Number of initial contour and number of param are inequal')
end

% warning for inequality param
numIter = param_cell{1}.numIter;
if length(unique(cellfun(@(x) x.numIter, param_cell))) > 1
    warning('numIter in params are inequal, using the first one')
end

% default not draw
if isfield(param_cell{1}, 'draw_step')
    draw_step = param_cell{1}.draw_step;
    if length(unique(cellfun(@(x) x.numIter, param_cell))) > 1
        warning('param1.draw_step ~= param2.draw_step using param1.draw_step')
    end
else
    draw_step = 0;
end

for aa = 1:length(param_cell)
    sigma = param_cell{aa}.sigma;
    param_cell{aa}.Ksigma = fspecial('gaussian', round(2 * sigma) * 2 + 1, sigma); % the Gaussian kernel
    param_cell{aa}.KI = conv2(Img, param_cell{aa}.Ksigma, 'same');
    param_cell{aa}.KONE = conv2(ones(size(Img)), param_cell{aa}.Ksigma, 'same');
end

% u = u0;
for aa = 1:numIter

    % RSF moethd
    for bb = 1:length(u_cell)
        u = u_cell{bb};
        prior = prior_cell{bb};
        param = param_cell{bb};
        u_other = [u_cell(1:(bb-1)), u_cell((bb+1):end)];
        if exist('callbacks', 'var')
            [u, ds] = mutual_rsf(u, u_other, prior, Img, param);
        else
            u = mutual_rsf(u, u_other, prior, Img, param);
        end
        u_cell{bb} = u;
    end

    % plot contour
    if draw_step > 0 && mod(aa, draw_step) == 0
        if ~exist('fhandle1', 'var')
            fhandle1 = figure;
        end
        figure(fhandle1);
        imagesc(Img, [min(Img(:)), max(Img(:))]);
        colormap(gray); hold on; axis off, axis equal
        for cc = 1:length(u_cell)
            contour(u_cell{cc}, [0 0], param_cell{cc}.color);
        end
        title([num2str(aa), ' iterations']);
        hold off;
    end

    % callbacks
    if exist('callbacks', 'var')
        for cb_cell = callbacks
            cb_fun = cb_cell{:};
            cb_fun(u, ds, aa);
        end
    end

end % ! END for

end % ! END RSF


function varargout = mutual_rsf(u, u_other, prior, Img, param)
%mutual_rsf - single step of RSF method with region mutual exclusion term
%
% Syntax: varargout = mutual_rsf(u, u_other, prior, Img, param)
%
% single step of RSF method with region mutual exclusion term

    lambda1   = param.lambda1;  
    lambda2   = param.lambda2;  
    nu        = param.nu;       % coefficient of the length term
    timestep  = param.timestep; % time step
    mu        = param.mu;       % coefficient of the level set (distance) regularization term P(\phi)
    epsilon   = param.epsilon;  % the papramater in the definition of smoothed Dirac function
    alpha     = param.alpha;    % coefficient of prior term
    eta       = param.eta;      % coefficient of mutual exclusion term
    KI        = param.KI;
    KONE      = param.KONE;
    Ksigma    = param.Ksigma;

    u = NeumannBoundCond(u);
    K = curvature_central(u);

    DrcU = (epsilon/pi) ./ (epsilon^2.+u.^2); % eq.(9)

    [f1, f2] = localBinaryFit(Img, u, KI, KONE, Ksigma, epsilon);


    % compute lambda1*e1-lambda2*e2
    s1 = lambda1 .* f1.^2 - lambda2 .* f2.^2;      % compute lambda1*e1-lambda2*e2 in the 1st term in eq. (15) in IEEE TIP 08
    s2 = lambda1 .* f1- lambda2.*f2;
    dataForce = (lambda1 -lambda2) * KONE .* Img .* Img + conv2(s1, Ksigma, 'same')-2.*Img.*conv2(s2,Ksigma,'same');
                                                   % eq.(15)
    A = -DrcU .* dataForce;                        % 1st term in eq. (15)
    P = mu * (4 * del2(u) - K);                    % 3rd term in eq. (15), where 4*del2(u) computes the laplacian (d^2u/dx^2 + d^2u/dy^2)
    L = nu .* DrcU .* K;                           % 2nd term in eq. (15)
    AT = -alpha * (u - prior);                     % prior constract term

    % mutual exclusion term
    MET = 0;
    for cc = 1:length(u_other)
        phi = u_other{cc};
        MET = MET + (-DrcU .* eta .* Heaviside(phi, epsilon)); % mutual exclusion term
    end

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
    c1 = conv2(Hu,Ksigma,'same');
    c2 = conv2(I,Ksigma,'same');          % the numerator of eq.(14) for i = 1
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
