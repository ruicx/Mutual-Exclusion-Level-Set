% name:       cycle_init.m
% usage:      --
% author:     Ruicheng
% date:       2020-08-05 09:19:17
% version:    1.0
% Env.:       MATLAB R2019b, WIN10


function masks = cycle_init(img_size, centers, radius)
%initial cycle contour
%
% Syntax: masks = cycle_init(img_size, centers)
%
% initial cycle contour

    h = img_size(1);
    w = img_size(2);
    assert(size(centers, 1) == length(radius),...
        "centers and radius must have the same number")

    n_radius = length(radius);

    masks = -2 * ones(h, w, n_radius);
    for aa = 1:n_radius
        c1 = centers(aa, 1);
        c2 = centers(aa, 2);
        r = radius(aa);
        temp = -2 * ones(h, w);
        % let the cycle (x-c1)^2 + (x-c2)^2 < r^2 to be 2
        temp(((1:h)' - c1).^2 + ((1:w) - c2).^2 < r^2) = 2;
        masks(:, :, aa) = temp;
    end
end
