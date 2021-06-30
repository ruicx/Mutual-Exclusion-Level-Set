% name:       Heaviside.m
% usage:      --
% date:       2020-08-19 14:46:02
% version:    1.0
% Env.:       MATLAB R2019b, WIN10


function h = Heaviside(x, epsilon)
    h = 0.5 * (1 + (2 / pi) .* atan(x ./ epsilon));
end
