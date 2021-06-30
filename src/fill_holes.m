% name:       fill_holes.m
% usage:      --
% author:     Ruicheng
% date:       2020-08-05 10:27:16
% version:    1.0
% Env.:       MATLAB R2019b, WIN10


function img_fill = fill_holes(img, threshold)
%fill_holes - fill hold and remove white point
%
% Syntax: fill_img = fill_holes(img, threshold)
%
% fill hold and remove white point
    if ~exist("threshold", "var")
        threshold = 10;
    end

    % fill holes
    img_fill = imfill(img, 'holes');

    % remove small white point
    img_fill = bwareaopen(img_fill > 0, threshold);

    % reset level set
    img_fill = double(img_fill);
    img_fill(img_fill == 1) = 2;
    img_fill(img_fill == 0) = -2;

end