% name:       draw_segmentation.m
% usage:      --
% author:     Ruicheng
% date:       2020-08-05 09:34:34
% version:    1.0
% Env.:       MATLAB R2019b, WIN10

function draw_segmentation(img, contours, colors, range)
%draw_segmentation - Description
%
% Syntax: draw_segmentation(img, contours, colors, range)
%
% Long description
    if ~exist('range', 'var')
        range = [min(img(:)), max(img(:))];
    end

    if islogical(contours)
        contour_range = [0.5, 0.5];
    else
        contour_range = [0, 0];
    end

    figure
    imagesc(img, range);
    colormap(gray); hold on; axis off, axis equal
    contours = double(contours);
    for aa = 1:size(contours, 3)
        contour(contours(:, :, aa), contour_range, colors(aa));
    end
end
