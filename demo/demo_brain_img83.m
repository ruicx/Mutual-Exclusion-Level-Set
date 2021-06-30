% name:       demo_brain_img83.m
% usage:      --
% author:     Ruicheng
% date:       2020-08-05 08:19:22
% version:    1.0
% Env.:       MATLAB R2019b, WIN10

clear
close all

data_name = "brain_img83";
addpath("../src")
load("../data/" + data_name + ".mat", "img")

param.numIter   = 50;
param.lambda1   = 1.3e-1;
param.lambda2   = 1.3e-1;
param.nu        = 0.008*255*255; % coefficient of the length term
param.timestep  = .05;           % time step
param.mu        = 1;             % coefficient of the level set (distance) regularization term P(\phi)
param.epsilon   = 1;             % the papramater in the definition of smoothed Dirac function
param.sigma     = 5;             % scale parameter in Gaussian kernel
param.alpha     = 0;             % coefficient of prior term
param.eta       = 2e3;           % coefficient of mutual exclusion term
param.draw_step = 1;             % draw step

init_centers = [94, 87; 97, 98];
init_raduis = [5, 5];

% zero prior to hold the place
prior1 = zeros(size(img));
prior2 = zeros(size(img));

% initial contour
init_contours = cycle_init(size(img), init_centers, init_raduis);
init1 = init_contours(:, :, 1);
init2 = init_contours(:, :, 2);

% segmentation
[u1, u2] = double_RSF(init1, init2, img, prior1, prior2, param);

% fill holes
seg1 = fill_holes(u1);
seg2 = fill_holes(u2);

% detail index
idx1 = 75:125;
idx2 = 65:115;

% show original image
figure
color_range = [min(img(:)), max(img(:))];
imagesc(img, color_range);
colormap(gray); axis off; axis equal
% print("../figures/" + data_name + "_original_image", '-depsc2', '-r600')
imagesc(img(idx1, idx2), [min(img(:)), max(img(:))])
colormap(gray); axis off; axis equal
% print("../figures/" + data_name + "_original_image_detail", '-depsc2', '-r600')

% show initial contour
draw_segmentation(img, init_contours, ['r', 'b'], color_range)
% print("../figures/" + data_name + "_initial_contour", '-depsc2', '-r600')
draw_segmentation(img(idx1, idx2), init_contours(idx1, idx2, :), ['r', 'b'], color_range)
% print("../figures/" + data_name + "_initial_contour_detail", '-depsc2', '-r600')

% show finial contour
draw_segmentation(img, cat(3, seg1, seg2), ['r', 'b'], color_range)
% print("../figures/" + data_name + "_result", '-depsc2', '-r600')
draw_segmentation(img(idx1, idx2), cat(3, seg1(idx1, idx2), seg2(idx1, idx2)), ['r', 'b'], color_range)
% print("../figures/" + data_name + "_result_detail", '-depsc2', '-r600')
