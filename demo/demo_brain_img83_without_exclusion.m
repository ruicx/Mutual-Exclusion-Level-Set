% name:       demo_brain_img83_without_exclusion.m
% usage:      --
% author:     Ruicheng
% date:       2020-08-21 09:07:21
% version:    1.0
% Env.:       MATLAB R2019b, WIN10

clear
close all

data_name = "brain_img83";
addpath("../src")
load("../data/" + data_name + ".mat", "img")
data_name = "brain_img83p";

param.numIter   = 50;
param.lambda1   = 1.5e-1;
param.lambda2   = 1.3e-1;
param.nu        = 0.010*255*255; % coefficient of the length term
param.timestep  = .05;           % time step
param.mu        = 1;             % coefficient of the level set (distance) regularization term P(\phi)
param.epsilon   = 1;             % the papramater in the definition of smoothed Dirac function
param.sigma     = 5;             % scale parameter in Gaussian kernel
param.alpha     = 0;             % coefficient of prior term
param.eta       = 0;             % coefficient of mutual exclusion term
param.draw_step = 1;             % draw step

init_centers = [94, 87; 97, 98];
init_raduis = [5, 5];

% zero prior to hold the place
prior1 = zeros(size(img));
prior2 = zeros(size(img));

% detail index
idx1 = 75:125;
idx2 = 65:115;

% initial contour
init_contours = cycle_init(size(img), init_centers, init_raduis);
init1 = init_contours(:, :, 1);
init2 = init_contours(:, :, 2);

% segmentation
param1 = param;
param1.numIter = 5;
[u1_process, u2_process] = double_RSF(init1, init2, img, prior1, prior2, param1);

% show process contour
draw_segmentation(img(idx1, idx2), cat(3, u1_process(idx1, idx2), u2_process(idx1, idx2)), ['r', 'b'], [min(img(:)), max(img(:))])
% print("../figures/"+ data_name + "_process_detail_without_exclusion.eps", '-depsc2', '-r600')

param.numIter = param.numIter - param1.numIter;
[u1, u2] = double_RSF(u1_process, u2_process, img, prior1, prior2, param);

% fill holes
seg1 = fill_holes(u1);
seg2 = fill_holes(u2);

% filter region
seg1 = filter_region(seg1, prior1);
seg2 = filter_region(seg2, prior2);

% show original image
figure
color_range = [min(img(:)), max(img(:))];
imagesc(img, color_range);
colormap(gray); axis off; axis equal
% print("../figures/" + data_name + "_original_image_without_exclusion", '-depsc2', '-r600')
imagesc(img(idx1, idx2), color_range)
colormap(gray); axis off; axis equal
% print("../figures/" + data_name + "_original_image_detail_without_exclusion", '-depsc2', '-r600')

% show initial contour
draw_segmentation(img, init_contours, ['r', 'b'], color_range)
% print("../figures/" + data_name + "_initial_contour_without_exclusion", '-depsc2', '-r600')
draw_segmentation(img(idx1, idx2), init_contours(idx1, idx2, :), ['r', 'b'], color_range)
% print("../figures/" + data_name + "_initial_contour_detail_without_exclusion", '-depsc2', '-r600')

% show finial contour
draw_segmentation(img, cat(3, seg1, seg2), ['r', 'b'], color_range)
% print("../figures/" + data_name + "_result_without_exclusion", '-depsc2', '-r600')
draw_segmentation(img(idx1, idx2), cat(3, seg1(idx1, idx2), seg2(idx1, idx2)), ['r', 'b'], color_range)
% print("../figures/" + data_name + "_result_detail_without_exclusion", '-depsc2', '-r600')
