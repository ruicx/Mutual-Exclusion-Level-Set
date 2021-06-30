% name:       demo_synthetic_inho1.m
% usage:      --
% author:     Ruicheng
% date:       2020-08-05 10:37:24
% version:    1.0
% Env.:       MATLAB R2019b, WIN10

clear
close all

data_name = "synthetic_inhomogeneous_img5";
addpath("../src")
img = load("../data/" + data_name + ".mat", "Img");
img = img.Img;

% segmentation parameter
param.numIter   = 150;
param.lambda1   = 1.15e-1;
param.lambda2   = 1.3e-1;
param.nu        = 0.022*255*255; % coefficient of the length term
param.timestep  = .08;           % time step
param.mu        = 1;             % coefficient of the level set (distance) regularization term P(\phi)
param.epsilon   = 1;             % the papramater in the definition of smoothed Dirac function
param.sigma     = 3;             % scale parameter in Gaussian kernel
param.alpha     = 1.6;           % coefficient of prior term
param.eta       = 1e2;           % coefficient of mutual exclusion term
param.draw_step = 1;             % draw step

init_centers = [85, 108; 168, 108; 168, 192; 85, 192];
init_raduis = [45, 45, 45, 45];

% initial contour
init_contours = cycle_init(size(img), init_centers, 0.8 * init_raduis);
init_contours = init_contours > 0;
init1 = init_contours(:, :, 1) | init_contours(:, :, 3);
init1 = bin_to_double(init1);
init2 = init_contours(:, :, 2) | init_contours(:, :, 4);
init2 = bin_to_double(init2);

% prior contour
prior_contours = cycle_init(size(img), init_centers, 1.3 * init_raduis);
prior_contours = prior_contours > 0;
prior1 = prior_contours(:, :, 1) | prior_contours(:, :, 3);
prior1 = bin_to_double(prior1);
prior2 = prior_contours(:, :, 2) | prior_contours(:, :, 4);
prior2 = bin_to_double(prior2);

% segmentation
[u1, u2] = double_RSF(init1, init2, img, prior1, prior2, param);

% compare segmented region
imshowpair(u1 > 0, u2 > 0, "montage");

% fill holes
seg1 = fill_holes(u1);
seg2 = fill_holes(u2);

% show original image
figure
imagesc(img, [min(img(:)), max(img(:))]);
colormap(gray); axis off; axis equal
% print("../figures/" + data_name + "_original_image", '-depsc2', '-r600')

% show initial contour
draw_segmentation(img, cat(3, init1, init2), ['r', 'b'])
% print("../figures/" + data_name + "_initial_contour", '-depsc2', '-r600')

% show prior contour
draw_segmentation(img, cat(3, prior1, prior2), ['r', 'b'])
% print("../figures/" + data_name + "_prior", '-depsc2', '-r600')

% show finial contour
draw_segmentation(img, cat(3, seg1, seg2), ['r', 'b'])
% print("../figures/" + data_name + "_result", '-depsc2', '-r600')
