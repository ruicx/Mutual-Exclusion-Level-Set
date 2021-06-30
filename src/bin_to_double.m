% name:       bin_to_double.m
% usage:      --
% author:     Ruicheng
% date:       2020-08-05 15:59:19
% version:    1.0
% Env.:       MATLAB R2019b, WIN10


function img = bin_to_double(img_bin, initial_value)
%bin_to_double - convert binary image to double level set
%
% Syntax: img = bin_to_double(img_bin, initial_value)
%
% convert binary image to double level set
    if ~exist('initial_value', 'var')
        initial_value = 2;
    end

    img = double(img_bin);
    img(img == 0) = -1;
    img = img * initial_value;
end