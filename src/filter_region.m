% name:       filter_region.m
% usage:      --
% author:     Ruicheng
% date:       2020-08-16 09:04:44
% version:    1.0
% Env.:       MATLAB R2019b, WIN10


function seg = filter_region(seg_raw, prior)
%filter_region - filter region accroding to prior
%
% Syntax: seg = filter_region(seg_raw, prior)
%
% filter region according to prior
    [seg_label, n_label] = bwlabel(seg_raw > 0);
    max_dice = 0;
    max_dice_label = 0;
    for aa = 1:n_label
        cur_dice = dice(seg_label == aa, prior > 0);
        if cur_dice > max_dice
            max_dice = cur_dice;
            max_dice_label = aa;
        end
    end
    seg = bin_to_double(seg_label == max_dice_label);
end
