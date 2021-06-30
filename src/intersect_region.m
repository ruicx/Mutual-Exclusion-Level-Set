% name:       intersect_region.m
% usage:      --
% author:     Ruicheng
% date:       2020-11-06 13:57:29
% version:    1.0
% Env.:       MATLAB R2019b, WIN10


function seg = intersect_region(u, prior, threshold)
%intersect_region - select region by prior
%
% Syntax: seg = intersect_region(u, prior, threshold)
%
% select region by prior
    seg = false(size(u));
    prior = prior > 0;
    [img_label, n_label] = bwlabel(u > 0);
    for aa = 1:n_label
        region = (img_label == aa);
        isc_region = region & prior;
        if sum(isc_region, 'all') / sum(region, 'all') > 0.5
            seg = seg | region;
        end
    end
    seg = fill_holes(seg, threshold);
end

 