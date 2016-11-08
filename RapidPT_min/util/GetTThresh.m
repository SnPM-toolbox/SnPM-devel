function [ t_threshold ] = GetTThresh(Distribution, alpha)
%getPVal Summary of this function goes here
%   Detailed explanation goes here

    t_threshold = prctile(Distribution, 100 - alpha);

end

