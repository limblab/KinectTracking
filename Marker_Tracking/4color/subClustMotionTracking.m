function [ cluster1] = subClustMotionTracking( cluster, influence, options )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    for i=1:length(cluster)
        if(~isempty(cluster{i}))
            cluster1{i} = subclust(cluster{i}, influence,[],  options);
        end
    end

end

