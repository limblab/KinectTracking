%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [md_smooth] = smoothMarkerData(marker_data)
%
%   This function smooths marker data with a median filter
%
% INPUTS:
%   marker_data  : struct containing marker data, in robot coordinates
%       Fields:
%           pos - position of markers, aligned spatially to match handle
%           movement
%           t - times for video frames of marker capture, aligned to
%           cerebus collection
%
% OUTPUTS:
%   md_smooth : struct of marker_data, smoothed
%       Fields:
%           pos - position of markers, smoothed
%           t - times for video frames of marker capture, aligned to
%           cerebus collection
%
%
% Written by Raeed Chowdhury and Joshua Glaser. Updated April 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function md_smooth = smoothMarkerData(marker_data)

md_smooth = marker_data;
md_smooth.pos=NaN(size(marker_data.pos));
for i=1:10
    for j=1:3
        temp=reshape(marker_data.pos(i,j,:),[1,size(marker_data.pos,3)]);
        md_smooth.pos(i,j,:)=medfilt1nan(temp,5);
    end
end