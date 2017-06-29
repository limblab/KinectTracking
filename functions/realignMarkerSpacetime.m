%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [kinect_pos,kinect_rotation] = realignMarkerSpace(cds,marker_data)
%
% OR
%
% function [kinect_pos,kinect_rotation] = realignMarkerSpace(cds,marker_data,alignment_settings)
%
%   This function realigns the data collected by the kinect to the data
%   taken by the Cerebus system. This function is meant to be used in
%   getTRCfromMarkers.
%
% INPUTS:
%   cds : cds with collected cerebus data
%   marker_data  : struct containing marker data, drawn directly from the
%   color tracking script
%   alignment_settings: optional argument that allows for using the same
%   spatial transformation to be applied to multiple files from the same
%   day
%
% OUTPUTS:
%   marker_data_aligned : struct of aligned marker data
%       Fields:
%           pos - position of markers, aligned spatially to match handle
%           movement
%           t - times for video frames of marker capture, aligned to
%           cerebus collection
%   alignment_settings : settings used for spatial transformation
%
%
% Written by Raeed Chowdhury and Joshua Glaser. Updated April 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [marker_data_aligned,alignment_settings] = realignMarkerSpacetime(cds,marker_data,varargin)

% figure out if rotation known
if nargin>3
    error('Too many arguments')
elseif nargin==3
    rotation_known = true;
    alignment_settings = varargin{1};
else
    rotation_known = false;
end

% First time align based on square wave
kinect_times = realignMarkerTimes(cds.analog{1},marker_data,true);

% 4a. Plot handle to determine some points to remove
%We want to remove the time points when the monkey has thrown away the
%handle, since then the hand won't be at the same position as the handle
if ~rotation_known
    figure; scatter(cds.kin.x(1:10:end),cds.kin.y(1:10:end))
    %Note- this plot can be removed if the limits (below) are always the same
    
    str2='Type in the handle limits on the x-axis \n';
    str3='e.g. [-10,10] (default is [-15,15]) \n';
    x_lim_handle=input([str2 str3]);
    
    if isempty(x_lim_handle)
        x_lim_handle=[-15,15]; %x limits (min and max)
    end
    
    str2='Type in the handle limits on the y-axis \n';
    str3='e.g. [-10,10] (default is [-50,-20]) \n';
    y_lim_handle=input([str2 str3]);

    if isempty(y_lim_handle)
        y_lim_handle=[-50,-20]; %y limits (min and max)
    end
    
    % find alignment
    plot_flag=1;
    [ affine_xform, times_good, pos_h, colors_xy ] = get_affine_xform( cds, kinect_times, marker_data.all_medians, x_lim_handle, y_lim_handle, plot_flag );
%     [ alignment_settings, times_good, pos_h, colors_xy ] = get_translation_rotation( cds, kinect_times, marker_data.all_medians, x_lim_handle, y_lim_handle, plot_flag );
    %Save a file w/ T and R, so it can be used for other files from the
    %same day
%     save([folder prefix '_kinect_rotation.mat'],'R','Tpre','Tpost')
    
% else
%     %Else load a file that has T and R
%     load([folder prefix '_kinect_rotation.mat']);
end

% 4d. Perform Translation and Rotation on the kinect data
% if ~rotation_known
%     plot_flag=1;
%     [ kinect_pos ] = do_translation_rotation( all_medians, R, Tpre, Tpost, plot_flag, times_good, pos_h, colors_xy );
% end

% if rotation_known
    [ kinect_pos ] = do_translation_rotation( marker_data.all_medians, alignment_settings);
% end

marker_data_aligned = struct('pos',kinect_pos,'t',kinect_times);