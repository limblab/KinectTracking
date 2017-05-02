%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [md_opensim] = transformForOpenSim(marker_data)
%
%   This function transforms marker data from lab coordinates into OpenSim
%   coordinates. This function is to be used in getTRCfromMarkers.
%
% INPUTS:
%   marker_data  : struct containing marker data, in robot coordinates
%       Fields:
%           pos - position of markers, aligned spatially to match handle
%           movement
%           t - times for video frames of marker capture, aligned to
%           cerebus collection
%   cds : CDS file to extract handle kinematics from
%
% OUTPUTS:
%   md_opensim : struct of marker_data in OpenSim coordinates
%       Fields:
%           pos - position of markers, aligned to OpenSim coordinates
%           t - times for video frames of marker capture, aligned to
%           cerebus collection
%   handle_opensim : position of the handle in OpenSim coordinates
%
%
% Written by Raeed Chowdhury and Joshua Glaser. Updated April 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [md_opensim,handle_opensim] = transformForOpenSim(marker_data,cds)
% Robot coordinates:
% Origin at shoulder joint center of robot, x is to right, y is towards screen, z is up
% OpenSim coordinates:
% Origin at shoulder joint center (marker 9), x is towards screen, y is up, z is to right

md_opensim = marker_data;

% extract and clean up shoulder jc data
shoulder_pos_lab = squeeze(marker_data.pos(9,:,:))';
marker_loss_points = find(diff(isnan(shoulder_pos_lab(:,1)))>0);
marker_reappear_points = find(diff(isnan(shoulder_pos_lab(:,1)))<0);
if(~isempty(marker_loss_points) || ~isempty(marker_reappear_points))
    if length(marker_loss_points)>length(marker_reappear_points) || marker_loss_points(end)>marker_reappear_points(end)
        %Means that a marker was lost but never found again at end of file
        %append last index
        marker_reappear_points(end+1) = length(shoulder_pos_lab);
    end
    if length(marker_loss_points)<length(marker_reappear_points) || marker_loss_points(1)>marker_reappear_points(1)
        %Means that a marker was lost to start
        %Dont know what to do here other than put a zero at start off loss
    %     marker_reappear_points(1) = [];
        marker_loss_points = [0;marker_loss_points];
        warning('Shoulder position marker not found at start of file. Replacing missing start points with first reappearance position')
    end
    for i=1:length(marker_loss_points)
        marker_loss = marker_loss_points(i);
        marker_reappear = marker_reappear_points(i);
        num_lost = marker_reappear-marker_loss;
        if marker_loss~=0
            rep_coord = repmat(shoulder_pos_lab(marker_loss,:),num_lost,1);
        else
            rep_coord = repmat(shoulder_pos_lab(marker_reappear,:),num_lost,1);
        end
        shoulder_pos_lab(marker_loss+1:marker_reappear,:) = rep_coord;
    end
end

% Recenter all markers on shoulder position
rep_shoulder_pos = repmat(shoulder_pos_lab,[1 1 11]);
rep_shoulder_pos = permute(rep_shoulder_pos,[3 2 1]);
kinect_pos_recenter = marker_data.pos-rep_shoulder_pos;

% switch axes of coordinate frame
% x->z, y->x, z->y
md_opensim.pos = kinect_pos_recenter(:,[2 3 1],:);
% md_opensim.pos(:,1,:) = kinect_pos_recenter(:,2,:); % new x=old y
% md_opensim.pos(:,2,:) = kinect_pos_recenter(:,3,:); % new y=old z
% md_opensim.pos(:,3,:) = kinect_pos_recenter(:,1,:); % new z=old x

% change from cm to meters
md_opensim.pos = md_opensim.pos/100;

% do same for handle position
handle_pos = [cds.kin.x cds.kin.y zeros(height(cds.kin),1)];
handle_pos=interp1(cds.kin.t,handle_pos,marker_data.t);
recenter_handle = handle_pos-shoulder_pos_lab;
handle_opensim =  recenter_handle(:,[2 3 1])/100;
