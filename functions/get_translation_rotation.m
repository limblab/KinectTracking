function [ R, Tpre, Tpost, times_good, pos_h, colors_xy ] = get_translation_rotation( bdf, kinect_times, all_medians, x_lim_handle, y_lim_handle, plot_flag )

%This function finds the translation and rotation that are necessary for
%converting the kinect markers into handle coordinates

%Inputs:
%bdf - cerebus file that contains handle information
%all_medians - kinect marker locations
%x_lim_handle, y_lim_handle - locations of handle to consider
%plot_flag - whether to plot the results

%Outputs:
%R - rotation matrix (to rotate the kinect positions)
%Tpre - translation of the kinect needed prior to rotation
%Tpost - translation of the kinect needed after rotation (to put back in
%handle coordinates)

%Note that the output times_good, pos_h, colors_xy, are all just
%used for making sure this works (via a plot). These can later be removed

%% 4a. Get handle information

handle_pos = bdf.pos(:,2:3);
handle_times = bdf.pos(:,1); %This should be the same as analog_ts

%% Get the handle times that correspond to the kinect times

%Note that there is a faster, but more complicated, way of doing this in the
%squarewave aligning function. Time isn't an issue here, so I haven't
%changed it.

n_times=size(all_medians,3); %Number of kinect time points
% 
% handle_time_idxs=zeros(length(n_times),1); %Initialize the handle time indexes
% %For each kinect time, find the closest handle time. handle_time_idxs has
% %the indexes of all the corresponding handle times.
% for i=1:n_times
%     [~,handle_time_idxs(i)]=min(abs(handle_times-kinect_times(i)));
% end
% 
% handle_times_ds=handle_times(handle_time_idxs); %These are the handle times that correspond to all the kinect times
% handle_pos_ds=handle_pos(handle_time_idxs,:); %The handle positions at all the kinect times

handle_pos_ds = interp1(handle_times,handle_pos,kinect_times);


%% Find time points that we need to remove (for several reasons)

%Remove kinect times from when the handle wasn't recorded
kin_early=kinect_times<1; %Find the kinect times from before the handle was recorded
kin_late=kinect_times>handle_times(length(handle_times)); %Find the kinect times from after the handle was recorded

marker=3; %This is the hand marker number that is about at the same point as the handle
kin_pos=reshape(all_medians(marker,:,:),[3,n_times]); %Positions of the above hand points (recorded by the kinect)

%Remove kinect times from when the hand is missing
missing=isnan(kin_pos(1,:)); %Times when hand points are missing

missing_smooth=smooth(single(missing),10); %Smooth when hand points are missing to find time areas when almost all hand points are there

%Find the good times that we should use for the alignment. Good times are
%those that:
%1. The hand is not missing
%2. The hand is not missing a lot around the time point (it's missing <20%
%of the time in the surrounding 10 frames)
%3. The handle was recorded (~kin_early & ~kin_late)
%4. The handle positions is within the specified limits
times_good=~missing' & missing_smooth<.2 & ~kin_early' & ~kin_late'...
    & handle_pos_ds(:,1)>x_lim_handle(1) & handle_pos_ds(:,1)<x_lim_handle(2) & handle_pos_ds(:,2)>y_lim_handle(1) & handle_pos_ds(:,2)<y_lim_handle(2);

%Get times and positions for the kinect and handle to be compared
kin_pos_good=kin_pos(:,times_good);
handle_pos_ds_good=handle_pos_ds(times_good,:);

%% Align

%Rename for simplicity
pos_k=kin_pos_good';
pos_h=handle_pos_ds_good;
pos_h(:,3)=0; %Set z-coordinate as 0 for handle

%X-coordinate of kinect is flipped: unflip it
pos_k(:,1)=-pos_k(:,1);

%Kinect is in meters, while handle is in cm: put everything in cm
pos_k=pos_k*100;

%Make mean of the kinect and handle positions zero (in order to be able to
%match them up via a rotation).
pos_k_shift=pos_k-repmat(mean(pos_k),[size(pos_k,1),1]);
pos_h_shift=pos_h-repmat(mean(pos_h),[size(pos_h,1),1]);

%Find the rotation that
[ pos_k_shift_rotate, R, q, fmin ] = rotate_match_func2( pos_k_shift, pos_h_shift );

Tpre=mean(pos_k); %The shift that was needed of the kinect positions prior to rotation
Tpost=mean(pos_h); %The shift that is needed after rotation to put back into handle coordinates
Tpost(3)=0; %Assume handle is at z=0

%% Plot to make sure it worked

if plot_flag
    
    %Make color map
    xVal=(pos_h(:,1)-min(pos_h(:,1)))/max(pos_h(:,1)-min(pos_h(:,1)));
    yVal=(pos_h(:,2)-min(pos_h(:,2)))/max(pos_h(:,2)-min(pos_h(:,2)));
    colors_xy = [.6*ones(size(xVal)),xVal,yVal];
    
    %Plot the kinect points that have been rotated (with coloring based on where
    %the handles initially were)
    figure; scatter3(pos_k_shift_rotate(:,1),pos_k_shift_rotate(:,2),pos_k_shift_rotate(:,3),[],colors_xy,'fill')
    title('Kinect')
    xlim([-15 15]);
    ylim([-15 15]);
    zlim([-15 15]);
    
    %Plot the handle points
    figure; scatter3(pos_h_shift(:,1),pos_h_shift(:,2),pos_h_shift(:,3),[],colors_xy,'fill')
    title('Handle')
    
end

end

