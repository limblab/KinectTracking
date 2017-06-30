function [ kinect_pos ] = do_affine_xform( all_medians, affine_xform, plot_flag )

%Puts the kinect markers into handle coordinates

%Inputs:
%R - rotation matrix (to rotate the kinect positions)
%Tpre - translation of the kinect needed prior to rotation
%Tpost - translation of the kinect needed after rotation (to put back in
%handle coordinates)
%all_medians - the kinect marker locations in kinect coordinates

%Outputs:
%kinect_pos - the kinect marker locations in handle coordinates

%Note that the inputs plot_flag, times_good, pos_h, colors_xy, are all just
%used for making sure this works (via a plot). These can later be removed

%% Deal with number of inputs

if nargin<3 %The final parameters are only for plotting
    plot_flag=0;
end


%% Put all kinect positions in handle coordinates

n_times = size(all_medians,3);

%Rename kinect positions
all_medians_v2=all_medians;

%%X-coordinate of kinect is flipped: unflip it
all_medians_v2(:,1,:)=-all_medians_v2(:,1,:);

%Kinect is in meters, while handle is in cm
all_medians_v2=100*all_medians_v2;

%These are the matrices with the kinect positions in handle coordinates
kinect_pos=NaN(size(all_medians_v2));

%Loop through the markers. For every marker, multiply the shited kinect position
%by the affine matrix.
for m=1:size(all_medians,1)
    m_pos=reshape(all_medians_v2(m,:,:),[3,n_times]);
    %switch to homogenous coordinates
    m_pos_prime = affine_xform * [m_pos; ones(1,size(m_pos,2))];
    %switch back to regular coordinates
    kinect_pos(m,:,:)=m_pos_prime(1:3,:);
end

%Possibly flip z-coordinate (since the rotation doesn't know up from down)
%Based on where we put the kinect, the z-coordinate (in handle coordinates) of marker 8 should be positive. Flip z's otherwise.
temp=nanmean(kinect_pos(8,:,:),3);
if temp(3)<0
    kinect_pos(:,3,:)=-kinect_pos(:,3,:);
end

%% Another plot test

if plot_flag
    colors = {'g','b','r','y','g','g','b','r','g','r'};
    figure;
    for ctr = 1:10
        scatter3(kinect_pos(ctr,1,:),kinect_pos(ctr,2,:),kinect_pos(ctr,3,:),[],colors{ctr});
        hold on
    end
    axis equal
end

end

