function [ sse ] = match_squarewave_func( x, led_vals, KinectSquare, times, handle_times )

%This function finds the sse between the shifted/scaled LED values from the
%kinect and from the cerebus.

shift=x(1);
scale=x(2);
% scale=1/(2*mean(led_vals));

%Scale LED values
led_vals_scale=scale*led_vals;
%Shift times
kinect_times=times+shift;

%Remove kinect times from when the handle wasn't recorded
kin_early=kinect_times<1; %Find the kinect times from before the handle was recorded
kin_late=kinect_times>handle_times(length(handle_times)); %Find the kinect times from after the handle was recorded
kinect_times(kin_early | kin_late)=[];
n_times=length(kinect_times); %Number of kinect time points

%Remove led vals from times when the handle wasn't recorded
led_vals_trim=led_vals_scale;
led_vals_trim(kin_early | kin_late)=[];

%Get the indexes (for the cerebus times) that match with the kinect times (in the below sections). 

combined_times=[kinect_times handle_times']; %Vector with times from both cerebus and kinect
[~,idx]=sort(combined_times); %Sort the times
[~,idx_rev]=sort(idx); %idx_rev says what place each time was in the sorted combined vector
a=idx_rev(1:n_times)-1; %Get the sorted location of each of the kinect times (the values from 1:n_times) in the combined vector.
%Then subtract 1 to get the indices of the handle_times with times closest
%(just below) the kinect times.

%Find where these values are in the combined time vector
c = ismember(idx_rev, a); 
indexes = find(c);
handle_time_idxs=indexes-n_times; %Subtract n_times since the first n_times values in the vector are from the kienct_times, and we want which one of the handle_times.

%Get the value of the voltage (corresponding to LED value) of the cerebus
%at the times found above
V=KinectSquare(handle_time_idxs);

%Get the SSE
sse=sum((V'-led_vals_trim).^2);

end

