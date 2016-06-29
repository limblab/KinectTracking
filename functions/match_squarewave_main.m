function [ kinect_times ] = match_squarewave_main( bdf, led_vals, times, kinect_start_guess, plot_flag)

%This function gets information about the cerebus squarewave, uses
%"match_squarewave_func" to find the timeshift between the kinect LED
%squarewave and the cerebus squarewave, and plots (if plot_flag=1).

%% Get information about Cerebus Squarewave

analog_ts = bdf.analog.ts';

% find KinectSyncPulse signal
pulse_idx = strcmp(bdf.analog.channel,'KinectSyncPulse');
KinectSyncPulse = bdf.analog.data(:,pulse_idx);

%Clean up the squarewave
KinectSquare=KinectSyncPulse>2000;

%Get the time (in Cerebus time) when the squarewave starts
start_ind=find(KinectSquare==1,1,'first');
cerebus_start=analog_ts(start_ind);


%% Find the time shift

[output,fval]=fminsearch(@(x) match_squarewave_func(x, led_vals, KinectSquare, times, analog_ts),[cerebus_start-kinect_start_guess 1/(2*mean(led_vals))]);

%% Update the kinect times

%Note, now the kinect times will now be called "kinect_times"

shift=output(1);
scale=output(2);
% scale=1/(2*mean(led_vals));

%Scale LED values
led_vals_scale=scale*led_vals;
%Shift times
kinect_times=times+shift;

%% Plot the resulting squarewave match (to make sure it was correct)

if plot_flag
    kinect_times_temp=kinect_times;
    
    %Remove kinect times from when the handle wasn't recorded
    kin_early=kinect_times_temp<1; %Find the kinect times from before the handle was recorded
    kin_late=kinect_times_temp>analog_ts(length(analog_ts)); %Find the kinect times from after the handle was recorded
    kinect_times_temp(kin_early | kin_late)=[];
    n_times=length(kinect_times_temp); %Number of kinect time points
    
    %Remove led vals from times when the handle wasn't recorded
    led_vals_trim=led_vals_scale;
    led_vals_trim(kin_early | kin_late)=[];
    
    %Get the indexes (for the cerebus times) that match with the kinect times (in the below sections)
    
    combined_times=[kinect_times_temp analog_ts']; %Vector with times from both cerebus and kinect
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
    
    %Plot kinect values in red and cerebus values in blue
    figure; plot(kinect_times_temp,led_vals_trim,'r');
    hold on;
    plot(analog_ts,KinectSquare);
    legend('Kinect','Cerebus')
    
end

end

