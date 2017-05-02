%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [kinect_times] = realignMarkerTimes(cerebus_vals,marker_data,plot_flag)
%
%   This function realigns the data collected by the kinect to the data
%   taken by the Cerebus system. This function is meant to be used in
%   realignMarkerSpacetime.
%
% INPUTS:
%   cerebus_vals : the analog cell from CDS containing KinectSyncPulse
%   marker_data  : struct containing marker data, drawn directly from the
%   color tracking script
%   plot_flag    : flag for whether to plot final alignment
%
% OUTPUTS:
%   kinect_times : vector of realigned times for marker data
%
%
% Written by Raeed Chowdhury and Joshua Glaser. Updated April 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kinect_times = realignMarkerTimes(cerebus_vals,marker_data,plot_flag)

times = marker_data.times;
led_vals = marker_data.led_vals;

h = figure;
plot(times,led_vals,'r');
title('Kinect LED vals')

%% 3b. Enter kinect start time estimate
str1='Enter a guess for kinect start time: \n';
kinect_start_guess=input(str1);
%Make sure entry is valid (either empty, or composed of markers)
while ~( isempty(kinect_start_guess) || ( all(kinect_start_guess>0) && all(kinect_start_guess<=times(end)) ) )
    kinect_start_guess=input('Enter valid time');
end
% kinect_start_guess=19.5;

close(h)

%% 3c. Align kinect led values with cerebus squarewave
kinect_times  = match_squarewave_main( cerebus_vals, led_vals, times, kinect_start_guess, plot_flag);