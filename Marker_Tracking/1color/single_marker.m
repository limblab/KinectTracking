%This script plots the original points
%and can be used for inputting the starting locations of the markers

%% All colors at once

n=length(color1);
xlims=[-.4 .5];
ylims=[-.4 .4]; 
zlims=[.8 1.5];

figure;

xlim(zlims)
ylim(xlims)
zlim(ylims)
set(gca,'NextPlot','replacechildren');


for i=300:2000
    temp=color3{i};
    x=temp(1:end/3);
    y=temp(end/3+1:2*end/3);
    z=temp(2*end/3+1:end);
    scatter3(z,x,y,'b')
    hold off
    title(i);
    xlabel('z')
    ylabel('x')
    zlabel('y')
    xlim(zlims)
    ylim(xlims)
    zlim(ylims)
%     pause(.03)
    pause;
end



%%

marker_inits= [1,0,-.08];
marker_inits_temp=marker_inits;
marker_inits(:,1)=marker_inits_temp(:,2);
marker_inits(:,2)=marker_inits_temp(:,3);
marker_inits(:,3)=marker_inits_temp(:,1);


%% Plotting Initializations

figure;
set(gca,'NextPlot','replacechildren');
xlims=[-.5 .5];
ylims=[-.5 .5];
zlims=[0.5 1.5];

pause_time=.03;

%% Blue Arm

%Initializations
start=300;
finish=length(color3);
n_times=finish-start+1;

plot_on=0; %Whether to plot while it's running
marker_ids=1; %Set the marker_ids specified in "Initializations"
color=color3; %Blue=1, Green=2, Red=3
prev_meds=marker_inits(marker_ids,:); %Set initial "previous marker locations" as the start locations input in "Initializations"
num_clust=length(marker_ids); %Number of clusters
within_clust_dist1=0.07; %How close points must be to the previous frame's marker to be considered
dist_min=0.07; %Minimum distance between markers (cluster medians aren't allowed w/ distance < min_dist)

medians=NaN(num_clust,3,n_times); %Has NaNs when a marker is missing
medians2=NaN(num_clust,3,n_times); %Has previous known positions when a marker is missing


% LOOP THROUGH TIME
t=0;
for i=start:finish
    
    t=t+1;
    
    %0. Get x,y,z positions
    temp=color{i};
    x=temp(1:end/3);
    y=temp(end/3+1:2*end/3);
    z=temp(2*end/3+1:end);
    loc=[x; y; z]';
    
    %1. Filter some bad points (those that are really far away)    
    %Get distances of all points to the marker in the previous frame
    if t==1
        D=pdist2(loc,prev_meds);
        prev_num_clust=num_clust;
    else
        D=pdist2(loc,medians2(:,:,t-1));
    end
    
    % Keep all the points close enough to the previous marker
    keep1=D(:,1)<within_clust_dist1;
    
    %Remove points (those we're not keeping)
    rmv=~(keep1);
    
    %Actually remove the points
    loc(rmv,:)=[];
        
    %2. Cluster and assign
    [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
    
    %3. Plot original image and cluster centers
    plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
    
end

%% Make matrices line up correctly with start time

marker=3;

%Make the final matrices begin at time 1 instead of time "start."

temp=medians;
temp2=medians2;

all_medians=NaN(11,3,finish);
all_medians2=NaN(11,3,finish);

all_medians(3,:,start:finish)=temp;
all_medians2(3,:,start:finish)=temp2;

%% Save

save('accelero_marker.mat','all_medians','all_medians2','led_vals','times');
