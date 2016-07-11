%Script to track the color markers over time

%A file needs to be loaded which contains the color pixels of each frame

%For the main tracking portions of the script (e.g. blue arm, red hand,
%etc.), there are more detailed comments in the "Blue Arm" section which
%comes first. For the subsequent sections, only unique aspects are
%commented on in detail.

%% 1. INITIALIZATIONS %%%%%%%%%%%%%%%%%%%%

%% Input File to Load

main_dir='/Users/jig289/Box Sync/Tracking_Data/';
%File to load
monkey='Chips';
date='12-10-15'; %mo-day-yr
exp='RW_PM';
num='002';

fname_load=ls([main_dir monkey '/Color_Tracking/' date '/Tracking/color_tracking ' exp '_' num '*']);
load(deblank(fname_load));


%% User Options / Initializations

first_time=1; %If this is the first file from a date, set equal to 1 (there are more initializations)

%Save the output?
savefile=0;

%Use default values for parameters?
use_defaults=0;

%Load all of the settings if it's not the first file 
if ~first_time
    date2=['20' num2str(date(7:8)) num2str(date(1:2)) num2str(date(4:5))];
    fname_load_settings=[main_dir monkey '/Color_Tracking/' date '/Markers/settings_' monkey '_' date2];
    load(fname_load_settings);
end

%TIME INITIALIZATIONS
start=536; %Time point we're starting at
n=length(color1);
finish=n; %Time point we're finishing at
n_times=finish-start+1; %Number of time points (frames)

n_times_prelim=n_times;

%Plot figure of schematic of marker locations?
plot_marker_locs=1;

%Manually type in initial marker locations? (Otherwise you get to click on
%them in the figure)
marker_init_manual=0;


%MARKER NUMBER INITIALIZATIONS
red_arm_marker_ids=[8,10];
blue_arm_marker_ids=[7];
green_shoulder_marker_ids=[9]; %Sometimes empty
green_elbow_marker_ids=[6];
red_hand_marker_ids=[3,4];
blue_hand_marker_ids=[2]; %Sometimes just 2, sometimes 2 and 11
green_hand_marker_ids=[1,5];

%% Plotting Initializations (Also user input)

plot_during=0; %If you'll be displaying the marker tracking while it's running at any point (only used during testing)

%The x and y limits matter for selecting the markers in the initial frame.
%The z limit will only matter if plot_during=1
xlims=[-.5 .5];
ylims=[-.5 .5];
zlims=[0.5 1.5];

pause_time=.03;

if plot_during
    figure;
    set(gca,'NextPlot','replacechildren');
end

%% Initializations of vectors/matrices

%Keeps track of all the cluster locations
all_medians=NaN(11,3,n_times); %Has NaNs when a marker is missing
all_medians2=NaN(11,3,n_times); %Estimates where markers are when they are missing

%Initialize some vectors that I use later for calculating the distance
%between points
dists=NaN(1,n_times_prelim);
dists1=NaN(1,n_times_prelim);
dists2=NaN(1,n_times_prelim);
dists3=NaN(1,n_times_prelim);
dists4=NaN(1,n_times_prelim);
dists5=NaN(1,n_times_prelim);

%% 2. SELECT THE INITIAL MARKER LOCATIONS IN THE "START" FRAME


%% Marker location schematic figure

marker_demo_locs=[0 0; 1 1; 1 -1; 2 1; 2 -1;...
    10 -1; 10 3; 10 6; 10 9; 9 0;...
    2 -3];
r=[1 0 0];
g=[0 1 0];
b=[0 1 1];
marker_demo_colors=[g; b; r; r; g; g; b; r; g; r; b];

if plot_marker_locs
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(1,2,1);
    scatter(marker_demo_locs(:,1),marker_demo_locs(:,2),200,marker_demo_colors,'filled');
    str={'1','2','3','4','5','6','7','8','9','10','11'};
    text(marker_demo_locs(:,1),marker_demo_locs(:,2),str)
    xlim([-5 15]);
    ylim([-5 15]);
end

%% Marker location initializations (interactive)

%We here set the locations of the markers in the "start" rame

if ~marker_init_manual
    
    marker_colors={'g','b','r','r','g','g','b','r','g','r'}; %The colors of each of our markers
    
    num_markers=10;
    marker_coords_xy=NaN(num_markers,2);
    
    if plot_marker_locs
        subplot(1,2,2);
    else
        figure
    end
    %Get x,y,z coordinates for points in all colors
    temp=color1{start};
    x1=temp(1:end/3);
    y1=temp(end/3+1:2*end/3);
    z1=temp(2*end/3+1:end);
    hold on;
    temp=color2{start};
    x2=temp(1:end/3);
    y2=temp(end/3+1:2*end/3);
    z2=temp(2*end/3+1:end);
    temp=color3{start};
    x3=temp(1:end/3);
    y3=temp(end/3+1:2*end/3);
    z3=temp(2*end/3+1:end);
    
    %Remove red cord in background
    keep=z3<2;   
    x3_keep=x3(keep);
    y3_keep=y3(keep);
    z3_keep=z3(keep);
    
    %Plot all the points in the x/y plane (z, which is depth, doesn't change
    %much between the points)
    scatter(x1,y1,'b')
    hold on;
    scatter(x2,y2,'g')
    scatter(x3_keep,y3_keep,'r')
    hold off
    xlabel('x')
    ylabel('y')
    xlim(xlims)
    ylim(ylims)
    
    %Have users select the markers
    for m=1:num_markers
        title(['Click marker ' num2str(m)])
        marker_coords_xy(m,:)=ginput(1);
    end
    
    %Get the 3d marker locations. To do so, we find the point (of the appropriate color) with the
    %closest x/y coordinate. We then get the x/y/z coordinates of that point.
    marker_inits=NaN(11,3); %Made large enough for an 11th marker (which we used at one point)
    for m=1:num_markers
        if marker_colors{1}=='r'
            closest_point=knnsearch([x3_keep' y3_keep'],marker_coords_xy(m,:),'k',1);
            marker_inits(m,:)=[x3_keep(closest_point) y3_keep(closest_point) z3_keep(closest_point)];
        end
        if marker_colors{1}=='g'
            closest_point=knnsearch([x2' y2'],marker_coords_xy(m,:),'k',1);
            marker_inits(m,:)=[x2(closest_point) y2(closest_point) z2(closest_point)];
        end
        if marker_colors{1}=='b'
            closest_point=knnsearch([x1' y1'],marker_coords_xy(m,:),'k',1);
            marker_inits(m,:)=[x1(closest_point) y1(closest_point) z1(closest_point)];
        end
    end
    
end

%% Marker location initializations (if you'd prefer to type it in)

%If you'd prefer to type in the initial marker locations ahead of time, you
%can do it here. You can get them from orioginal_colos_plot_4colors

if marker_init_manual
    marker_inits=NaN(11,3);
    marker_inits(1,:)=[.9,-.06,-.15];
    marker_inits(2,:)=[.9,-.04,-.15];
    marker_inits(3,:)=[.9,-.05,-.16];
    marker_inits(4,:)=[.9,-.03,-.14];
    marker_inits(5,:)=[.9,-.02,-.16];
    marker_inits(6,:)=[1,.17,-.16];
    marker_inits(7,:)=[1,.16,-.11];
    marker_inits(8,:)=[1,.15,-.06];
    marker_inits(9,:)=[1,.15,0.0];
    marker_inits(10,:)=[1,.15,-.15];
    
    %I plot z,x,y (instead of x,y,z), so I input z,x,y above. Here, switch to x,y,z
    marker_inits_temp=marker_inits;
    marker_inits(:,1)=marker_inits_temp(:,2);
    marker_inits(:,2)=marker_inits_temp(:,3);
    marker_inits(:,3)=marker_inits_temp(:,1);
end

%% 3. TRACK ARM MARKERS (GET THE LOCATION OF ARM MARKERS)

%% Blue Arm

%Initializations
plot_on=0; %Whether to plot while it's running
marker_ids=blue_arm_marker_ids; %Set the marker_ids specified in "Initializations"
color=color1; %Blue=1, Green=2, Red=3
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

%Put the markers found here in the matrix of all markers
all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;

%% Red Arm (Preliminary)

if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Initializations
    plot_on=0;
    marker_ids=red_arm_marker_ids;
    color=color3;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids); 
    within_clust_dist1=.07; %How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
    within_clust_dist2=.07; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered
    dist_min=0.07; 
    
    medians=NaN(num_clust,3,n_times);
    medians2=NaN(num_clust,3,n_times);
        
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
        
        % Keep all the points close enough to either of the previous markers
        keep1=D(:,1)<within_clust_dist1;
        keep2=D(:,2)<within_clust_dist2;       
        
        %Remove points (those we're not keeping)
        rmv=~(keep1 | keep2);
        
        %Actually remove the points
        loc(rmv,:)=[];
        
        
        %2. Cluster and assign
        %Note that this uses "cluster_func" instead of "cluster_func2"
        %which is slightly faster but less accurate. This is because we
        %will be redoing this later with cluster_func2. This current run is
        %only to determine the distances from the red arm to blue arm
        %markers (which will help in the next run)
        [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func(t, loc, num_clust, prev_num_clust, dist_min, .05, prev_meds, medians, medians2 );
        
        %3. Plot original image and cluster centers
        plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
        
    end
    
    %Put the markers found here in the matrix of all markers
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
    
end

%% Set limits on red arm to blue arm distances

%PLOT RED ARM TO BLUE ARM DISTANCES
if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Calculate distances for each time point
    for i=1:n_times_prelim
        dists(i)=pdist2(all_medians(10,:,i),all_medians(7,:,i)); %Distance between markers 7 and 10 (blue arm and red elbow)
        dists2(i)=pdist2(all_medians(8,:,i),all_medians(7,:,i)); %Distance between markers 7 and 8 (blue arm and red arm)
    end
    
    if use_defaults %User defaults
        red_elbow_dist_from_blue=nanmean(dists)+4*nanstd(dists);
        red_blue_arm_dist_max=nanmean(dists2)+4*nanstd(dists2);
        
    else %If not using defaults
        
        %Plot
        figure; plot(dists);
        hold on;
        plot(dists2)
        legend('7-10','7-8')
        title('Distance between blue arm marker (7) and red arm markers (8 and 10)');
        
        %VISUALIZE FRAMES
        user_input=1; %A value so that it enters the while loop below
        while ~isempty(user_input)
            str1='Next, you will set maximum values of distances between the blue and red arm markers \n';
            str2='First, to help make this decision, enter time point you want to visualize (or just press enter to continue) \n';
            user_input=input([str1 str2]);
            %Make sure the input was valid (an integer between start and finish)
            if ~isempty(user_input)
                while ~(isnumeric(user_input) && mod(user_input,1)==0 && user_input>=start && user_input<=finish)
                    user_input=input('Re-enter valid time point \n');
                end
            end
            if ~isempty(user_input)
                plot_together_3colors_func(user_input, [7 8 10], [1:10], all_medians, color1, color2, color3,  start, finish, 1)
            end
        end
        
        % SET RED ARM TO BLUE ARM DISTANCES
        str1='Input red_elbow_dist_from_blue \n';
        str2='The blue values in the above plot should be generally be below this value (the red elbow should be within this distance of the blue arm)\n';
        str3='The purpose of this is to keep all points w/in this distance of the blue as marker candidates (useful if the red elbow marker was gone the previous frame) \n';
        str4='Value is generally ~ .05-.1 \n';
        red_elbow_dist_from_blue=input([str1 str2 str3 str4]);
        %Make sure it's a valid entry
        while ~(isnumeric(red_elbow_dist_from_blue))
            red_elbow_dist_from_blue=input('Re-enter valid value');
        end
        
        
        str1='Input red_blue_arm_dist_max \n';
        str2='All values in above plot should be below this value (Maximum distance from a red arm point to the blue)\n';
        str3='The purpose of this is to remove all points farther than this from the blue marker (to get rid of noise)\n';
        str4='Value is generally ~ .05-.1 \n';
        red_blue_arm_dist_max=input([str1 str2 str3 str4]);
        %Make sure it's a valid entry
        while ~(isnumeric(red_blue_arm_dist_max))
            red_blue_arm_dist_max=input('Re-enter valid value');
        end
        
    end
end


%% Red Arm (Redo)
%Note that this is different from the previous version of "Red Arm" because
%now there are constraints involving distance from the blue arm marker

%Initializations
plot_on=0;
marker_ids=red_arm_marker_ids;
color=color3;
prev_meds=marker_inits(marker_ids,:);
num_clust=length(marker_ids); %Number of clusters
within_clust_dist1=.07; %How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
within_clust_dist2=.07; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered   
dist_min=0.06; %Minimum distance between markers (cluster medians aren't allowed w/ distance < min_dist)

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
    
    %Get distance to blue arm marker from the current frame
    D2=pdist2(loc,all_medians2(7,:,t)); 
    
    % Keep all the points close enough to either of the previous markers
    keep1=D(:,1)<within_clust_dist1;
    keep2=D(:,2)<within_clust_dist2;
    
    %Also keep if the it's near the blue arm marker (in case one of the
    %others disappears for a while)    
    keep3=D2<red_elbow_dist_from_blue;
    
    %Remove points that are too far from the blue marker
    rmv0=D2>red_blue_arm_dist_max;
    
    %Remove points (those we're not keeping, or those we're removing)
    rmv=~(keep1 | keep2 | keep3) | rmv0;
  
    %Actually remove the points
    loc(rmv,:)=[];
       
    %2. Cluster and assign
    [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
    
    %3. Plot original image and cluster centers
    plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
    
end

all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;

%% Remove faulty red elbow points

%This calculates (and plots) the angle made by points 7,8,10
%Problems with the red elbow marker (point 10) will make this angle wrong
%We will remove those points

angle=NaN(1,n_times); %Initialize vector of angles for each frame

for i=1:n_times
    
    if all(~isnan(all_medians([7 8 10],1,i))) %Only find angle for frames when all markers 7/8/10 are present
        a=all_medians(10,:,i);
        b=all_medians(7,:,i);
        c=all_medians2(8,:,i);
        
        u=a-b; %Vector from 10 to 7
        v=c-b; %Vector from 8 to 7
        
        angle(i)=acos(dot(u,v)/norm(u)/norm(v)); %Angle made by 10,7,8
    end
end

%Plot
figure; plot(angle)
red_elbow_angle_thresh=nanmean(angle)-4*nanstd(angle); %Frames with an angle below this will have marker 10 removed (default threshold)
title(['Red Elbow Angles: Default Threshold=' num2str(red_elbow_angle_thresh)]);


if ~use_defaults %Default was set above: nanmean(angle)-4*nanstd(angle). If you're not using the default:
    
    %VISUALIZE FRAMES
    user_input=1; %A value so that it enters the while loop below
    while ~isempty(user_input)
        str1='Next, you will enter a minimum threshold for the red elbow angle \n';
        str2='First, to decide how to set this threshold, you can enter a time point you want to visualize (or just press enter to continue) \n';
        user_input=input([str1 str2]);
        %Make sure the input was valid (an integer between start and finish)
        if ~isempty(user_input)
            while ~(isnumeric(user_input) && mod(user_input,1)==0 && user_input>=start && user_input<=finish)
                user_input=input('Re-enter valid time point \n');
            end
        end
        if ~isempty(user_input)
            plot_together_3colors_func(user_input, [7 8 10], [1:10], all_medians, color1, color2, color3, start, finish, 1)
        end
    end
    
    %SET ANGLE THRESHOLD FOR RED ELBOW REMOVAL
    str1='Enter angle threshold for red elbow removal. Press enter for default. \n';
    temp=input(str1);
    %Make sure it's a valid entry
    while ~(isnumeric(temp) || isempty(temp))
        temp=input('Re-enter valid value');
    end
    
    if ~isempty(temp)
        red_elbow_angle_thresh=temp;
    end
    
end

%Remove red elbow points (based on angle)
rmv10=angle<red_elbow_angle_thresh;
all_medians(10,:,rmv10)=NaN;

%% Green Shoulder

if ~isempty(green_shoulder_marker_ids) %Only do this if it is a file with a green shoulder marker
    %Initializations
    plot_on=0;
    marker_ids=green_shoulder_marker_ids;
    color=color2;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids);
    within_clust_dist1=.07;
    dist_min=0.07; 
    
    medians=NaN(num_clust,3,n_times); 
    medians2=NaN(num_clust,3,n_times);
    
    
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
        
        if t==1
            D=pdist2(loc,prev_meds);
            prev_num_clust=num_clust;
        else
            D=pdist2(loc,medians2(:,:,t-1));
        end
        
        % Keep all the points close enough to the previous marker
        keep1=D(:,1)<within_clust_dist1;
        
        % Remove
        rmv=~(keep1);
        loc(rmv,:)=[];
        
        
        %2. Cluster and assign
        [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
        
        %3. Plot original image and cluster centers
        plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
        
    end
    
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
end


%% Green Elbow

%Initializations
plot_on=0;
marker_ids=green_elbow_marker_ids;
color=color2;
prev_meds=marker_inits(marker_ids,:);
num_clust=length(marker_ids);
within_clust_dist1=.07;
dist_min=0.07;

medians=NaN(num_clust,3,n_times);
medians2=NaN(num_clust,3,n_times);


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
    
    if t==1
        D=pdist2(loc,prev_meds);
        prev_num_clust=num_clust;
    else
        D=pdist2(loc,medians2(:,:,t-1));
    end
    
    % Keep all the points close enough to the previous marker
    keep1=D(:,1)<within_clust_dist1;
    
    %Also use distance from red elbow marker
    D2=pdist2(loc,all_medians(10,:,t));   
    keep2=D2<.03; %Keep points that are close to red elbow marker        
    rmv0=D2>.03; %Remove points that are too far from red elbow marker
   
    rmv=~(keep1|keep2) | rmv0; %Keep points that are either close enough to 
    %the previous marker or the red elbow marker. Additionally, remove 
    %points that are too far from the red elbow marker (even if they're 
    %close enough to the previous marker)

    %Remove
    loc(rmv,:)=[];
        
    %2. Cluster and assign
    [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
    
    %3. Plot original image and cluster centers
    plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
    
end

all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;

%% Remove faulty green elbow points

%This calculates (and plots) the angle made by points 7,8,6
%Problems with the green elbow marker (point 6) will make this angle wrong
%We will remove those points

angle=NaN(1,n_times); %Initialize vector of angles for each frame

for i=1:n_times
    
    if all(~isnan(all_medians([7 8 6],1,i))) %Only find angle for frames when all markers 6/7/8 are present
        
        a=all_medians(6,:,i);
        b=all_medians(7,:,i);
        c=all_medians2(8,:,i);
        
        u=a-b; %Vector from 6 to 7
        v=c-b; %Vector from 8 to 7
        
        angle(i)=acos(dot(u,v)/norm(u)/norm(v)); %Angle made by 6,7,8
    end
end

%Plot
figure; plot(angle)
green_elbow_angle_thresh=nanmean(angle)-4*nanstd(angle); %Frames with an angle below this will have marker 6 removed
title(['Green Elbow Angles: Default Threshold=' num2str(green_elbow_angle_thresh)]);

if ~use_defaults %Default was set above: nanmean(angle)-4*nanstd(angle). If you're not using the default:
    
    %VISUALIZE FRAMES
    user_input=1; %A value so that it enters the while loop below
    while ~isempty(user_input)
        str1='Next, you will enter a minimum threshold for the green elbow angle \n';
        str2='First, to decide how to set this threshold, you can enter a time point you want to visualize (or just press enter to continue) \n';
        user_input=input([str1 str2]);
        %Make sure the input was valid (an integer between start and finish)
        if ~isempty(user_input)
            while ~(isnumeric(user_input) && mod(user_input,1)==0 && user_input>=start && user_input<=finish)
                user_input=input('Re-enter valid time point \n');
            end
        end
        if ~isempty(user_input)
            plot_together_3colors_func(user_input, [6 7 8], [1:10], all_medians, color1, color2, color3, start, finish, 1)
        end
    end
    
    % SET ANGLE THRESHOLD FOR GREEN ELBOW REMOVAL
    str1='Enter angle threshold for green elbow removal. Press enter for default. \n';
    temp=input(str1);
    while ~(isnumeric(temp) || isempty(temp))
        temp=input('Re-enter valid value');
    end
    if ~isempty(temp)
        green_elbow_angle_thresh=temp;
    end
    
end

% Remove green elbow points (based on angle)
rmv6=angle<green_elbow_angle_thresh;
all_medians(6,:,rmv6)=NaN;

%% 4. PRELIMINARY TRACKING OF HAND MARKERS (IN ORDER TO SET DISTANCE CONSTRAINTS TO ARM MARKERS)

%% Red Hand (Preliminary)
if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Initializations
    plot_on=0;
    marker_ids=red_hand_marker_ids;
    color=color3;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids); 
    within_clust_dist1=.05;  %How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
    within_clust_dist2=.05; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered
    dist_min=0.02;
    
    medians=NaN(num_clust,3,n_times); 
    medians2=NaN(num_clust,3,n_times); 
    
    
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
        
        % Keep all the points close enough to either of the previous markers
        keep1=D(:,1)<within_clust_dist1;
        keep2=D(:,2)<within_clust_dist2;
        
        %Remove points (those we're not keeping)
        rmv=~(keep1 | keep2);
        
        %Actually remove the points
        loc(rmv,:)=[];
        
        
        %2. Cluster and assign
        %Note that this uses "cluster_func" instead of "cluster_func2"
        %which is slightly faster but less accurate. This is because we
        %will be redoing this later with cluster_func2. This current run is
        %only to determine the distances from the red hand markers to arm
        %markers (which will help in the next run)
        [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func(t, loc, num_clust, prev_num_clust, dist_min, .05, prev_meds, medians, medians2 );
               
        %3. Plot original image and cluster centers
        plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
        
    end
    
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
    
end

%% Green Hand (Preliminary)
if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Initializations
    plot_on=0;
    marker_ids=green_hand_marker_ids;
    color=color2;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids);
    within_clust_dist1=.07; %How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
    within_clust_dist2=.07; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered
    dist_min=0.03; 
    
    medians=NaN(num_clust,3,n_times); 
    medians2=NaN(num_clust,3,n_times); 
    
    
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
        % Keep all the points close enough to either of the previous markers
        keep1=D(:,1)<within_clust_dist1;
        keep2=D(:,2)<within_clust_dist2;
        
        %Remove points (those we're not keeping)
        rmv=~(keep1 | keep2);
        
        %Actually remove the points
        loc(rmv,:)=[];
        
        
        %2. Cluster and assign
        %Note that this uses "cluster_func" instead of "cluster_func2"
        %which is slightly faster but less accurate. This is because we
        %will be redoing this later with cluster_func2. This current run is
        %only to determine the distances from the green hand to arm
        %markers (which will help in the next run)
        [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func(t, loc, num_clust, prev_num_clust, dist_min, .05, prev_meds, medians, medians2 );
        
        %3. Plot original image and cluster centers
        plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
        
    end
    
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
    
end

%% Blue Hand (Preliminary)
if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Initializations
    plot_on=0;
    marker_ids=blue_hand_marker_ids;
    color=color1;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids);
    within_clust_dist1=.07;
    dist_min=0.07;
    
    medians=NaN(num_clust,3,n_times); 
    medians2=NaN(num_clust,3,n_times); 
    
    
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
        
        %Keep points close enough to previous marker
        keep1=D(:,1)<within_clust_dist1;        
        rmv=~(keep1);
        
        %Actually remove
        loc(rmv,:)=[];
        
        
        %2. Cluster and assign
        %Note that this uses "cluster_func" instead of "cluster_func2"
        %which is slightly faster but less accurate. This is because we
        %will be redoing this later with cluster_func2. This current run is
        %only to determine the distances from the blue hand to arm
        %markers (which will help in the next run)
        [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func(t, loc, num_clust, prev_num_clust, dist_min, .05, prev_meds, medians, medians2 );
        
        %3. Plot original image and cluster centers
        plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
        
    end
    
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
end


%% 5. SET DISTANCE CONSTRAINTS OF HAND MARKERS TO ARM MARKERS

%% Set distance limits of hand markers to red elbow marker (after plotting the distances first)

%PLOT
%Plots the distances of every hand marker to to red elbow marker in order
%to determine what distances are allowed (for rerunning the hand marker
%tracking)
if first_time
    
    %Calculate distances
    for i=1:n_times_prelim
        dists1(i)=pdist2(all_medians(10,:,i),all_medians(1,:,i)); %Distance from point 10 (red elbow to point 1)
        dists2(i)=pdist2(all_medians(10,:,i),all_medians(2,:,i)); %etc...
        dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i));
        dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i));
        dists5(i)=pdist2(all_medians(10,:,i),all_medians(5,:,i));
    end
    
    %Plot
    figure; plot(dists1,'g'); hold on;
    plot(dists2,'b');
    plot(dists3,'r');
    plot(dists4,'r');
    plot(dists5,'g');
    title('Distance between red elbow marker (10) and hand markers');
    legend('10-1','10-2','10-3','10-4','10-5');
    
    
    if use_defaults %If you use the defaults
        
        str1='Enter vector of times to calculate default distances (e.g. 1:20000) \n';
        str2='Or just press enter to use all times \n';
        times_inc=input([str1 str2]);
        if isempty(times_inc)
            times_inc=1:n_times_prelim;
        end 
        
        %Get the distances at the times you include
        dists1=dists1(times_inc);
        dists2=dists2(times_inc);
        dists3=dists2(times_inc);
        dists4=dists2(times_inc);
        dists5=dists2(times_inc);
        
        green_hand_dists_elbow=[nanmean(dists5)-5*nanstd(dists5) nanmean(dists1)+5*nanstd(dists1)];
        red_hand_dists_elbow=[nanmean(dists3)-5*nanstd(dists3) nanmean(dists3)+5*nanstd(dists3)];
        blue_hand_dists_elbow=[nanmean(dists2)-5*nanstd(dists2) nanmean(dists2)+5*nanstd(dists2)];
        green_separator=(nanmean(dists5)+nanmean(dists1))/2; 
        
    else %If you don't use the defaults
        
  
        %Visualize Frames
        user_input=1; %A value so that it enters the while loop below
        while ~isempty(user_input)
            str1='Below you will enter limits for the distances between the red elbow marker and hand markers \n';
            str2='First, to help set these limits, you can enter a time point you want to visualize (or just press enter to continue) \n';
            user_input=input([str1 str2]);
            %Make sure the input was valid (an integer between start and finish)
            if ~isempty(user_input)
                while ~(isnumeric(user_input) && mod(user_input,1)==0 && user_input>=start && user_input<=finish)
                    user_input=input('Re-enter valid time point \n');
                end
            end
            if ~isempty(user_input)
                plot_together_3colors_func(user_input, [1:5], [1:10], all_medians, color1, color2, color3, start, finish, 1)
            end
        end
        
        
        
        % Enter hand distance limits from red elbow
        
        str1='Input green_hand_dists_elbow \n';
        str2='Lower and upper limits of distances of the green hand markers to the red elbow marker\n';
        str3='Value is generally ~ [.15,.26] (around 2 or 3 cm from most points) \n';
        green_hand_dists_elbow=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(green_hand_dists_elbow)==2)
            green_hand_dists_elbow=input('Re-enter valid values');
        end
        
        str1='Input red_hand_dists_elbow \n';
        str2='Lower and upper limits of distances of the red hand markers to the red elbow marker\n';
        str3='Value is generally ~ [.17,.23] (around 2 or 3 cm from most points) \n';
        red_hand_dists_elbow=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(red_hand_dists_elbow)==2)
            red_hand_dists_elbow=input('Re-enter valid values');
        end
        
        str1='Input blue_hand_dists_elbow \n';
        str2='Lower and upper limits of distances of the blue hand markers to the red elbow marker\n';
        str3='Value is generally ~ [.17,.23] (around 2 or 3 cm from most points) \n';
        blue_hand_dists_elbow=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(blue_hand_dists_elbow)==2)
            blue_hand_dists_elbow=input('Re-enter valid values');
        end
        
        str1='Input green_separator \n';
        str2='Distance that separates the green hand points (marker 1 and 5) \n';
        str3='Value is generally ~ .2\n';
        green_separator=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(isnumeric(green_separator))
            green_separator=input('Re-enter valid values');
        end
        
    end
    
end

%% Set distance limits of hand markers to blue arm marker (after plotting the distances first)

%PLOT
%Plots the distances of every hand marker to to blue arm marker in order
%to determine what distances are allowed (for rerunning the hand marker
%tracking)

%Note that using the distance from the hand to the blue arm marker only is helpful for a task with holding the handle

if first_time
    
    %Calculate distances
    for i=1:n_times_prelim
        dists1(i)=pdist2(all_medians(7,:,i),all_medians(1,:,i)); %Distance from point 7 (blue arm) to point 1
        dists2(i)=pdist2(all_medians(7,:,i),all_medians(2,:,i));
        dists3(i)=pdist2(all_medians(7,:,i),all_medians(3,:,i));
        dists4(i)=pdist2(all_medians(7,:,i),all_medians(4,:,i));
        dists5(i)=pdist2(all_medians(7,:,i),all_medians(5,:,i));
    end
    
    
    if use_defaults %If you use the defaults
        
        %Get the distances at the times you include (specified when setting
        %distances to red elbow marker)
        dists1=dists1(times_inc);
        dists2=dists2(times_inc);
        dists3=dists2(times_inc);
        dists4=dists2(times_inc);
        dists5=dists2(times_inc);
        
        green_hand_dists_bluearm=[nanmean(dists5)-6*nanstd(dists5) nanmean(dists1)+6*nanstd(dists1)];
        red_hand_dists_bluearm=[nanmean(dists3)-6*nanstd(dists3) nanmean(dists3)+6*nanstd(dists3)];
        blue_hand_dists_bluearm=[nanmean(dists2)-6*nanstd(dists2) nanmean(dists2)+6*nanstd(dists2)];
        
    else %If you don't use the defaults
        
        
        %Plot
        figure; %Blue
        plot(dists2,'b-x');
        title('Distance from blue arm marker to blue hand marker');
        figure; %Red
        plot(dists3,'r-x');
        hold on;
        plot(dists4,'m-x');
        title('Distance from blue arm marker to red hand marker');
        legend('Dist to pt3', 'Dist to pt4')
        figure; %Green
        plot(dists1,'g-x');
        hold on;
        plot(dists5,'c-x');
        title('Distance from blue arm marker to green hand markers');
        legend('Dist to pt1','Dist to pt5');
        
        
        
        % Enter hand distance limits from blue arm
        
        str1='Input green_hand_dists_bluearm \n';
        str2='Lower and upper limits of distances of the green hand markers (green and cyan above) to the blue arm marker\n';
        str3='Value is generally ~ [.15,.30] (around 4 cm from most points) \n';
        green_hand_dists_bluearm=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(green_hand_dists_bluearm)==2)
            green_hand_dists_bluearm=input('Re-enter valid values');
        end
        
        str1='Input red_hand_dists_bluearm \n';
        str2='Lower and upper limits of distances of the red hand markers (red above) to the blue arm marker\n';
        str3='Value is generally ~ [.16,.28] (around 4 cm from most points) \n';
        red_hand_dists_bluearm=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(red_hand_dists_bluearm)==2)
            red_hand_dists_bluearm=input('Re-enter valid values');
        end
        
        str1='Input blue_hand_dists_bluearm \n';
        str2='Lower and upper limits of distances of the blue hand markers (blue above) to the blue arm marker\n';
        str3='Value is generally ~ [.16,.28] (around 4 cm from most points) \n';
        blue_hand_dists_bluearm=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(blue_hand_dists_bluearm)==2)
            blue_hand_dists_bluearm=input('Re-enter valid values');
        end
        
    end
end

%% Set distance limits of hand markers to red arm marker (after plotting the distances first)

%Plots the distances of every hand marker to to red arm marker in order
%to determine what distances are allowed (for rerunning the hand marker
%tracking)

%Note that using the distance from the hand to the red arm marker only is helpful for a task with holding the handle


if first_time
    %Calculate distances
    for i=1:n_times_prelim
        dists1(i)=pdist2(all_medians(8,:,i),all_medians(1,:,i)); %Distance from point 8 (red arm to point 1)
        dists2(i)=pdist2(all_medians(8,:,i),all_medians(2,:,i));
        dists3(i)=pdist2(all_medians(8,:,i),all_medians(3,:,i));
        dists4(i)=pdist2(all_medians(8,:,i),all_medians(4,:,i));
        dists5(i)=pdist2(all_medians(8,:,i),all_medians(5,:,i));
    end
    
    
    if use_defaults %If you use the defaults
        
        %Get the distances at the times you include (specified when setting
        %distances to red elbow marker)
        dists1=dists1(times_inc);
        dists2=dists2(times_inc);
        dists3=dists2(times_inc);
        dists4=dists2(times_inc);
        dists5=dists2(times_inc);
        
        green_hand_dists_redarm=[nanmean(dists5)-6*nanstd(dists5) nanmean(dists1)+6*nanstd(dists1)];
        red_hand_dists_redarm=[nanmean(dists3)-6*nanstd(dists3) nanmean(dists3)+6*nanstd(dists3)];
        blue_hand_dists_redarm=[nanmean(dists2)-6*nanstd(dists2) nanmean(dists2)+6*nanstd(dists2)];
        
    else %If you don't use the defaults
        
        
        %Plot
        figure; %Blue
        plot(dists2,'b-x');
        title('Distance from red arm marker to blue hand markers');
        figure; %Red
        plot(dists3,'r-x');
        hold on;
        plot(dists4,'m-x');
        title('Distance from blue arm marker to red hand marker');
        legend('Dist to pt3', 'Dist to pt4')
        figure; %Green
        plot(dists1,'g-x');
        hold on;
        plot(dists5,'c-x');
        title('Distance from red arm marker to green hand markers');
        legend('Dist to pt1','Dist to pt5');
        
        
        
        % Enter hand distance limits from red arm
        
        str1='Input green_hand_dists_redarm \n';
        str2='Lower and upper limits of distances of the green hand markers (green and cyan above) to the red arm marker\n';
        str3='Value is generally ~ [.15,.35] (around 5 cm from most points) \n';
        green_hand_dists_redarm=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(green_hand_dists_redarm)==2)
            green_hand_dists_redarm=input('Re-enter valid values');
        end
        
        str1='Input red_hand_dists_redarm \n';
        str2='Lower and upper limits of distances of the red hand markers (red above) to the red arm marker\n';
        str3='Value is generally ~ [.15,.34] (around 5 cm from most points) \n';
        red_hand_dists_redarm=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(red_hand_dists_redarm)==2)
            red_hand_dists_redarm=input('Re-enter valid values');
        end
        
        str1='Input blue_hand_dists_redarm \n';
        str2='Lower and upper limits of distances of the blue hand markers (blue above) to the red arm marker\n';
        str3='Value is generally ~ [.15,.34] (around 5 cm from most points) \n';
        blue_hand_dists_redarm=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(length(blue_hand_dists_redarm)==2)
            blue_hand_dists_redarm=input('Re-enter valid values');
        end
        
    end
end

%% Set distance limits between green hand markers

%Plots the distances of the green hand markers to each other, in order to
%determine what distances are allowed (for rerunning the hand marker tracking)

if first_time
    %Calculate distances
    for i=1:n_times_prelim
        dists1(i)=pdist2(all_medians(5,:,i),all_medians(1,:,i)); %Distances between green hand markers
    end
    
    
    if use_defaults %If you use the defaults
        
        %Get the distances at the times you include (specified when setting
        %distances to red elbow marker)
        dists1=dists1(times_inc);
        
        green_dist_min=nanmean(dists1)-5*nanstd(dists1);
        
        if green_dist_min<.02 %Should never be less than this value
            green_dist_min=.02;
        end
        
    else %If you don't use the defaults
        
     
        %Plot
        figure; plot(dists1,'g');
        title('Distance between green hand markers');
        
        
        % Enter minimum hand distances from each other
        
        str1='Input green_dist_min \n';
        str2='Minimum distance allowed between green hand markers (green above)\n';
        str3='Value is generally ~ .03 \n';
        green_dist_min=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(isnumeric(green_dist_min))
            green_dist_min=input('Re-enter valid values');
        end
        
    end
    
end

%% Set distance limits between red hand markers

%Plots the distances of the red hand markers to each other, in order to
%determine what distances are allowed (for rerunning the hand marker tracking)

if first_time
    %Calculate distances
    for i=1:n_times_prelim
        dists2(i)=pdist2(all_medians(3,:,i),all_medians(4,:,i)); %Distances between green hand markers
    end
    
    
    if use_defaults %If you use the defaults
        
        %Get the distances at the times you include (specified when setting
        %distances to red elbow marker)
        dists2=dists2(times_inc);
        
        red_dist_min=nanmean(dists2)-5*nanstd(dists2);
        
%         if red_dist_min<.01 %Should never be less than this value
%             red_dist_min=.01;
%         end
        
    else %If you don't use the defaults
        
     
        %Plot
        figure; plot(dists2,'r');
        title('Distance between red hand markers');
        
        
        % Enter minimum hand distances from each other
        
        str1='Input red_dist_min \n';
        str2='Minimum distance allowed between red hand markers (red above)\n';
        str3='Value is generally ~ .02 \n';
        red_dist_min=input([str1 str2 str3]);
        %Make sure entry was valid
        while ~(isnumeric(red_dist_min))
            red_dist_min=input('Re-enter valid values');
        end
        
    end
    
end




%% Red Hand (Redo)

%Initializations
plot_on=0;
marker_ids=red_hand_marker_ids;
color=color3;
prev_meds=marker_inits(marker_ids,:);
num_clust=length(marker_ids); 
within_clust_dist1=.05; %How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
within_clust_dist2=.05; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered
dist_min=red_dist_min; 

medians=NaN(num_clust,3,n_times);
medians2=NaN(num_clust,3,n_times); 

num_gone1=0; %The number of frames the first marker has been gone
num_gone2=0; %The number of frames the second marker has been gone

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
    
    % Keep all the points close enough to either of the previous markers
    keep1=D(:,1)<within_clust_dist1;
    keep2=D(:,2)<within_clust_dist2;
    
    %Remove points that are too close or far from the red elbow marker
    D2=pdist2(loc,all_medians(10,:,t));
    rmv0=D2<red_hand_dists_elbow(1) | D2>red_hand_dists_elbow(2);
    
    %Remove points that are too close or far from the blue arm marker
    D3=pdist2(loc,all_medians(7,:,t));
    rmv1=D3<red_hand_dists_bluearm(1) | D3>red_hand_dists_bluearm(2);
    
    %Remove points that are too close or far from the red arm marker
    D4=pdist2(loc,all_medians(8,:,t));
    rmv2=D4<red_hand_dists_redarm(1) | D4>red_hand_dists_redarm(2);
    
    %Use above criteria to set points for removal
    %We will always remove points that are too close or far from the arm markers (rmv0, rmv1, rmv2).
    %Depending on how many frames the markers have been missing, we additionally use different criteria for removing.
    %If both markers have been missing for <=4 frames, keep points close
    %enough to the marker's locations in the previous frame
    if num_gone1<=4 & num_gone2<=4
        rmv=~(keep1 | keep2) | rmv0 | rmv1 | rmv2; 
    %If the second marker has been missing for >4 frames, only keep points close
    %enough to the first marker's location in the previous frame
    else if num_gone1<=4 & num_gone2>4
            rmv=~(keep1) | rmv0 | rmv1 | rmv2;
    %If the first marker has been missing for >4 frames, only keep points close
    %enough to the second marker's location in the previous frame
        else if num_gone1>4 & num_gone2<=4
                rmv=~(keep2) | rmv0 | rmv1 | rmv2;
   %If both markers have been missing for >4 frames, don't keep any points
   %based on distance to the markers' locations in the previous frame
            else
                rmv=rmv0 | rmv1 | rmv2;
            end
        end
    end
    
    
    %Actually remove
    loc(rmv,:)=[];
    
    
    %2. Cluster and assign
    [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
    
    %Update how many frames markers have been missing
    %If marker1 is missing, add 1 to num_gone1. Otherwise set num_gone1=0
    %(since it's been missing 0 frames)
    if isnan(medians(1,1,t))
        num_gone1=num_gone1+1;
    else
        num_gone1=0;
    end
    %If marker2 is missing, add 1 to num_gone2. Otherwise set num_gone2=0
    %(since it's been missing 0 frames)
    if isnan(medians(2,1,t))
        num_gone2=num_gone2+1;
    else
        num_gone2=0;
    end
    
    %If the red elbow marker is not missing, make sure the the first marker
    %(marker # 3) is farther from the red elbow marker than the second
    %marker (marker #4). If not, switch their assignment.
    %If the red elbow marker is missing, but the blue arm marker is not
    %missing, then do the same as above w/ the blue arm marker.
    if t>1
        %     if isnan(medians(1,1,t-1)) || isnan(medians(2,1,t-1))
        if ~isnan(all_medians(10,1,t))
            if pdist2(medians(1,:,t),all_medians(10,:,t))<pdist2(medians(2,:,t),all_medians(10,:,t))
                temp=medians(1,:,t);
                temp2=medians2(1,:,t);
                medians(1,:,t)=medians(2,:,t);
                medians(2,:,t)=temp;
                medians2(1,:,t)=medians2(2,:,t);
                medians2(2,:,t)=temp2;
            end
        else if ~isnan(all_medians(7,1,t))
                if pdist2(medians(1,:,t),all_medians(7,:,t))<pdist2(medians(2,:,t),all_medians(7,:,t))
                    temp=medians(1,:,t);
                    temp2=medians2(1,:,t);
                    medians(1,:,t)=medians(2,:,t);
                    medians(2,:,t)=temp;
                    medians2(1,:,t)=medians2(2,:,t);
                    medians2(2,:,t)=temp2;
                end
            end
        end
        %     end
        
        
        %If the markers are too far away from the markers in the prevoius
        %frame, remove them.
        %For the first marker
        if pdist2(medians(1,:,t),medians(1,:,t-1))>within_clust_dist1
            medians(1,:,t)=NaN;
            medians2(1,:,t)=medians2(1,:,t-1);
        end
        %For the second marker
        if pdist2(medians(2,:,t),medians(2,:,t-1))>within_clust_dist2
            medians(2,:,t)=NaN;
            medians2(2,:,t)=medians2(2,:,t-1);
        end
        
        
    end
       
    %3. Plot original image and cluster centers
    plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
    
end

all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;

%% Green Hand (Redo)

%Initializations
plot_on=0;
marker_ids=green_hand_marker_ids;
color=color2;
prev_meds=marker_inits(marker_ids,:);
num_clust=length(marker_ids); %Number of clusters
within_clust_dist1=.07;%How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
within_clust_dist2=.07; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered
dist_min=green_dist_min; %Minimum distance between markers (cluster medians aren't allowed w/ distance < min_dist)

medians=NaN(num_clust,3,n_times); %Has NaNs when a marker is missing
medians2=NaN(num_clust,3,n_times); %Has previous known positions when a marker is missing

num_gone1=0;%The number of frames the first marker has been gone
num_gone2=0; %The number of frames the second marker has been gone

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
    
    % Keep all the points close enough to either of the previous markers
    keep1=D(:,1)<within_clust_dist1;
    keep2=D(:,2)<within_clust_dist2;
    
    %Remove points that are too close or far from the red elbow marker
    D2=pdist2(loc,all_medians(10,:,t));
    rmv0=D2<green_hand_dists_elbow(1) | D2>green_hand_dists_elbow(2);
    
    %Remove points that are too close or far from the blue arm marker
    D3=pdist2(loc,all_medians(7,:,t));
    rmv1=D3<green_hand_dists_bluearm(1) | D3>green_hand_dists_bluearm(2);

    %Remove points that are too close or far from the red arm marker    
    D4=pdist2(loc,all_medians(8,:,t));
    rmv2=D4<green_hand_dists_redarm(1) | D4>green_hand_dists_redarm(2);
    
    %Use above criteria to set points for removal
    %We will always remove points that are too close or far from the arm markers (rmv0, rmv1, rmv2).
    %Depending on how many frames the markers have been missing, we additionally use different criteria for removing.
    
    %If both markers have been missing for <=4 frames, keep points close
    %enough to the marker's locations in the previous frame
    if num_gone1<=4 & num_gone2<=4
        rmv=~(keep1 | keep2) | rmv0 | rmv1 | rmv2;
    %If the second marker has been missing for >4 frames, only keep points close
    %enough to the first marker's location in the previous frame
    else if num_gone1<=4 & num_gone2>4
            rmv=~(keep1) | rmv0 | rmv1 | rmv2;
    %If the first marker has been missing for >4 frames, only keep points close
    %enough to the second marker's location in the previous frame
        else if num_gone1>4 & num_gone2<=4
                rmv=~(keep2) | rmv0 | rmv1 | rmv2;
    %If both markers have been missing for >4 frames, don't keep any points
    %based on distance to the markers' locations in the previous frame
            else
                rmv=rmv0 | rmv1 | rmv2;
            end
        end
    end
    
    
    %Actually remove
    loc(rmv,:)=[];
        
    %2. Cluster and assign
    [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
    
    %Update how many frames markers have been missing
    %If marker1 is missing, add 1 to num_gone1. Otherwise set num_gone1=0
    %(since it's been missing 0 frames)
    if isnan(medians(1,1,t))
        num_gone1=num_gone1+1;
    else
        num_gone1=0;
    end
    %If marker2 is missing, add 1 to num_gone2. Otherwise set num_gone2=0
    %(since it's been missing 0 frames)
    if isnan(medians(2,1,t))
        num_gone2=num_gone2+1;
    else
        num_gone2=0;
    end
    
    %If the red elbow marker is not missing, make sure the the first marker
    %(marker # 1) is farther from the red elbow marker than the second
    %marker (marker #5). If not, switch their assignment.
    %If the red elbow marker is missing, but the blue arm marker is not
    %missing, then do the same as above w/ the blue arm marker.
    if t>1
        %     if isnan(medians(1,1,t-1)) || isnan(medians(2,1,t-1))
        if ~isnan(all_medians(10,1,t))
            if pdist2(medians(1,:,t),all_medians(10,:,t))<pdist2(medians(2,:,t),all_medians(10,:,t))
                temp=medians(1,:,t);
                temp2=medians2(1,:,t);
                medians(1,:,t)=medians(2,:,t);
                medians(2,:,t)=temp;
                medians2(1,:,t)=medians2(2,:,t);
                medians2(2,:,t)=temp2;
            end
        else if ~isnan(all_medians(7,1,t))
                if pdist2(medians(1,:,t),all_medians(7,:,t))<pdist2(medians(2,:,t),all_medians(7,:,t))
                    temp=medians(1,:,t);
                    temp2=medians2(1,:,t);
                    medians(1,:,t)=medians(2,:,t);
                    medians(2,:,t)=temp;
                    medians2(1,:,t)=medians2(2,:,t);
                    medians2(2,:,t)=temp2;
                end
            end
        end
        %     end
        
        %If both markers were gone previous frame, and now there's one, assume
        %it's not the first marker (marker #1  by the fingers)
        if isnan(medians(1,1,t-1)) && isnan(medians(2,1,t-1))
            if ~isnan(medians(1,1,t)) && isnan(medians(2,1,t)) %If only the first marker shows up, then flip the assignment
                temp=medians(1,:,t);
                medians(1,:,t)=medians(2,:,t);
                medians(2,:,t)=temp;
            end
        end
        
        %If the markers are too far away from the markers in the prevoius
        %frame, remove them.
        %For the first marker
        if abs(pdist2(medians(1,:,t),medians(1,:,t-1)))>within_clust_dist1
            medians(1,:,t)=NaN;
            medians2(1,:,t)=medians2(1,:,t-1);
        end
        %For the second marker
        if abs(pdist2(medians(2,:,t),medians(2,:,t-1)))>within_clust_dist2
            medians(2,:,t)=NaN;
            medians2(2,:,t)=medians2(2,:,t-1);
        end
        
        
        %If there's only a single marker (one is missing this frame), 
        %determine its label based on the distance from the red elbow 
        
        %The distance from the second marker (marker #5) to the red elbow
        %should be less than "green_separator." If it's more, then
        %change the label of this marker to the first marker (marker #1).
        if isnan(medians(1,1,t)) && ~isnan(medians(2,1,t))
            if pdist2(medians(2,:,t),all_medians(10,:,t))>green_separator
                medians(1,:,t)=medians(2,:,t);
                medians(2,:,t)=NaN;
            end
        end
        %The distance from the first marker (marker #1) to the red elbow
        %should be greater than "green_separator." If it's less, then
        %change the label of this marker to the second marker (marker #5).
        if ~isnan(medians(1,1,t)) && isnan(medians(2,1,t))
            if pdist2(medians(1,:,t),all_medians(10,:,t))<green_separator
                medians(2,:,t)=medians(1,:,t);
                medians(1,:,t)=NaN;
            end
        end
        
        
    end
    %3. Plot original image and cluster centers
    plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
    
end

all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;

%% Blue Hand (Redo)

%Initializations
plot_on=0;
marker_ids=blue_hand_marker_ids;
color=color1;
prev_meds=marker_inits(marker_ids,:);
num_clust=length(marker_ids); 
within_clust_dist1=.07;
dist_min=0.07; 

medians=NaN(num_clust,3,n_times); 
medians2=NaN(num_clust,3,n_times); 

num_gone=0; %The number of frames the marker has been gone

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
    
    %Remove points that are too close or far from the red elbow marker
    D2=pdist2(loc,all_medians(10,:,t));
    rmv0=D2<blue_hand_dists_elbow(1) | D2>blue_hand_dists_elbow(2);
    
    %Remove points that are too close or far from the blue arm marker
    D3=pdist2(loc,all_medians(7,:,t));
    rmv1=D3<blue_hand_dists_bluearm(1) | D3>blue_hand_dists_bluearm(2);
    
    %Remove points that are too close or far from the red arm marker
    D4=pdist2(loc,all_medians(8,:,t));
    rmv2=D4<blue_hand_dists_redarm(1) | D4>blue_hand_dists_redarm(2);
    
    %Use above criteria to set points for removal
    %We will always remove points that are too close or far from the arm markers (rmv0, rmv1, rmv2).    
    %If the marker has been missing for <=4 frames, keep points close
    %enough to the marker's location in the previous frame
    if num_gone<=4
        rmv=~(keep1)| rmv0 | rmv1 | rmv2;
    else
    %If the marker has been missing for >4 frames, don't keep any points
    %based on distance to the marker's location in the previous frame
        rmv=rmv0 | rmv1 | rmv2;
    end
    
    %Actually remove
    loc(rmv,:)=[];
    
    
    %2. Cluster and assign
    [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
    
    
    %Update how many frames markers have been missing
    %If the marker is missing, add 1 to num_gone. Otherwise set num_gone=0
    %(since it's been missing 0 frames)
    if isnan(medians(1,1,t))
        num_gone=num_gone+1;
    else
        num_gone=0;
    end
    
    %If marker 2 is too far from the wrist
    
    %     if pdist2(medians(1,:,t),all_medians(4,:,t))>dist_from_wrist
    %         if pdist2(medians(1,:,t),all_medians(5,:,t))>dist_from_wrist
    %             medians(1,:,t)=NaN;
    %             medians2(1,:,t)=medians2(1,:,t-1);
    %         end
    %     end
    
    
    %3. Plot original image and cluster centers
    plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
    
end

all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;

%% REMOVE/SWITCH HAND MARKERS BELOW

%% Remove marker 3 that is too far away from other hand markers - Plot

%Calculate distance from marker 3 to other hand markers
for i=1:n_times
    dists4(i)=pdist2(all_medians(3,:,i),all_medians(4,:,i)); %Distance between marker 3 and marker 4...
    dists5(i)=pdist2(all_medians(3,:,i),all_medians(5,:,i));
    dists2(i)=pdist2(all_medians(3,:,i),all_medians(2,:,i));
end

%Plot
if first_time   
    figure;
    hold on;
    plot(dists4,'-x')
    plot(dists5,'-x')
    plot(dists2,'-x')
end
%% Remove marker 3 that is too far away from other hand markers

%Determine frames when each of the distances is greater than expected
d2=dists2>nanmean(dists2)+.02;
d4=dists4>nanmean(dists4)+.02;
d5=dists5>nanmean(dists5)+.02;

%Remove when distance from all points is too large (including times when
%marker 2 is missing)
rmv=(d4 & d5 & isnan(dists2)) | (d4 & d5 & d2);
all_medians(3,:,rmv)=NaN;

%% Remove marker 2 that is too far away from other hand markers - Plot

%Calculate distance from marker 2 to other hand markers
for i=1:n_times
    dists3(i)=pdist2(all_medians(2,:,i),all_medians(3,:,i)); %Distances from marker 2 to 3...
    dists4(i)=pdist2(all_medians(2,:,i),all_medians(4,:,i));
    dists5(i)=pdist2(all_medians(2,:,i),all_medians(5,:,i));
end

%Plot
if first_time    
    figure;
    hold on;
    plot(dists3,'-x')
    plot(dists4,'-x')
    plot(dists5,'-x')    
end
%% Remove marker 2 that is too far away from other hand markers

%Determine frames when each of the distances is greater than expected
d3=dists3>nanmean(dists3)+.03;
d4=dists4>nanmean(dists4)+.03;
d5=dists5>nanmean(dists5)+.03;

%Remove when distance from all points is too large (including times when
%marker 3 is missing)
rmv=(d4 & d5 & isnan(dists3)) | (d4 & d5 & d3);
all_medians(2,:,rmv)=NaN;

%% Calculate/Plot hand distances from elbow

%Calculate distances from red elbow to points on the hand
for i=1:n_times
    dists1(i)=pdist2(all_medians(10,:,i),all_medians(1,:,i));
    dists2(i)=pdist2(all_medians(10,:,i),all_medians(2,:,i));
    dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i));
    dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i));
    dists5(i)=pdist2(all_medians(10,:,i),all_medians(5,:,i));
    
    % dists1(i)=pdist2(all_medians(6,:,i),all_medians(1,:,i));
    % dists2(i)=pdist2(all_medians(6,:,i),all_medians(2,:,i));
    % dists3(i)=pdist2(all_medians(6,:,i),all_medians(3,:,i));
    % dists4(i)=pdist2(all_medians(6,:,i),all_medians(4,:,i));
    % dists5(i)=pdist2(all_medians(6,:,i),all_medians(5,:,i));
    
    % dists1(i)=pdist2(all_medians(7,:,i),all_medians(1,:,i));
    % dists2(i)=pdist2(all_medians(7,:,i),all_medians(2,:,i));
    % dists3(i)=pdist2(all_medians(7,:,i),all_medians(3,:,i));
    % dists4(i)=pdist2(all_medians(7,:,i),all_medians(4,:,i));
    % dists5(i)=pdist2(all_medians(7,:,i),all_medians(5,:,i));
end

%Plot
figure;
plot(dists1,'g-x');
hold on;
plot(dists2,'b-x');
plot(dists3,'r-x');
plot(dists4,'m-x');
plot(dists5,'c-x');

%% SET 9. Determine red hand point 3 to remove (based on having a larger distance to the elbow than point 1)

wdw=20; %Time window for plots below

%Find times (idxs) where the distance from the elbow to point 3 are
%greater than the distance from the elbow to point 1 (which shouldn't
%happen)
idxs=find(dists3>dists1); 

%Initialize points to be removed
rmv3_cell=cell(1,length(idxs));
rmv1_cell=cell(1,length(idxs));
switch1_cell=cell(1,length(idxs));

%Look at a time window around each of these points to determine what point
%to remove (which is specified in the next section).
for i=1:length(idxs)
    if idxs(i)+wdw<=n_times & idxs(i)-wdw>0
        figure;
        plot(idxs(i)-wdw:idxs(i)+wdw,dists1(idxs(i)-wdw:idxs(i)+wdw),'g-x');
        hold on;
        plot(idxs(i)-wdw:idxs(i)+wdw,dists2(idxs(i)-wdw:idxs(i)+wdw),'b-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists3(idxs(i)-wdw:idxs(i)+wdw),'r-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists4(idxs(i)-wdw:idxs(i)+wdw),'m-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists5(idxs(i)-wdw:idxs(i)+wdw),'c-x');        
        title(num2str(idxs(i)));
        
        rmv3_cell{i}=input('9. Times to remove hand point 3 (red) \n');
        rmv1_cell{i}=input('9. Times to remove hand point 1 (green) \n');
        switch1_cell{i}=input('9. Times to switch markers 1 and 5 (green and cyan) \n');
    end
end

%% Set red hand point 3 and green hand point 1 to remove (and green switch), based on above section

rmv3=[rmv3_cell{:}];  %Times to remove of red hand point 3
rmv1=[rmv1_cell{:}];  %Times to remove of green hand point 1
switch1=[switch1_cell{:}]; %Times to switch markers 1 and 5 (green and cyan)

%% Remove and switch (based on above section)

%Remove (based on above section)
all_medians(3,:,rmv3)=NaN;
all_medians(1,:,rmv1)=NaN;

%Switch Markers 1 and 5 (based on above section)
temp=all_medians(1,:,switch1);
temp2=all_medians2(1,:,switch1);
all_medians(1,:,switch1)=all_medians(5,:,switch1);
all_medians2(1,:,switch1)=all_medians2(5,:,switch1);
all_medians(5,:,switch1)=temp;
all_medians2(5,:,switch1)=temp2;

%% SET 10. Determine blue hand point 2 to remove (based on having a larger distance to the elbow than point 1)

wdw=20; %Time window for plots below

%Find times (idxs) where the distance from the elbow to point 2 are
%greater than the distance from the elbow to point 1 (which shouldn't
%happen)
idxs=find(dists2>dists1);

%Initialize points to be removed
rmv2_cell=cell(1,length(idxs));
rmv1_cell=cell(1,length(idxs));

%Look at a time window around each of these points to determine what point
%to remove (which is specified in the next section).
for i=1:length(idxs)
    if idxs(i)+wdw<=n_times & idxs(i)-wdw>0
        figure;
        plot(idxs(i)-wdw:idxs(i)+wdw,dists1(idxs(i)-wdw:idxs(i)+wdw),'g-x');
        hold on;
        plot(idxs(i)-wdw:idxs(i)+wdw,dists2(idxs(i)-wdw:idxs(i)+wdw),'b-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists3(idxs(i)-wdw:idxs(i)+wdw),'r-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists4(idxs(i)-wdw:idxs(i)+wdw),'m-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists5(idxs(i)-wdw:idxs(i)+wdw),'c-x');
        
        title(num2str(idxs(i)));
        
        rmv2_cell{i}=input('10. Times to remove hand point 2 (blue) \n');
        rmv1_cell{i}=input('10. Times to remove hand point 1 (green) \n');
    end
end

%% Set blue hand point 2 and green hand point 1 to remove, based on above section

%When blue is greater than green in above plot

rmv2=[rmv2_cell{:}]; %Times to remove of blue hand point 2
rmv1=[rmv1_cell{:}]; %Times to remove of green hand point 1

all_medians(2,:,rmv2)=NaN;
all_medians(1,:,rmv1)=NaN;

%% SET 11. Determine red hand points to remove (based on having similar distance from elbow)

%Calculate distances of red hand markers 3 and 4, and the red elbow
for i=1:n_times    
    dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i));
    dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i));   
end


wdw=20; %Time window for plots below

%Find times when marker 3 and marker 4 are a similar distance to the elbow
%marker (which is a problem)
idxs=find(abs(dists3-dists4)<.01);

%Initialize points to be removed
red_test_cell=cell(1,length(idxs));
pink_test_cell=cell(1,length(idxs));


%Plot those times
for i=1:length(idxs)
    figure;
    plot(idxs(i)-wdw:idxs(i)+wdw,dists1(idxs(i)-wdw:idxs(i)+wdw),'g-x');
    hold on;
    plot(idxs(i)-wdw:idxs(i)+wdw,dists2(idxs(i)-wdw:idxs(i)+wdw),'b-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists3(idxs(i)-wdw:idxs(i)+wdw),'r-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists4(idxs(i)-wdw:idxs(i)+wdw),'m-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists5(idxs(i)-wdw:idxs(i)+wdw),'c-x');
    
    title(num2str(idxs(i)));
    
    red_test_cell{i}=input('11. Times pink is in red area \n');
    pink_test_cell{i}=input('11. Times red is in pink area \n');
    
end

%% Set times when red hand points should be removed (based on having similar distance from elbow)

red_test=[red_test_cell{:}]; %Pink is in red area
pink_test=[pink_test_cell{:}]; %Red is in pink area

%% Determine which red hand point (3 or 4) to remove, at the above times

%Above we specified times when the two red markers had the same distance from
%the elbow. 

%Sometimes both were had the expected distance for marker 3 (red
%points). In this case, we want to determine which of these markers really
%should be marker 3, and which should be removed. We do this by determining
%which is closer to the previous frame's marker 3.
red_points=sort(red_test);
for i=1:length(red_points)
    t=red_points(i);
    %Compare the distances from the current markers 3 and 4 to the
    %previous marker 3. Make the one that's closer the new marker 3.
    if pdist2(all_medians2(3,:,t-1),all_medians(3,:,t))>pdist2(all_medians2(3,:,t-1),all_medians(4,:,t))
        all_medians(3,:,t)=all_medians(4,:,t);
    end
    all_medians(4,:,t)=NaN;
end

%Sometimes both were had the expected distance for marker 4 (pink
%points). In this case, we want to determine which of these markers really
%should be marker 4, and which should be removed. We do this by determining
%which is closer to the previous frame's marker 4.
pink_points=sort(pink_test);
for i=1:length(pink_points)
    t=pink_points(i);
    %Compare the distances from the current markers 3 and 4 to the
    %previous marker 4. Make the one that's closer the new marker 4.
    if pdist2(all_medians2(4,:,t-1),all_medians(4,:,t))>pdist2(all_medians2(4,:,t-1),all_medians(3,:,t))
        all_medians(4,:,t)=all_medians(3,:,t);
    end
    all_medians(3,:,t)=NaN;
end

%% Calculate hand distances from red elbow marker

%Redo, since some points have been updated in the previous sections
for i=1:n_times
    dists1(i)=pdist2(all_medians(10,:,i),all_medians(1,:,i));
    dists2(i)=pdist2(all_medians(10,:,i),all_medians(2,:,i));
    dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i));
    dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i));
    dists5(i)=pdist2(all_medians(10,:,i),all_medians(5,:,i));
end

%Plot

% figure;
% plot(dists1,'g-x');
% hold on;
% plot(dists2,'b-x');
% plot(dists3,'r-x');
% plot(dists4,'m-x');
% plot(dists5,'c-x');

%% Switch red hand points when one is missing (the assignment was wrong)

%Here, we'll determine what time points to switch the red hand markers (3
%and 4) when one of those markers was missing. That is, there was a single
%marker that was incorrectly assigned.
%We will determine this based on having a expected/unexpected distance from
%the red elbow marker. We will go both forward and backward in time, to see
%whether switching the marker assignment would lead to distances from the
%elbow that are more consistent with the previous time.

%Initializations for forward run
dists3_forward=NaN(1,n_times);
dists4_forward=NaN(1,n_times);
dists33_forward=NaN(1,n_times);
dists44_forward=NaN(1,n_times);
dists34_forward=NaN(1,n_times);
dists43_forward=NaN(1,n_times);
change_forward=zeros(1,n_times);

%Calculate the distances from the red hand markers (3 and 4) to the red
%elbow marker (10)
for i=1:n_times    
    dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i)); %Distances from marker 3 to the elbow (marker 10)
    dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i)); %Distances from marker 4 to the elbow (marker 10)
    % dists3(i)=pdist2(all_medians(6,:,i),all_medians_temp(3,:,i));
    % dists4(i)=pdist2(all_medians(6,:,i),all_medians_temp(4,:,i));
end


for i=1:n_times
    %dists3_forward is the most recently known distance from marker 3 to
    %the elbow (the distance from the last frame where neither marker is missing).
    if ~isnan(dists3(i))
        dists3_forward(i)=dists3(i);
    else
        dists3_forward(i)=dists3_forward(i-1);
    end
    %dists3_forward is the most recently known distance from marker 3 to
    %the elbow (the distance from the last frame where neither marker is missing).
    if ~isnan(dists4(i))
        dists4_forward(i)=dists4(i);
    else
        dists4_forward(i)=dists4_forward(i-1);
    end
    %Calculate the differences in distances between successive timepoints.
    if i>1
        dists33_forward(i)=abs(dists3_forward(i)-dists3_forward(i-1)); %The difference in the distance of marker 3 from in current frame and marker 3 from the prevous frame
        dists44_forward(i)=abs(dists4_forward(i)-dists4_forward(i-1)); %The difference in the distance of marker 4 from in current frame and marker 4 from the prevous frame
        dists34_forward(i)=abs(dists4_forward(i)-dists3_forward(i-1)); %The difference in the distance of marker 4 from in current frame and marker 3 from the prevous frame
        dists43_forward(i)=abs(dists3_forward(i)-dists4_forward(i-1)); %The difference in the distance of marker 3 from in current frame and marker 4 from the prevous frame
    end
  
    %According to the forward direction, we want to switch markers 3 and 4
    %(correct the incorrect assignment) when:
    %The marker was identified as marker 4, and the current marker 4
    %distance (to the elbow) is closer to the previous marker 3 distance
    %than the previous marker 4 distance.
    if ((dists34_forward(i)<dists44_forward(i)) & isnan(dists3(i)) & ~isnan(dists4(i)))
        change_forward(i)=1;
        dists3_forward(i)=dists4_forward(i);
        dists4_forward(i)=dists4_forward(i-1);
    end
    %The marker was identified as marker 3, and the current marker 3
    %distance (to the elbow) is closer to the previous marker 4 distance
    %than the previous marker 3 distance.
    if ((dists43_forward(i)<dists33_forward(i)) & ~isnan(dists3(i)) & isnan(dists4(i)))
        change_forward(i)=1;
        dists4_forward(i)=dists3_forward(i);
        dists3_forward(i)=dists3_forward(i-1);
    end
end

%Next we do the same as above, but backwards in time.

%Initializations for backwards run
dists3_backward=NaN(1,n_times);
dists4_backward=NaN(1,n_times);
dists33_backward=NaN(1,n_times);
dists44_backward=NaN(1,n_times);
dists34_backward=NaN(1,n_times);
dists43_backward=NaN(1,n_times);
change_backward=zeros(1,n_times);

for i=n_times-1:-1:1
    %dists3_backward is the next known distance from marker 3 to
    %the elbow (the distance from the next frame where neither marker is missing).
    if ~isnan(dists3(i))
        dists3_backward(i)=dists3(i);
    else
        dists3_backward(i)=dists3_backward(i+1);
    end
    %dists4_backward is the next known distance from marker 4 to
    %the elbow (the distance from the next frame where neither marker is missing).
    if ~isnan(dists4(i))
        dists4_backward(i)=dists4(i);
    else
        dists4_backward(i)=dists4_backward(i+1);
    end
    
    %Calculate the differences in distances between successive timepoints.
    if i<n_times
        dists33_backward(i)=abs(dists3_backward(i+1)-dists3_backward(i)); %The difference in the distance of marker 3 from in current frame and marker 3 from the next frame
        dists44_backward(i)=abs(dists4_backward(i+1)-dists4_backward(i)); %The difference in the distance of marker 4 from in current frame and marker 4 from the next frame
        dists34_backward(i)=abs(dists4_backward(i+1)-dists3_backward(i)); %The difference in the distance of marker 3 from in current frame and marker 4 from the next frame
        dists43_backward(i)=abs(dists3_backward(i+1)-dists4_backward(i)); %The difference in the distance of marker 4 from in current frame and marker 3 from the next frame
    end
    
    %According to the backward direction, we want to switch markers 3 and 4
    %(correct the incorrect assignment) when:
    %The marker was identified as marker 4, and the current marker 4
    %distance (to the elbow) is closer to the next marker 3 distance
    %than the next marker 4 distance.
    if ((dists43_backward(i)<dists44_backward(i)) & isnan(dists3(i)) & ~isnan(dists4(i)))
        change_backward(i)=1;
        dists3_backward(i)=dists4_backward(i);
        dists4_backward(i)=dists4_backward(i+1);
    end
    %The marker was identified as marker 3, and the current marker 3
    %distance (to the elbow) is closer to the next marker 4 distance
    %than the next marker 3 distance.
    if ((dists34_backward(i)<dists33_forward(i)) & ~isnan(dists3(i)) & isnan(dists4(i)))
        change_backward(i)=1;
        dists4_backward(i)=dists3_backward(i);
        dists3_backward(i)=dists3_backward(i+1);
    end
    
end

%Automatically make the switch if there should be a switch (according to
%above) in both the foward and backward directions. In the next section,
%we'll manually look at the cases of disagreement between the forward and
%backward directions.
switch34=change_forward & change_backward;
temp=all_medians(3,:,switch34);
temp2=all_medians2(3,:,switch34);
all_medians(3,:,switch34)=all_medians(4,:,switch34);
all_medians2(3,:,switch34)=all_medians2(4,:,switch34);
all_medians(4,:,switch34)=temp;
all_medians2(4,:,switch34)=temp2;

%% SET 12. Plot questionable ones

%Above, we switched the red markers when we were sure about it. That is
%when going forward and backward through time, both agreed that the markers
%should be switched. There are also questionable cases, when going froward
%or backward in time (but not both), it appears they should be swiched. We
%plot those cases here to determine which should be switched.

%Recalculate distances after above switches
for i=1:n_times    
    dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i));
    dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i));
end

wdw=20; %Time window for plots below

%Find times when the markers should be switched according to forward time
%or backward time (see previous section), but not both.
idxs=find((change_forward & ~change_backward) | (~change_forward & change_backward));

%Initialize points to switch
switch34_cell=cell(1,length(idxs));

%Plot those times
for i=1:length(idxs)
    figure;
    plot(idxs(i)-wdw:idxs(i)+wdw,dists1(idxs(i)-wdw:idxs(i)+wdw),'g-x');
    hold on;
    plot(idxs(i)-wdw:idxs(i)+wdw,dists2(idxs(i)-wdw:idxs(i)+wdw),'b-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists3(idxs(i)-wdw:idxs(i)+wdw),'r-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists4(idxs(i)-wdw:idxs(i)+wdw),'m-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists5(idxs(i)-wdw:idxs(i)+wdw),'c-x');
    title(num2str(idxs(i)));
    
    switch34_cell{i}=input('12. Time points to switch (red and pink) \n');
    
end

%% Set more red hand markers to switch (the questionable ones from above)

switch34=[switch34_cell{:}]; %The time points to switch
temp=all_medians(3,:,switch34);
temp2=all_medians(3,:,switch34);
all_medians(3,:,switch34)=all_medians(4,:,switch34);
all_medians2(3,:,switch34)=all_medians2(4,:,switch34);
all_medians(4,:,switch34)=temp;
all_medians2(4,:,switch34)=temp2;

%% SET 13. Determine green hand points to remove (based on having similar distance from elbow)

wdw=20; %Time window for plots below

%Find times when marker 1 and marker 5 are a similar distance to the elbow
%marker (which is a problem)
idxs=find(abs(dists1-dists5)<.01);

%Initialize points to be removed
green_test_cell=cell(1,length(idxs));
cyan_test_cell=cell(1,length(idxs));

%Plot those times
for i=1:length(idxs)
    figure;
    plot(idxs(i)-wdw:idxs(i)+wdw,dists1(idxs(i)-wdw:idxs(i)+wdw),'g-x');
    hold on;
    plot(idxs(i)-wdw:idxs(i)+wdw,dists2(idxs(i)-wdw:idxs(i)+wdw),'b-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists3(idxs(i)-wdw:idxs(i)+wdw),'r-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists4(idxs(i)-wdw:idxs(i)+wdw),'m-x');
    plot(idxs(i)-wdw:idxs(i)+wdw,dists5(idxs(i)-wdw:idxs(i)+wdw),'c-x');
    
    title(num2str(idxs(i)));
    
    green_test_cell{i}=input('13. Times cyan is in green area \n');
    cyan_test_cell{i}=input('13. Times green is in cyan area \n');
end

%% Set green hand points to remove

green_test=[green_test_cell{:}]; %Cyan is in green area
cyan_test=[cyan_test_cell{:}]; %Green is in cyan area

%% Determine which green hand point (1 or 5) to remove, at the above times

%Above we specified times when the two green markers had the same distance from
%the elbow.

%Sometimes both were had the expected distance for marker 1 (green
%points). In this case, we want to determine which of these markers really
%should be marker 1, and which should be removed. We do this by determining
%which is closer to the previous frame's marker 1.
if ~isempty(green_test)
    green_points=sort(green_test);
    for i=1:length(green_points)
        t=green_points(i);
        %Compare the distances from the current markers 1 and 5 to the
        %previous marker 1. Make the one that's closer the new marker 1.
        if pdist2(all_medians2(1,:,t-1),all_medians(1,:,t))>pdist2(all_medians2(1,:,t-1),all_medians(5,:,t))
            all_medians(1,:,t)=all_medians(5,:,t);
        end
        all_medians(5,:,t)=NaN;
    end
end

%Sometimes both were had the expected distance for marker 5 (cyan
%points). In this case, we want to determine which of these markers really
%should be marker 5, and which should be removed. We do this by determining
%which is closer to the previous frame's marker 5.
if ~isempty(cyan_test)
    cyan_points=sort(cyan_test);
    for i=1:length(cyan_points)
        t=cyan_points(i);
        %Compare the distances from the current markers 1 and 5 to the
        %previous marker 5. Make the one that's closer the new marker 5.
        if pdist2(all_medians2(5,:,t-1),all_medians(5,:,t))>pdist2(all_medians2(5,:,t-1),all_medians(1,:,t))
            all_medians(5,:,t)=all_medians(1,:,t);
        end
        all_medians(1,:,t)=NaN;
    end
end

%% Compare hand distances to shoulder (for times other markers are missing)

%Calculate the distances from marker 9 (green shoulder) to all the hand
%points. Note that all_medians2 is used for marker 9, since the marker will
%not move significantly across frames, and will allow us to have a value
%for every frame.
for i=1:n_times
    dists1(i)=pdist2(all_medians2(9,:,i),all_medians(1,:,i));
    dists2(i)=pdist2(all_medians2(9,:,i),all_medians(2,:,i));
    dists3(i)=pdist2(all_medians2(9,:,i),all_medians(3,:,i));
    dists4(i)=pdist2(all_medians2(9,:,i),all_medians(4,:,i));
    dists5(i)=pdist2(all_medians2(9,:,i),all_medians(5,:,i));
end

%Plot
if first_time
    figure; %Blue
    plot(dists2,'b-x');
    figure; %Red
    plot(dists3,'r-x');
    hold on;
    plot(dists4,'m-x');
    figure; %Green
    plot(dists1,'g-x');
    hold on;
    plot(dists5,'c-x');
end

%% SET 14. hand distance limits from the shoulder

if first_time
    
    str1='14A. Input green_keep \n';
    str2='Lower and upper limits of distances of the green hand markers (green and cyan above) to the green shoulder marker\n';
    str3='Value is generally ~ [.15,.45] \n';    
    green_keep=input([str1 str2 str3]);
    
    str1='14B. Input red_keep \n';
    str2='Lower and upper limits of distances of the red hand markers (red and magenta above) to the green shoulder marker\n';
    str3='Value is generally ~ [.15,.45] \n';    
    red_keep=input([str1 str2 str3]);
    
    str1='14C. Input blue_keep \n';
    str2='Lower and upper limits of distances of the blue hand markers (blue above) to the green shoulder marker\n';
    str3='Value is generally ~ [.15,.45] \n';    
    blue_keep=input([str1 str2 str3]);
    
end
%% Remove hand points (because they're too close or far from shoulder)

rmv1=dists1<green_keep(1) | dists1>green_keep(2);
all_medians(1,:,rmv1)=NaN;
rmv2=dists2<blue_keep(1) | dists2>blue_keep(2);
all_medians(2,:,rmv2)=NaN;
rmv3=dists3<red_keep(1) | dists3>red_keep(2);
all_medians(3,:,rmv3)=NaN;
rmv4=dists4<red_keep(1) | dists4>red_keep(2);
all_medians(4,:,rmv4)=NaN;
rmv5=dists5<green_keep(1) | dists5>green_keep(2);
all_medians(5,:,rmv5)=NaN;

%% Recalculate hand distances from the red elbow marker (in order to make comparisons of hand points)

for i=1:n_times
    dists1(i)=pdist2(all_medians(10,:,i),all_medians(1,:,i));
    dists2(i)=pdist2(all_medians(10,:,i),all_medians(2,:,i));
    dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i));
    dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i));
    dists5(i)=pdist2(all_medians(10,:,i),all_medians(5,:,i));
end

%% SET 15. Check whether elbow-5 distance is greater than elbow-3 distance
%Which it shouldn't be

wdw=20; %Time window for plots below

%Find times when marker 5 has a greater distance to the elbow
%than marker 3 (which is a problem)
idxs=find(dists5>dists3);

%Initialize points to be removed
rmv3_cell=cell(1,length(idxs));
rmv5_cell=cell(1,length(idxs));

%Plot those times
for i=1:length(idxs)
    if idxs(i)+wdw<=n_times & idxs(i)-wdw>0
        figure;
        plot(idxs(i)-wdw:idxs(i)+wdw,dists1(idxs(i)-wdw:idxs(i)+wdw),'g-x');
        hold on;
        plot(idxs(i)-wdw:idxs(i)+wdw,dists2(idxs(i)-wdw:idxs(i)+wdw),'b-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists3(idxs(i)-wdw:idxs(i)+wdw),'r-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists4(idxs(i)-wdw:idxs(i)+wdw),'m-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists5(idxs(i)-wdw:idxs(i)+wdw),'c-x');
        
        title(num2str(idxs(i)));
        
        rmv3_cell{i}=input('15. Times to remove hand point 3 (red) \n');
        rmv5_cell{i}=input('15. Times to remove hand point 5 (cyan) \n');
    end
end

%% Set marker 3 or 5 to remove based on above

%If cyan greater than red

rmv3=[rmv3_cell{:}]; %Times to remove marker of red hand point 3
rmv5=[rmv5_cell{:}]; %Times to remove marker of green hand point 5

all_medians(3,:,rmv3)=NaN;
all_medians(5,:,rmv5)=NaN;

%% SET 16. Check whether elbow-4 distance is greater than elbow-2 distance
%Which it shouldn't be

wdw=20; %Time window for plots below

%Find times when marker 4 has a greater distance to the elbow
%than marker2 (which is a problem)
idxs=find(dists4>dists2);

%Initialize points to be removed
rmv2_cell=cell(1,length(idxs));
rmv4_cell=cell(1,length(idxs));


%Plot those times
for i=1:length(idxs)
    if idxs(i)+wdw<=n_times & idxs(i)-wdw>0
        figure;
        plot(idxs(i)-wdw:idxs(i)+wdw,dists1(idxs(i)-wdw:idxs(i)+wdw),'g-x');
        hold on;
        plot(idxs(i)-wdw:idxs(i)+wdw,dists2(idxs(i)-wdw:idxs(i)+wdw),'b-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists3(idxs(i)-wdw:idxs(i)+wdw),'r-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists4(idxs(i)-wdw:idxs(i)+wdw),'m-x');
        plot(idxs(i)-wdw:idxs(i)+wdw,dists5(idxs(i)-wdw:idxs(i)+wdw),'c-x');
        title(num2str(idxs(i)));
        
        rmv2_cell{i}=input('16. Times to remove hand point 2 (blue) \n');
        rmv4_cell{i}=input('16. Times to remove hand point 4 (magenta) \n');
    end
end

%% Set marker 2 or 4 to remove based on above

%If magenta greater than blue

rmv2=[rmv2_cell{:}]; %Times to remove marker of blue hand point 2
rmv4=[rmv4_cell{:}]; %Times to remove marker of red hand point 4

all_medians(2,:,rmv2)=NaN;
all_medians(4,:,rmv4)=NaN;


%% Update all_medians2 to deal with removals

%If the marker is present at a given time, set all_medians2=all_medians
%If the marker isn't present at a given time, set all_medians 2 as
%all_medians from the previous frame

for j=1:10 %Loop through markers
    for t=1:n_times %Loop through times
        if ~isnan(all_medians(j,1,t))
            all_medians2(j,:,t)=all_medians(j,:,t);
        else
            all_medians2(j,:,t)=all_medians2(j,:,t-1);
        end
    end
end

%% Make matrices line up correctly with start time

%Make the final matrices begin at time 1 instead of time "start."

temp=all_medians;
temp2=all_medians2;

all_medians=NaN(11,3,finish);
all_medians2=NaN(11,3,finish);

all_medians(:,:,start:finish)=temp;
all_medians2(:,:,start:finish)=temp2;


%% Get smoothed time points

all_medians_smooth=NaN(size(all_medians));
for i=1:10
    for j=1:3
        temp=reshape(all_medians(i,j,:),[1,size(all_medians,3)]);
        all_medians_smooth(i,j,:)=medfilt1nan(temp,5);
    end
end


%% Hand Removals when there is no elbow

% for i=1:n_times
%     dists(i)=pdist2(all_medians(4,:,i),all_medians(5,:,i));
%     dists2(i)=pdist2(all_medians(3,:,i),all_medians(5,:,i));
%     dists3(i)=pdist2(all_medians(3,:,i),all_medians(4,:,i));
% end
% figure; plot(dists)
% hold on;
% plot(dists2)
% plot(dists3)


%% Update the hand markers

% num_clusts=6; %5 hand and 1 elbow
% num_dists=num_clusts*(num_clusts-1)/2;
%
% hand_dists=NaN(n_times,num_dists);
% hand_dists_matrix=NaN(num_clusts,num_clusts,n_times);
% for i=1:n_times
%     hand_dists(i,:)=pdist(all_medians([1:5 10],:,i));
%     %12 13 14 15 23 24 25 34 35 45 (when it was just 1-5)
%     hand_dists_matrix(:,:,i)=squareform(hand_dists(i,:));
% end
% med_hand_dists=nanmedian(hand_dists);
% med_hand_dists_matrix=squareform(med_hand_dists);
% med_hand_dists_matrix_reshape=repmat(med_hand_dists_matrix,[1,1,n_times]);
%
%
% perm1=[5 2 3 4 1 6]; %1 and 5 switched
% perm2=[1 2 4 3 5 6]; %3 and 4 switched
% perm3=[5 2 4 3 1 6]; %Both switched
%
% %Permutation of hand distances
% hand_dists_matrix1=hand_dists_matrix(perm1,perm1,:);
% hand_dists_matrix2=hand_dists_matrix(perm2,perm2,:);
% hand_dists_matrix3=hand_dists_matrix(perm3,perm3,:);
%
% num_hand_marks=reshape(sum(~isnan(all_medians(1:5,1,:))),[n_times,1]);
% elbow_weight=1./num_hand_marks;
% elbow_weight(elbow_weight==Inf)=1;
% weight_matrix=ones(size(hand_dists_matrix));
% % weight_matrix(6,1,:)=elbow_weight;
% % weight_matrix(6,2,:)=elbow_weight;
% % weight_matrix(6,3,:)=elbow_weight;
% % weight_matrix(6,4,:)=elbow_weight;
% % weight_matrix(6,5,:)=elbow_weight;
% % weight_matrix(6,6,:)=elbow_weight;
% % weight_matrix(1,6,:)=elbow_weight;
% % weight_matrix(2,6,:)=elbow_weight;
% % weight_matrix(3,6,:)=elbow_weight;
% % weight_matrix(4,6,:)=elbow_weight;
% % weight_matrix(5,6,:)=elbow_weight;
% % weight_matrix(6,6,:)=elbow_weight;
%
% %Differences of hand distances from expected (for actual and permutations)
%
% sse=reshape(nansum(nansum(weight_matrix.*(hand_dists_matrix-med_hand_dists_matrix_reshape).^2,1),2),[1,n_times]);
% sse1=reshape(nansum(nansum(weight_matrix.*(hand_dists_matrix1-med_hand_dists_matrix_reshape).^2,1),2),[1,n_times]);
% sse2=reshape(nansum(nansum(weight_matrix.*(hand_dists_matrix2-med_hand_dists_matrix_reshape).^2,1),2),[1,n_times]);
% sse3=reshape(nansum(nansum(weight_matrix.*(hand_dists_matrix3-med_hand_dists_matrix_reshape).^2,1),2),[1,n_times]);
%
% [~,m]=min([sse; sse1; sse2; sse3]);
%
%
%
%
%
%
%
% %     if m(i)==2 || m(i)==4
% %         temp=all_medians(5,:,i);
% %         all_medians(5,:,i)=all_medians(8,:,i);
% %         all_medians(8,:,i)=temp;
% %     end
% %     if m(i)==3 || m(i)==4
% %         temp=all_medians(6,:,i);
% %         all_medians(6,:,i)=all_medians(9,:,i);
% %         all_medians(9,:,i)=temp;
% %     end


%% save
if savefile
    date2=['20' num2str(date(7:8)) num2str(date(1:2)) num2str(date(4:5))];
    fname_save=[main_dir monkey '/Color_Tracking/' date '/Markers/markers_' monkey '_' date2 '_' exp '_' num];
    save(fname_save,'all_medians','all_medians2','led_vals','times');
    
    if first_time
        fname_save_settings=[main_dir monkey '/Color_Tracking/' date '/Markers/settings_' monkey '_' date2];
        save(fname_save_settings,'red_elbow_dist_from_blue','red_blue_arm_dist_max','red_arm_thresh1','red_arm_thresh2',...
        'green_hand_dists_elbow','red_hand_dists_elbow','blue_hand_dists_elbow','green_separator',...
        'green_hand_dists_bluearm','red_hand_dists_bluearm','blue_hand_dists_bluearm',...
        'green_hand_dists_redarm', 'red_hand_dists_redarm', 'blue_hand_dists_redarm',...
        'red_dist_min','green_dist_min','red_keep','green_keep','blue_keep','marker_inits');     
    end
end