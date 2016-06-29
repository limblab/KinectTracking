%Get arm positions (normally)
%Maybe I don't have to correct for the elbow??

%Get approximate distances to hand points from elbow
%--based on a few examples I type in

%Get hand markers based on these approximate distances (assignments will be wrong)

%Determine better hand marker distance metrics from elbow

%Get hand markers again
%--Only do it when we have the elbow markers
%--Start a segment when there are at least 3 markers
%--Stop the segment when there is 1 marker (or 0)
%--Do assignments like usual (but w/ extra blue marker)


%% Load ColorTracking File

main_dir='/Users/jig289/Box Sync/Tracking_Data/';


%File to load
monkey='Whistlepig';
date='11-24-15'; %mo-day-yr
exp='reach';
num='002';

fname_load=ls([main_dir monkey '/Color_Tracking/' date '/Tracking/color_tracking ' exp '_' num '*']);
load(deblank(fname_load));

% load(ls([main_dir monkey '/Color_Tracking/' date '/Tracking/color_tracking ' exp '_' num '*']))

%% Marker location reminder figure

marker_demo_locs=[0 0; 1 1; 1 -1; 2 1; 2 -1;...
    10 -1; 10 3; 10 6; 10 9; 9 0;...
    2 -3];
r=[1 0 0];
g=[0 1 0];
b=[0 1 1];
marker_demo_colors=[g; b; r; r; g; g; b; r; g; r; b];

figure; scatter(marker_demo_locs(:,1),marker_demo_locs(:,2),200,marker_demo_colors,'filled');
str={'1','2','3','4','5','6','7','8','9','10','11'};
text(marker_demo_locs(:,1),marker_demo_locs(:,2),str)
xlim([-5 15]);
ylim([-5 15]);

%% SET 0. Initializations


first_time=1; %If this is the first file from a date, set equal to 1 (there are more initializations)

%Load all of the settings if it's not the first file 
if ~first_time
    date2=['20' num2str(date(7:8)) num2str(date(1:2)) num2str(date(4:5))];
    fname_load_settings=[main_dir monkey '/Color_Tracking/' date '/Markers/settings_' monkey '_' date2];
    load(fname_load_settings);
end

%TIME INITIALIZATIONS
start=190; %Time point we're starting at
n=length(color1);
finish=n; %Time point we're finishing at
n_times=finish-start+1; %Number of time points (frames)

%MARKER NUMBER INITIALIZATIONS
red_arm_marker_ids=[8,10];
blue_arm_marker_ids=[7];
green_shoulder_marker_ids=[9]; %Sometimes empty
green_elbow_marker_ids=[6];
red_hand_marker_ids=[3,4];
blue_hand_marker_ids=[2,11]; %Sometimes just 2, sometimes 2 and 11
green_hand_marker_ids=[1,5];

%MARKER LOCATION INITIALIZATIONS
marker_inits=NaN(11,3);
marker_inits(1,:)=[0 0 0];
marker_inits(2,:)=[0 0 0];
marker_inits(3,:)=[0 0 0];
marker_inits(4,:)=[0 0 0];
marker_inits(5,:)=[0 0 0];
marker_inits(6,:)=[.86,.08,-.05];
marker_inits(7,:)=[.88,.04,-.04];
marker_inits(8,:)=[.89,-.02,-.04];
marker_inits(9,:)=[.91,-.07,-.05];
marker_inits(10,:)=[.87,.08,-.04];
marker_inits(11,:)=[0 0 0];

%I plot z,x,y (instead of x,y,z), so I input z,x,y above. Here, switch to x,y,z
marker_inits_temp=marker_inits;
marker_inits(:,1)=marker_inits_temp(:,2);
marker_inits(:,2)=marker_inits_temp(:,3);
marker_inits(:,3)=marker_inits_temp(:,1);

%Keeps track of all the cluster locations
all_medians=NaN(11,3,n_times); %Has NaNs when a marker is missing
all_medians2=NaN(11,3,n_times); %Estimates where markers are when they are missing

%Initialize some vectors that I use later for calculating the distance
%between points
dists=NaN(1,n_times);
dists1=NaN(1,n_times);
dists2=NaN(1,n_times);
dists3=NaN(1,n_times);
dists4=NaN(1,n_times);
dists5=NaN(1,n_times);

%% Plotting Initializations

fig=figure;
dcm_obj = datacursormode(fig);
set(gca,'NextPlot','replacechildren');
xlims=[-.5 .5];
ylims=[-.5 .5];
zlims=[0.5 1.5];

pause_time=.03;

%% Blue Arm

%Initializations
plot_on=1; %Whether to plot while it's running
marker_ids=blue_arm_marker_ids; %Set the marker_ids specified in "Initializations"
color=color1; %Blue=1, Green=2, Red=3
prev_meds=marker_inits(marker_ids,:); %Set initial "previous marker locations" as the start locations input in "Initializations"
num_clust=length(marker_ids); %Number of clusters
within_clust_dist1=0.08; %How close points must be to the previous frame's marker to be considered
dist_min=0.08; %Minimum distance between markers (cluster medians aren't allowed w/ distance < min_dist)

medians=NaN(num_clust,3,n_times); %Has NaNs when a marker is missing
medians2=NaN(num_clust,3,n_times); %Has previous known positions when a marker is missing

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
    
    %Remove points (those we're not keeping)
    rmv=~(keep1);
    
    %Actually remove the points
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
    
    
    %3. Plot original image and cluster centers
    if num_gone>2
    plot_colors=[1 0 0];
        temp=color1{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'b')
        hold on;
        temp=color2{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'g')
        temp=color3{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'r')
    scatter3(medians(:,1,t),medians(:,2,t),medians(:,3,t),200,plot_colors,'filled')
    hold off;
    xlim(xlims)
    ylim(ylims)
    zlim(zlims)
    title(i)
    pause(pause_time);
    
    
    k=waitforbuttonpress;
    if ~k
        f = getCursorInfo(dcm_obj);
        medians(1,:,t)=f.Position;
        medians2(1,:,t)=f.Position;
        prev_meds=f.Position;
        prev_num_clust=1;
        num_gone=0;
    end
    end
    
    
end

%Put the markers found here in the matrix of all markers
all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;


%% Green Shoulder

if ~isempty(green_shoulder_marker_ids) %Only do this if it is a file with a green shoulder marker
    %Initializations
    plot_on=0;
    marker_ids=green_shoulder_marker_ids;
    color=color2;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids);
    within_clust_dist1=.1;
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
        if plot_on
        plot_colors=[1 0 0];
        temp=color1{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'b')
        hold on;
        temp=color2{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'g')
        temp=color3{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'r')
        scatter3(medians(:,1,t),medians(:,2,t),medians(:,3,t),200,plot_colors,'filled')
        hold off;
        xlim(xlims)
        ylim(ylims)
        zlim(zlims)
        title(i)
        pause(pause_time);
        end
    end
    
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
end


%% Red Arm

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
%% Plot Red elbow to Blue Arm Distance

if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Calculate distances for each time point
    for i=1:n_times
        dists(i)=pdist2(all_medians(10,:,i),all_medians(7,:,i)); %Distance between markers 7 and 10 (blue arm and red elbow)
        dists2(i)=pdist2(all_medians(8,:,i),all_medians(7,:,i)); %Distance between markers 7 and 9 (blue arm and red arm)
    end
    
    %Plot
    figure; plot(dists); 
    hold on;
    plot(dists2)
    legend('7-10','7-9')
end
%% SET 1: Red Elbow/Arm to Blue Arm Distance

if first_time
    
    str1='1A. Input red_elbow_dist_from_blue \n';
    str2='The blue values in the above plot should be below this value (the red elbow should always be within this distance of the blue arm)\n';
    str3='Value is generally ~ .05 \n';    
    red_elbow_dist_from_blue=input([str1 str2 str3]);
    
    str1='1B. Input red_blue_arm_dist_max \n';
    str2='%All values in above plot should be below this value (Maximum distance from a red arm point to the blue)\n';
    str3='Value is generally ~ .08 \n';   
    red_blue_arm_dist_max=input([str1 str2 str3]);
        
end


%% Plot Red elbow/arm to Green Shoulder Distance

if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Calculate distances for each time point
    for i=1:n_times
        dists(i)=pdist2(all_medians(10,:,i),all_medians2(9,:,i)); %Distance between markers 9 and 10 (green arm and red elbow)
        dists2(i)=pdist2(all_medians(8,:,i),all_medians2(9,:,i)); %Distance between markers 9 and 8 (green arm and red arm)
    end
    
    %Plot
    figure; plot(dists); 
    hold on;
    plot(dists2)
    legend('9-10','9-8')
end
%% Set threshold of red arm points from green shoulder

red_green_thresh=.11;

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
    
    %Correct assignments
    
    %If one missing
    if (isnan(medians(1,1,t)) && ~isnan(medians(2,1,t)))
        if pdist2(medians(2,:,t),all_medians2(9,:,t))<red_green_thresh
            temp=medians(1,:,t);
            temp2=medians2(1,:,t);
            medians(1,:,t)=medians(2,:,t);
            medians(2,:,t)=temp;
            medians2(1,:,t)=medians2(2,:,t);
            medians2(2,:,t)=temp2;
        end
    end
    if (~isnan(medians(1,1,t)) && isnan(medians(2,1,t)))
        if pdist2(medians(1,:,t),all_medians2(9,:,t))>red_green_thresh
            temp=medians(1,:,t);
            temp2=medians2(1,:,t);
            medians(1,:,t)=medians(2,:,t);
            medians(2,:,t)=temp;
            medians2(1,:,t)=medians2(2,:,t);
            medians2(2,:,t)=temp2;
        end
    end
    
    %If both there
    if pdist2(medians(1,:,t),all_medians2(9,:,t))>pdist2(medians(2,:,t),all_medians2(9,:,t))
        temp=medians(1,:,t);
        temp2=medians2(1,:,t);
        medians(1,:,t)=medians(2,:,t);
        medians(2,:,t)=temp;
        medians2(1,:,t)=medians2(2,:,t);
        medians2(2,:,t)=temp2;
    end
    
    
%     if ~isnan(all_medians(9,1,t))
%             if pdist2(medians(1,:,t),all_medians(9,:,t))<pdist2(medians(2,:,t),all_medians(9,:,t))
%                 temp=medians(1,:,t);
%                 temp2=medians2(1,:,t);
%                 medians(1,:,t)=medians(2,:,t);
%                 medians(2,:,t)=temp;
%                 medians2(1,:,t)=medians2(2,:,t);
%                 medians2(2,:,t)=temp2;
%             end
%     else if ~isnan(all_medians(7,1,t))
%             if pdist2(medians(1,:,t),all_medians(7,:,t))>pdist2(medians(2,:,t),all_medians(7,:,t))
%                 temp=medians(1,:,t);
%                 temp2=medians2(1,:,t);
%                 medians(1,:,t)=medians(2,:,t);
%                 medians(2,:,t)=temp;
%                 medians2(1,:,t)=medians2(2,:,t);
%                 medians2(2,:,t)=temp2;
%             end
%         end
%     end
    
    
    %3. Plot original image and cluster centers
    plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
    
end

all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;



%% REMOVE RED ARM POINTS


%% Plot red elbow angle relative to arm (based on angle)

%This calculates (and plots) the angle made by points 7,8,10
%Problems with the red elbow marker (point 10) will make this angle wrong

%Note- look back at version 10 for code that has to do with distance from
%the line, instead of the angle

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
if first_time   
    figure; plot(angle)   
end
%% SET 3. red elbow points to remove (based on angle)

% red_elbow_angle_thresh=2.55;
red_elbow_angle_thresh=nanmean(angle)-4*nanstd(angle); %Frames with an angle below this will have marker 10 removed
%% Remove red elbow points (based on angle)

rmv10=angle<red_elbow_angle_thresh;
all_medians(10,:,rmv10)=NaN;

%% Green Elbow

%Initializations
plot_on=0;
marker_ids=green_elbow_marker_ids;
color=color2;
prev_meds=marker_inits(marker_ids,:);
num_clust=length(marker_ids);
within_clust_dist1=.12;
dist_min=0.1;

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
%     plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
    if plot_on
        plot_colors=[0.5 0.5 0];
        temp=color1{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'b')
        hold on;
        temp=color2{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'g')
        temp=color3{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(x,y,z,'r')
        scatter3(medians(:,1,t),medians(:,2,t),medians(:,3,t),200,plot_colors,'filled')
        hold off;
        xlim(xlims)
        ylim(ylims)
        zlim(zlims)
        title(i)
        pause(pause_time);
        end
    
end

all_medians(marker_ids,:,:)=medians;
all_medians2(marker_ids,:,:)=medians2;

%% Plot green elbow angle relative to arm (based on angle)

%This calculates (and plots) the angle made by points 7,8,6
%Problems with the green elbow marker (point 6) will make this angle wrong

%Note- look back at version 10 for code that has to do with distance from
%the line, instead of the angle


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
if first_time   
    figure; plot(angle)   
end
%% SET 4. green elbow points to remove (based on angle)

% green_elbow_angle_thresh=2.55;
green_elbow_angle_thresh=nanmean(angle)-4*nanstd(angle); %Frames with an angle below this will have marker 6 removed
%% Remove green elbow points

rmv6=angle<green_elbow_angle_thresh;
all_medians(6,:,rmv6)=NaN;

%Note that in version 10 and previous, there were other attempted methods
%for removing green elbow problems (based on distance to the other arm
%points)



%% HAND

%% Very Approximate distances from elbow
    
green_hand_dists_elbow=[.15 .29];   
red_hand_dists_elbow=[.15 .27];      
blue_hand_dists_elbow=[.15 .27];
    

%% Green Hand
if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Initializations
    plot_on=0;
    marker_ids=green_hand_marker_ids;
    color=color2;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids); 
    within_clust_dist1=.05;  %How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
    within_clust_dist2=.05; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered
    dist_min=0.04;
    
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
        
        %Remove points that are too close or far from the red elbow marker
        D2=pdist2(loc,all_medians(10,:,t));
        rmv0=D2<green_hand_dists_elbow(1) | D2>green_hand_dists_elbow(2);
        
        D3=pdist2(loc,all_medians(6,:,t));
        rmv1=D3<green_hand_dists_elbow(1) | D3>green_hand_dists_elbow(2);
        
        %Remove points that are close to the shoulder marker
        D4=pdist2(loc,all_medians(9,:,t));
        rmv2=D4<.05;

        %Remove "reflection" points
        rmv3=((z<all_medians(10,3,t)-.05) | (z<all_medians(6,3,t)-.05)) & z<.75;
        
        
        %Remove bad areas (write y,z,x in terms of other plot)
        bad_points1=[-.032 -.234 .7];
        bad_points_min1=bad_points1-[.025 .01 .015];
        bad_points_max1=bad_points1+[.025 .01 .015];       
        rmv_bad_points1=y>bad_points_min1(2) & y<bad_points_max1(2) & z>bad_points_min1(3) & z<bad_points_max1(3); %I don't include the x variable
            
        
        %Remove points (those we're not keeping)
        rmv= rmv0 | rmv1 | rmv2 | rmv3' | rmv_bad_points1';
        
        
        if isnan(all_medians(6,1,t)) && isnan(all_medians(10,1,t))
            rmv=1:length(x);
        end
        
        %Actually remove the points
        loc(rmv,:)=[];
        
        
        %2. Cluster and assign
        %Note that this uses "cluster_func" instead of "cluster_func2"
        %which is slightly faster but less accurate. This is because we
        %will be redoing this later with cluster_func2. This current run is
        %only to determine the distances from the red hand markers to arm
        %markers (which will help in the next run)
        [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func_reach2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
               
        %3. Plot original image and cluster centers
        plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
        
    end
    
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
    
end


%% Calculate green hand distances from the elbow

dists1=NaN(1,n_times);
dists5=NaN(1,n_times);

%Calculate distances from red elbow to points on the hand
for i=1:n_times
    dists1(i)=pdist2(all_medians(10,:,i),all_medians(1,:,i));
    dists5(i)=pdist2(all_medians(10,:,i),all_medians(5,:,i));

end

%%
%Plot
figure;
plot(dists1,'g-x');
hold on;
plot(dists5,'c-x');


%% Switch green hand markers

    green_thresh=0.2;

    %If one missing
    
    for t=1:n_times
        
        if (isnan(all_medians(1,1,t)) && ~isnan(all_medians(5,1,t)))
            if dists5(t)>green_thresh
                temp=all_medians(1,:,t);
                all_medians(1,:,t)=all_medians(5,:,t);
                all_medians(5,:,t)=temp;
            end
        end
        if (~isnan(all_medians(1,1,t)) && isnan(all_medians(5,1,t)))
            if dists1(t)<green_thresh
                temp=all_medians(1,:,t);
                all_medians(1,:,t)=all_medians(5,:,t);
                all_medians(5,:,t)=temp;
            end
        end
        
        %If both there
        if dists1(t)<dists5(t)
            temp=all_medians(1,:,t);
            all_medians(1,:,t)=all_medians(5,:,t);
            all_medians(5,:,t)=temp;
        end
        
    end

    %% Redo w/ distance from other elbow marker
    
dists1=NaN(1,n_times);
dists5=NaN(1,n_times);

%Calculate distances from red elbow to points on the hand
for i=1:n_times
    dists1(i)=pdist2(all_medians(6,:,i),all_medians(1,:,i));
    dists5(i)=pdist2(all_medians(6,:,i),all_medians(5,:,i));
end   
    
figure;
plot(dists1,'g-x');
hold on;
plot(dists5,'c-x');    
    


%% Switch green hand markers

    green_thresh=0.215;

    %If one missing
    
    for t=1:n_times
        
        if (isnan(all_medians(1,1,t)) && ~isnan(all_medians(5,1,t)))
            if dists5(t)>green_thresh
                temp=all_medians(1,:,t);
                all_medians(1,:,t)=all_medians(5,:,t);
                all_medians(5,:,t)=temp;
            end
        end
        if (~isnan(all_medians(1,1,t)) && isnan(all_medians(5,1,t)))
            if dists1(t)<green_thresh
                temp=all_medians(1,:,t);
                all_medians(1,:,t)=all_medians(5,:,t);
                all_medians(5,:,t)=temp;
            end
        end
        
        %If both there
        if dists1(t)<dists5(t)
            temp=all_medians(1,:,t);
            all_medians(1,:,t)=all_medians(5,:,t);
            all_medians(5,:,t)=temp;
        end
        
    end

 



%% Red Hand
if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Initializations
    plot_on=0;
    marker_ids=red_hand_marker_ids;
    color=color3;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids); 
    within_clust_dist1=.05;  %How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
    within_clust_dist2=.05; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered
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
        
        %Remove points that are too close or far from the red elbow marker
        D2=pdist2(loc,all_medians(10,:,t));
        D4=pdist2(all_medians(1,:,t),all_medians(10,:,t)); %And Remove points farther than the green end-of-hand marker
        rmv0=D2<red_hand_dists_elbow(1) | D2>red_hand_dists_elbow(2) | D2>D4;
        
        D3=pdist2(loc,all_medians(6,:,t));
        D5=pdist2(all_medians(1,:,t),all_medians(6,:,t)); %And Remove points farther than the green end-of-hand marker
        rmv1=D3<red_hand_dists_elbow(1) | D3>red_hand_dists_elbow(2) | D3>D5;
        
        %Remove "reflection" points
        rmv2=z<.8 & z<all_medians(1,3,t)-.04;
        
        %Remove "reflection" points
        rmv3=((z<all_medians(10,3,t)-.05) | (z<all_medians(6,3,t)-.05)) & z<.75;
        
        
        %Remove bad areas (write y,z,x in terms of other plot)
        bad_points1=[-.006 0.032 1.02];        
        bad_points_min1=bad_points1-.01;
        bad_points_max1=bad_points1+.01;       
        rmv_bad_points1=x>bad_points_min1(1) & x<bad_points_max1(1) & y>bad_points_min1(2) & y<bad_points_max1(2) & z>bad_points_min1(3) & z<bad_points_max1(3);
            
        bad_points2=[-.032 -.234 .7];        
        bad_points_min2=bad_points2-[.025 .01 .015];
        bad_points_max2=bad_points2+[.025 .01 .015];       
        rmv_bad_points2=y>bad_points_min2(2) & y<bad_points_max2(2) & z>bad_points_min2(3) & z<bad_points_max2(3); %I don't include the x variable
       
        bad_points3=[.043 -.254 .765];        
        bad_points_min3=bad_points3-.01;
        bad_points_max3=bad_points3+.01;       
        rmv_bad_points3=x>bad_points_min3(1) & x<bad_points_max3(1) & y>bad_points_min3(2) & y<bad_points_max3(2) & z>bad_points_min3(3) & z<bad_points_max3(3);
           
        
        %Remove points (those we're not keeping)
        rmv= rmv0 | rmv1 | rmv2' | rmv3' | rmv_bad_points1' | rmv_bad_points2';
        
        if isnan(all_medians(6,1,t)) && isnan(all_medians(10,1,t))
            rmv=1:length(x);
        end
        
        %Actually remove the points
        loc(rmv,:)=[];
        
        
        %2. Cluster and assign
        %Note that this uses "cluster_func" instead of "cluster_func2"
        %which is slightly faster but less accurate. This is because we
        %will be redoing this later with cluster_func2. This current run is
        %only to determine the distances from the red hand markers to arm
        %markers (which will help in the next run)
        [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func_reach2(t, loc, num_clust, prev_num_clust, dist_min, prev_meds, medians, medians2 );
               
        %3. Plot original image and cluster centers
        plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
        
        
        %Correct if both there (switch based on distance from the elbow)
    if pdist2(medians(1,:,t),all_medians2(10,:,t))<pdist2(medians(2,:,t),all_medians2(10,:,t))
        temp=medians(1,:,t);
        temp2=medians2(1,:,t);
        medians(1,:,t)=medians(2,:,t);
        medians(2,:,t)=temp;
        medians2(1,:,t)=medians2(2,:,t);
        medians2(2,:,t)=temp2;
    end
    if pdist2(medians(1,:,t),all_medians2(6,:,t))<pdist2(medians(2,:,t),all_medians2(6,:,t))
        temp=medians(1,:,t);
        temp2=medians2(1,:,t);
        medians(1,:,t)=medians(2,:,t);
        medians(2,:,t)=temp;
        medians2(1,:,t)=medians2(2,:,t);
        medians2(2,:,t)=temp2;
    end    
        
    end
    
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
    
end



%% Blue Hand
if first_time %If this is not the first file from a date, we don't need to run this.
    
    %Initializations
    plot_on=0;
    marker_ids=blue_hand_marker_ids;
    color=color1;
    prev_meds=marker_inits(marker_ids,:);
    num_clust=length(marker_ids); 
    within_clust_dist1=.05;  %How close points must be to the previous frame's first marker, # marker_ids(1), to be considered
    within_clust_dist2=.05; %How close points must be to the previous frame's second marker, # marker_ids(2), to be considered
    dist_min=0.04;
    
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
        
        %Remove points that are too close or far from the red elbow marker
        D2=pdist2(loc,all_medians(10,:,t));
        D4=pdist2(all_medians(1,:,t),all_medians(10,:,t)); %And Remove points farther than the green end-of-hand marker
        rmv0=D2<blue_hand_dists_elbow(1) | D2>blue_hand_dists_elbow(2) | D2>D4;
        
        D3=pdist2(loc,all_medians(6,:,t));
        D5=pdist2(all_medians(1,:,t),all_medians(6,:,t)); %And Remove points farther than the green end-of-hand marker
        rmv1=D3<blue_hand_dists_elbow(1) | D3>blue_hand_dists_elbow(2) | D3>D5;
        
        
        
        
        %Remove points (those we're not keeping)
        rmv= rmv0 | rmv1;
        
        if isnan(all_medians(6,1,t)) && isnan(all_medians(10,1,t))
            rmv=1:length(x);
        end
        
        %Actually remove the points
        loc(rmv,:)=[];
        
        
        %2. Cluster and assign
        %Note that this uses "cluster_func" instead of "cluster_func2"
        %which is slightly faster but less accurate. This is because we
        %will be redoing this later with cluster_func2. This current run is
        %only to determine the distances from the red hand markers to arm
        %markers (which will help in the next run)
        [ prev_num_clust, prev_meds, medians, medians2  ] = cluster_func_reach(t, loc, num_clust, prev_num_clust, dist_min, .05, prev_meds, medians, medians2 );
               
        %3. Plot original image and cluster centers
        plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
        
    end
    
    all_medians(marker_ids,:,:)=medians;
    all_medians2(marker_ids,:,:)=medians2;
    
end













%% Plot Red elbow/arm to Green Shoulder Distance

    
    num_present=nansum(~isnan(all_medians([1 2 3 4 5 11],1,1:n_times)));
    




%% Determine times to get hand points

%Do medfilt to determine times to even consider...
all_medians_smooth=NaN(size(all_medians));
for i=1:11
    for j=1:3
        temp=reshape(all_medians(i,j,:),[1,size(all_medians,3)]);
        all_medians_smooth(i,j,:)=medfilt1nan(temp,5);
    end
end

num_present_smooth=nansum(~isnan(all_medians_smooth([1 2 3 4 5 11],1,1:n_times)));

good_frames_temp=find(reshape(num_present_smooth>=3,[1 n_times]));
% good_frames_temp=unique([good_frames-2 good_frames-1 good_frames good_frames+1 good_frames+2]);
% good_frames(good_frames<1 | good_frames>n_times)=[];

good_frames_temp_logical=zeros(1,n_times);
good_frames_temp_logical(good_frames_temp)=1;

%Find when there's 5 good frames in a row
temp=round(5*(smooth(good_frames_temp_logical,5)))';
good_frames_original=find(temp==5);
good_frames=unique([good_frames_original-2 good_frames_original-1 good_frames_original good_frames_original+1 good_frames_original+2]);
good_frames(good_frames<1 | good_frames>n_times)=[];

%% Calculate hand distances from the elbow

dists1=NaN(1,n_times);
dists2=NaN(1,n_times);
dists3=NaN(1,n_times);
dists4=NaN(1,n_times);
dists5=NaN(1,n_times);
dists11=NaN(1,n_times);

%Calculate distances from red elbow to points on the hand
for i=good_frames
    dists1(i)=pdist2(all_medians(10,:,i),all_medians(1,:,i));
    dists2(i)=pdist2(all_medians(10,:,i),all_medians(2,:,i));
    dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i));
    dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i));
    dists5(i)=pdist2(all_medians(10,:,i),all_medians(5,:,i));
    dists11(i)=pdist2(all_medians(10,:,i),all_medians(11,:,i));
    
%     dists1(i)=pdist2(all_medians(6,:,i),all_medians(1,:,i));
%     dists2(i)=pdist2(all_medians(6,:,i),all_medians(2,:,i));
%     dists3(i)=pdist2(all_medians(6,:,i),all_medians(3,:,i));
%     dists4(i)=pdist2(all_medians(6,:,i),all_medians(4,:,i));
%     dists5(i)=pdist2(all_medians(6,:,i),all_medians(5,:,i));
%     dists11(i)=pdist2(all_medians(6,:,i),all_medians(11,:,i));

end

%%
%Plot
figure;
plot(dists1,'g-x');
hold on;
plot(dists5,'c-x');

figure;
plot(dists2,'b-x');
hold on;
plot(dists11,'k-x');

figure;
plot(dists3,'r-x');
hold on;
plot(dists4,'m-x');
    
      

%% Manual removal

all_medians(4,:,5160:5240)=NaN;



%% SET 11. Determine red hand points to remove (based on having similar distance from elbow)

%Calculate distances of red hand markers 3 and 4, and the red elbow

dists3=NaN(1,n_times);
dists4=NaN(1,n_times);

for i=good_frames
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


%% Calculate hand distances from the elbow

dists=NaN(1,n_times);
dists1=NaN(1,n_times);
dists2=NaN(1,n_times);
dists3=NaN(1,n_times);
dists4=NaN(1,n_times);
dists5=NaN(1,n_times);
dists11=NaN(1,n_times);

%Calculate distances from red elbow to points on the hand
for i=good_frames
    dists1(i)=pdist2(all_medians(10,:,i),all_medians(1,:,i));
    dists2(i)=pdist2(all_medians(10,:,i),all_medians(2,:,i));
    dists3(i)=pdist2(all_medians(10,:,i),all_medians(3,:,i));
    dists4(i)=pdist2(all_medians(10,:,i),all_medians(4,:,i));
    dists5(i)=pdist2(all_medians(10,:,i),all_medians(5,:,i));
    dists11(i)=pdist2(all_medians(10,:,i),all_medians(11,:,i));
    
    % dists1(i)=pdist2(all_medians(6,:,i),all_medians(1,:,i));
    % dists2(i)=pdist2(all_medians(6,:,i),all_medians(2,:,i));
    % dists3(i)=pdist2(all_medians(6,:,i),all_medians(3,:,i));
    % dists4(i)=pdist2(all_medians(6,:,i),all_medians(4,:,i));
    % dists5(i)=pdist2(all_medians(6,:,i),all_medians(5,:,i));

end    
    
    
%% Switch red hand points when one is missing (the assignment was wrong)

%Here, we'll determine what time points to switch the red hand markers (3
%and 4) when one of those markers was missing. That is, there was a single
%marker that was incorrectly assigned.
%We will determine this based on having a expected/unexpected distance from
%the red elbow marker. We will go both forward and backward in time, to see
%whether switching the marker assignment would lead to distances from the
%elbow that are more consistent with the previous time.

start_frame=find(~isnan(dists3) & ~isnan(dists4),1);
end_frame=find(~isnan(dists3) & ~isnan(dists4),1,'last');



%Initializations for forward run
dists3_forward=NaN(1,n_times);
dists4_forward=NaN(1,n_times);
dists33_forward=NaN(1,n_times);
dists44_forward=NaN(1,n_times);
dists34_forward=NaN(1,n_times);
dists43_forward=NaN(1,n_times);
change_forward=zeros(1,n_times);

for i=start_frame:n_times
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

for i=end_frame-1:-1:1
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
idxs=find(((change_forward & ~change_backward) | (~change_forward & change_backward)) & good_frames_logical);

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
    
    
       
    
 %% Assign blue hand markers
 
 
%% Distances from red elbow 
dists2=NaN(1,n_times);
dists11=NaN(1,n_times);

%Calculate distances from red elbow to points on the hand
for i=good_frames
    dists2(i)=pdist2(all_medians(2,:,i),all_medians(10,:,i));
    dists11(i)=pdist2(all_medians(11,:,i),all_medians(10,:,i));
end   
    
figure;
plot(dists2,'b-x');
hold on;
plot(dists11,'k-x');  
  
%% Switch blue hand markers

    blue_thresh=0.18;

    %If one missing
    
    for t=good_frames
        
        if (isnan(all_medians(2,1,t)) && ~isnan(all_medians(11,1,t)))
            if dists11(t)>blue_thresh
                temp=all_medians(2,:,t);
                all_medians(2,:,t)=all_medians(11,:,t);
                all_medians(11,:,t)=temp;
            end
        end
        if (~isnan(all_medians(2,1,t)) && isnan(all_medians(11,1,t)))
            if dists2(t)<blue_thresh
                temp=all_medians(2,:,t);
                all_medians(2,:,t)=all_medians(11,:,t);
                all_medians(11,:,t)=temp;
            end
        end
        
        %If both there
        if dists2(t)<dists11(t)
            temp=all_medians(2,:,t);
            all_medians(2,:,t)=all_medians(11,:,t);
            all_medians(11,:,t)=temp;
        end
        
    end

%% Distances from green elbow 
dists2=NaN(1,n_times);
dists11=NaN(1,n_times);

%Calculate distances from red elbow to points on the hand
for i=good_frames
    dists2(i)=pdist2(all_medians(2,:,i),all_medians(6,:,i));
    dists11(i)=pdist2(all_medians(11,:,i),all_medians(6,:,i));
end   
    
figure;
plot(dists2,'b-x');
hold on;
plot(dists11,'k-x');  
  
%% Switch blue hand markers

    blue_thresh=0.19;

    %If one missing
    
    for t=good_frames
        
        if (isnan(all_medians(2,1,t)) && ~isnan(all_medians(11,1,t)))
            if dists11(t)>blue_thresh
                temp=all_medians(2,:,t);
                all_medians(2,:,t)=all_medians(11,:,t);
                all_medians(11,:,t)=temp;
            end
        end
        if (~isnan(all_medians(2,1,t)) && isnan(all_medians(11,1,t)))
            if dists2(t)<blue_thresh
                temp=all_medians(2,:,t);
                all_medians(2,:,t)=all_medians(11,:,t);
                all_medians(11,:,t)=temp;
            end
        end
        
        %If both there
        if dists2(t)<dists11(t)
            temp=all_medians(2,:,t);
            all_medians(2,:,t)=all_medians(11,:,t);
            all_medians(11,:,t)=temp;
        end
        
    end

    
    
    
    
    
    
    
    %%
%     
%     temp1=reshape(~isnan(all_medians(1,1,:)),[1 n_times]);
%     temp2=reshape(~isnan(all_medians(2,1,:)),[1 n_times]);
%     temp3=reshape(~isnan(all_medians(3,1,:)),[1 n_times]);
%     temp4=reshape(~isnan(all_medians(4,1,:)),[1 n_times]);
%     temp5=reshape(~isnan(all_medians(5,1,:)),[1 n_times]);
%     temp11=reshape(~isnan(all_medians(11,1,:)),[1 n_times]);
%     
%     idx=find(temp11(good_frames) & temp5(good_frames))
    
    
  
%% Switch blue hand markers (manually enter)

    switch_idx=[572:576 2187:2189 3803:3807 3100 3102 3103 4570 4571 4573 4576 4967:4969 4971:4973 4976 8401:8406 11790:11791];

    %If one missing
    
    for t=switch_idx
            temp=all_medians(2,:,t);
            all_medians(2,:,t)=all_medians(11,:,t);
            all_medians(11,:,t)=temp;
    end
        
    
 
    
    
%% Remove blue hand markers (manually enter)

rmv_idx=[370 495 513 2718 1018:1021 7995 16805];
    
for t=rmv_idx
            
            all_medians(2,:,t)=NaN;
    end    
    
    
    
    
    %% Update all_medians2 to deal with removals

%If the marker is present at a given time, set all_medians2=all_medians
%If the marker isn't present at a given time, set all_medians 2 as
%all_medians from the previous frame

for j=1:11 %Loop through markers
    for t=1:n_times %Loop through times
        if ~isnan(all_medians(j,1,t))
            all_medians2(j,:,t)=all_medians(j,:,t);
        else
            if t==1
                all_medians2(j,:,t)=0;
            else
                all_medians2(j,:,t)=all_medians2(j,:,t-1);
            end
        end
    end
end


%% redo goodframes


%Do medfilt to determine times to even consider...
all_medians_smooth=NaN(size(all_medians));
for i=1:11
    for j=1:3
        temp=reshape(all_medians(i,j,:),[1,size(all_medians,3)]);
        all_medians_smooth(i,j,:)=medfilt1nan(temp,5);
    end
end

num_present_smooth=nansum(~isnan(all_medians_smooth([1 2 3 4 5 11],1,1:n_times)));

good_frames_temp=find(reshape(num_present_smooth>=3,[1 n_times]));
% good_frames_temp=unique([good_frames-2 good_frames-1 good_frames good_frames+1 good_frames+2]);
% good_frames(good_frames<1 | good_frames>n_times)=[];

good_frames_temp_logical=zeros(1,n_times);
good_frames_temp_logical(good_frames_temp)=1;

%Find when there's 5 good frames in a row
temp=round(5*(smooth(good_frames_temp_logical,5)))';
good_frames_original=find(temp==5);
good_frames=unique([good_frames_original-2 good_frames_original-1 good_frames_original good_frames_original+1 good_frames_original+2]);
good_frames(good_frames<1 | good_frames>n_times)=[];


    
%% Make matrices line up correctly with start time

%Make the final matrices begin at time 1 instead of time "start."

all_medians_orig=all_medians;
all_medians_orig2=all_medians2;

all_medians_good=NaN(11,3,n_times);
all_medians_good2=NaN(11,3,n_times);

all_medians_good(:,:,good_frames)=all_medians(:,:,good_frames);
all_medians_good2(:,:,good_frames)=all_medians2(:,:,good_frames);

all_medians=NaN(11,3,finish);
all_medians2=NaN(11,3,finish);

all_medians(:,:,start:finish)=all_medians_good;
all_medians2(:,:,start:finish)=all_medians_good2;








    %% Save
    

    
savefile=1;
if savefile
    date2=['20' num2str(date(7:8)) num2str(date(1:2)) num2str(date(4:5))];
    fname_save=[main_dir monkey '/Color_Tracking/' date '/Markers/markers_' monkey '_' date2 '_' exp '_' num 'original'];
    save(fname_save,'all_medians_orig','all_medians_orig2','led_vals','times');    
    
    
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% Testing things

fig=figure;
dcm_obj = datacursormode(fig);
for i=1:5
scatter3(rand(1,10),rand(1,10),rand(1,10));
k=waitforbuttonpress;
if ~k
f = getCursorInfo(dcm_obj);
f.Position
end
end
% for i=start:5000
%     figure;
%     h=scatter3(rand(1,10),rand(1,10),rand(1,10));
%     points=get(h,'Children');
% %     pos=get(gca,'CurrentPoint');
% 
% % end
% %%
% 
% 
% h=scatter(x,y,S,C); points=get(h,'Children');
% for i=1:numel(children)
% set(points(i),'HitTest','on','ButtonDownFcn',{'myFunction',i});
% end