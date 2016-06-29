% This script will show a video of the marker locations


%% Align

all_medians_aligned=NaN(11,3,finish);
all_medians_aligned(:,:,start:finish)=all_medians;

%% User initializations
markers=1:10; %The markers you want to plot
pause_manually=0; %Whether you want it to pause until you click after each frame
plot_original=0; %Whether the original points (pixels) should be plotted
start_frame=1000+start-1; %Frame you want to start at
finish_frame=2000+start-1; %Frame you want to end at

%% Plot

medians=all_medians_aligned; %Renamed (to make compatible with my previous plotting script)

%Set plotting limits as minimums and maximums
xlims=[min(min(medians(markers,1,:))) max(max(medians(markers,1,:)))];
ylims=[min(min(medians(markers,2,:))) max(max(medians(markers,2,:)))];
zlims=[min(min(medians(markers,3,:))) max(max(medians(markers,3,:)))];

%     xlims=[-.5 .5];
%     ylims=[-.4 .4];
%     zlims=[.7 1.3];

%Initialize figure
figure;
set(gca,'NextPlot','replacechildren');
xlim(zlims)
ylim(xlims)
zlim(ylims)

%Colors of each marker
plot_colors_all=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 0 0 0; 1 .5 0; .5 0 1; 0 1 .5; 0 0 0];
plot_colors=plot_colors_all(markers,:);


for i=start_frame:finish_frame %Loop through time
    
    %Plot Original Points
    if plot_original
        temp=color1{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(z,x,y,'b')
        hold on;
        temp=color2{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(z,x,y,'g')
        temp=color3{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(z,x,y,'r')
        temp=color4{i};
        x=temp(1:end/3);
        y=temp(end/3+1:2*end/3);
        z=temp(2*end/3+1:end);
        scatter3(z,x,y,'y')
    end
    
    %Plot Markers
    scatter3(medians(markers,3,i),medians(markers,1,i),medians(markers,2,i),200,plot_colors,'filled')
    hold off
    
    title(i); %Give time as title
    xlabel('z')
    ylabel('x')
    zlabel('y')
    
    %Set limits
    xlim(zlims)
    ylim(xlims)
    zlim(ylims)
    
    %Pause between frames
    if pause_manually
        pause;
    else
        pause(.03)
    end
end