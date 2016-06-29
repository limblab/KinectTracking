%% Align

all_medians_aligned=NaN(11,3,finish);

all_medians_aligned(:,:,start:finish)=all_medians;


%% Set things
markers=[7 9 10];%1:11;
pause_manually=1;
plot_original=1;
start_frame=2625+start-1;
finish=n;

%% Plot

medians=all_medians_aligned; %Renamed (to make compatible with my previous plotting script)

%Set plotting limits as minimums and maximums
xlims=[min(min(medians(markers,1,:))) max(max(medians(markers,1,:)))];
ylims=[min(min(medians(markers,2,:))) max(max(medians(markers,2,:)))];
zlims=[min(min(medians(markers,3,:))) max(max(medians(markers,3,:)))];

%Initialize figure
figure;
set(gca,'NextPlot','replacechildren');
xlim(zlims)
ylim(xlims)
zlim(ylims)

%Colors of each marker
plot_colors_all=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 0 0 0; 1 .5 0; .5 0 1; 0 1 .5; 0 0 0];
plot_colors=plot_colors_all(markers,:);


for i=start_frame:finish %Loop through time
    
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