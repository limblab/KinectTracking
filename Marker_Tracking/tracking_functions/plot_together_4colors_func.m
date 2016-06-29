function [] = plot_together_4colors_func(frame,text_markers,markers,all_medians, color1, color2, color3, color4, start, finish, plot_original)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Align

all_medians_aligned=NaN(11,3,finish);
all_medians_aligned(:,:,start:finish)=all_medians;

%% Plot

medians=all_medians_aligned; %Renamed (to make compatible with my previous plotting script)

%Set plotting limits as minimums and maximums
xlims=[min(min(medians(markers,1,:))) max(max(medians(markers,1,:)))];
ylims=[min(min(medians(markers,2,:))) max(max(medians(markers,2,:)))];
zlims=[min(min(medians(markers,3,:))) max(max(medians(markers,3,:)))];

xlims=[-.4 .5];
ylims=[-.4 .4];
zlims=[.8 1.4];

%Initialize figure
figure;
set(gca,'NextPlot','replacechildren');
xlim(zlims)
ylim(xlims)
zlim(ylims)

%Colors of each marker
plot_colors_all=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 0 0 0; 1 .5 0; .5 0 1; 0 1 .5; 0 0 0];
plot_colors=plot_colors_all(markers,:);


i=frame;

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

%Label markers
text_str=cell(1,length(text_markers));
for j=1:length(text_markers)
    text_str{j}=num2str(text_markers(j));
end
text(medians(text_markers,3,i),medians(text_markers,1,i),medians(text_markers,2,i),text_str,'FontSize',30)

% text(medians(text_markers,3,i),medians(text_markers,1,i),medians(text_markers,2,i),{'1','2','3','4','5','6','7','8','9','10'},'FontSize',20)





title(i); %Give time as title
xlabel('z')
ylabel('x')
zlabel('y')

%Set limits
xlim(zlims)
ylim(xlims)
zlim(ylims)


end


%%


