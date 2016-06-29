medians=kinect_pos_smooth;

xlims=[min(min(medians(:,1,:))) max(max(medians(:,1,:)))];
ylims=[min(min(medians(:,2,:))) max(max(medians(:,2,:)))];
zlims=[min(min(medians(:,3,:))) max(max(medians(:,3,:)))];

figure;
set(gca,'NextPlot','replacechildren');
xlim(xlims)
ylim(ylims)
zlim(zlims)

plot_colors=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 0 0 0; 1 .5 0; .5 0 1; 0 1 .5; 0 0 0];

t = kinect_times;
handle_pos = interp1(bdf.pos(:,1),bdf.pos(:,2:3),t);

% vid=[];    
for i=3500:4000
    scatter3(medians(:,1,i),medians(:,2,i),medians(:,3,i),20,plot_colors,'filled')
    
%     pause;
    title(t(i))
    
    hold on
    start_ind = max(1,i-60);
    plot3(handle_pos(start_ind:i,1),handle_pos(start_ind:i,2),zeros(i-start_ind+1,1),'-b')
    hold off
    
%     view([58.1250 36.9294])
    view([0 90])
    xlim(xlims)
    ylim(ylims)
    zlim(zlims)

    pause(.03)
    
end

%%
medians=kinect_pos;

xlims=[min(min(medians(:,1,:))) max(max(medians(:,1,:)))];
ylims=[min(min(medians(:,2,:))) max(max(medians(:,2,:)))];

figure;
set(gca,'NextPlot','replacechildren');
xlim(xlims)
ylim(ylims)

plot_colors=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 0 0 0; 1 .5 0; .5 0 1; 0 1 .5; 0 0 0];

for i=1:2000
    scatter(medians(:,1,i),medians(:,2,i),200,plot_colors,'filled')
%     view([0.8, 0.1, -0.1 ])
    pause(.03)
%     pause;
    title(i)

    
end

%%
medians=kinect_pos;

% samp_rate = 1/mean(diff(kinect_times));
% [b,a] = butter(4,10/(samp_rate/2));
% for i = 1:size(medians,1)
%     median_i = squeeze(medians(i,:,:))';
%     medians(i,:,:) = filtfilt(b,a,median_i)';
% end

xlims=[min(min(medians(:,1,:))) max(max(medians(:,1,:)))];
ylims=[min(min(medians(:,2,:))) max(max(medians(:,2,:)))];
zlims=[min(min(medians(:,3,:))) max(max(medians(:,3,:)))];

fhandle = figure('position',[100,50,1000,800]);
subplot(221)
set(gca,'NextPlot','replacechildren');
xlim(xlims)
ylim(ylims)
zlim(zlims)

subplot(222)
set(gca,'NextPlot','replacechildren');
xlim(xlims)
ylim(ylims)
title 'Top down'

subplot(223)
set(gca,'NextPlot','replacechildren');
xlim(xlims)
ylim(ylims)
title 'Behind'

subplot(224)
set(gca,'NextPlot','replacechildren');
xlim(xlims)
ylim(ylims)
title 'Side'

plot_colors=[1 0 0; 0 1 0; 0 0 1; 1 1 0; 0 1 1; 1 0 1; 0 0 0; 1 .5 0; .5 0 1; 0 1 .5; 0 0 0];

t = kinect_times;
handle_pos = interp1(bdf.pos(:,1),bdf.pos(:,2:3),t);

% vid=[];    
for i=100:2000
    subplot(221)
    scatter3(medians(:,1,i),medians(:,2,i),medians(:,3,i),200,plot_colors,'filled')
    
    hold on
    start_ind = max(1,i-60);
    plot3(handle_pos(start_ind:i,1),handle_pos(start_ind:i,2),zeros(i-start_ind+1,1),'-b')
    hold off
    
    view([58.1250 36.9294])
    xlim(xlims)
    ylim(ylims)
    zlim(zlims)
    title(t(i))
    
    subplot(222)
    scatter(medians(:,1,i),medians(:,2,i),200,plot_colors,'filled')
    hold on
    plot(handle_pos(start_ind:i,1),handle_pos(start_ind:i,2),'-b')
    hold off
    xlim(xlims)
    ylim(ylims)
    
%     subplot(223)
%     scatter(medians(:,2,i),medians(:,3,i),200,plot_colors,'filled')
%     hold on
%     plot(handle_pos(start_ind:i,2),zeros(i-start_ind+1,1),'-b')
%     hold off
%     xlim(ylims)
%     ylim(zlims)
%     
%     subplot(224)
%     scatter(medians(:,1,i),medians(:,3,i),200,plot_colors,'filled')
%     hold on
%     plot(handle_pos(start_ind:i,1),zeros(i-start_ind+1,1),'-b')
%     hold off
%     xlim(xlims)
%     ylim(zlims)
    
    pause(.001)
end