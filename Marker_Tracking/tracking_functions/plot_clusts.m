function [ ] = plot_clusts( plot_on, num_clust, x, y, z, medians, i, t, pause_time, xlims, ylims, zlims )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if num_clust==3    
    plot_colors=[1 0 0; 0 1 0; 0 0 1];
end
if num_clust==2
    plot_colors=[1 0 0; 0 1 0];
end
if num_clust==1
    plot_colors=[1 0 0];
end

    if plot_on
    
    scatter3(x,y,z)
    hold on;
    scatter3(medians(:,1,t),medians(:,2,t),medians(:,3,t),200,plot_colors,'filled')
    hold off;
    xlim(xlims)
    ylim(ylims)
    zlim(zlims)
    title(i)
    pause(pause_time);
    end




end

