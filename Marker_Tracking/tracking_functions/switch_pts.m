function [ all_medians, all_medians2 ] = switch_pts(switch_pts,t,all_medians,all_medians2)

pt1=switch_pts(1);
pt2=switch_pts(2);

temp=all_medians(pt1,:,t);
temp2=all_medians2(pt1,:,t);
all_medians(pt1,:,t)=all_medians(pt2,:,t);
all_medians2(pt1,:,t)=all_medians2(pt2,:,t);
all_medians(pt2,:,t)=temp;
all_medians2(pt2,:,t)=temp2;

