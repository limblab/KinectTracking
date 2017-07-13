function affine_xform = getTRCfromMarkers(cds,marker_data,saveFolder,varargin)
% Given a cds, marker data from color tracking, and a save location, will
% spatiotemporally align the markers to the data in the CDS, smooth the
% markers, transform the coordinates for use in OpenSim, and write out a
% TRC file for the marker locations and a MOT file for the GRF from the
% handle interaction force into the save folder.


%% 4. PUT KINECT MARKER LOCATIONS IN HANDLE COORDINATES
% rotation_known=0; %Whether the rotation matrix is already known (from another file from that day)
% figure out if rotation known
if nargin>4
    error('Too many arguments')
elseif nargin==4
    affine_xform = varargin{1};
    [md,affine_xform] = realignMarkerSpacetime(cds,marker_data,affine_xform);
else
    % first file of the day, affine xform unknown
    [md,affine_xform] = realignMarkerSpacetime(cds,marker_data);
end


%% 5. SMOOTH OUT MARKERS

md = smoothMarkerData(md);

%% 5. FIND TIMES TO EXCLUDE (BECAUSE THE MONKEY THREW AWAY THE HANDLE)

%% 5a. Calculate the distances of the hand marker to the handle (and plot)
%This can be used to determine times when the monkey has thrown away the
%handle

% n_times = size(md.pos,3);
% k=reshape(md.pos(3,:,:),[3,n_times]);
% h=interp1(cds.kin.t,[cds.kin.x cds.kin.y],md.t);
% h(:,3)=0;
% 
% err=NaN(1,n_times);
% for i=1:n_times    
%     err(i)=pdist2(k(:,i)',h(i,:));
% end
% 
% figure; plot(md.t,err)
% figure; scatter3(md.pos(3,1,:),md.pos(3,2,:),md.pos(3,3,:))
% axis equal
% %This can be used in combination w/ the z-force
% 
% clear err
% clear n_times
% clear h
% clear k

%% 6. PUT KINECT DATA INTO OPENSIM COORDINATES

[md,handle_opensim] = transformForOpenSim(md,cds);

%% 7. Save TRC and GRF files

prefix=cds.meta.rawFileName;
if ~iscell(prefix)
    prefix={prefix};
end
if saveFolder(end) ~= filesep
    saveFolder = [saveFolder filesep];
end
writeToTRC(md,[saveFolder prefix{1} '_markerData.trc'])
writeHandleForceOpenSim(md,cds,handle_opensim,[saveFolder prefix{1} '_grf.mot'])

%% 9. PUT TARGET DATA INTO TRC FORMAT
% % find meta data
% frame_rate = 1/mean(diff(kinect_times));
% start_idx = find(kinect_times>=0,1,'first');
% num_frames = length(kinect_times)-start_idx+1;
% % marker_names = {'Marker_1','Marker_2','Marker_3','Marker_4','Marker_5','Marker_6','Marker_7','Marker_8','Shoulder JC','Pronation Pt1'};
% marker_names = {'Handle', 'Target'};
% num_markers = length(marker_names); % ONLY USED 10 MARKERS FOR CHIPS DATA
% 
% % open file and write header
% fid = fopen([folder prefix '_targets.trc'],'w');
% fprintf(fid,['PathFileType\t4\tX/Y/Z\t' prefix '.trc\n']);
% fprintf(fid,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
% fprintf(fid,'%5.2f\t%5.2f\t%d\t%d\tm\t%5.2f\t%d\t%d\n',[frame_rate frame_rate num_frames num_markers frame_rate 1 num_frames]);
% 
% % write out data header
% fprintf(fid,'Frame#\tTime\t');
% for i = 1:num_markers
%     fprintf(fid,'%s\t\t\t',marker_names{i});
% end
% fprintf(fid,'\n');
% fprintf(fid,'\t\t');
% for i = 1:num_markers
%     fprintf(fid,'X%d\tY%d\tZ%d\t',[i,i,i]);
% end
% fprintf(fid,'\n\n');
% 
% % write out data
% for j=1:num_frames
%     frame_idx = j-1+start_idx;
%     fprintf(fid,'%d\t%f\t',[j kinect_times(frame_idx)]);
%     
%     % print handle position
%     if isnan(handle_opensim(frame_idx,1))
%         fprintf(fid,'\t\t\t');
%     else
%         fprintf(fid,'%f\t%f\t%f\t',handle_opensim(frame_idx,:));
%     end
%     
%     %print target position
%     if isnan(target_opensim(frame_idx,1))
%         fprintf(fid,'\t\t\t');
%     else
%         fprintf(fid,'%f\t%f\t%f\t',target_opensim(frame_idx,:));
%     end
%     fprintf(fid,'\n');
% end
% 
% % close file
% fclose(fid);
% clear fid
