%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function writeToTRC(marker_data,savefile_name)
%
%   This function writes out marker data into a TRC format to be read by
%   OpenSim
%
% INPUTS:
%   marker_data  : struct containing marker data, in OppenSim coordinates
%       Fields:
%           pos - position of markers, in OpenSim coordinates
%           t - times for video frames of marker capture, aligned to
%           cerebus collection
%   savefile_name : name of TRC to write out (including folder)
%
% OUTPUTS:
%   md_opensim : struct of marker_data in OpenSim coordinates
%       Fields:
%           pos - position of markers, aligned to OpenSim coordinates
%           t - times for video frames of marker capture, aligned to
%           cerebus collection
%
%
% Written by Raeed Chowdhury and Joshua Glaser. Updated April 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeToTRC(marker_data,savefile_name)

% throw away all frames that have 1 or fewer hand markers
good_frame = true(size(marker_data.pos,3),1);
for i=1:length(good_frame)
    if sum(isnan(marker_data.pos(1:5,1,i)))>3
        good_frame(i) = false;
    end
    
    if sum(isnan(marker_data.pos(1:3,1,i)))>1
        good_frame(i) = false;
    end
end

%% 8. PUT KINECT MOTION TRACKING DATA INTO TRC FORMAT
% find meta data
frame_rate = 1/mean(diff(marker_data.t));
num_markers = 10; % ONLY USED 10 MARKERS FOR CHIPS DATA
start_idx = find(marker_data.t>=0,1,'first');
num_frames = length(marker_data.t)-start_idx+1;
marker_names = {'Marker_1','Marker_2','Marker_3','Marker_4','Marker_5','Marker_6','Marker_7','Marker_8','Shoulder JC','Pronation Pt1'};

% open file and write header
[~,prefix,ext] = fileparts(savefile_name);
if ~strcmp(ext,'.trc')
    warning('Not writing with .trc extension. You might want to rename that...')
end
fid = fopen(savefile_name,'w');
fprintf(fid,['PathFileType\t4\tX/Y/Z\t' prefix '.trc\n']);
fprintf(fid,'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid,'%5.2f\t%5.2f\t%d\t%d\tm\t%5.2f\t%d\t%d\n',[frame_rate frame_rate num_frames num_markers frame_rate 1 num_frames]);

% write out data header
fprintf(fid,'Frame#\tTime\t');
for i = 1:num_markers
    fprintf(fid,'%s\t\t\t',marker_names{i});
end
fprintf(fid,'\n');
fprintf(fid,'\t\t');
for i = 1:num_markers
    fprintf(fid,'X%d\tY%d\tZ%d\t',[i,i,i]);
end
fprintf(fid,'\n\n');

% write out data
for j=1:num_frames
    frame_idx = j-1+start_idx;
    if(good_frame(frame_idx))
        fprintf(fid,'%d\t%f\t',[j marker_data.t(frame_idx)]);
        marker_pos = marker_data.pos(1:num_markers,:,frame_idx);
        for i = 1:num_markers
            if isnan(marker_pos(i,1))
                fprintf(fid,'\t\t\t');
            else
                fprintf(fid,'%f\t%f\t%f\t',marker_pos(i,:));
            end
        end
        fprintf(fid,'\n');
    end
end

% close file
fclose(fid);
clear fid