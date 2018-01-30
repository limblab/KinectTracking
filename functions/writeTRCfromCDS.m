%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function writeToTRC(cds,savefile_name)
%
%   This function writes out marker data into a TRC format to be read by
%   OpenSim
%
% INPUTS:
%   cds  : cds object containing marker data, in OpenSim coordinates
%   savefile_name : name of TRC to write out (including folder)
%
% OUTPUTS:
%   none
%
%
% Written by Raeed Chowdhury
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeTRCfromCDS(cds,savefile_name)

%% Figure out if we have marker data and load it from cds analog
marker_analog_idx = 0;
for i=1:length(cds.analog)
    header = cds.analog{i}.Properties.VariableNames;
    if any(contains(header,'Frame')) || any(contains(header,'Marker'))
        marker_analog_idx = i;
        break
    end
end
if marker_analog_idx > 0
    marker_data=cds.analog{marker_analog_idx};
else
    error('No marker information in CDS')
end


%% throw away all frames that have 1 or fewer hand markers
missing_wristhand = sum(ismissing(marker_data(:,3:7)),2) > 3;
missing_hand = sum(ismissing(marker_data(:,3:5)),2) > 1;
good_frame = ~(missing_wristhand | missing_hand);
% for i=1:length(good_frame)
%     if sum(isnan(marker_data{i,3:7}))>9
%         good_frame(i) = false;
%     end
%     
%     if sum(isnan(marker_data{i,3:5}))>3
%         good_frame(i) = false;
%     end
% end

%% 8. PUT KINECT MOTION TRACKING DATA INTO TRC FORMAT
% find meta data
frame_rate = 1/mean(diff(marker_data.t));
num_markers = width(marker_data)-2;
num_frames = height(marker_data);
marker_names = marker_data.Properties.VariableNames(3:end);

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
for frame_idx=1:num_frames
    if(good_frame(frame_idx))
        fprintf(fid,'%d\t',marker_data.Frame(frame_idx));
        
        for i = 2:width(marker_data)
            if isnan(marker_data{frame_idx,i})
                fprintf(fid,'\t\t\t');
            else
                fprintf(fid,'%f\t',marker_data{frame_idx,i});
            end
        end
        
        fprintf(fid,'\n');
    end
end

% close file
fclose(fid);
clear fid