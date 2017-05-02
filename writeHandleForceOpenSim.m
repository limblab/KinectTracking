function writeHandleForceOpenSim(marker_data,cds,handle_opensim,savefile_name)

% Convert forces to OpenSim coordinates
% Robot coordinates:
% Origin at shoulder joint center of robot, x is to right, y is towards screen, z is up
% OpenSim coordinates:
% Origin at shoulder joint center (marker 9), x is towards screen, y is up, z is to right

% First take care of forces
force = -cds.force{:,[3 4 2 6 7 5]}; % permute forces to match OpenSim coordinates and flips sign for force applied to arm
% force = interp1(cds.force.t,force,marker_data.t);

% find meta data
version = 1;
nRows = size(force,1);
nColumns = 10; % 6 DOF for load cell and 1 for time plus 3 for contact point
header_names = {'time','fx','fy','fz','mx','my','mz','px','py','pz'};

% open file and write header
[~,prefix,ext] = fileparts(savefile_name);
if ~strcmp(ext,'.mot')
    warning('Not writing with .mot extension. You might want to rename that...')
end
fid = fopen(savefile_name,'w');
fprintf(fid,[prefix '.mot\n']);
fprintf(fid,'version=%d\n',version);
fprintf(fid,'nRows=%d\n',nRows);
fprintf(fid,'nColumns=%d\n',nColumns);
fprintf(fid,'inDegrees=yes\n');
fprintf(fid,'endheader\n');

% write out data header
for i = 1:nColumns
    fprintf(fid,'%s\t',header_names{i});
end
fprintf(fid,'\n');

% write out data
for j=1:nRows
    
%     % skip if shoulder marker didn't exist
%     if ~isnan(handle_opensim(j,1))
%         % print times
%         fprintf(fid,'%f\t',marker_data.t(j));
%         
%         % print forces
%         for i = 2:7
%             fprintf(fid,'%f\t',force(j,i-1));
%         end
%         
%         % print application points
%         fprintf(fid,'%f\t%f\t%f\t',handle_opensim(j,:));
%         
%         fprintf(fid,'\n');
%     end

        fprintf(fid,'%f\t',cds.force.t(j));
        
        % print forces
        for i = 2:7
            fprintf(fid,'%f\t',force(j,i-1));
        end
        
        % print application points in hand coordinates for OpenSim
        fprintf(fid,'%f\t%f\t%f\t',[0 -0.02 0]);
        
        fprintf(fid,'\n');
end

% close file
fclose(fid);
clear fid