function [sse]=rotate_func(x,rec_pos_trans,true_pos_trans)

%This function rotates the reconstructed position and gives the sum squared error
%between this rotation and the true position. 

%The inputs are the vector of angles (x), centered reconstructed positions (rec_pos_trans), and
%centered true positions (true_pos_trans). The output is the sum squared
%error (sse)

%Note that for our kinect/handle alignment applications, the "true
%position" is the handle locations. And the "reconstructed position" is the
% kinect position.

theta=x(1);
psi=x(2);
phi=x(3);

%This is the rotation matrix in cartesian coordinates
R=[cos(theta)*cos(psi) cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi);...
        -cos(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);...
        sin(theta) -sin(phi)*cos(theta) cos(phi)*cos(theta)];

%Get the sum squared error
diff=rec_pos_trans*R-true_pos_trans;    
sse=norm(diff);  