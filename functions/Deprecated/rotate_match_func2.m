function [ rec_pos_trans2, R, q, fval] = rotate_match_func2( rec_pos_trans, true_pos_trans )
%This function finds the rotation to best match
%the reconstructed position to the true position.

%The input are centered reconstructed positions (rec_pos_trans) and
%centered true positions (true_pos_trans). The output is the rotated
%reconstructed positions (rec_pos_trans2), the rotation matrix (R), and
%the rotation angles (q).
%This function relies on the function "rotate_func"

%Note that for our kinect/handle alignment applications, the "true
%position" is the handle locations. And the "reconstructed position" is the
% kinect position.

%First do rough grid search of possible rotations
num_angles=16;

%The angles we'll search through
dthetas=linspace(0,2*pi,num_angles+1);
dsis=linspace(0,2*pi,num_angles+1);
dphis=linspace(0,2*pi,num_angles+1);

sse_min=Inf;

temp_rec_pos_trans=rec_pos_trans;
%Test all rotation combinations
for j=1:num_angles
    for k=1:num_angles
        for l=1:num_angles
            dtheta=dthetas(j);
            dsi=dsis(k);
            dphi=dphis(l);
            sse=rotate_func([dtheta,dsi,dphi],temp_rec_pos_trans,true_pos_trans); %Get the sum-squared error or the rotation
            if sse<sse_min %If this is the best rotation, we will use these parameters as initial parameters for a more thorough search
                sse_min=sse;
                init=[dtheta, dsi, dphi];
            end
        end
    end
end

%Find the best rotation (using fmincon instead of a gridsearch),
%using the gridsearch optimal parameters as the initial parameters
[q,fval]=fmincon(@(x) rotate_func(x,temp_rec_pos_trans,true_pos_trans),init,[],[],[],[],[0 0 0],[2*pi 2*pi 2*pi]);

% break out angles
theta=q(1);
psi=q(2);
phi=q(3);

%This is the best rotation matrix in cartesian coordinates
R=[cos(theta)*cos(psi) cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi) sin(phi)*sin(psi)-cos(phi)*sin(theta)*cos(psi);...
    -cos(theta)*sin(psi) cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi) sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);...
    sin(theta) -sin(phi)*cos(theta) cos(phi)*cos(theta)];

%Do the rotation
rec_pos_trans2=rec_pos_trans*R;

end

