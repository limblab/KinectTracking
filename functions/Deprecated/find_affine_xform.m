function affine_xform = find_affine_xform(pos_k,pos_h)

% Affine transformation is as follows:
% p_prime = A*p
% where p and p_prime are in homogenous coordinates
% [p1; p2; p3; 1]
% and A is a 3D (4x4)
% affine transformation matrix
% [a11 a12 a13 a14;
%  a21 a22 a23 a24;
%  a31 a32 a33 a34;
%  0   0   0   1]
%
% Rearrange equation to solve for elements of A in least squares:
% P = [p1 p2 p3 1 0  0  0  0 0  0  0  0;
%      0  0  0  0 p1 p2 p3 1 0  0  0  0;
%      0  0  0  0 0  0  0  0 p1 p2 p3 1]
%
% a = [a11; a12; a13; a14; a21; a22; a23; a24; a31; a32; a33; a34]
% 
% such that
% p_hat = [p1_prime; p2_prime; p3_prime] = P*a
%
% Stack up all points and solve for a with simple least squares, i.e.
% a = P \ p_hat;
%
% In this case, kinect pos (pos_k) is p and handle pos (pos_h) is p_prime

% Convert kinect pos to matrix form
P_kinect = kron( eye(3), [pos_k ones(size(pos_k,1),1)] );

% Stack handle pos in the same way
% x-coordinate first, then y, then z (to match kron output)
p_hat_handle = pos_h(:);

% solve for affine vector
a_vect = P_kinect \ p_hat_handle;

% reshape into affine matrix
affine_xform = reshape(a_vect,4,3)';

% tack on last row
affine_xform(4,:) = [0 0 0 1];