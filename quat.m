function [eta,ep1,ep2,ep3] = quat(phi,rotmat)
%QUAT This function returns the quaternion for a given phi and rotation
%matrix
eta = cos(phi/2);
ep1 = (rotmat(2,3)-rotmat(3,2))/(4*eta);
ep2 = (rotmat(3,1)-rotmat(1,3))/(4*eta);
ep3 = (rotmat(1,2)-rotmat(2,1))/(4*eta);
end