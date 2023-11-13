function [phi,theta,psi] = princi(rotmat)
%PRINCI This function returns yaw, pitch, roll in radians
%   Detailed explanation goes here
phi = atan(rotmat(2,3)/rotmat(3,3));
theta = -1*asin(rotmat(1,3));
psi = atan(rotmat(1,2)/rotmat(1,1));
end