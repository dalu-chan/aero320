function [phi,a] = panglerot(rotmat)
%PANGLEROT This function returns the principle angle and axis of rotation
%for a given rotation matrix
%   Detailed explanation goes here
phi = acos((trace(rotmat)-1)/2);
a1 = (rotmat(2,3)-rotmat(3,2))/(2*sin(phi));
a2 = (rotmat(3,1)-rotmat(1,3))/(2*sin(phi));
a3 = (rotmat(1,2)-rotmat(2,1))/(2*sin(phi));
a = [a1 a2 a3];
end