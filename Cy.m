function [Cy] = Cy(theta,type)
%CY This function gives the y principle rotation matrix for input theta (rad) 
%   Detailed explanation goes here
if type=="rad"
Cy = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
elseif type=="deg"
Cy = [cosd(theta) 0 -sind(theta); 0 1 0; sind(theta) 0 cosd(theta)];
end
end