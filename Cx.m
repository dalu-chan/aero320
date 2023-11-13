function [Cx] = Cx(phi,type)
%CX This function gives the x principle rotation matrix for input phi (rad)
%   Detailed explanation goes here

if type=="rad"
Cx = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
elseif type=="deg"
Cx = [1 0 0; 0 cosd(phi) sind(phi); 0 -sind(phi) cosd(phi)];
end
end