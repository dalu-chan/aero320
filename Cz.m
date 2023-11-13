function [Cz] = Cz(psi,type)
%CZ This function gives the z principle rotation matrix for input psi (rad)
%   Detailed explanation goes here
if type=="rad"
Cz = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
elseif type=="deg"
Cz = [cosd(psi) sind(psi) 0; -sind(psi) cosd(psi) 0; 0 0 1];
end
end