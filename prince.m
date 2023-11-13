function [Cx,Cy,Cz] = prince(x,y,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Cx =  [1 0 0; 0 cosd(x) sind(x); 0 -sind(x) cosd(x)];
Cy =  [cosd(y) 0 -sind(y); 0 1 0; sind(y) 0 cosd(y)];
Cz =  [cosd(z) sind(z) 0; -sind(z) cosd(z) 0; 0 0 1];
end