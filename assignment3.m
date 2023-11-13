clc
clear
close all

%% Question 10
syms x y z dx dy dz
W21 = [dx;0;0]+Cx(x)*[0;0;dz]+Cx(x)*Cz(z)*[0;dy;0];

%% Question 12
syms alpha thetai thetao
% Part a
Coi = Cy(alpha);
disp("Rotation matrix Coi is")
disp(Coi)
% Part b
C2i = Cx(thetai);
disp("Rotation matrix C2i is")
disp(C2i)
Co3 = Cx(thetao);
C3o = transpose(Co3);
disp("Rotation matrix C3o is")
disp(C3o)
% Part c
Cc3 = transpose(Coi);
disp("Rotation matrix Cc3 is")
disp(Cc3)
C3c = Coi;
% Part d
Zco = Co3*C3c*[0;0;1];
disp("The form given by sequence one is")
disp(Zco)
Ci2 = transpose(C2i);
C2c = eye(3);
Zco = Coi*Ci2*C2c*[0;0;1];
disp("The form given by sequence two is")
disp(Zco)
%% Problem 5

r = [6738;3391;1953];
v = [-3.5;4.39;4.44];

zlvlh = -r/norm(r);
ylvlh = -cross(r,v)/norm(cross(r,v));
xlvlh = cross(ylvlh,zlvlh);
lvlh = [xlvlh,ylvlh,zlvlh];

xeci = [1;0;0];
yeci = [0;1;0];
zeci = [0;0;1];
eci = [xeci,yeci,zeci];

% Rotation matrix
rotmat=rotmatrix(lvlh,eci);
disp("The rotation matrix is ")
disp(rotmat)

% Principle angle and rotation
phi = acos((trace(rotmat)-1)/2)
a1 = (rotmat(2,3)-rotmat(3,2))/(2*sin(phi));
a2 = (rotmat(3,1)-rotmat(1,3))/(2*sin(phi));
a3 = (rotmat(1,2)-rotmat(2,1))/(2*sin(phi));
a = [a1 a2 a3]

% Quaternion
eta = cos(phi/2)
ep1 = (rotmat(2,3)-rotmat(3,2))/(4*eta)
ep2 = (rotmat(3,1)-rotmat(1,3))/(4*eta)
ep3 = (rotmat(1,2)-rotmat(2,1))/(4*eta)

% Principle Rotations
phi = atan(rotmat(2,3)/rotmat(3,3))
theta = -1*asin(rotmat(1,3))
psi = atan(rotmat(1,2)/rotmat(1,1))

%% Functions used
function [matrix] = rotmatrix(basis2,basis1)
%rotmatrix This function finds the rotation matrix between two input basis'
%   Detailed explanation goes here

matrix(1,1) = dot(basis1(:,1),basis2(:,1));
matrix(1,2) = dot(basis1(:,2),basis2(:,1));
matrix(1,3) = dot(basis1(:,3),basis2(:,1));
matrix(2,1) = dot(basis1(:,2),basis2(:,1));
matrix(2,2) = dot(basis1(:,2),basis2(:,2));
matrix(2,3) = dot(basis1(:,2),basis2(:,3));
matrix(3,1) = dot(basis1(:,3),basis2(:,1));
matrix(3,2) = dot(basis1(:,3),basis2(:,2));
matrix(3,3) = dot(basis1(:,3),basis2(:,3));

end
function [Cx] = Cx(phi)
%CX This function gives the x principle rotation matrix for input phi (rad)
%   Detailed explanation goes here
Cx = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
end
function [Cy] = Cy(theta)
%CY This function gives the y principle rotation matrix for input theta (rad) 
%   Detailed explanation goes here
Cy = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
end
function [Cz] = Cz(psi)
%CZ This function gives the z principle rotation matrix for input psi (rad)
%   Detailed explanation goes here
Cz = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];
end