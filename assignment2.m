clc
clear
close

%% Question 10
syms x y z
Cx =  [1 0 0; 0 cos(x) sin(x); 0 -sin(x) cos(x)];
Cy =  [cos(y) 0 -sin(y); 0 1 0; sin(y) 0 cos(y)];
Cz =  [cos(z) sin(z) 0; -sin(z) cos(z) 0; 0 0 1];
C = Cx*Cz*Cy;
disp(C)

%% Question 13
a = [1; 0; 0]; %column unit vector where aT*a = 1
across = crossm(a); 
disp(across*across*across)
disp(-1.*across)
disp("They're the same!")
%% Problem 2
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
rotmat=rotmatrix(lvlh,eci);
disp("The rotation matrix is ")
disp(rotmat)

%% functions
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

function crossProductMatrix = crossm(vector)
    % Check if the input vector is of the correct size
    if numel(vector) ~= 3
        error('Input vector must have exactly 3 elements');
    end

    % Extract components of the vector
    x = vector(1);
    y = vector(2);
    z = vector(3);

    % Create the cross product matrix
    crossProductMatrix = [0, -z, y;
                          z, 0, -x;
                         -y, x, 0];
end