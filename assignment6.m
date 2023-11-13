%% Problem 3
clc
clear
close all

omega = 10000; %rpm about +z
omega = 2*pi*omega/60; %omega converted to rad/s 
omega = [0;0;omega];
m = 4; %mass of disk kg
r = 7/100; %radius of disk, cm to m
th = 1/100; %thickness of disk, cm to m
dist = 2/100;
dist = [0;0;dist];

dgamma = 2; %rad/s about +x
dgamma = [dgamma;0;0];

Iz = .5*m*r^2; %kg*m^2
Iy = (.25*m*r^2)+(1/12*m*th^2);
Ix = Iy;
I = diag([Ix, Iy, Iz]);
J = I - m*crossm(dist)*crossm(dist);

tau = crossm(dgamma)*J*omega;
force = tau/2/(2/100);
disp(force)