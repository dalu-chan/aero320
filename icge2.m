clc
clear
close

z = 20;
y = -5;
x = 10;

Cx =  [1 0 0; 0 cosd(x) sind(x); 0 -sind(x) cosd(x)];
Cy =  [cosd(y) 0 -sind(y); 0 1 0; sind(y) 0 cosd(y)];
Cz =  [cosd(z) sind(z) 0; -sind(z) cosd(z) 0; 0 0 1];

Ctfto = Cx*Cy*Cz;
disp(Ctfto)

x = 190;
Cx =  [1 0 0; 0 cosd(x) sind(x); 0 -sind(x) cosd(x)];

Ctfenu = Cx*Cy*Cz;
disp(Ctfenu)

north = [1, 0, 0; 0 -1 0; 0 0 -1]*Ctfenu*[0;1;0];
roots(@ (x) 133319*x^2+348142*x-535319)
