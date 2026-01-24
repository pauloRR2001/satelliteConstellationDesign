clc
clear
close all

mu = 398600.4415;
e = 0.65;
R_E = 6378.1363;

a = 3*R_E;
p = a*(1-e^2);

figure
plotOrbit3(0, 0, 0, p, e, linspace(0,2*pi,1000), 'k', 1, 1, [0,0,0],0,1)
hold on
plotOrbit3(0, 0, 0, p, e, deg2rad(-55), 'r*', 1, 1, [0,0,0],0,1)
plot([0 0],[-p p], 'k--')
plot([-a*(1+e),a*(1-e)], [0,0],'k--')
%xline(0,'k')
%yline(0,'k')
view(0,90)
axis equal
xlim([-8*R_E,3*R_E])

%%

mu = 398600.4415;
e = 1.1;
R_E = 6378.1363;

a = -2*R_E;
p = a*(1-e^2);

figure
plotOrbit3(0, 0, 0, p, e, linspace(deg2rad(-95),deg2rad(95),1000), 'k', 1, 1, [0,0,0],0,1)
hold on
plotOrbit3(0, 0, 0, p, e, deg2rad(-55), 'r*', 1, 1, [0,0,0],0,1)
plot([0 0],[-p p], 'k--')
plot([-p,p], [0,0],'k--')
%xline(0,'k')
%yline(0,'k')
view(0,90)
axis equal
%xlim([-8*R_E,3*R_E])