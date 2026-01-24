clc; clear; close all

%% Constants
mu = 398600.4418;      % km^3/s^2
R_E = 6378.137;        % km
t = 66;                % number of satellites
a = 7158;              % km
inc = deg2rad(86.4);   % radians
p = 6; f = 2;
T = 2*pi*sqrt(a^3/mu); % orbital period [s]
time = linspace(0,T,10000);

%% RAAN and Mean Anomaly setup (radians)
raan = zeros(t, 1);
M = zeros(t, 1);
idx = 1;
for i = 1:p
    for j = 1:t/p
        raan(idx) = mod((2 * pi * (i - 1) / p) / 2, 2*pi);
        M(idx) = mod(2 * pi * (p / t) * (j - 1) + 2 * pi * (f / t) * (i - 1), 2*pi);
        idx = idx + 1;
    end
end

%% Propagate using Newtonian gravity (Cartesian, no perturbations)
Z_all = zeros(length(time),6,t);
options = odeset('RelTol',1e-12, 'AbsTol',1e-12);

for j = 1:t
    % Initial Cartesian state
    [r0, v0] = kep2car(a, 0, inc, raan(j), 0, M(j), mu);
    Y0 = [r0; v0];  % Initial state [x y z vx vy vz]

    % Propagate using two-body dynamics
    [~, Y] = ode45(@(t,y) twoBodyODE(t,y,mu), time, Y0, options);
    Z_all(:,:,j) = Y;
end

%% Compute closest approach between any two satellites
d_min = 100*R_E;
s1 = 0; s2 = 0; t_min = 0;

for i = 1:65
    for j = i+1:66
        pos_i = squeeze(Z_all(:,1:3,i));
        pos_j = squeeze(Z_all(:,1:3,j));
        d = vecnorm(pos_i - pos_j, 2, 2);
        [d_ij_min, t_ij] = min(d);

        if d_ij_min < d_min
            d_min = d_ij_min;
            s1 = i;
            s2 = j;
            t_min = t_ij;
        end
    end
end

fprintf('Part d\n')
fprintf('Minimum distance: %.6f km\n', d_min);
fprintf('Between satellites %d and %d at time index %d (t = %.2f sec)\n', ...
    s1, s2, t_min, time(t_min));

fprintf('Part e\n')
fprintf('It happened at X = [%.2f,%.2f,%.2f,%.2f,%.2f,%.2f]\n',Z_all(t_min,:,s1));
fprintf('and X = [%.2f,%.2f,%.2f,%.2f,%.2f,%.2f] for each satellite respectivelly\n', Z_all(t_min,:,s2));

fprintf('Sat 1 Keplerian:')
[a,e,inc,raan,argp,f] = cart2kep(Z_all(t_min,1,s1),Z_all(t_min,2,s1),Z_all(t_min,3,s1),Z_all(t_min,4,s1),Z_all(t_min,5,s1),...
    Z_all(t_min,6,s1),mu);
X1 = real([a,e,rad2deg(inc),rad2deg(raan),rad2deg(argp),rad2deg(f)]);
fprintf(' X = [%.2f,%.2f,%.2f,%.2f,%.2f,%.2f]\n',X1);

fprintf('Sat 2 Keplerian:')
[a,e,inc,raan,argp,f] = cart2kep(Z_all(t_min,1,s2),Z_all(t_min,2,s2),Z_all(t_min,3,s2),Z_all(t_min,4,s2),Z_all(t_min,5,s2),...
    Z_all(t_min,6,s2),mu);
X2 = [a,e,rad2deg(inc),rad2deg(raan),rad2deg(argp),rad2deg(f)];
fprintf(' X = [%.2f,%.2f,%.2f,%.2f,%.2f,%.2f] for each satellite respectivelly\n', X2);

%% Plot all satellite orbits
figure('Name', 'All Satellite Orbits');
hold on; axis equal; grid on
xlabel('$X$ [km]', 'Interpreter', 'latex')
ylabel('$Y$ [km]', 'Interpreter', 'latex')
zlabel('$Z$ [km]', 'Interpreter', 'latex')
title('3D Trajectories of All Satellites', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex')
colors = lines(t);

for k = 1:t
    pos = squeeze(Z_all(:,1:3,k));
    plot3(pos(:,1), pos(:,2), pos(:,3), 'Color', colors(k,:), 'LineWidth', 1.2)
end
pos1 = squeeze(Z_all(:,1:3,s1));
plot3(pos1(t_min,1), pos1(t_min,2), pos1(t_min,3), 'r*')

pos2 = squeeze(Z_all(:,1:3,s2));
plot3(pos2(t_min,1), pos2(t_min,2), pos2(t_min,3), 'r*')

plot3([pos1(t_min,1),pos2(t_min,1)], [pos1(t_min,2),pos2(t_min,2)],...
    [pos1(t_min,3),pos2(t_min,3)], 'k','LineWidth',1.5)

earthy(R_E, 'Earth', 0.9, [0; 0; 0]); hold on

function dYdt = twoBodyODE(~, Y, mu)
    r = Y(1:3);
    v = Y(4:6);
    a = -mu * r / norm(r)^3;
    dYdt = [v; a];
end

function [r_eci, v_eci] = kep2car(a,e,i,RAAN,argp,M,mu)
    % From orbital elements to Cartesian state
    E = solveKepler(M,e);  % Solve Kepler’s Equation
    nu = 2 * atan2(sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2));  % true anomaly

    r_pqw = a*(1 - e*cos(E)) * [cos(nu); sin(nu); 0];
    v_pqw = sqrt(mu*a) / (a*(1 - e*cos(E))) * [-sin(E); sqrt(1-e^2)*cos(E); 0];

    R = rotationMatrix(i, RAAN, argp);
    r_eci = R * r_pqw;
    v_eci = R * v_pqw;
end

function R = rotationMatrix(i, RAAN, argp)
    % Classical orbital rotation matrix: PQW -> ECI
    R3_W = [cos(RAAN) -sin(RAAN) 0; sin(RAAN) cos(RAAN) 0; 0 0 1];
    R1_i = [1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
    R3_w = [cos(argp) -sin(argp) 0; sin(argp) cos(argp) 0; 0 0 1];
    R = R3_W * R1_i * R3_w;
end

function E = solveKepler(M, e)
    % Solve Kepler’s Equation using Newton-Raphson
    E = M;
    for j = 1:100
        f = E - e*sin(E) - M;
        fprime = 1 - e*cos(E);
        dE = -f / fprime;
        E = E + dE;
        if abs(dE) < 1e-12
            break;
        end
    end
end
