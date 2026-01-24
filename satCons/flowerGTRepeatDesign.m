%% Part (a): Compute Semi-Major Axis from Ground Track Repeat Condition

clc; clear; close all

% Author: Paulo Ramirez
% Institution: Purdue University
% Course: AAE 590 – Satellite Constellations

% === Constants ===
mu = 398600.4418;             % Earth's gravitational parameter [km^3/s^2]
TE = 86164;                   % Sidereal day [s]
Np = 17;                      % Earth rotations in repeat cycle
Nd = 10;                      % Satellite revolutions in repeat cycle
omega_E = 7.2921151467e-5;    % Earth rotation rate [rad/s]
R_E = 6378.137;               % Earth radius [km]

% === Compute Orbital Period from Ground Track Condition ===
T_orbit = Nd * 2*pi / omega_E;  % orbital period based on ground track repeat [s]

% === Compute Semi-Major Axis using modified Kepler’s Third Law ===
a = (mu / omega_E^2 * (Nd / Np)^2)^(1/3);  % [km]

fprintf('[a] Semi-major axis for 2D Flower Ground Track Repeat: %.3f km\n', a);

%% Part (b): RAAN–Mean Anomaly Distribution (2D Lattice Flower)

% === Constellation Definition ===
N_o = 3;        % Number of orbital planes
N_so = 9;       % Satellites per orbit
N_c = 2;        % Petal count
N = N_o * N_so; % Total number of satellites

raan = zeros(N,1);     % RAAN [rad]
M = zeros(N,1);        % Mean Anomaly [rad]

% === Satellite Orbital Element Storage Structure ===
sats_kep(N).a = [];  % Preallocate struct array

idx = 1;
for i = 1:N_o
    for j = 1:N_so
        raan_val = 2*pi * (i - 1)/N_o;
        M_val = 2*pi * ((j - 1)/N_so - (N_c/N_so)*(i - 1)/N_o);
        M_val = mod(M_val, 2*pi);  % wrap to [0, 2π]

        raan(idx) = raan_val;
        M(idx) = M_val;

        % Store orbital elements (a, e, i, RAAN, ω, M)
        sats_kep(idx).a = a;
        sats_kep(idx).e = 0;
        sats_kep(idx).i = deg2rad(56);   % Inclination [rad]
        sats_kep(idx).RAAN = raan_val;
        sats_kep(idx).argp = 0;          % Argument of perigee
        sats_kep(idx).M = M_val;

        idx = idx + 1;
    end
end

% === Plot RAAN vs Mean Anomaly Distribution ===
figure
scatter(rad2deg(raan), rad2deg(M), 50, 'filled', 'b')
xlabel('RAAN $\Omega$ [deg]', 'Interpreter', 'latex')
ylabel('Mean Anomaly $M$ [deg]', 'Interpreter', 'latex')
title('2D Lattice Flower Constellation: RAAN vs Mean Anomaly', 'Interpreter', 'latex')
grid on; axis equal
xlim([0 360]); ylim([0 360])
xticks(0:60:360); yticks(0:60:360)
set(gca, 'TickLabelInterpreter', 'latex')

%% Part (c): Inertial 3D Representation of the Constellation at t = 0

% === Preallocate ECI position matrix ===
Y = zeros(N, 3);  % Satellite positions [x, y, z] in km

for k = 1:N
    % Extract stored orbital elements
    a_k     = sats_kep(k).a;
    e_k     = sats_kep(k).e;
    i_k     = sats_kep(k).i;
    RAAN_k  = sats_kep(k).RAAN;
    argp_k  = sats_kep(k).argp;
    M_k     = sats_kep(k).M;

    % Compute ECI position at t = 0 (true anomaly = M since e = 0)
    [x, y, z, ~, ~, ~] = kep2cart(a_k, e_k, i_k, RAAN_k, argp_k, M_k, mu);
    Y(k, :) = [x, y, z];
end

% === Plot: Polar View (Top-Down Z-axis) ===

figure('Name', 'ECI Polar View')
earthy(R_E, 'Earth', 0.9, [0; 0; 0]); hold on
scatter3(Y(:,1), Y(:,2), Y(:,3), 50, 'filled', 'b')
xlabel('$X$ [km]', 'Interpreter', 'latex')
ylabel('$Y$ [km]', 'Interpreter', 'latex')
zlabel('$Z$ [km]', 'Interpreter', 'latex')
title('Initial Galileo Constellation in ECI (Polar View)', 'Interpreter', 'latex')
axis equal; grid on
view(0, 90)
set(gca, 'TickLabelInterpreter', 'latex')

% Optional: Plot orbit paths
for j = 1:N
    plotOrbit3(sats_kep(j).RAAN, sats_kep(j).i, sats_kep(j).argp, ...
               sats_kep(j).a, sats_kep(j).e, linspace(0, 2*pi, 100), ...
               '', 1, 1, [0,0,0], 0, 1);
end

% === Plot: Isometric View ===
figure('Name', 'ECI Isometric View')
earthy(R_E, 'Earth', 0.9, [0; 0; 0]); hold on
scatter3(Y(:,1), Y(:,2), Y(:,3), 50, 'filled', 'b')
xlabel('$X$ [km]', 'Interpreter', 'latex')
ylabel('$Y$ [km]', 'Interpreter', 'latex')
zlabel('$Z$ [km]', 'Interpreter', 'latex')
title('Initial Galileo Constellation in ECI (Isometric View)', 'Interpreter', 'latex')
axis equal; grid on
view(45, 30)
set(gca, 'TickLabelInterpreter', 'latex')

% Optional: Plot orbit paths
for j = 1:N
    plotOrbit3(sats_kep(j).RAAN, sats_kep(j).i, sats_kep(j).argp, ...
               sats_kep(j).a, sats_kep(j).e, linspace(0, 2*pi, 100), ...
               '', 1, 1, [0,0,0], 0, 1);
end

%% Part (d): Minimum Time for Full Relative Motion in the Constellation

% === Time to Repeat Full Configuration ===
N_repeat = lcm(N_o, N_so);          % LCM of number of planes and sats/orbit

% Nodal period (sidereal day)
T_nodal = 2*pi / omega_E;           % [s] Earth's sidereal rotation period

% Total time to repeat constellation geometry
T_full = Nd * T_nodal;              % [s] same as 10 sidereal days

fprintf('[d] Time to repeat full relative geometry: %.5f days (%.5f hours)\n', ...
        T_full / 86400, T_full / 3600);

%% Part (e): Visibility from ARMS Building Using Stored Satellite States

% === Observer Location (ARMS Building) ===
lat_obs = deg2rad(40.43094);
lon_obs = deg2rad(-86.915798);
alt_obs = 0;                % [km]
el_mask = 15;               % Elevation mask [deg]

% === Time Vector for Simulation ===
dt = 60;                    % [s] time step
t_vec = 0:dt:T_full;        % full cycle
nt = length(t_vec);

sats_ECI = zeros(length(t_vec),3,N);

% === Ground Station ECEF at t = 0 ===
xc0 = R_E * cos(lat_obs) * cos(lon_obs);
yc0 = R_E * cos(lat_obs) * sin(lon_obs);
zc  = R_E * sin(lat_obs);

% === Rotation Matrix to ENU Frame at Observer ===
g = [1, 0, 0;
     0, cos(lat_obs), sin(lat_obs);
     0, -sin(lat_obs), cos(lat_obs)] * ...
    [sin(lon_obs), -cos(lon_obs), 0;
     cos(lon_obs),  sin(lon_obs), 0;
     0,             0,            1];

% === Preallocate Visibility Matrix: [sat x time] ===
visibility = zeros(N, nt);

% === Loop Over All Satellites and Times ===
for k = 1:N
    % Unpack satellite orbital elements
    a_k     = sats_kep(k).a;
    e_k     = sats_kep(k).e;
    i_k     = sats_kep(k).i;
    RAAN_k  = sats_kep(k).RAAN;
    argp_k  = sats_kep(k).argp;
    M0_k    = sats_kep(k).M;

    for i = 1:nt
        t = t_vec(i);

        % === Propagate Mean Anomaly ===
        Mk = mod(M0_k + sqrt(mu / a_k^3) * t, 2*pi);

        % === ECI Satellite Position ===
        [x, y, z, ~, ~, ~] = kep2cart(a_k, e_k, i_k, RAAN_k, argp_k, Mk, mu);

        sats_ECI(i, :, k) = [x, y, z];

        % === Rotate to ECEF ===
        theta = omega_E * t;
        xrel =  x * cos(theta) + y * sin(theta);
        yrel = -x * sin(theta) + y * cos(theta);
        zrel =  z;

        % === Topocentric Vector in ECEF ===
        xrr = xrel - xc0;
        yrr = yrel - yc0;
        zrr = zrel - zc;

        % === ENU Frame Conversion ===
        rv = g * [xrr; yrr; zrr];

        % === Elevation Calculation ===
        if rv(2) > 0
            elevation = asin(rv(2) / norm(rv)) * 180/pi;
            if elevation > el_mask
                visibility(k, i) = 1;
            end
        end
    end
end

% === Compute Total Number of Visible Satellites at Each Time ===
visible_satellites = sum(visibility, 1);  % [1 x nt]

% === Plot Result ===
figure
scatter(t_vec / 86400, visible_satellites,'.', 'LineWidth',0.5)
xlabel('Time [days]', 'Interpreter', 'latex')
ylabel('Number of Visible Satellites', 'Interpreter', 'latex')
title('Galileo Constellation Visibility from ARMS Building', 'Interpreter', 'latex')
grid on
set(gca, 'TickLabelInterpreter', 'latex')

%% Part (f): GDOP Over Time from ARMS Building

inc = deg2rad(56); % Galileo inclination [deg]

% Elevation mask
el_mask = deg2rad(15);

% Time setup
dt = 60;  % seconds
t_vec = 0:dt:T_full;
gdop_vec = NaN(size(t_vec));

% Loop over time
for ti = 1:length(t_vec)
    t = t_vec(ti);
    G = [];  % Geometry matrix for visible satellites
    % Inside the time loop:
theta = omega_E * t;
g = [1, 0, 0; ...
     0, cos(lat_obs), sin(lat_obs); ...
     0, -sin(lat_obs), cos(lat_obs)] * ...
    [sin(lon_obs + theta), -cos(lon_obs + theta), 0; ...
     cos(lon_obs + theta),  sin(lon_obs + theta), 0; ...
     0,                    0,                     1];
    % Observer ECEF position (fixed on Earth's surface)
    theta = omega_E * t;  % Earth rotation angle
    x_obs = R_E * cos(lat_obs) * cos(lon_obs + theta);
    y_obs = R_E * cos(lat_obs) * sin(lon_obs + theta);
    z_obs = R_E * sin(lat_obs);
    r_obs = [x_obs; y_obs; z_obs];
    
    for k = 1:N
        % Propagate satellite M (true anomaly = M since e = 0)
        Mk = mod(M(k) + sqrt(mu / a^3) * t, 2*pi);
        [x, y, z, ~, ~, ~] = kep2cart(a, 0, inc, raan(k), 0, Mk, mu);
        r_sat = [x; y; z];
        
        % Line-of-sight vector
        rho = r_sat - r_obs;
        los_enu = g * (rho / norm(rho));  % ECEF→ENU and normalize
        
        % Elevation calculation
        enu = g * rho;  % Using same ENU matrix from part (e)
        el = asin(enu(2) / norm(rho));
        
        if el >= el_mask
            G = [G; los_enu', 1];  % Append LOS vector + clock row
        end
    end
    
    % Compute GDOP if at least 4 satellites visible
    if size(G,1) >= 4
        Q = inv(G' * G);
        gdop_vec(ti) = sqrt(trace(Q));
    end
end

% Plot GDOP
figure
scatter(t_vec / 86400, gdop_vec, 5, 'filled')
xlabel('Time [days]', 'Interpreter', 'latex')
ylabel('Geometric Dilution of Precision (GDOP)', 'Interpreter', 'latex')
title('Geometric Dilution of Precision from ARMS Building Over Time', 'Interpreter', 'latex')
grid on
set(gca, 'TickLabelInterpreter', 'latex')

