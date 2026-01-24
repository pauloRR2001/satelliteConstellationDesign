clc; clear; close all

%% === Constants ===
mu = 398600.4418;            % km^3/s^2
TE = 86164.0905;             % sidereal day [s]
Np = 14;                     % Earth rotations
Nd = 1;                      % Satellite revolutions
omega_E = 7.2921151467e-5;   % rad/s
R_E = 6378.137;              % km

% === Observer (ARMS building) ===
lat_obs = deg2rad(40.43094);
lon_obs = deg2rad(-86.915798);
alt_obs = 0;                         % km

el_mask = deg2rad(15);              % radians

% === Orbit Setup ===
T_orbit = Nd * TE / Np;                          % [s]
a = ((T_orbit^2 * mu) / (4 * pi^2))^(1/3);       % km
n = sqrt(mu / a^3);                              % [rad/s]
T_sim = TE;                                      % simulate 1 sidereal day
dt = 60;                                         % time step
t_vec = 0:dt:T_sim;

inc = deg2rad(49.4);
ecc = 0;
raan0 = deg2rad(19.5);

%% === Constellation Search (No, Nso, Nc) ===
Ns_total = 30;              % total number of satellites
best_coverage = 0;          
best_config = [];           % [No, Nso, Nc]

for No = 1:Ns_total
    if mod(Ns_total, No) ~= 0
        continue
    end
    Nso = Ns_total / No;

    for Nc = 1:No-1
        N = Nso * No;
        raan = zeros(N,1);
        M = zeros(N,1);

        % === RAAN–Mean Anomaly Distribution (2D Lattice) ===
        idx = 1;
        for i = 1:No
            for j = 1:Nso
                raan(idx) = raan0+2*pi * (i - 1)/No;
                M(idx) = 2*pi * ((j - 1)/Nso - (Nc/Nso)*(i - 1)/No);
                idx = idx + 1;
            end
        end
        M = mod(M, 2*pi);

        % === Visibility Simulation ===
        visible_total = 0;
        for ti = 1:length(t_vec)
            t = t_vec(ti);
            r_obs = spherical2ECF(lat_obs, lon_obs, alt_obs, R_E);

            visible = false;
            for k = 1:N
                Mk = mod(M(k) + n * t, 2*pi);
                [x_ECI, y_ECI, z_ECI, ~, ~, ~] = kep2cart(a, 0, deg2rad(49.4), raan(k), 0, Mk, mu);
                r_sat = ECI2ECEF([x_ECI, y_ECI, z_ECI], t, omega_E);

                rho = r_sat' - r_obs;

                % Convert to ENU for elevation
                R_ENU = [-sin(lon_obs),              cos(lon_obs),            0;
                         -sin(lat_obs)*cos(lon_obs), -sin(lat_obs)*sin(lon_obs), cos(lat_obs);
                          cos(lat_obs)*cos(lon_obs),  cos(lat_obs)*sin(lon_obs), sin(lat_obs)];
                enu = R_ENU * rho;
                %el = asin(enu(3) / norm(enu));

                zen_hat = r_obs / norm(r_obs);
                el = asin(dot(rho, zen_hat) / norm(rho));

                if el >= el_mask
                    visible = true;
                    break
                end
            end
            if visible
                visible_total = visible_total + dt;
            end
        end

        % === Track Best Configuration ===
        if visible_total > best_coverage
            best_coverage = visible_total;
            best_config = [No, Nso, Nc];  % correct order
        end
    end
end

% === Output Results ===
fprintf('[a] Optimal Constellation Configuration:\n');
fprintf('    No  = %d\n', best_config(1));   % number of orbital planes
fprintf('    Nso = %d\n', best_config(2));   % satellites per orbit
fprintf('    Nc  = %d\n', best_config(3));   % phasing
fprintf('    Total visible time: %.5f hours (%.5f%% of sidereal day)\n', ...
        best_coverage / 3600, 100 * best_coverage / TE);

No = best_config(1);
Nso = best_config(2);
Nc = best_config(3);
N = No*Nso;

%% Part (b): RAAN–Mean Anomaly Distribution Plot

% Allocate arrays
raan = zeros(N, 1);  % Right Ascension of the Ascending Node [rad]
M    = zeros(N, 1);  % Mean Anomaly [rad]

% Fill in using 2D Lattice Flower formulation
idx = 1;
for i = 1:No
    for j = 1:Nso
        raan(idx) = raan0+2*pi * (i - 1) / No;
        M(idx)    = 2*pi * ((j - 1)/Nso - (Nc/Nso) * (i - 1)/No);
        idx = idx + 1;
    end
end

% Wrap into [0, 2pi]
M = mod(M, 2*pi);
raan = mod(raan, 2*pi);

% Plot
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

% === Compute ECI positions of all satellites at t = 0 ===
for k = 1:N
    Mk = M(k);              % mean anomaly (true anomaly since ecc = 0)
    RAAN_k = raan(k);
    [x, y, z, ~, ~, ~] = kep2cart(a, ecc, inc, RAAN_k, 0, Mk, mu);
    Y(k, :) = [x, y, z];
end

% === Polar View (Z-axis top-down) ===
figure('Name', 'ECI Polar View')
earthy(R_E, 'Earth', 0.9, [0; 0; 0]); hold on
scatter3(Y(:,1), Y(:,2), Y(:,3), 50, 'filled', 'b')
xlabel('$X$ [km]', 'Interpreter', 'latex')
ylabel('$Y$ [km]', 'Interpreter', 'latex')
zlabel('$Z$ [km]', 'Interpreter', 'latex')
title('Initial Constellation in ECI (Polar View)', 'Interpreter', 'latex')
axis equal; grid on
view(0, 90)
set(gca, 'TickLabelInterpreter', 'latex')

% Plot full orbits for each satellite
for k = 1:N
    plotOrbit3(raan(k), inc, 0, a, ecc, linspace(0, 2*pi, 100), ...
               '', 1, 1, [0,0,0], 0, 1);
end

% === Isometric View ===
figure('Name', 'ECI Isometric View')
earthy(R_E, 'Earth', 0.9, [0; 0; 0]); hold on
scatter3(Y(:,1), Y(:,2), Y(:,3), 50, 'filled', 'b')
xlabel('$X$ [km]', 'Interpreter', 'latex')
ylabel('$Y$ [km]', 'Interpreter', 'latex')
zlabel('$Z$ [km]', 'Interpreter', 'latex')
title('Initial Constellation in ECI (Isometric View)', 'Interpreter', 'latex')
axis equal; grid on
view(45, 30)
set(gca, 'TickLabelInterpreter', 'latex')

% Plot full orbits again
for k = 1:N
    plotOrbit3(raan(k), inc, 0, a, ecc, linspace(0, 2*pi, 100), ...
               '', 1, 1, [0,0,0], 0, 1);
end

%% Part D

t_track = linspace(0, T_orbit, 1000);

figure
hold on
for j=1:N
    plotGroundTrackLag(a,ecc,inc,raan(j),0,t_track,0,mu,omega_E)
end
plot(rad2deg(lon_obs),rad2deg(lat_obs), 'b*')

%% psrt F

% Compute number of unique ground-tracks
num_tracks = gcd(No, Nc);

fprintf('[f] Number of unique ground-tracks: %d\n', num_tracks);

%% Part g: new GS

%% Part (g): Coverage from Red Pyramid in Egypt

% === Red Pyramid coordinates ===
lat_egypt = deg2rad(29 + 48/60 + 31.4/3600);   % [rad]
lon_egypt = deg2rad(31 + 12/60 + 21.3/3600);   % [rad]
alt_egypt = 0;                                % [km]

% === Initialize coverage counter ===
visible_total_egypt = 0;

for ti = 1:length(t_vec)
    t = t_vec(ti);
    r_obs = spherical2ECF(lat_egypt, lon_egypt, alt_egypt, R_E);

    visible = false;
    for k = 1:N
        Mk = mod(M(k) + n * t, 2*pi);
        [x_ECI, y_ECI, z_ECI, ~, ~, ~] = kep2cart(a, 0, inc, raan(k), 0, Mk, mu);
        r_sat = ECI2ECEF([x_ECI, y_ECI, z_ECI], t, omega_E);

        rho = r_sat' - r_obs;

        % Convert to ENU for elevation using dot product method
        zen_hat = r_obs / norm(r_obs);
        el = asin(dot(rho, zen_hat) / norm(rho));

        if el >= el_mask
            visible = true;
            break
        end
    end

    if visible
        visible_total_egypt = visible_total_egypt + dt;
    end
end

%% === Output Coverage Comparison ===
fprintf('[g] Coverage for Red Pyramid (Egypt): %.5f hours (%.5f%% of sidereal day)\n', ...
    visible_total_egypt / 3600, 100 * visible_total_egypt / TE);

fprintf('[g] Difference from ARMS: %.5f hours\n', ...
    (best_coverage - visible_total_egypt) / 3600);

figure
hold on
for j=1:N
    plotGroundTrackLag(a,ecc,inc,raan(j),0,t_track,0,mu,omega_E)
end
plot(rad2deg(lon_obs),rad2deg(lat_obs), 'b*')
plot(rad2deg(lon_egypt),rad2deg(lat_egypt), 'k*')