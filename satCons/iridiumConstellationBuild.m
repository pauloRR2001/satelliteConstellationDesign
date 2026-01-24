clc
clear
close all

%%  Given and Constants 

% Parameters (shared between parts)
t = 66;           % Total satellites
p = 6;            % Number of orbital planes
f = 2;            % Phasing parameter
inc = 86.4;       % Inclination [deg]
a = 7158;         % Semi-major axis [km]
earth_radius = 6378.1363; % Earth radius [km]
mu = 3.986005e5;     % Earth's gravitational parameter [km³/s²]

%% Part a) Satellite Distribution in Raan-M Space
% Initialize arrays
omega = zeros(t, 1);
mean_anomaly = zeros(t, 1);

% Calculate RAAN and Mean Anomaly for each satellite
idx = 1;
for i = 1:p
    for j = 1:t/p
        delta_omega = (2 * pi * (i - 1) / p)/2;
        delta_M = 2 * pi * (p / t) * (j - 1) + 2 * pi * (f / t) * (i - 1);
        
        omega_deg = rad2deg(delta_omega);
        M_deg = rad2deg(delta_M);
        
        % Normalize to [0, 360) degrees
        omega_deg = mod(omega_deg, 360);
        M_deg = mod(M_deg, 360);
        
        omega(idx) = omega_deg;
        mean_anomaly(idx) = M_deg;
        idx = idx + 1;
    end
end

% Plot Ω-M distribution
figure;
scatter(omega, mean_anomaly, 20, 'filled');
xlabel('Right Ascension of the Ascending Node (Ω) [degrees]');
ylabel('Mean Anomaly (M) [degrees]');
title('Distribution of Satellites in the Iridium Constellation');
grid on;
xticks(0:60:360);
yticks(0:60:360);
xlim([0 360]);
ylim([0 360]);

%% Part b) 3D Visualization 
% Initialize arrays
sat_positions = zeros(t, 3); % ECEF positions [km]

% Reuse omega values from Part a (with /2 term)
idx = 1;
for i = 1:p
    for j = 1:t/p
        raan = deg2rad(omega(idx)); % Uses RAAN from Part a (with /2 term)
        M = deg2rad(mean_anomaly(idx)); % Uses mean anomaly from Part a
        
        % Convert Keplerian to ECEF (circular orbit: ecc=0, argp=0, f=M)
        [x, y, z, ~, ~, ~] = kep2cart(a, 0, deg2rad(inc), raan, 0, M, mu);
        sat_positions(idx, :) = [x, y, z];
        idx = idx + 1;
    end
end

% Generate orbit paths with matching RAAN (/2 term)
num_points = 100; % Points per orbit
orbit_points = zeros(p, num_points, 3);
for i = 1:p
    raan = deg2rad(360 * (i-1) / p / 2); % Added /2 to match Part a
    for k = 1:num_points
        M_k = deg2rad(360 * (k-1) / num_points);
        [x, y, z, ~, ~, ~] = kep2cart(a, 0, deg2rad(inc), raan, 0, M_k, mu);
        orbit_points(i, k, :) = [x, y, z];
    end
end

% Polar view
figure;
earthy(earth_radius, 'Earth', 0.9, [0; 0; 0]);
hold on;

% Plot orbits (6 planes, 30° RAAN spacing due to /2 term in part a, streets of coverage)
colors = lines(p);
for i = 1:p
    plot3(squeeze(orbit_points(i, :, 1)), squeeze(orbit_points(i, :, 2)), squeeze(orbit_points(i, :, 3)), ...
         'Color', colors(i,:), 'LineWidth', 1.5);
end

% Plot satellites
scatter3(sat_positions(:, 1), sat_positions(:, 2), sat_positions(:, 3), 30, 'k', 'filled');

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Iridium Constellation - Polar View');
axis equal;
grid on;
view(0, 90); % Top-down view aka polar 

% Isometric view
figure;
earthy(earth_radius, 'Earth', 0.9, [0; 0; 0]);
hold on;
for i = 1:p
    plot3(squeeze(orbit_points(i, :, 1)), squeeze(orbit_points(i, :, 2)), squeeze(orbit_points(i, :, 3)), ...
         'Color', colors(i,:), 'LineWidth', 1.5);
end
scatter3(sat_positions(:, 1), sat_positions(:, 2), sat_positions(:, 3), 30, 'k', 'filled');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Iridium Constellation - Isometric View');
axis equal;
grid on;
view(45, 30); % Isometric view

%% Part c) Minimum distance between two co-orbital satellites
delta_M_deg = 360 / (t/p); % Angular separation in mean anomaly [deg]
delta_M_rad = deg2rad(delta_M_deg);

% Minimum distance formula (chord length)
min_dis_same_orbit = a * abs(sin(delta_M_rad / 2))