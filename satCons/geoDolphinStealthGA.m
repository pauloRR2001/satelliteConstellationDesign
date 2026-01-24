clc; clear; close all

% GEO Dolphin Stealth Constellation via Integer GA
% Goal: Place 10 GEO satellites (circular, equatorial by default) to
% maximize the average Euclidean distance from all satellites to a group
% of 100 dolphins moving NE over 6 months.
% Integer GA optimizes satellite longitudes (degrees, integer).

addpath('..\geneticAlgorithm');
addpath('..\functions');
addpath('..\functions\converters');
addpath('..\functions\plotters');

%% Constants
mu      = 398600.4418;      % km^3/s^2
R_E     = 6378.137;         % km
TE      = 86164.0905;       % sidereal day [s]
omega_E = 7.2921151467e-5;  % rad/s

% Geostationary radius (from Earth's center)
a_geo = (mu / omega_E^2)^(1/3);  % ~42164 km
alt_geo = a_geo - R_E;            % altitude above Earth [km]

%% Dolphin group setup
Nd = 100;                         % number of dolphins
rng(42);                          % deterministic randomness
% Initial lat/lon near equator within +/- 2 deg
lat0_deg = (rand(Nd,1) * 4 - 2);  % [-2, 2]
lon0_deg = (rand(Nd,1) * 4 - 2);  % [-2, 2]

% Monthly NE drift
months = 6;                       % study duration (months)
lat_drift_per_month = 1.0;        % deg/month
lon_drift_per_month = 1.5;        % deg/month

% Build monthly dolphin lat/lon (degrees)
lat_month_deg = zeros(Nd, months+1);
lon_month_deg = zeros(Nd, months+1);
for m = 0:months
    lat_month_deg(:, m+1) = lat0_deg + m * lat_drift_per_month;
    lon_month_deg(:, m+1) = lon0_deg + m * lon_drift_per_month;
    % clamp latitude to [-90, 90] and wrap longitude to [-180, 180]
    lat_month_deg(:, m+1) = max(min(lat_month_deg(:, m+1), 90), -90);
    lon_month_deg(:, m+1) = mod(lon_month_deg(:, m+1) + 180, 360) - 180;
end

%% GA: satellites at GEO
Nsats = 10;
% Chromosome per satellite: [RAAN_deg, INC_deg, F0_deg] (length 3*Nsats)
% Bounds: RAAN ∈ [0, 360], INC ∈ [0, 90], F0 ∈ [0, 360]
lb = zeros(1, 3*Nsats);
ub = zeros(1, 3*Nsats);
for s = 1:Nsats
    lb(3*(s-1)+1) = 0;    ub(3*(s-1)+1) = 360; % RAAN
    lb(3*(s-1)+2) = 0;    ub(3*(s-1)+2) = 90;  % INC
    lb(3*(s-1)+3) = 0;    ub(3*(s-1)+3) = 360; % F0
end

% Initial population
Npop = 60; maxGens = 100;
initPop = zeros(Npop, 3*Nsats);
for i = 1:Npop
    chrom = zeros(1, 3*Nsats);
    for s = 1:Nsats
        chrom(3*(s-1)+1) = randi([0, 360]);  % RAAN
        chrom(3*(s-1)+2) = randi([0, 90]);   % INC
        chrom(3*(s-1)+3) = randi([0, 360]);  % F0
    end
    initPop(i, :) = chrom;
end

% Fitness (minimize): negative average distance (so GA maximizes distance)
objectiveFcn = @(chrom) fitnessGEO(chrom, Nsats, a_geo, R_E, lat_month_deg, lon_month_deg, mu, omega_E, TE);

% GA options
opts.mutationRate   = 0.2;
opts.crossoverRate  = 0.85;
opts.eliteCount     = 2;
opts.tournamentSize = 3;
opts.mutationStep   = 1;           % 1 degree steps
opts.lowerBounds    = lb;
opts.upperBounds    = ub;

[bestChrom, bestFitness, hist] = runIntegerGA(objectiveFcn, initPop, maxGens, opts);

fprintf('Best RAAN/INC/F0 per satellite (deg):\n');
for s = 1:Nsats
    fprintf(' Sat %2d: RAAN=%6.2f, INC=%6.2f, F0=%6.2f\n', s, bestChrom(3*(s-1)+1), bestChrom(3*(s-1)+2), bestChrom(3*(s-1)+3));
end
fprintf('Objective (negative avg distance): %.6f\n', bestFitness);

%% Plot: ECI/ECEF snapshot at final month (m = 6)
% Build satellites ECI then ECEF at t corresponding to month
% Use monthly steps: 1 sidereal month ~ 30*TE seconds
m_final = months;
t_month_seconds = m_final * 30 * TE;  % coarse mapping

% Satellite ECI/ECEF positions
sat_ecef = zeros(Nsats, 3);
for s = 1:Nsats
    inc = deg2rad(bestChrom(3*(s-1)+2));
    raan = deg2rad(bestChrom(3*(s-1)+1));
    argp = 0;
    f0   = deg2rad(bestChrom(3*(s-1)+3));
    n    = sqrt(mu / a_geo^3);
    f_t  = f0 + n * t_month_seconds;  % ECI angle advance
    [x_ECI, y_ECI, z_ECI, ~, ~, ~] = kep2cart(a_geo, 0, inc, raan, argp, f_t, mu);
    r_ecef = ECI2ECEF([x_ECI, y_ECI, z_ECI], t_month_seconds, omega_E);
    sat_ecef(s, :) = r_ecef;
end

% Dolphins ECEF at final month
Dol_ecef = zeros(Nd, 3);
for i = 1:Nd
    latr = deg2rad(lat_month_deg(i, m_final+1));
    lonr = deg2rad(lon_month_deg(i, m_final+1));
    Dol_ecef(i, :) = spherical2ECF(latr, lonr, 0, R_E).';
end

figure('Name','GEO Satellites and Dolphins at Month 6 (ECEF)'); hold on
earthy(R_E, 'Earth', 0.9, [0;0;0]);
scatter3(sat_ecef(:,1), sat_ecef(:,2), sat_ecef(:,3), 40, 'filled', 'MarkerFaceColor', [1,0,0]);
scatter3(Dol_ecef(:,1), Dol_ecef(:,2), Dol_ecef(:,3), 15, 'filled', 'MarkerFaceColor', [0,0.4,1]);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]'); grid on; axis equal
legend({'Earth','Satellites','Dolphins'}, 'Location','bestoutside');
view(45, 20)

% Also plot initial dolphin locations (month 0) in ECEF
Dol0_ecef = zeros(Nd, 3);
for i = 1:Nd
    latr0 = deg2rad(lat_month_deg(i, 1));
    lonr0 = deg2rad(lon_month_deg(i, 1));
    Dol0_ecef(i, :) = spherical2ECF(latr0, lonr0, 0, R_E).';
end
scatter3(Dol0_ecef(:,1), Dol0_ecef(:,2), Dol0_ecef(:,3), 12, 'filled', 'MarkerFaceColor', [0.2,0.8,0.2])
legend({'Earth','Satellites','Dolphins M6','Dolphins M0'}, 'Location','bestoutside');

%% Inertial-frame (ECI) orbits and dolphins start/end
nu_vec = linspace(0, 2*pi, 240);
figure('Name','Sat Orbits (ECI) and Dolphins Start/End (ECI)'); hold on
earthy(R_E, 'Earth', 0.6, [0;0;0]);
% Plot satellite orbits in ECI
for s = 1:Nsats
    raan = deg2rad(bestChrom(3*(s-1)+1));
    inc  = deg2rad(bestChrom(3*(s-1)+2));
    argp = 0;
    XYZ = zeros(numel(nu_vec), 3);
    for m = 1:numel(nu_vec)
        [x,y,z,~,~,~] = kep2cart(a_geo, 0, inc, raan, argp, nu_vec(m), mu);
        XYZ(m,:) = [x,y,z];
    end
    plot3(XYZ(:,1), XYZ(:,2), XYZ(:,3), 'r-', 'LineWidth', 0.8);
end
% Dolphins start (ECI at t=0 equals ECEF)
Dol0_eci = Dol0_ecef;
scatter3(Dol0_eci(:,1), Dol0_eci(:,2), Dol0_eci(:,3), 12, 'filled', 'MarkerFaceColor', [0.2,0.8,0.2]);
% Dolphins final (convert ECEF->ECI at t_month_seconds)
lag = omega_E * t_month_seconds;
NCF = [cos(lag), sin(lag), 0; -sin(lag), cos(lag), 0; 0, 0, 1];
Dol_final_eci = (NCF' * Dol_ecef')';
scatter3(Dol_final_eci(:,1), Dol_final_eci(:,2), Dol_final_eci(:,3), 15, 'filled', 'MarkerFaceColor', [0,0.4,1]);
xlabel('X_{ECI} [km]'); ylabel('Y_{ECI} [km]'); zlabel('Z_{ECI} [km]'); grid on; axis equal
legend({'Earth','Sat orbits','Dolphins M0 (ECI)','Dolphins M6 (ECI)'}, 'Location','bestoutside');
view(45, 20)

%% Rotating-frame (ECEF) satellite trajectories over 6 months
plot_dt = 6*3600;  % 6 hours
t_plot = 0:plot_dt:(months*30*TE);
Nt_plot = numel(t_plot);
figure('Name','Sat Trajectories over 6 Months (ECEF)'); hold on
earthy(R_E, 'Earth', 0.8, [0;0;0]);
for s = 1:Nsats
    raan = deg2rad(bestChrom(3*(s-1)+1));
    inc  = deg2rad(bestChrom(3*(s-1)+2));
    f0   = deg2rad(bestChrom(3*(s-1)+3));
    argp = 0; n = sqrt(mu / a_geo^3);
    X = zeros(Nt_plot,1); Y = zeros(Nt_plot,1); Z = zeros(Nt_plot,1);
    for ti = 1:Nt_plot
        t = t_plot(ti);
        f_t = f0 + n*t;
        [xECI,yECI,zECI,~,~,~] = kep2cart(a_geo, 0, inc, raan, argp, f_t, mu);
        rECEF = ECI2ECEF([xECI,yECI,zECI], t, omega_E);
        X(ti) = rECEF(1); Y(ti) = rECEF(2); Z(ti) = rECEF(3);
    end
    plot3(X, Y, Z, 'r-', 'LineWidth', 0.8);
end
% Dolphins M0 and M6 in ECEF
scatter3(Dol0_ecef(:,1), Dol0_ecef(:,2), Dol0_ecef(:,3), 12, 'filled', 'MarkerFaceColor', [0.2,0.8,0.2]);
scatter3(Dol_ecef(:,1), Dol_ecef(:,2), Dol_ecef(:,3), 15, 'filled', 'MarkerFaceColor', [0,0.4,1]);
xlabel('X_{ECEF} [km]'); ylabel('Y_{ECEF} [km]'); zlabel('Z_{ECEF} [km]'); grid on; axis equal
legend({'Earth','Sat trajectories','Dolphins M0 (ECEF)','Dolphins M6 (ECEF)'}, 'Location','bestoutside');
view(45, 20)

%% Ground tracks (lat–lon) over 6 months (ECEF)
% Expect near-constant longitude for geostationary (inc≈0) and analemmas for inclined GEO
figure('Name','Sat Ground Tracks and Dolphins (Lat–Lon)'); hold on
% Plot world coastlines for context
load('topo.mat','topo'); topoplot = [topo(:,181:360), topo(:,1:180)];
contour(-180:179, -90:89, topoplot, [0, 0], 'k');
for s = 1:Nsats
    raan = deg2rad(bestChrom(3*(s-1)+1));
    inc  = deg2rad(bestChrom(3*(s-1)+2));
    f0   = deg2rad(bestChrom(3*(s-1)+3));
    argp = 0; n = sqrt(mu / a_geo^3);
    lat_deg = zeros(Nt_plot,1); lon_deg = zeros(Nt_plot,1);
    for ti = 1:Nt_plot
        t = t_plot(ti);
        f_t = f0 + n*t;
        [xECI,yECI,zECI,~,~,~] = kep2cart(a_geo, 0, inc, raan, argp, f_t, mu);
        rECEF = ECI2ECEF([xECI,yECI,zECI], t, omega_E);
        rmag = norm(rECEF);
        lat = asin(rECEF(3)/rmag);
        lon = atan2(rECEF(2), rECEF(1));
        lat_deg(ti) = rad2deg(lat);
        lon_deg(ti) = rad2deg(lon);
    end
    % Handle wrap-around for plotting continuity
    jumps = find(abs(diff(lon_deg)) > 180);
    lon_deg(jumps+1) = NaN;
    plot(lon_deg, lat_deg, 'r-', 'LineWidth', 0.8);
end
% Dolphins start/end in lat–lon
Dol0_lat = lat_month_deg(:,1); Dol0_lon = lon_month_deg(:,1);
Dol6_lat = lat_month_deg(:,months+1); Dol6_lon = lon_month_deg(:,months+1);
scatter(Dol0_lon, Dol0_lat, 12, [0.2,0.8,0.2], 'filled');
scatter(Dol6_lon, Dol6_lat, 15, [0,0.4,1], 'filled');
xlabel('Longitude [deg]'); ylabel('Latitude [deg]'); grid on
xlim([-180 180]); ylim([-90 90]); xticks(-180:30:180); yticks(-90:15:90);
legend({'Coastlines','Sat ground tracks','Dolphins M0','Dolphins M6'}, 'Location','bestoutside');

figure('Name','GA Fitness per Generation')
plot(0:maxGens, hist.fitnessTrace, '-o'); grid on
xlabel('Generation'); ylabel('Best Fitness (neg avg distance)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitness function: negative average distance over 6 months
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = fitnessGEO(chrom, Nsats, a_geo, R_E, lat_month_deg, lon_month_deg, mu, omega_E, TE)
    Nd = size(lat_month_deg,1);
    months = size(lat_month_deg,2) - 1;

    n = sqrt(mu / a_geo^3);

    avg_dist_accum = 0; total_counts = 0;
    for m = 0:months
        t_m = m * 30 * TE;  % coarse month mapping
        sat_ecef = zeros(Nsats, 3);
        for s = 1:Nsats
            raan = deg2rad(chrom(3*(s-1)+1));
            inc  = deg2rad(chrom(3*(s-1)+2));
            f0   = deg2rad(chrom(3*(s-1)+3));
            argp = 0;
            f_t  = f0 + n * t_m;
            [x_ECI, y_ECI, z_ECI, ~, ~, ~] = kep2cart(a_geo, 0, inc, raan, argp, f_t, mu);
            sat_ecef(s, :) = ECI2ECEF([x_ECI, y_ECI, z_ECI], t_m, omega_E);
        end
        % Dolphins ECEF at month m
        Dol_ecef = zeros(Nd, 3);
        for i = 1:Nd
            latr = deg2rad(lat_month_deg(i, m+1));
            lonr = deg2rad(lon_month_deg(i, m+1));
            Dol_ecef(i, :) = spherical2ECF(latr, lonr, 0, R_E).';
        end
        % Average Euclidean distance across all pairs
        for s = 1:Nsats
            for i = 1:Nd
                d = norm(sat_ecef(s,:) - Dol_ecef(i,:));
                avg_dist_accum = avg_dist_accum + d;
            end
        end
        total_counts = total_counts + Nsats * Nd;
    end

    avg_distance = avg_dist_accum / total_counts;
    f = -avg_distance;  % maximize average distance by minimizing negative
end
