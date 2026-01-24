clc; clear; close all

%% Spy Flower Constellation Optimization (minimize Iridium coverage)
% Objective: Choose Flower integer parameters (No, Nso, Nc) to MINIMIZE
% unobstructed satellite-to-satellite visibility between the Spy constellation
% and the Iridium constellation over one sidereal month.
% Visibility metric: Unobstructed line-of-sight (LOS) between any Spy–Iridium
% pair at a time step (i.e., line segment does not intersect Earth sphere).
% Coverage is the total time with any visible pair.

%% Constants
mu      = 398600.4418;      % km^3/s^2
R_E     = 6378.137;         % km
TE      = 86164.0905;       % sidereal day [s]
omega_E = 7.2921151467e-5;  % Earth rotation rate [rad/s]

%% Iridium Constellation (fixed)
t_ir   = 66;                % total satellites
p_ir   = 6;                 % planes
f_ir   = 2;                 % phasing
inc_ir = deg2rad(86.4);     % radians
a_ir   = 7158;              % km (from problem4)
ecc_ir = 0;                 % circular
argp_ir= 0;                 % argument of perigee
n_ir   = sqrt(mu/a_ir^3);   % mean motion [rad/s]

[raan_ir, M0_ir] = buildIridium(p_ir, t_ir, f_ir);

%% Spy Flower Constellation parameters (to search)
Ns_total   = 30;            % total spy satellites (fixed for fair comparison)
a_spy      = a_ir;          % km (match altitude for geometry relevance)
inc_spy    = deg2rad(49.4); % radians (typical Flower incl.)
ecc_spy    = 0;
argp_spy   = 0;
n_spy      = sqrt(mu/a_spy^3);
raan0_spy  = 0;             % base RAAN offset [rad]

%% Simulation setup
num_days   = 30;            % one sidereal month (~30 days)
T_sim      = num_days * TE; % [s]
dt         = 600;           % step [s] (10 min; adjust for speed/accuracy)
t_vec      = 0:dt:T_sim;
Nt         = numel(t_vec);

%% Precompute Iridium positions for all time steps (performance)
posIr = zeros(Nt, t_ir, 3);
for ti = 1:Nt
    t = t_vec(ti);
    Mk = mod(M0_ir + n_ir * t, 2*pi);
    for k = 1:t_ir
        [x, y, z, ~, ~, ~] = kep2cart(a_ir, ecc_ir, inc_ir, raan_ir(k), argp_ir, Mk(k), mu);
        posIr(ti, k, :) = [x, y, z];
    end
end

%% Search over Flower integer parameters (No, Nso, Nc)
best_coverage = inf;
best_config   = [NaN, NaN, NaN];  % [No, Nso, Nc]

for No = 1:Ns_total
    if mod(Ns_total, No) ~= 0
        continue
    end
    Nso = Ns_total / No;

    for Nc = 1:No-1
        % Build Spy Flower RAAN–M distribution
        [raan_spy, M0_spy] = buildFlower(No, Nso, Nc, raan0_spy);

        % Accumulate coverage time (seconds) when ANY Spy–Iridium pair is visible
        coverage_seconds = 0;

        for ti = 1:Nt
            t = t_vec(ti);
            Mk_spy = mod(M0_spy + n_spy * t, 2*pi);

            % Spy ECI positions at time t
            Spos = zeros(Ns_total, 3);
            for s = 1:Ns_total
                [xs, ys, zs, ~, ~, ~] = kep2cart(a_spy, ecc_spy, inc_spy, raan_spy(s), argp_spy, Mk_spy(s), mu);
                Spos(s, :) = [xs, ys, zs];
            end

            % Iridium ECI positions at time t (precomputed)
            Ipos = squeeze(posIr(ti, :, :));  % [66 x 3]

            % If any pair visible at time t, count dt
            if anyPairVisible(Spos, Ipos, R_E)
                coverage_seconds = coverage_seconds + dt;
            end
        end

        % Track best (minimal) coverage
        if coverage_seconds < best_coverage
            best_coverage = coverage_seconds;
            best_config   = [No, Nso, Nc];
        end
    end
end

%% Report results
fprintf('Spy Flower optimization complete over %d days, dt = %d s\n', num_days, dt);
fprintf('Best config: No = %d, Nso = %d, Nc = %d (Ns_total = %d)\n', best_config(1), best_config(2), best_config(3), Ns_total);
fprintf('Minimum Iridium coverage: %.3f hours (%.3f%% of month)\n', best_coverage/3600, 100*best_coverage/T_sim);

%% Optional: Plot RAAN–M distribution for best Spy config
[raan_best, M0_best] = buildFlower(best_config(1), best_config(2), best_config(3), raan0_spy);
figure('Name','Spy RAAN–M Distribution');
scatter(rad2deg(mod(raan_best,2*pi)), rad2deg(mod(M0_best,2*pi)), 40, 'filled');
xlabel('RAAN \Omega [deg]', 'Interpreter','latex');
ylabel('Mean Anomaly M [deg]', 'Interpreter','latex');
title('Best Spy Flower: RAAN–Mean Anomaly', 'Interpreter','latex');
grid on; axis([0 360 0 360]); xticks(0:60:360); yticks(0:60:360);
set(gca,'TickLabelInterpreter','latex');

%% Plot 3D orbits for both constellations (Spy=red, Iridium=blue)
nu_vec = linspace(0, 2*pi, 200);

figure('Name','Spy (red) vs Iridium (blue) Orbits')
hold on

% Iridium: plot orbit lines (blue)
for k = 1:t_ir
    XYZ = zeros(numel(nu_vec), 3);
    for m = 1:numel(nu_vec)
        [x, y, z, ~, ~, ~] = kep2cart(a_ir, ecc_ir, inc_ir, raan_ir(k), argp_ir, nu_vec(m), mu);
        XYZ(m, :) = [x, y, z];
    end
    plot3(XYZ(:,1), XYZ(:,2), XYZ(:,3), 'Color', [0, 0, 1], 'LineWidth', 1.0)
end

% Spy: plot orbit lines (red)
for s = 1:Ns_total
    XYZ = zeros(numel(nu_vec), 3);
    for m = 1:numel(nu_vec)
        [x, y, z, ~, ~, ~] = kep2cart(a_spy, ecc_spy, inc_spy, raan_best(s), argp_spy, nu_vec(m), mu);
        XYZ(m, :) = [x, y, z];
    end
    plot3(XYZ(:,1), XYZ(:,2), XYZ(:,3), 'Color', [1, 0, 0], 'LineWidth', 1.0)
end
%earthy(R_E, 'Earth', 0.9, [0; 0; 0]); hold on

% Dots at satellite positions at t = 0
% Iridium dots (blue)
I0 = zeros(t_ir, 3);
for k = 1:t_ir
    [x, y, z, ~, ~, ~] = kep2cart(a_ir, ecc_ir, inc_ir, raan_ir(k), argp_ir, M0_ir(k), mu);
    I0(k, :) = [x, y, z];
end
scatter3(I0(:,1), I0(:,2), I0(:,3), 20, 'filled', 'MarkerFaceColor', [0, 0, 1], 'MarkerEdgeColor', 'none')

% Spy dots (red)
S0 = zeros(Ns_total, 3);
for s = 1:Ns_total
    [x, y, z, ~, ~, ~] = kep2cart(a_spy, ecc_spy, inc_spy, raan_best(s), argp_spy, M0_best(s), mu);
    S0(s, :) = [x, y, z];
end
scatter3(S0(:,1), S0(:,2), S0(:,3), 20, 'filled', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', 'none')

xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('3D Orbits: Spy (red) and Iridium (blue)');
axis equal; grid on
view(45, 30)
legend({'Iridium orbits','Spy orbits','Iridium sats','Spy sats'}, 'Location', 'bestoutside')

%% Rotating-frame (ECEF) trajectories over full 30 days
% Downsample plotting to keep performance reasonable
plot_dt = 1800;                        % 30 min sampling for plotting
t_plot  = 0:plot_dt:T_sim;
Nt_plot = numel(t_plot);

figure('Name','ECEF Trajectories over 30 Days (Spy=red, Iridium=blue)')
hold on
earthy(R_E, 'Earth', 0.8, [0; 0; 0]);

% Iridium ECEF trajectories (blue)
for k = 1:t_ir
    X = zeros(Nt_plot,1); Y = zeros(Nt_plot,1); Z = zeros(Nt_plot,1);
    for ti = 1:Nt_plot
        t = t_plot(ti);
        Mk = mod(M0_ir(k) + n_ir * t, 2*pi);
        [x, y, z, ~, ~, ~] = kep2cart(a_ir, ecc_ir, inc_ir, raan_ir(k), argp_ir, Mk, mu);
        r_ecef = ECI2ECEF([x, y, z], t, omega_E);
        X(ti) = r_ecef(1); Y(ti) = r_ecef(2); Z(ti) = r_ecef(3);
    end
    plot3(X, Y, Z, 'Color', [0, 0, 1], 'LineWidth', 0.8)
end

% Spy ECEF trajectories (red) using best config
for s = 1:Ns_total
    X = zeros(Nt_plot,1); Y = zeros(Nt_plot,1); Z = zeros(Nt_plot,1);
    for ti = 1:Nt_plot
        t = t_plot(ti);
        Mk = mod(M0_best(s) + n_spy * t, 2*pi);
        [x, y, z, ~, ~, ~] = kep2cart(a_spy, ecc_spy, inc_spy, raan_best(s), argp_spy, Mk, mu);
        r_ecef = ECI2ECEF([x, y, z], t, omega_E);
        X(ti) = r_ecef(1); Y(ti) = r_ecef(2); Z(ti) = r_ecef(3);
    end
    plot3(X, Y, Z, 'Color', [1, 0, 0], 'LineWidth', 0.8)
end

% Dots at t=0 positions in ECEF
I0_ecef = zeros(t_ir, 3);
for k = 1:t_ir
    [x, y, z, ~, ~, ~] = kep2cart(a_ir, ecc_ir, inc_ir, raan_ir(k), argp_ir, M0_ir(k), mu);
    I0_ecef(k, :) = ECI2ECEF([x, y, z], 0, omega_E);
end
scatter3(I0_ecef(:,1), I0_ecef(:,2), I0_ecef(:,3), 10, 'filled', 'MarkerFaceColor', [0, 0, 1], 'MarkerEdgeColor', 'none')

S0_ecef = zeros(Ns_total, 3);
for s = 1:Ns_total
    [x, y, z, ~, ~, ~] = kep2cart(a_spy, ecc_spy, inc_spy, raan_best(s), argp_spy, M0_best(s), mu);
    S0_ecef(s, :) = ECI2ECEF([x, y, z], 0, omega_E);
end
scatter3(S0_ecef(:,1), S0_ecef(:,2), S0_ecef(:,3), 10, 'filled', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', 'none')

xlabel('X_{ECEF} [km]'); ylabel('Y_{ECEF} [km]'); zlabel('Z_{ECEF} [km]');
title('Earth-Fixed Trajectories over 30 Days');
axis equal; grid on
view(45, 30)
legend({'Iridium ECEF traj','Spy ECEF traj','Iridium t=0','Spy t=0'}, 'Location', 'bestoutside')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions (local)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [raan, M] = buildIridium(p, t, f)
    raan = zeros(t,1);
    M    = zeros(t,1);
    idx = 1;
    for i = 1:p
        for j = 1:t/p
            raan(idx) = mod((2*pi*(i-1)/p)/2, 2*pi); % half-spacing per problem4
            M(idx)    = mod(2*pi*(p/t)*(j-1) + 2*pi*(f/t)*(i-1), 2*pi);
            idx = idx + 1;
        end
    end
end

function [raan, M] = buildFlower(No, Nso, Nc, raan0)
    N = No * Nso;
    raan = zeros(N,1);
    M    = zeros(N,1);
    idx = 1;
    for i = 1:No
        for j = 1:Nso
            raan(idx) = mod(raan0 + 2*pi*(i-1)/No, 2*pi);
            M(idx)    = 2*pi * ((j-1)/Nso - (Nc/Nso)*(i-1)/No);
            idx = idx + 1;
        end
    end
    M = mod(M, 2*pi);
end

function vis = anyPairVisible(Spos, Ipos, R_E)
    % Returns true if ANY Spy–Iridium pair has unobstructed LOS (Earth not blocking).
    Ns = size(Spos,1);
    Ni = size(Ipos,1);
    vis = false;
    for s = 1:Ns
        r1 = Spos(s, :).';
        for i = 1:Ni
            r2 = Ipos(i, :).';
            if unobstructedLOS(r1, r2, R_E)
                vis = true;
                return
            end
        end
    end
end

function ok = unobstructedLOS(r1, r2, R_E)
    % Earth occlusion test via closest approach from origin to line segment r1->r2.
    d  = r2 - r1;
    d2 = dot(d, d);
    if d2 == 0
        ok = (norm(r1) > R_E); % coincident points edge-case
        return
    end
    s = -dot(r1, d) / d2;     % param along segment
    if s < 0
        p = r1;                % closest at r1
    elseif s > 1
        p = r2;                % closest at r2
    else
        p = r1 + s * d;        % closest point on segment
    end
    ok = (norm(p) > R_E);      % LOS unobstructed if segment stays outside Earth
end
