% AAE 590 SC Homework 6 - Problem 3
% Author: Paulo Ramirez
% Institution: Purdue University
% Course: AAE 590 - Satellite Constellations

clc; clear; close all;

%% Constants and Parameters
mu = 398600.4418;              % Gravitational parameter [km^3/s^2]
R_E = 6378.137;                % Earth radius [km]
omega_E = 7.2921151467e-5;     % Earth's rotation rate [rad/s]
h = 705;                       % Orbit altitude [km]
an = R_E + h;                  % Nominal semi-major axis [km]

% Aqua Satellite
Cd_aqua = 2.2;
S_aqua = 38.4 * 1e-6;          % [km^2]
m_aqua = 2934;                 % [kg]
S_m_aqua = S_aqua / m_aqua;    % Ballistic coefficient [km^2/kg]
rho_aqua = atmDensity(h);      % [kg/km^3]

% GCOM-W1 Satellite
Cd_gcom = 2.2;
S_gcom = 26.95 * 1e-6;         % [km^2]
m_gcom = 1991;                 % [kg]
S_m_gcom = S_gcom / m_gcom;
rho_gcom = atmDensity(h);

% Max cross-track deviation Aqua [km]
Delta_d_aqua = 20;
Delta_lambda_aqua = Delta_d_aqua / R_E; % [rad]

%% (a) Time between maneuvers (equal for both satellites)
t_m = 4 * sqrt(Delta_lambda_aqua / (3 * omega_E * rho_aqua * Cd_aqua * S_m_aqua) * sqrt(an / mu));
td = 2 * 86400; % Time delay between maneuvers [s]

%% (b) Maneuver size for GCOM-W1 under equal TM
dadt_gcom = rho_gcom * Cd_gcom * S_m_gcom * sqrt(mu * an);  % [km/s]
Delta_a_gcom = dadt_gcom * t_m;  % [km]
del_a0_gcom = Delta_a_gcom / 2;

%% (c) Impulse required for GCOM-W1
vc = sqrt(mu / an);  % Circular orbital speed [km/s]
Delta_v_gcom = 0.5 * vc * (Delta_a_gcom / an);  % [km/s]

fprintf('[a] Time between maneuvers (days): %.4f\n', t_m / 86400);
fprintf('[b] Maneuver size GCOM-W1 (m): %.4f\n', Delta_a_gcom * 1e3);
fprintf('[c] Required delta-v GCOM-W1 (m/s): %.4f\n', Delta_v_gcom * 1e3);

%% (d) Trajectory evolution
% Use max along-track deviation from ballistic ratio
delta_tau_aqua = Delta_lambda_aqua / omega_E;
delta_tau_gcom = delta_tau_aqua * (S_m_gcom / S_m_aqua);
delta_lambda_gcom = delta_tau_gcom * omega_E;

% Time span
tspan = linspace(0, t_m, 1000);

% Aqua trajectory
del_a0_aqua = (4 * sqrt(Delta_lambda_aqua / (3 * omega_E) * rho_aqua * Cd_aqua * S_m_aqua * sqrt(mu * an^3))) / 2;
aqua_tau = -1.5 * del_a0_aqua / an * tspan + 0.75 * rho_aqua * Cd_aqua * S_m_aqua * sqrt(mu / an) * tspan.^2;
aqua_lambda = omega_E * aqua_tau;

% GCOM-W1 trajectory
gcom_tau = -1.5 * del_a0_gcom / an * tspan + 0.75 * rho_gcom * Cd_gcom * S_m_gcom * sqrt(mu / an) * tspan.^2;
gcom_lambda = omega_E * gcom_tau;

% Altitudes from semi-major axis decay (via ode45)
aqua_struct = ode45(@(t, y) SMArate(t, y, an, mu, S_m_aqua, Cd_aqua, R_E, omega_E), tspan, an + del_a0_aqua);
gcom_struct = ode45(@(t, y) SMArate(t, y, an, mu, S_m_gcom, Cd_gcom, R_E, omega_E), tspan, an + del_a0_gcom);
aqua_alt = deval(aqua_struct, tspan) - an;
gcom_alt = deval(gcom_struct, tspan) - an;

% Separation computation
tau_diff = 85.5;  % [s]
sep = tau_diff + delta_tau_aqua/2 + delta_tau_gcom/2;
sep = tau_diff+86;

%% Plot (d)
figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1)
hold on
plot(rad2deg(aqua_lambda+(Delta_lambda_aqua)/2), aqua_alt, 'b', 'LineWidth', 1.5)
plot(rad2deg(gcom_lambda+(delta_lambda_gcom)/2), gcom_alt, 'r', 'LineWidth', 1.5)
plot([rad2deg(aqua_lambda(1)+(Delta_lambda_aqua)/2),...
    rad2deg(aqua_lambda(end)+(Delta_lambda_aqua)/2)], [aqua_alt(1),aqua_alt(end)], 'b', 'LineWidth', 1.5)
plot([rad2deg(gcom_lambda(1)+(delta_lambda_gcom)/2),...
    rad2deg(gcom_lambda(end)+(delta_lambda_gcom)/2)], [gcom_alt(1),gcom_alt(end)], 'r', 'LineWidth', 1.5)
xlabel('Cross-Track [deg]','interpreter','latex');
ylabel('Altitude Variation [km]','interpreter','latex');
title('Cross-Track vs Altitude','interpreter','latex'); grid minor
legend('Aqua', 'GCOM-W1','interpreter','latex','Location','bestoutside'); axis tight

subplot(1,2,2)
hold on
plot(aqua_tau, aqua_alt, 'b', 'LineWidth', 1.5)
plot(gcom_tau + sep, gcom_alt, 'r', 'LineWidth', 1.5)
plot([aqua_tau(1),aqua_tau(end)], [aqua_alt(1),aqua_alt(end)], 'b', 'LineWidth', 1.5)
plot([gcom_tau(1),gcom_tau(end)] + sep, [gcom_alt(1),gcom_alt(end)], 'r', 'LineWidth', 1.5)
xlabel('Along-Track [s]','interpreter','latex');
ylabel('Altitude Variation [km]','interpreter','latex');
title('Along-Track vs Altitude','interpreter','latex'); grid minor
%legend('Aqua', 'GCOM-W1','interpreter','latex');
sgtitle('Variations to Nominal Orbit','interpreter','latex')

%% (e) Relative Along-Track
rel_tau = (gcom_tau - aqua_tau) + sep;
figure
plot(tspan / 86400, rel_tau, 'k', 'LineWidth', 1.5)
xlabel('Time [days]','interpreter','latex');
ylabel('Relative Along-Track [s]','interpreter','latex');
title('Relative Along-Track: GCOM-W1 vs Aqua','interpreter','latex'); grid minor

%% (f) Tumbling case
tumble_time = t_m - 2*24*3600;
[~, tumble_idx] = min(abs(tspan - tumble_time));

tspan_pretumble = tspan(1:tumble_idx);
tspan_posttumble = tspan(tumble_idx+1:end);

S_gcom_tumble = 2 * S_gcom;
S_m_gcom_tumble = S_gcom_tumble / m_gcom;

gcom_tau_pretumble = -1.5 * del_a0_gcom / an * tspan_pretumble + ...
    0.75 * rho_gcom * Cd_gcom * S_m_gcom_tumble * sqrt(mu / an) * tspan_pretumble.^2;

tail = (tspan_posttumble-tspan_posttumble(1))*0.0001;

gcom_tau_posttumble = -1.5 * del_a0_gcom / an * tspan_posttumble + ...
    0.75 * rho_gcom * Cd_gcom * S_m_gcom_tumble * sqrt(mu / an) * tspan_posttumble.^2+tail;

gcom_tau_tumble = [gcom_tau_pretumble,gcom_tau_posttumble];

rel_tau_tumble = (gcom_tau_tumble - aqua_tau) + tau_diff +...
    delta_tau_aqua/2 + delta_tau_gcom/2;

figure
plot(tspan / 86400, rel_tau_tumble, 'm', 'LineWidth', 1.5)
xlabel('Time [days]','interpreter','latex');
ylabel('Relative Along-Track [s]','interpreter','latex');
title('Tumbling Case: Relative Along-Track','interpreter','latex'); grid minor
legend('GCOM-W1 (Tumbling)','interpreter','latex')

%% Functions
function adot = SMArate(t,a,an,mu,S_m,Cd,R_E, omega_E)
    acc_D = simpleDrag(a, mu, S_m, Cd, R_E, omega_E);
    adot = 2 * sqrt((a^3)/mu) * acc_D;
end

function acc_D = simpleDrag(a, mu, S_m, Cd, R_E, omega_E)

    r_vec = a*[1, 0, 0];
    r = norm(r_vec);
    
    v_c = sqrt(mu/r);
    v_vec = v_c*[0, 1, 0];

    alt = r-R_E;
    rho = atmDensity(alt);

    %v_rel_vec = v_vec - cross(omega_E*[0,0,1], r_vec);
    v_rel_vec = v_vec; % ignore earth roatation
    v_rel = norm(v_rel_vec);

    % Compute drag perturbation acceleration vector (km/s²)
    acc_D = -0.5*rho*Cd*S_m*(v_rel^2)*v_rel_vec/norm(v_rel_vec);

    acc_D = acc_D(2);
end
