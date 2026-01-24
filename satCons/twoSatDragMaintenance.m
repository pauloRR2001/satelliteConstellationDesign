clc; clear; close all

%% Constants
mu = 398600.4418;              % km^3/s^2
R_E = 6378.137;                 % km
omega_E = 7.2921151467e-5;      % rad/s

h = 705;                       % km
an = R_E + h;                  % nominal SMA

%% Aqua
Cd_aqua = 2.2;
S_aqua = 38.4 * (1/1000^2);                      % km^2
m_aqua = 2934;                      % kg
S_m_aqua = S_aqua/m_aqua;              % km^2/kg

rho_aqua = atmDensity(h);

Delta_d_aqua = 20;                  % max cross-track deviation in km
Delta_lambda_aqua = Delta_d_aqua / R_E; % convert arc to angle [rad]

t_m_aqua = 4*sqrt(Delta_lambda_aqua/(3*omega_E*rho_aqua*Cd_aqua*S_m_aqua)*sqrt(an/mu));

Delta_a_aqua = 4*sqrt(Delta_lambda_aqua/(3*omega_E)*rho_aqua*Cd_aqua*S_m_aqua*sqrt(mu*an^3));

del_a0_aqua = Delta_a_aqua/2;

a0_aqua = an + del_a0_aqua;

t_aqua = linspace(0, t_m_aqua, 1000);
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

a_struct_aqua = ode45(@(t_aqua,Y) SMArate(t_aqua,Y,an,mu,S_m_aqua,Cd_aqua,R_E,omega_E),t_aqua,a0_aqua,opt);
a_aqua = deval(a_struct_aqua, t_aqua);

del_lambda_aqua = -3/2*omega_E*del_a0_aqua/an*t_aqua+...
    3/4*omega_E*rho_aqua*Cd_aqua*S_m_aqua*sqrt(mu/an)*t_aqua.^2;

del_tau_aqua = -3/2*del_a0_aqua/an*t_aqua + 3/4*rho_aqua*S_m_aqua*Cd_aqua*sqrt(mu/an)*t_aqua.^2;

Delta_tau_aqua = del_tau_aqua(end);

%% GCOM-W1 Parameters
Cd_gcom = 2.2;
S_gcom = 26.95 * (1/1000^2);       % km^2
m_gcom = 1991;                     % kg
S_m_gcom = S_gcom / m_gcom;                  % km^2/kg

rho_gcom = atmDensity(h);        % kg/km^3

% Constraints
Delta_tau_gcom = 86;              % s

% Maneuver Timing
t_m_gcom = 4 * sqrt(Delta_tau_gcom / (3 * rho_gcom * Cd_gcom * S_m_gcom) * sqrt(an / mu));
%t_m_aqua = 4 * sqrt(Delta_tau_aqua / (3 * rho_aqua * Cd_gcom * S_m_aqua) * sqrt(an / mu));
fprintf('[a] Time between maneuvers:\n GCOM-W1: %.4f days\n Aqua: %.4f days\n', ...
    t_m_gcom/86400, t_m_aqua/86400)

%% Maneuver Size
Delta_a_gcom = 4 * sqrt(Delta_tau_gcom / 3 * rho_gcom * Cd_gcom * S_m_gcom * sqrt(mu * an^3));
del_a0_gcom = Delta_a_gcom / 2;
fprintf('[b] Maneuver size:\n GCOM-W1: %.4f m\n Aqua: %.4f m\n', ...
    Delta_a_gcom*1e3, Delta_a_aqua*1e3)

%% Impulse
v0 = sqrt(mu / (an - del_a0_gcom));
v1 = sqrt(mu / (an + del_a0_gcom));
Delta_v_gcom = abs(v1 - v0);
fprintf('[c] Total impulse for GCOM-W1: %.6f m/s\n', Delta_v_gcom * 1e3)

%% Time vectors
t_gcom = linspace(0, t_m_gcom, 1000);

%% Drift computation

del_lambda_gcom = -1.5 * omega_E * del_a0_gcom / an * t_gcom + ...
    0.75 * omega_E * rho_gcom * Cd_gcom * S_m_gcom * sqrt(mu / an) * t_gcom.^2;
del_tau_gcom = -1.5 * del_a0_gcom / an * t_gcom + ...
    0.75 * rho_gcom * Cd_gcom * S_m_gcom * sqrt(mu / an) * t_gcom.^2;

tau_diff = 85.5;
sep = tau_diff + del_tau_aqua/2 + del_tau_gcom/2;

%% Propagate altitude
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);
%aqua_struct = ode45(@(t,Y) SMArate(t,Y,an,mu,S_m_aqua,Cd_gcom,R_E,omega_E), t_aqua, an + del_a0_aqua, opt);
a_gcom_struct = ode45(@(t,Y) SMArate(t,Y,an,mu,S_m_gcom,Cd_gcom,R_E,omega_E), t_gcom, an + del_a0_gcom, opt);
%aqua_a = deval(aqua_struct, t_aqua);
a_gcom = deval(a_gcom_struct, t_gcom);

%% Plot: Cross-track vs altitude
figure('Position',[100 100 1000 400])
subplot(1,2,1)
hold on
plot(rad2deg(del_lambda_aqua), a_aqua - an, 'b', 'LineWidth', 1.5)
plot(rad2deg(del_lambda_gcom), a_gcom - an, 'r', 'LineWidth', 1.5)
grid on; grid minor
xlabel('$\Delta \lambda$ [deg]', 'Interpreter','latex')
ylabel('$\delta a$ [km]', 'Interpreter','latex')
title('Cross-Track vs Altitude', 'Interpreter','latex')
legend('Aqua','GCOM-W1','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

%% Plot: Along-track vs altitude
subplot(1,2,2)
hold on
plot(del_tau_aqua-sep/2, a_aqua - an, 'b', 'LineWidth', 1.5)
plot(del_tau_gcom-sep/2, a_gcom - an, 'r', 'LineWidth', 1.5)
grid on; grid minor
xlabel('$\Delta \tau$ [s]', 'Interpreter','latex')
ylabel('$\delta a$ [km]', 'Interpreter','latex')
title('Along-Track vs Altitude', 'Interpreter','latex')
legend('Aqua','GCOM-W1','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

%% Relative along-track plot (Part e)
t_common = linspace(0, min(t_m_aqua, t_m_gcom), 1000);
tau_gcom_shared = -1.5 * del_a0_gcom / an * t_common + ...
    0.75 * rho_gcom * Cd_gcom * S_m_gcom * sqrt(mu / an) * t_common.^2;
tau_aqua_shared = -1.5 * del_a0_aqua / an * t_common + ...
    0.75 * rho_aqua * Cd_gcom * S_m_aqua * sqrt(mu / an) * t_common.^2;
dtau_rel = tau_gcom_shared - tau_aqua_shared;

figure
plot(t_common/3600, dtau_rel, 'k', 'LineWidth', 1.5)
grid on; grid minor
xlabel('Time [hours]', 'Interpreter','latex')
ylabel('$\Delta \tau_{\text{rel}}$ [km]', 'Interpreter','latex')
title('Relative Along-Track: GCOM-W1 w.r.t. Aqua', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')


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


