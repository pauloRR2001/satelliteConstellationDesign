clc; clear; close all

%% Constants
mu = 398600.4418;              % km^3/s^2
R_E = 6378.137;                 % km
omega_E = 7.2921151467e-5;      % rad/s

Cd = 2.2;
S = 38.4 * (1/1000^2);                      % km^2
m = 2934;                      % kg
S_m = S/m;              % km^2/kg

%% Orbit Parameters
h = 705;                       % km
an = R_E + h;                  % nominal SMA
Delta_d = 20;                  % max cross-track deviation in km
Delta_lambda = Delta_d / R_E; % convert arc to angle [rad]

%% part a

t_mcross = 4*sqrt(Delta_lambda/(3*omega_E*atmDensity(an-R_E)*Cd*S_m)*sqrt(an/mu))

%% part b
Delta_a = 4*sqrt(Delta_lambda/(3*omega_E)*atmDensity(an-R_E)*Cd*S_m*sqrt(mu*an^3))

del_a0 = Delta_a/2;
%% part c
vc0 = sqrt(mu/(an-del_a0));
vc1 = sqrt(mu/(an+del_a0));

DV = abs(vc1 - vc0)

%% part d

a0 = an + del_a0;

t = linspace(0, t_mcross, 10000);
opt = odeset('RelTol',1e-12,'AbsTol',1e-12);

a_struct = ode45(@(t,Y) SMArate(t,Y,an,mu,S_m,Cd,R_E,omega_E),t,a0,opt);
a = deval(a_struct, t);

del_lambda = -3/2*omega_E*del_a0/an*t+...
    3/4*omega_E*atmDensity(an-R_E)*Cd*S_m*sqrt(mu/an)*t.^2;

del_tau = -3/2*del_a0/an*t + 3/4*atmDensity(an-R_E)*S_m*Cd*sqrt(mu/an)*t.^2;

del_tau_max = min(del_tau)

%% plott
figure('Position',[100 100 1000 400])

% Define y-values (delta a)
del_a = a - an;

% Determine dynamic x-limits
lambda_deg = rad2deg(del_lambda);
lambda_min = min(lambda_deg);
lambda_max = max(lambda_deg);
dx_lambda = abs(lambda_min)/6;
lambda_xlim = [lambda_min-dx_lambda, dx_lambda];

tau_min = min(del_tau);
tau_max = max(del_tau);
dx_tau = abs(tau_min)/6;
tau_xlim = [tau_min-dx_tau, dx_tau];

% Left Plot: Cross-track
subplot(1,2,1)
hold on
plot(lambda_deg, del_a, 'b', 'LineWidth', 1.5)
plot([lambda_deg(1) lambda_deg(end)], ...
     [del_a(1) del_a(end)], 'b', 'LineWidth', 1.5)
grid on
grid minor
xlabel('$\Delta \lambda$ [deg]', 'Interpreter', 'latex')
ylabel('$\delta a = a - a_n$ [km]', 'Interpreter', 'latex')
title('Cross-track vs Semi-major Axis Variations', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex')
xlim(lambda_xlim)

% Right Plot: Along-track
subplot(1,2,2)
hold on
plot(del_tau, del_a, 'r', 'LineWidth', 1.5)
plot([del_tau(1) del_tau(end)], ...
     [del_a(1) del_a(end)], 'r', 'LineWidth', 1.5)
grid on
grid minor
xlabel('$\Delta \tau$ [s]', 'Interpreter', 'latex')
ylabel('$\delta a = a - a_n$ [km]', 'Interpreter', 'latex')
title('Along-track (Lag) vs Semi-major Axis Variations', 'Interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex')
xlim(tau_xlim)
sgtitle('Variation from Nominal Orbit','interpreter','latex')



%% functions

function adot = SMArate(t,a,an,mu,S_m,Cd,R_E,omega_E)
    
    acc_D = simpleDrag(a, mu, S_m, Cd, R_E, omega_E);
    adot = 2*sqrt((a^3)/mu)*acc_D;
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


