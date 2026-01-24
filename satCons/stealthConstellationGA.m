clc; clear; close all

% Stealth Flower Constellation via Integer GA
% Minimize Spy–Iridium unobstructed LOS over a short horizon, subject to
% integer parameters and basic uniform constellation constraints.

addpath('..\geneticAlgorithm');
addpath('..\functions');
addpath('..\functions\converters');
addpath('..\functions\plotters');

%% Constants
mu      = 398600.4418;      % km^3/s^2
R_E     = 6378.137;         % km
TE      = 86164.0905;       % sidereal day [s]
omega_E = 7.2921151467e-5;  % rad/s

%% Fixed Iridium Constellation
t_ir   = 66;                % total satellites
p_ir   = 6;                 % planes
f_ir   = 2;                 % phasing
inc_ir = deg2rad(86.4);
a_ir   = 7158;              % km
ecc_ir = 0; argp_ir = 0;
n_ir   = sqrt(mu/a_ir^3);
[raan_ir, M0_ir] = buildIridium(p_ir, t_ir, f_ir);

%% Spy catalog and GA chromosome
Ns_total   = 30;              % total spy satellites
alist      = [6900, 7158, 7300];      % discrete altitudes (km)
inclist    = deg2rad([49.4, 56, 70]); % discrete inclinations

% GA chromosome: [No, Nso, Nc, a_idx, inc_idx]
% Bounds
lb = [1, 1, 1, 1, 1];
ub = [Ns_total, Ns_total, Ns_total, numel(alist), numel(inclist)];

%% Evaluation horizon for GA (keep small for speed)
evalDays = 3;                % 3 sidereal days
T_eval   = evalDays * TE;
dt       = 3600;             % 1-hour sampling
% Precompute Iridium positions over t_vec for reuse in fitness

% Time vector
% Note: Keep this global to reuse inside fitness via nested function
% to avoid recomputing Iridium
% (MATLAB nested functions inherit workspace variables)
t_vec = 0:dt:T_eval;
Nt    = numel(t_vec);

posIr = zeros(Nt, t_ir, 3);
for ti = 1:Nt
    t = t_vec(ti);
    Mk = mod(M0_ir + n_ir * t, 2*pi);
    for k = 1:t_ir
        [x, y, z, ~, ~, ~] = kep2cart(a_ir, ecc_ir, inc_ir, raan_ir(k), argp_ir, Mk(k), mu);
        posIr(ti, k, :) = [x, y, z];
    end
end

%% Fitness function
% Minimize f = coverage_seconds + penalties
pen_NoNso    = 1e5;  % penalty weight for No*Nso != Ns_total
pen_Nc       = 1e4;  % penalty for invalid Nc
pen_bounds   = 1e4;  % penalty for out-of-range indices
pen_sep      = 1e3;  % penalty weight per km below d_min
minSepKm     = 50;   % minimum Spy–Spy separation threshold [km]

objectiveFcn = @(chrom) fitnessSpy(chrom, Ns_total, alist, inclist, t_vec, Nt, dt, posIr, mu, R_E, pen_bounds, pen_NoNso, pen_Nc, pen_sep, minSepKm);

%% Initial population
Npop   = 40;
initPop = zeros(Npop, 5);
initPop(:,1) = randi([lb(1), ub(1)], Npop, 1);   % No
initPop(:,2) = randi([lb(2), ub(2)], Npop, 1);   % Nso
initPop(:,3) = randi([lb(3), ub(3)], Npop, 1);   % Nc
initPop(:,4) = randi([lb(4), ub(4)], Npop, 1);   % a_idx
initPop(:,5) = randi([lb(5), ub(5)], Npop, 1);   % inc_idx

%% GA options
opts.mutationRate   = 0.2;
opts.crossoverRate  = 0.85;
opts.eliteCount     = 2;
opts.tournamentSize = 3;
opts.mutationStep   = 1;
opts.lowerBounds    = lb;
opts.upperBounds    = ub;

maxGens = 60;
[bestChrom, bestFitness, hist] = runIntegerGA(objectiveFcn, initPop, maxGens, opts);

% Decode best
No_best   = bestChrom(1);
Nso_best  = bestChrom(2);
Nc_best   = bestChrom(3);
a_best    = alist(bestChrom(4));
inc_best  = inclist(bestChrom(5));

fprintf('Best Spy chromosome: No=%d, Nso=%d, Nc=%d, a=%.1f km, inc=%.2f deg\n', ...
    No_best, Nso_best, Nc_best, a_best, rad2deg(inc_best));
fprintf('Best fitness: %.4f\n', bestFitness);

%% Quick plots for best config
[raan_spy_best, M0_spy_best] = buildFlower(No_best, Nso_best, Nc_best, 0);
figure('Name','Best Spy RAAN–M');
scatter(rad2deg(mod(raan_spy_best,2*pi)), rad2deg(mod(M0_spy_best,2*pi)), 40, 'filled');
xlabel('RAAN [deg]'); ylabel('Mean Anomaly [deg]'); grid on; axis([0 360 0 360]);

ecc_spy = 0; argp_spy = 0; n_spy = sqrt(mu/a_best^3);
nu_vec = linspace(0, 2*pi, 150);
figure('Name','Best Spy vs Iridium (ECI)'); hold on
earthy(R_E, 'Earth', 0.9, [0; 0; 0]);
for s = 1:No_best*Nso_best
    XYZ = zeros(numel(nu_vec),3);
    for m = 1:numel(nu_vec)
        [x,y,z,~,~,~] = kep2cart(a_best, ecc_spy, inc_best, raan_spy_best(s), argp_spy, nu_vec(m), mu);
        XYZ(m,:) = [x,y,z];
    end
    plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'r-');
end
for k = 1:t_ir
    XYZ = zeros(numel(nu_vec),3);
    for m = 1:numel(nu_vec)
        [x,y,z,~,~,~] = kep2cart(a_ir, ecc_ir, inc_ir, raan_ir(k), argp_ir, nu_vec(m), mu);
        XYZ(m,:) = [x,y,z];
    end
    plot3(XYZ(:,1),XYZ(:,2),XYZ(:,3),'b-');
end

% Dots at satellite positions at t = 0 (ECI)
I0 = zeros(t_ir, 3);
for k = 1:t_ir
    [x, y, z, ~, ~, ~] = kep2cart(a_ir, ecc_ir, inc_ir, raan_ir(k), argp_ir, M0_ir(k), mu);
    I0(k, :) = [x, y, z];
end
scatter3(I0(:,1), I0(:,2), I0(:,3), 20, 'filled', 'MarkerFaceColor', [0, 0, 1], 'MarkerEdgeColor', 'none')

S0 = zeros(No_best*Nso_best, 3);
for s = 1:No_best*Nso_best
    [x, y, z, ~, ~, ~] = kep2cart(a_best, ecc_spy, inc_best, raan_spy_best(s), argp_spy, M0_spy_best(s), mu);
    S0(s, :) = [x, y, z];
end
scatter3(S0(:,1), S0(:,2), S0(:,3), 20, 'filled', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', 'none')
axis equal; grid on; xlabel('X'); ylabel('Y'); zlabel('Z');
legend({'Spy orbits','Iridium orbits','Iridium sats','Spy sats'}, 'Location', 'bestoutside');

figure('Name','GA Fitness per Generation');
plot(0:maxGens, hist.fitnessTrace, '-o'); grid on;
xlabel('Generation'); ylabel('Best Fitness (lower is better)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested/local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = fitnessSpy(chrom, Ns_total, alist, inclist, t_vec, Nt, dt, posIr, mu, R_E, pen_bounds, pen_NoNso, pen_Nc, pen_sep, minSepKm)
    % Decode
    No   = chrom(1); Nso = chrom(2); Nc = chrom(3);
    a_i  = chrom(4); inc_i = chrom(5);

    % Penalties for bounds and constraints
    f = 0;
    if (No < 1) || (Nso < 1) || (Nc < 1) || (a_i < 1) || (inc_i < 1)
        f = f + pen_bounds;
    end
    if (No > Ns_total) || (Nso > Ns_total) || (Nc > Ns_total) || (a_i > numel(alist)) || (inc_i > numel(inclist))
        f = f + pen_bounds;
    end

    % Divisibility / size constraint
    if No*Nso ~= Ns_total
        f = f + pen_NoNso * abs(No*Nso - Ns_total);
    end
    % Nc range depends on No
    if (Nc < 1) || (Nc >= No)
        f = f + pen_Nc * abs(Nc - max(1,min(No-1,Nc)));
    end

    % If obviously invalid, short-circuit
    if f > pen_NoNso
        f = f + 1e6; return; 
    end

    % Build Spy constellation
    a_spy   = alist(a_i);
    inc_spy = inclist(inc_i);
    [raan_spy, M0_spy] = buildFlower(No, Nso, Nc, 0);

    % Coverage metric over t_vec
    coverage_seconds = 0;
    minSep = inf;
    for ti = 1:Nt
        t = t_vec(ti);
        Mk_spy = mod(M0_spy + sqrt(mu/a_spy^3)*t, 2*pi);
        Spos = zeros(Ns_total,3);
        for s = 1:Ns_total
            [xs,ys,zs,~,~,~] = kep2cart(a_spy, 0, inc_spy, raan_spy(s), 0, Mk_spy(s), mu);
            Spos(s,:) = [xs,ys,zs];
        end
        Ipos = squeeze(posIr(ti,:,:));
        if anyPairVisible(Spos, Ipos, R_E)
            coverage_seconds = coverage_seconds + dt;
        end
        % Spy–Spy min separation (coarse)
        for i = 1:Ns_total-1
            for j = i+1:Ns_total
                d = norm(Spos(i,:) - Spos(j,:));
                if d < minSep
                    minSep = d;
                end
            end
        end
    end
    % Penalty for min separation below threshold
    if minSep < minSepKm
        f = f + pen_sep * (minSepKm - minSep);
    end
    % Objective
    f = f + coverage_seconds;
end

function [raan, M] = buildIridium(p, t, f)
    raan = zeros(t,1);
    M    = zeros(t,1);
    idx = 1;
    for i = 1:p
        for j = 1:t/p
            raan(idx) = mod((2*pi*(i-1)/p)/2, 2*pi);
            M(idx)    = mod(2*pi*(p/t)*(j-1) + 2*pi*(f/t)*(i-1), 2*pi);
            idx = idx + 1;
        end
    end
end

function [raan, M] = buildFlower(No, Nso, Nc, raan0)
    N = No*Nso;
    raan = zeros(N,1);
    M    = zeros(N,1);
    idx = 1;
    for i = 1:No
        for j = 1:Nso
            raan(idx) = mod(raan0 + 2*pi*(i-1)/No, 2*pi);
            M(idx)    = 2*pi*((j-1)/Nso - (Nc/Nso)*(i-1)/No);
            idx = idx + 1;
        end
    end
    M = mod(M, 2*pi);
end

function vis = anyPairVisible(Spos, Ipos, R_E)
    Ns = size(Spos,1); Ni = size(Ipos,1);
    vis = false;
    for s = 1:Ns
        r1 = Spos(s,:).';
        for i = 1:Ni
            r2 = Ipos(i,:).';
            if unobstructedLOS(r1, r2, R_E)
                vis = true; return
            end
        end
    end
end

function ok = unobstructedLOS(r1, r2, R_E)
    d  = r2 - r1; d2 = dot(d,d);
    if d2 == 0, ok = (norm(r1) > R_E); return; end
    s = -dot(r1,d)/d2;
    if s < 0
        p = r1;
    elseif s > 1
        p = r2;
    else
        p = r1 + s*d;
    end
    ok = (norm(p) > R_E);
end
