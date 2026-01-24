clc; clear; close all

% Test: Integer GA to maximize rectangular fence area with W+L = 42
% Chromosome: [W, L] (integers)
% Fitness (minimize): f = -(W*L) + penalty*abs(W+L-42)

% Add GA path
addpath('..\geneticAlgorithm');

S = 42;                   % total blocks for W + L
Npop = 60;                % population size
maxGens = 80;             % generations

% Initial population: integers in [1, S-1]
initPop = randi([1, S-1], Npop, 2);

% Objective function (penalty-based)
penaltyWeight = 1000;     % strong weight to enforce W+L=42
obj = @(x) (-(x(1)*x(2)) + penaltyWeight*abs(x(1)+x(2)-S));

% GA options
opts.mutationRate   = 0.2;
opts.crossoverRate  = 0.9;
opts.eliteCount     = 2;
opts.tournamentSize = 3;
opts.mutationStep   = 1;
opts.lowerBounds    = [1, 1];
opts.upperBounds    = [S-1, S-1];

% Run GA
[bestChrom, bestFitness, hist] = runIntegerGA(obj, initPop, maxGens, opts);
W = bestChrom(1);
L = bestChrom(2);
area = W * L;

fprintf('Best W = %d, L = %d (W+L = %d), Area = %d\n', W, L, W+L, area);

% Compare to theoretical optimum under W+L=S (unconstrained integers): W=L=S/2
Wopt = floor(S/2);
Lopt = S - Wopt;
areaOpt = Wopt * Lopt;
fprintf('Theoretical (integer) best: W=%d, L=%d, Area=%d\n', Wopt, Lopt, areaOpt);

% Plot fitness trace
figure('Name','Fence GA Fitness per Generation');
plot(0:maxGens, hist.fitnessTrace, '-o'); grid on;
xlabel('Generation'); ylabel('Best Fitness (lower is better)');

% Plot the fence in 2D
figure('Name','Rectangular Fence (Best Solution)'); hold on; axis equal; grid on;
% Draw rectangle corners
plot([0 W W 0 0], [0 0 L L 0], 'k-', 'LineWidth', 2);
% Mark blocks along edges (conceptual placement)
% Using unit spacing markers along width and length edges
for i = 0:W
    plot(i, 0, 'bs', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
    plot(i, L, 'bs', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
end
for j = 0:L
    plot(0, j, 'rs', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    plot(W, j, 'rs', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
end

xlabel('Width (blocks)'); ylabel('Length (blocks)');
title(sprintf('Fence: W=%d, L=%d, Area=%d (Constraint W+L=%d)', W, L, area, S));
legend({'Fence outline','Width blocks','Length blocks'}, 'Location','bestoutside');
