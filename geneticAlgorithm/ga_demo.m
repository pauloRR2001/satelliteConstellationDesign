clc; clear; close all

% Demo: Integer GA minimizing distance to a target integer vector
% Target integer vector
TARGET = [3 7 2 5 9 1];

% Objective: sum of squared differences
obj = @(x) sum((x - TARGET).^2);

% Initial population: random integers in [0, 12]
Npop  = 40;
Nvars = numel(TARGET);
lb    = zeros(1,Nvars);
ub    = 12*ones(1,Nvars);
initPop = randi([0 12], Npop, Nvars);

% GA options
opts.mutationRate   = 0.15;
opts.crossoverRate  = 0.85;
opts.eliteCount     = 2;
opts.tournamentSize = 3;
opts.mutationStep   = 1;
opts.lowerBounds    = lb;
opts.upperBounds    = ub;

maxGens = 60;

[bestChrom, bestFitness, hist] = runIntegerGA(obj, initPop, maxGens, opts);

fprintf('Best chromosome: [%s]\n', num2str(bestChrom));
fprintf('Best fitness: %.4f\n', bestFitness);

figure('Name','GA Fitness per Generation')
plot(0:maxGens, hist.fitnessTrace, '-o')
xlabel('Generation')
ylabel('Best Fitness in Generation')
grid on
