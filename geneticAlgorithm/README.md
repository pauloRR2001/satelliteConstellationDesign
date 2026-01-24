# Integer-Only Genetic Algorithm (MATLAB)

Minimal functions to evolve integer chromosomes for constellation design or general integer optimization.

## Files
- `runIntegerGA.m`: GA driver (selection, crossover, mutation, elitism).
- `selectTournament.m`: Tournament selection (minimization).
- `crossoverOnePoint.m`: One-point crossover.
- `mutateRandomInt.m`: Random per-gene ±step mutation with optional bounds.
- `evaluatePopulation.m`: Batch fitness evaluation and sorting.
- `ga_demo.m`: Tiny demo minimizing distance to a target integer vector.

## Usage
```matlab
cd('c:/Users/paulo/MATLAB Drive/2BP/geneticAlgorithm');
% Define objective
obj = @(x) myFitness(x);       % x is integer row vector
% Build initial integer population
initPop = randi([0 10], 50, 6);   % 50 individuals, 6 genes
% Options
opts.mutationRate   = 0.15;
opts.crossoverRate  = 0.8;
opts.eliteCount     = 2;
opts.tournamentSize = 3;
opts.mutationStep   = 1;
opts.lowerBounds    = zeros(1,6);
opts.upperBounds    = 10*ones(1,6);
% Run
[bestChrom, bestFitness, hist] = runIntegerGA(obj, initPop, 100, opts);
```

Run the demo:
```matlab
cd('c:/Users/paulo/MATLAB Drive/2BP/geneticAlgorithm');
ga_demo
```
