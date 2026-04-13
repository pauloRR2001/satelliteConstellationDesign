# Satellite Constellation Design (MATLAB)

This repository contains MATLAB scripts and functions used for designing and analysing satellite constellations. It includes an integer-only genetic algorithm implementation, several constellation example scripts, orbital converters/propagators, Lambert solvers and plotting utilities.

This README is based on the files present in the repository and does not assume additional scripts or toolboxes beyond what is included here.

---

## Quick overview

Key folders and files (present in this repository):

- geneticAlgorithm/
  - runIntegerGA.m               % Integer-only GA driver
  - selectTournament.m           % Tournament selection
  - crossoverOnePoint.m          % One-point crossover
  - mutateRandomInt.m            % Integer mutation operator
  - evaluatePopulation.m         % Batch fitness evaluation
  - ga_demo.m                    % Small demo (see geneticAlgorithm/README.md)
  - README.md                    % Notes for the GA utilities

- satCons/
  - flowerGTRepeatDesign.m       % Flower ground-track / plotting example
  - iridiumConstellationBuild.m  % Iridium-like constellation build & plots
  - optimalFlowerGSCoverage.m    % (constellation coverage example)
  - stealthConstellationGA.m     % GA example for stealth constellation (script)
  - geoDolphinStealthGA.m       % Another GA-based example
  - adversialConstellationOptimization.m
  - bruteMinInterSatDis.m
  - countUniqueFlowerCombinations.m
  - satDragMaintenance.m
  - twoSatDragMaintenance.m
  - twoSatTumblingDrag.m
  - uniqueNumberOfWalkerDelta.m

- functions/
  - converters/                   % kep2cart, cart2kep, ECI2ECEF, etc.
  - propagators/                  % keplerianPropUnperturbed.m
  - lambert solvers (functions/lambertSolver*.m)
  - plotters/                     % plotOrbit3, plotGroundTrackLag, earthy (planet plotting)
  - genTools/                     % helper utilities (atmDensity, rotateXYZrth, nonuniqueAngle)

- tests/
  - fence_ga_test.m              % a test script referencing GA utilities

Use these scripts as starting points for experiments and to see example usage patterns.

---

## Getting started (example)

1. Open MATLAB and add the repository to the path. From the repository root in MATLAB:

   ```matlab
   addpath(genpath(pwd))
   ```

2. Run the GA demo (example):

   ```matlab
   cd('geneticAlgorithm')
   ga_demo
   ```

   The GA utilities are implemented in `runIntegerGA.m` and companion files. See `geneticAlgorithm/README.md` for details and a short usage example showing how to call `runIntegerGA` with an objective function and an initial integer population.

3. Run an example constellation script (examples are in `satCons/`):

   - flower ground-track / constellation plotting:
     ```matlab
     cd('satCons')
     flowerGTRepeatDesign
     ```

   - Iridium-like constellation build and plots:
     ```matlab
     cd('satCons')
     iridiumConstellationBuild
     ```

These scripts are full MATLAB scripts (not packaged functions). Inspect them to adapt parameters for your experiments.

---

## How the GA is organized

The GA implementation is intentionally minimal and operates on integer chromosomes. The main entry point is `runIntegerGA.m` which expects:

- objectiveFcn: function handle f(x) returning a scalar fitness to minimize
- initialPopulation: integer matrix (rows are individuals)
- maxGenerations: number of generations to run
- options: struct with fields used by the GA (mutationRate, crossoverRate, eliteCount, tournamentSize, mutationStep, lowerBounds, upperBounds)

See `geneticAlgorithm/README.md` for a short example of building a population and calling `runIntegerGA`.

---

## Notes about the codebase

- Many example scripts in `satCons/` are structured as scripts (they run when executed) and produce plots. Inspect and adapt them directly if you want to run parts of the analysis or reuse functions.
- Utility functions for coordinate conversions, orbital propagation and plotting are located under `functions/`. Useful files include `functions/converters/kep2cart.m`, `functions/propagators/keplerianPropUnperturbed.m`, and plotting helpers in `functions/plotters/`.
- There is no assumption in this README about external toolboxes. Some scripts may use MATLAB toolboxes (for example plotting or parallelization) — check the top of each script for toolbox-specific calls.

---

## Running and adapting examples

- To reuse the GA, write a fitness function that accepts an integer row vector and returns a scalar fitness. Build an initial integer population (rows) and call `runIntegerGA` as shown in `geneticAlgorithm/README.md`.
- Example scripts in `satCons/` show common constellation design workflows (ground-track repeat calculation, plotting, simple coverage/visibility checks). They can be edited directly to change parameters and save results.

---

## Tests

A simple test script is available at `tests/fence_ga_test.m`. Use it as a reference for how tests or small experiments are structured in this repository.

---

## Contributing

To contribute updates or fixes:
- Open an issue describing the change or bug.
- Fork the repository and create a branch with your changes.
- Submit a pull request with a clear description of the change.

---

## License

No license file is present in the repository root. If you want to apply a license, add a LICENSE file (for example MIT) and update this README accordingly.

---

If you want the README to include additional runnable examples or to expose a specific script as the canonical entry point, tell me which file(s) you want promoted and I will update the README to show exact commands and example input/outputs.
