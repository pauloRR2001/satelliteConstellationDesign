# Satellite Constellation Design Toolkit

Design and evaluate satellite constellations using Genetic Algorithms (GA) and exhaustive brute-force sweeps. This MATLAB-based toolkit helps explore constellation parameter spaces (planes, satellites per plane, inclination, phasing, RAAN spacing) and rank candidate designs by mission-specific metrics such as coverage, revisit time and cost.

---

## Features

- Genetic Algorithm optimizer for large search spaces and multi-objective trade-offs.
- Brute-force sweep mode for exhaustive evaluation of small discrete domains.
- Modular evaluation pipeline so you can plug in custom metrics and constraints.
- CSV/JSON export of results and basic visualization utilities.
- Designed to work with plain MATLAB; optionally compatible with MATLAB Global Optimization Toolbox.

## Requirements

- MATLAB (recent release recommended).
- Optional: MATLAB Global Optimization Toolbox (only required for toolbox-specific GA functions). The repository contains pure-MATLAB implementations where possible.

## Quick start

1. Clone the repository:

   ```bash
   git clone https://github.com/pauloRR2001/satelliteConstellationDesign.git
   ```

2. Open MATLAB and add the repository to your path:

   ```matlab
   addpath(genpath('path/to/satelliteConstellationDesign'))
   ```

3. Configure an experiment

   - Edit the configuration script or file in `config/` or `examples/` to set GA parameters (population size, generations, mutation rate) or sweep ranges (parameter bounds and resolution).
   - Specify evaluation metrics and output directory.

4. Run an experiment

   - Run the genetic algorithm entry script (look for files named like `main_ga.m`, `run_ga.m` or similar in the repository). Example:

     ```matlab
     run_ga(config)
     ```

   - Run the brute-force sweep entry script (look for `main_sweep.m`, `run_sweep.m` or similar). Example:

     ```matlab
     run_sweep(config)
     ```

Note: Replace the example function names above with the actual entry points present in this repository. See the `examples/` folder for concrete scripts.

## Configuration options (typical)

- Mission geometry: target region (lat/lon bounds), minimum elevation angle, ground station locations.
- Constellation parameters: numberOfPlanes, satsPerPlane, inclination, RAAN spacing, phasing.
- GA settings: populationSize, generations, crossoverRate, mutationRate, elitism.
- Sweep settings: parameter ranges and discretization steps.
- Objective weighting: single objective or weighted combination of coverage, revisit time, cost.

## Output

- Ranked candidate constellations with evaluation metrics.
- CSV/JSON summary files suitable for offline analysis.
- Visualization figures (ground tracks, coverage heatmaps) if enabled.
- Convergence diagnostics for the GA run.

## Project layout (typical)

- `/src` or `/matlab` — core functions and modules
- `/examples` — example configs and run scripts
- `/config` — default configuration templates
- `/data` — input files (ground stations, target region polygons)
- `/results` — output produced by runs
- `/docs` — additional documentation and references

Adjust locations above to match the layout in this repository.

## Tips

- Evaluations are independent and typically expensive — use MATLAB's Parallel Computing Toolbox (`parfor`, `parpool`) to speed up sweeps and fitness evaluations.
- Use brute-force sweep for small discrete searches or to validate GA results.
- Keep evaluation functions simple and vectorized where possible for performance.

## Contributing

- Open an issue to report bugs or request features.
- Fork the repo, create a feature branch, add tests/examples if applicable, and submit a pull request.

## License

This project is provided under the MIT License. Update the LICENSE file if you prefer a different license.

## Citation / References

If you use this code in research, please cite the repository and list any relevant academic references for the algorithms you used.

---
