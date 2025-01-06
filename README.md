# Bubbles on String: Viscoelastic Fluid Dynamics Simulation

[![License](https://img.shields.io/github/license/comphy-lab/BubblesOnString?style=flat)](LICENSE)
![GitHub repo size](https://img.shields.io/github/repo-size/comphy-lab/BubblesOnString)
![GitHub last commit](https://img.shields.io/github/last-commit/comphy-lab/BubblesOnString)
[![GitHub issues](https://img.shields.io/github/issues/comphy-lab/BubblesOnString)](https://github.com/comphy-lab/BubblesOnString/issues)

This repository contains simulation code for studying the dynamics of bubbles that are formed when a jet hits a liquid pool. The code uses the Basilisk framework with custom implementations for viscoelastic fluid dynamics using the log-conformation method as described in [Viscoelastic3D repository](https://github.com/comphy-lab/Viscoelastic3D).

## Overview

The project investigates how viscoelastic behavior influences bubble dynamics that are formed when a jet hits a liquid pool.
- Two-phase flow interactions with viscoelastic fluids
- Log-conformation method for numerical stability at high Weissenberg numbers
- Axisymmetric implementations
- Surface tension and interface tracking using Volume of Fluid (VoF) method

## Key Features

- **Two-Phase Viscoelastic Solver**: Implementation of viscoelastic fluid dynamics with density/viscosity/elastic property interpolation between phases
- **Log-Conformation Method**: Advanced numerical approach for handling high Weissenberg number problems
- **Scalar-Based Implementation**: Both 2D (+axisymmetric) and 3D versions available
- **Surface Tension**: Brackbill method implementation for interface forces
- **Adaptive Mesh Refinement**: Efficient grid management for interface tracking

## Installation and Setup

### Prerequisites
- Basilisk C (automatically installed by setup script)
- C compiler (gcc/clang)
- Make build system
- MPI for parallel execution (optional)

### Basic Installation

1. Clone the repository:
```bash
git clone https://github.com/comphy-lab/BubblesOnString.git
cd BubblesOnString
```

2. Install Basilisk and set up the environment:
```bash
./reset_install_requirements.sh
```

### Environment Setup

The installation script will:
- Check for and install Basilisk if not present
- Configure environment variables
- Create a `.project_config` file with proper paths

For a clean reinstall:
```bash
./reset_install_requirements.sh --hard
```

## Project Structure

- `basilisk/src/`: Core Basilisk codebase with CFD solvers
  - `navier-stokes/`: Navier-Stokes equation solvers
  - `vof.h`: Volume of Fluid implementation
  - Other core modules for CFD

- `src-local/`: Custom implementations for this project
  - `two-phaseVE.h`: Two-phase viscoelastic flow solver
  - `log-conform-viscoelastic-scalar-2D.h`: 2D log-conformation solver
  - `log-conform-viscoelastic-scalar-3D.h`: 3D log-conformation solver
  - `eigen_decomposition.h`: Matrix operations for stress tensor calculations

- `testCases/`: Example simulations and validation cases
  - `jetOnPool.c`: Jet impact on pool surface simulation
  - Various test configurations

## Running Simulations

### Using Make (Recommended)

1. Navigate to test case directory:
```bash
cd testCases
```

2. Run with visualization:
```bash
CFLAGS=-DDISPLAY=-1 make jetOnPool.tst
```

3. Run without visualization:
```bash
make jetOnPool.tst
```

### Direct Compilation

```bash
qcc -O2 -Wall -disable-dimensions -I$(PWD)/src-local jetOnPool.c -o jetOnPool -lm
./jetOnPool
```

### Parallel Execution (MPI)

For cluster environments:

1. Compile with MPI support:
```bash
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions jetOnPool.c -o jetOnPool -lm
```

2. Run with MPI:
```bash
mpirun -np 4 ./jetOnPool
```

## Contributing

We welcome contributions! The project uses standardized issue templates to streamline communication and problem-solving. When contributing, please:

1. Fork the repository
2. Create a feature branch
3. Submit a Pull Request

### Issue Templates

We provide three types of issue templates to help you contribute effectively:

1. **Bug Reports** (`bug_report.md`):
   - Detailed reproduction steps
   - System information
   - Simulation parameters
   - Expected vs actual behavior
   - Relevant outputs and visualizations

2. **Feature Requests** (`feature_request.md`):
   - Problem description
   - Proposed solution
   - Implementation details
   - Impact assessment
   - Relevant literature

3. **Questions/Discussions** (`question.md`):
   - Clear topic statement
   - Context and background
   - Related resources
   - System information
   - Supporting materials

### External Resources

For additional help and discussions:
- [Basilisk Forum](http://groups.google.com/d/forum/basilisk-fr): For Basilisk-specific questions
- [Viscoelastic3D Repository](https://github.com/comphy-lab/Viscoelastic3D): For viscoelastic solver implementation
- [Physics of Fluids Group](https://pof.tnw.utwente.nl/): For academic collaborations

## Authors

- Vatsal Sanjay (University of Twente)
  - Email: vatsalsanjay@gmail.com
  - Physics of Fluids Group

## License

This project is licensed under standard academic terms. Please cite the associated papers if you use this code in your research.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{BubblesOnString,
  author = {Sanjay, V.},
  title = {BubblesOnString: Viscoelastic Fluid Dynamics Simulation},
  year = {2024},
  publisher = {GitHub},
  url = {https://github.com/comphy-lab/BubblesOnString}
}
```

