# Integer Quantum Hall Effect (IQHE) Simulation

This repository contains Python code to simulate and analyze the Integer Quantum Hall Effect (IQHE) on a square lattice. The simulation involves constructing lattice models under a perpendicular magnetic field, generating Hamiltonian matrices, and analyzing the energy levels and degeneracies. This was used in the first part of my thesis: [Integer Quantum Hall Effect on Lattices with Cylindrical and Toroidal Geometries](https://scholarworks.calstate.edu/concern/theses/7d279044v).

## Repository Structure

### Files

1. **`lattice.py`**
   - Defines the `Lattice` and `DoppedLattice` classes to model square lattices with and without impurities.
   - Key features:
     - Define lattice dimensions (`Nx`, `Ny`).
     - Specify onsite potential (`eps`), hopping energy (`t`), and magnetic flux (`1/Q`).
     - Optionally, add impurities with a specified width (`w`) in `DoppedLattice`.
     - Calculate impurity potentials and total sites on the lattice.

2. **`hamiltonian.py`**
   - Defines the `Hamiltonian` class to construct the Hamiltonian matrix for the lattice under a magnetic field.
   - Key features:
     - Incorporates onsite energy, hopping energies in `x` and `y` directions, and periodic boundary conditions.
     - Applies the substitution to introduce magnetic flux.
     - Supports both clean and doped lattices.

3. **`analysis.py`**
   - Defines the `Analysis` class to compute and visualize energy levels of the Hamiltonian matrix.
   - Key features:
     - Calculate eigenvalues and their degeneracies.
     - Visualize energy levels and their degeneracies using histograms.
     - Analyze the effect of impurities on energy levels.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/tiffanyxchen/iqhe
   cd iqhe
   ```

2. Install dependencies:
   ```bash
   pip install numpy matplotlib
   ```

## Usage

### Example Workflow

1. **Define a Lattice:**
   ```python
   from lattice import DoppedLattice

   lattice = DoppedLattice(Nx=10, Ny=10, eps=0, t=1, Q=10, w=1)
   print(lattice)
   ```

2. **Generate a Hamiltonian:**
   ```python
   from hamiltonian import Hamiltonian

   hamiltonian = Hamiltonian(lattice)
   print(hamiltonian.getMatrix())
   ```

3. **Analyze and Visualize Results:**
   ```python
   from analysis import Analysis
   
   lattice       = Lattice(10, 10, 0, 2, 10)
   doped_lattice = DoppedLattice(10, 10, 0, 2, 10, 2)
   hamiltonian   = Hamiltonian(lattice)
   hamiltonian_d = Hamiltonian(doped_lattice)

   result = Analysis(hamiltonian.getMatrix(), hamiltonian_d.getMatrix())
   result.graphEnergyLevels(doped_lattice)
   ```

### Outputs
- Hamiltonian matrix.
- Eigenvalues of the Hamiltonian.
- Histogram showing energy levels and degeneracies.

## Example Visualization
The histogram visualizes the energy levels and their degeneracies for clean and doped lattices.

![Energy Levels Histogram](https://github.com/tiffanyxchen/iqhe/blob/main/10x10%20Square%20Lattice%3A%20eps%3D0%2C%20t%3D2%2C%20Q%3D10%2C%20w%3D2.png)

## Authors
- **Xing [Tiffany] Chen** (Author of the scripts)

## License

## Acknowledgments
- The implementation leverages Python libraries like `numpy` and `matplotlib` for numerical computations and visualizations.

---
For further inquiries or contributions, feel free to open an issue or submit a pull request!
