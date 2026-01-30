# Laminate_Solver

Laminate_Solver is a standalone, CLI-based engineering tool for the analysis of composite laminates using **Classical Laminate Theory (CLT)**.

Designed for ease of distribution, the entire engine is contained in a **single Python file** (`Master_Laminate.py`), making it perfect for students, researchers, and engineers who need quick, reliable stiffness and deformation analysis without managing complex dependencies. The solver adheres to the sign conventions and integration methods found in textbooks by *R.M. Jones* and *J.N. Reddy*.

## Key Features

### 1. Advanced Constitutive Analysis
- **Full ABD Matrix Calculation:** Computes Extensional [A], Coupling [B], and Bending [D] stiffness matrices.
- **Rigidly Aligned Reporting:** Generates professional text reports with perfectly aligned matrix equations ($\{N\} = [ABD] \times \{\epsilon\}$), regardless of number magnitude.
- **Coupling Detection:** Automatically flags physical phenomena:
  - *Extension-Shear Coupling* ($A_{16}, A_{26} \neq 0$)
  - *Bending-Extension Coupling* ($B_{ij} \neq 0$ - Warping Risk)
  - *Bend-Twist Coupling* ($D_{16}, D_{26} \neq 0$)
- **Z-Axis Toggle:** Switch between **Z-Down** (Aerospace/Jones) and **Z-Up** (Physics/Reddy) conventions instantly.

### 2. Flexible Material Definition
- **Multiple Input Modes:**
  - **Engineering Constants:** ($E_1, E_2, \nu_{12}, G_{12}$)
  - **Direct Q-Matrix:** Input plane-stress stiffness directly.
  - **3D Stiffness (C-Matrix):** Inputs full $6 \times 6$ constitutive matrix; the tool automatically performs **Plane Stress Reduction**.
- **Pre-Rotated Matrix Support:** Specifically handles scenarios where input matrices are already transformed (e.g., $Q_{45}$), ensuring the math remains accurate while preserving the correct stacking notation (e.g., `[45/-45]`).

### 3. Solver & Optimization
- **Deformation Solver:** Calculates global strains $\{\epsilon^0\}$ and curvatures $\{\kappa\}$ for given Load/Moment vectors $\{N\}$ and $\{M\}$.
- **Monte Carlo Optimization:** Iteratively shuffles stacking sequences to minimize specific coupling terms (e.g., minimizing $D_{16}/D_{11}$ ratio) while maintaining symmetry.

## Installation

1. Clone the repository:
   ```bash
   git clone [https://github.com/naftalk/Laminate_Solver.git](https://github.com/naftalk/Laminate_Solver.git)
   
2. Install the required dependency:

    ```pip install numpy```

## Usage

Simply run the python script

```python Master_Laminate.py```

## Workflow
   
- Material Wizard: Define your materials.

- Stack Builder: Add layers, define angles and thickness.

- Analysis: Generate a Stiffness Report or run the Deformation Solver.

- Export: Reports are automatically saved as Report_Stiff.txt and Report_Deform.txt.

## Verification

This tool has been verified against standard textbook problems. Please see ```Examples.md``` for step-by-step verification cases, including:

    Symmetric Angle-Ply: Proving zero coupling for pre-rotated matrices.

    Unsymmetric Hybrid: Verifying B-matrix coupling logic.

    Standard Carbon: Verifying engineering constant transformation.

## Theory & Conventions

The tool implements the standard integration of transformed reduced stiffness components Qˉ​ij​:

$$ A_{ij} = \sum_{k=1}^{N} (\bar{Q}{ij})k (z_k - z{k-1}) $$ 
$$ B{ij} = \frac{1}{2} \sum_{k=1}^{N} (\bar{Q}{ij})k (z_k^2 - z{k-1}^2) $$ 
$$ D{ij} = \frac{1}{3} \sum_{k=1}^{N} (\bar{Q}_{ij})k (z_k^3 - z{k-1}^3) $$

## License

Distributed under the MIT License. See ```LICENSE``` for more information.
