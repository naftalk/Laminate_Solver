Verification Examples

Use these examples to verify the accuracy of the Laminate_Solver tool.

Example 1: Symmetric Angle-Ply (Pre-Rotated Matrices)
Configuration: [45 / -45 / -45 / 45]
Total Thickness: 12.0 mm (4 layers * 3 mm each)

1. Define Materials
You need to define two separate materials because the matrices are already rotated.

Material A: MAT_POS (+45)
- Mode: B (Plane Stress Q Matrix)
- Input:
[
  [6.55  5.15  4.50]
  [5.15  6.55  4.50]
  [4.50  4.50  5.15]
]
Material B: MAT_NEG (-45)
- Mode: B (Plane Stress Q Matrix)
- Input:
[
 [6.55  5.15 -4.50]
 [5.15  6.55 -4.50]
 [-4.50 -4.50  5.15]
]

2. Build Stack
Important: Since these matrices already contain the 45° transformation, we tell the code they are Pre-Rotated.
-Add Layer 1: Material MAT_POS, Angle 45, Thick 3.
-Prompt: "Is this Q-matrix PRE-ROTATED?" → YES (y)
-Add Layer 2: Material MAT_NEG, Angle -45, Thick 3.
-Prompt: "Is this Q-matrix PRE-ROTATED?" → YES (y)
-Make symmetry to have a final stacking of [45/-45/-45/45]

3. Expected Results
  [B] Matrix: All Zeros (Symmetric Laminate).
  [A] Matrix: The shear coupling terms (A16​,A26​) must be 0.

Example 2
Unsymmetric Hybrid [0 / 45]
Configuration: [0 (5mm) / 45 (3mm)] Total Thickness: 8.0 mm
1. Define Material
-Material: MAT_BASE
-Mode: B (Plane Stress Q Matrix)
-Input: (This is the 0-degree property)
[
[20.0  0.7  0.0]
[0.7   2.0  0.0]
[0.0   0.0  0.7]
]
2. Build Stack
-Add Layer 1 (Bottom): Material MAT_BASE, Angle 0, Thick 5.
-Prompt: "Is this Q-matrix PRE-ROTATED?" → NO (n)
                            (This uses the matrix exactly as typed)
-Add Layer 2 (Top): Material MAT_BASE, Angle 45, Thick 3.
-Prompt: "Is this Q-matrix PRE-ROTATED?" → NO (n)
                            (The code will automatically rotate the matrix by 45°)
3. Expected Results
-[B] Matrix: Non-Zero.
-The laminate is unsymmetric (different thicknesses, different angles), so it will show coupling.
-Physics Check: Look for the warning [!] BENDING-EXTENSION (B != 0).


Example 3
Engineering Constants (Standard Carbon)
Configuration: Quasi-Isotropic [0 / 90 / 45 / -45]s Material: Standard Carbon/Epoxy

1. Define Material
-Material: CARBON_STD
-Mode: A (Engineering Constants)
Input:
-E1: 181 GPa
-E2: 10.3 GPa
-ν12: 0.28
-G12: 7.17 GPa

2. Build Stack
-Add Layer: 0 deg, 0.125 mm
-Add Layer: 90 deg, 0.125 mm
-Add Layer: 45 deg, 0.125 mm
-Add Layer: -45 deg, 0.125 mm
Action: Select 3. Make Symmetric

3. Expected Results
-[A] Matrix: Should be isotropic in the plane (A11​≈A22​).
-[A] Shear: A16​ and A26​ should be 0.
-[D] Matrix: D16​ and D26​ will be Non-Zero (Bend-Twist coupling exists).
