import numpy as np
import json
import os
import random
from datetime import datetime


# ==============================================================================
# 1. CORE MATH & PHYSICS ENGINE
# ==============================================================================

class UnitSystem:
    def __init__(self):
        self.stress = "MPa"
        self.length = "mm"
        self.force = "N/mm"
        self.moment = "N"

    def set_units(self, choice):
        if choice == '1':  # SI (MPa, mm)
            self.stress, self.length, self.force, self.moment = "MPa", "mm", "N/mm", "N"
        elif choice == '2':  # SI (GPa, mm)
            self.stress, self.length, self.force, self.moment = "GPa", "mm", "kN/mm", "kN"
        elif choice == '3':  # SI (Pa, m)
            self.stress, self.length, self.force, self.moment = "Pa", "m", "N/m", "N"
        elif choice == '4':  # Imperial (psi, in)
            self.stress, self.length, self.force, self.moment = "psi", "in", "lb/in", "lb"
        elif choice == '5':  # Imperial (Msi, in)
            self.stress, self.length, self.force, self.moment = "Msi", "in", "kips/in", "kips"


class Material:
    def __init__(self, name, mode, props=None, matrix=None):
        self.name = name.upper()
        self.mode = mode
        self.data = props if props else {}
        self.raw_matrix = np.array(matrix) if matrix else None

        if mode == 'props_2d':
            self.Q_0 = self._calc_Q_from_2d_props(props)
        elif mode == 'direct_q':
            self.Q_0 = self.raw_matrix
        elif mode == 'stiffness_c':
            self.Q_0 = self._calc_Q_from_C_matrix(self.raw_matrix)

    def _calc_Q_from_2d_props(self, p):
        nu21 = p['nu12'] * (p['E2'] / p['E1'])
        denom = 1 - p['nu12'] * nu21
        if denom <= 0: raise ValueError("Invalid Properties (1 - nu12*nu21 <= 0)")
        return np.array([
            [p['E1'] / denom, (p['nu12'] * p['E2']) / denom, 0],
            [(p['nu12'] * p['E2']) / denom, p['E2'] / denom, 0],
            [0, 0, p['G12']]
        ])

    def _calc_Q_from_C_matrix(self, C):
        if C.shape != (6, 6): raise ValueError("Mode C requires 6x6 Matrix.")
        c33 = C[2, 2]
        if c33 == 0: raise ValueError("C33 cannot be zero.")
        q11 = C[0, 0] - (C[0, 2] ** 2 / c33)
        q12 = C[0, 1] - (C[0, 2] * C[1, 2] / c33)
        q22 = C[1, 1] - (C[1, 2] ** 2 / c33)
        q66 = C[5, 5]
        return np.array([[q11, q12, 0], [q12, q22, 0], [0, 0, q66]])

    def get_transformed_Q(self, angle_deg):
        if angle_deg == 0: return self.Q_0
        theta = np.radians(angle_deg)
        c = np.cos(theta);
        s = np.sin(theta)
        c2, s2 = c * c, s * s
        c4, s4 = c2 * c2, s2 * s2

        Q11, Q12, Q22 = self.Q_0[0, 0], self.Q_0[0, 1], self.Q_0[1, 1]
        Q66 = self.Q_0[2, 2]

        Q_bar = np.zeros((3, 3))
        Q_bar[0, 0] = Q11 * c4 + Q22 * s4 + 2 * (Q12 + 2 * Q66) * s2 * c2
        Q_bar[1, 1] = Q11 * s4 + Q22 * c4 + 2 * (Q12 + 2 * Q66) * s2 * c2
        Q_bar[0, 1] = (Q11 + Q22 - 4 * Q66) * s2 * c2 + Q12 * (c4 + s4)
        Q_bar[1, 0] = Q_bar[0, 1]
        Q_bar[2, 2] = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * s2 * c2 + Q66 * (c4 + s4)
        Q_bar[0, 2] = (Q11 - Q12 - 2 * Q66) * s * c2 * c + (Q12 - Q22 + 2 * Q66) * s * s2 * c
        Q_bar[2, 0] = Q_bar[0, 2]
        Q_bar[1, 2] = (Q11 - Q12 - 2 * Q66) * s2 * s * c + (Q12 - Q22 + 2 * Q66) * s * c2 * c
        Q_bar[2, 1] = Q_bar[1, 2]
        return Q_bar

    def to_dict(self):
        return {
            'name': self.name, 'mode': self.mode, 'data': self.data,
            'matrix': self.raw_matrix.tolist() if self.raw_matrix is not None else None
        }


class LaminateSolver:
    def __init__(self, unit_sys):
        self.layers = []
        self.units = unit_sys
        self.z_sign = -1.0

    def add_layer(self, material, display_angle, thickness, is_pre_rotated=False):
        math_angle = 0.0 if is_pre_rotated else float(display_angle)
        self.layers.append({
            'mat': material, 'disp_angle': float(display_angle),
            'calc_angle': math_angle, 'thick': float(thickness)
        })

    def remove_last_layer(self):
        if self.layers: self.layers.pop()

    def make_symmetric(self):
        mirror = [l.copy() for l in reversed(self.layers)]
        self.layers.extend(mirror)

    def get_notation(self, layer_list=None):
        target = layer_list if layer_list is not None else self.layers
        if not target: return "[]"
        angles = [f"{int(l['disp_angle'])}" if l['disp_angle'].is_integer() else f"{l['disp_angle']:.1f}" for l in
                  target]
        return "[" + "/".join(angles) + "]"

    def classify_laminate(self, A, B, D):
        norm_B = np.linalg.norm(B)
        A16, A26 = A[0, 2], A[1, 2]
        is_sym = norm_B < 1e-4
        is_balanced = (abs(A16) < 1e-4) and (abs(A26) < 1e-4)
        angles = [l['disp_angle'] for l in self.layers]
        unique_angles = set(abs(a) for a in angles)
        desc = []
        if is_sym:
            desc.append("SYMMETRIC")
        else:
            desc.append("UNSYMMETRIC")
        if len(unique_angles) == 1 and 0 in unique_angles:
            desc.append("UNIDIRECTIONAL")
        elif len(unique_angles) == 2 and 0 in unique_angles and 90 in unique_angles:
            desc.append("CROSS-PLY")
        elif is_balanced:
            desc.append("BALANCED")
        else:
            desc.append("GENERAL ANGLE-PLY")
        return " ".join(desc)

    def calculate_stiffness(self, specific_layers=None):
        # Optimized to calculate once and return tuple
        target = specific_layers if specific_layers is not None else self.layers
        if not target: return np.zeros((3, 3)), np.zeros((3, 3)), np.zeros((3, 3)), 0

        total_h = sum(L['thick'] for L in target)
        z = -total_h / 2.0
        A = np.zeros((3, 3));
        B = np.zeros((3, 3));
        D = np.zeros((3, 3))

        for ply in target:
            Q_bar = ply['mat'].get_transformed_Q(ply['calc_angle'])
            z_top = z + ply['thick']
            A += Q_bar * (z_top - z)
            B += 0.5 * Q_bar * (z_top ** 2 - z ** 2)
            D += (1.0 / 3.0) * Q_bar * (z_top ** 3 - z ** 3)
            z = z_top
        B = B * self.z_sign
        return A, B, D, total_h

    def get_coupling_ratio(self, specific_layers=None):
        _, _, D, _ = self.calculate_stiffness(specific_layers)
        if D[0, 0] == 0: return 0.0
        return abs(D[0, 2] / D[0, 0])

    def solve_deformations(self, A, B, D, loads):
        # Uses pre-calculated A, B, D
        ABD = np.block([[A, B], [B, D]])
        try:
            ABD_inv = np.linalg.inv(ABD)
            deformations = ABD_inv @ loads
            return deformations
        except np.linalg.LinAlgError:
            return None

    def optimize_stack_iterative(self, iterations=5000):
        if not self.layers: return None, 0.0, 0.0
        current_ratio = self.get_coupling_ratio()
        best_layers = [l.copy() for l in self.layers]
        best_ratio = current_ratio
        _, B, _, _ = self.calculate_stiffness()
        is_symmetric = np.allclose(B, 0, atol=1e-5)

        if is_symmetric:
            search_pool = [l.copy() for l in self.layers[:len(self.layers) // 2]]
        else:
            search_pool = [l.copy() for l in self.layers]

        print(f"[OPT] Running {iterations} permutations...")
        for _ in range(iterations):
            candidate_pool = list(search_pool)
            random.shuffle(candidate_pool)
            if is_symmetric:
                if len(self.layers) % 2 != 0:
                    mid = self.layers[len(self.layers) // 2].copy()
                    mirror = [l.copy() for l in reversed(candidate_pool)]
                    candidate_stack = candidate_pool + [mid] + mirror
                else:
                    mirror = [l.copy() for l in reversed(candidate_pool)]
                    candidate_stack = candidate_pool + mirror
            else:
                candidate_stack = candidate_pool

            ratio = self.get_coupling_ratio(candidate_stack)
            if ratio < best_ratio:
                best_ratio = ratio
                best_layers = candidate_stack
        return best_layers, current_ratio, best_ratio


# ==============================================================================
# 2. REPORT GENERATION
# ==============================================================================

def generate_report(solver, A, B, D, h, deformations=None, load_vals=None):
    # Optimization: Inputs A, B, D are passed in, not re-calculated.
    laminate_type = solver.classify_laminate(A, B, D)
    u = solver.units
    layer_count = len(solver.layers) #Layer count

    lines = []
    lines.append("=" * 105)
    lines.append(f"                    LAMINATE ANALYSIS REPORT")
    lines.append(f"                    {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("=" * 105)
    lines.append(f"Class:    {laminate_type}")
    lines.append(f"Notation: {solver.get_notation()}")
    lines.append(f"Total Layers: {layer_count}")
    lines.append(f"Total H:  {h:.4f} {u.length}")
    lines.append(f"Z-Axis:   {'Z-Down' if solver.z_sign == -1 else 'Z-Up'}")

    def fmt(val):
        if abs(val) > 1e5 or (abs(val) < 1e-2 and val != 0): return f"{val:11.3e}"
        return f"{val:11.3f}"

    def print_3x3(M, title, unit):
        # Center the title over 38 characters (approx width of 3 columns)
        header = f"{title} ({unit})"
        lines.append(f"\n{header:^40}")
        lines.append("-" * 40)
        for i in range(3):
            row = " | " + " ".join([fmt(M[i, j]) for j in range(3)]) + " |"
            lines.append(row)
        lines.append("-" * 40)

    # --- PART 1: INDIVIDUAL MATRICES ---
    lines.append("\n" + "=" * 105)
    lines.append("SECTION 1: CONSTITUTIVE MATRICES")
    lines.append("=" * 105)
    print_3x3(A, "[A] EXTENSIONAL STIFFNESS", u.force)
    print_3x3(B, "[B] COUPLING STIFFNESS", u.force)
    print_3x3(D, "[D] BENDING STIFFNESS", u.moment)

    # --- PART 2: THEORETICAL BACKGROUND ---
    lines.append("\n" + "=" * 105)
    lines.append("SECTION 2: THEORY & PHYSICS")
    lines.append("=" * 105)
    lines.append("Integration Formulas (Jones / Reddy):")
    lines.append("  A_ij = Σ [ Q_bar * (z_k - z_k-1) ]")
    lines.append("  B_ij = 1/2 * Σ [ Q_bar * (z_k^2 - z_k-1^2) ]")
    lines.append("  D_ij = 1/3 * Σ [ Q_bar * (z_k^3 - z_k-1^3) ]")

    # Physics Analysis
    lines.append("\nPhysics Coupling Analysis:")
    if abs(A[0, 2]) > 1 or abs(A[1, 2]) > 1:
        lines.append("  [!] EXTENSION-SHEAR (A16, A26 != 0): Pulling causes Shearing.")
    else:
        lines.append("  [OK] No Extension-Shear Coupling.")

    if np.linalg.norm(B) > 1e-2:
        lines.append("  [!] BENDING-EXTENSION (B != 0): UNSYMMETRIC. Warping risk high.")
    else:
        lines.append("  [OK] Symmetric (B ~ 0). No Warping.")

    if abs(D[0, 2]) > 1 or abs(D[1, 2]) > 1:
        lines.append("  [!] BEND-TWIST (D16, D26 != 0): Bending causes Twisting.")
    else:
        lines.append("  [OK] No Bend-Twist Coupling.")

    # --- PART 3: FULL EQUATION VISUALIZATION ---
    lines.append("\n" + "=" * 105)
    lines.append("SECTION 3: CONSTITUTIVE EQUATION")
    lines.append("=" * 105)

    ABD = np.block([[A, B], [B, D]])

    #Centering
    # Each number column is ~12 chars. 3 cols = 36 chars.
    # We force the headers to be centered over exactly 38 characters.
    t_A = f"[ A ] Matrix ({u.force})"
    t_B = f"[ B ] Matrix ({u.force})"
    t_D = f"[ D ] Matrix ({u.moment})"

    header_top = (f" {{ N }}       {t_A:^39}            {t_B:^39}           {{ ε^0 }}")
    lines.append(header_top)
    lines.append("-" * 105)

    force_labels = ["Nx ", "Ny ", "Nxy"]
    moment_labels = ["Mx ", "My ", "Mxy"]
    strain_labels = ["ε_x^0", "ε_y^0", "γ_xy^0"]
    curv_labels = ["κ_x  ", "κ_y  ", "κ_xy "]

    for i in range(6):
        # LEFT VECTOR (LOADS)
        if i < 3:
            label = force_labels[i]
        else:
            label = moment_labels[i - 3]

        row_str = f" {{ {label:<3} }}  =  | "

        # MATRIX
        for j in range(6):
            row_str += fmt(ABD[i, j]) + " "
            if j == 2: row_str += " | "  # Mid line separator

        row_str += " | "

        # RIGHT VECTOR (STRAINS)
        # If deformations exist, show values. If not, show symbols.
        if deformations is not None:
            val_str = f"{deformations[i]:.4e}"
            row_str += f" * {{ {val_str:<9} }}"
        else:
            if i < 3:
                s_lab = strain_labels[i]
            else:
                s_lab = curv_labels[i - 3]
            row_str += f" * {{ {s_lab:<8} }}"

        lines.append(row_str)

        # Mid Separator
        if i == 2:
            lines.append("-" * 105)
            header_bot = f" {{ M }}       {t_B:^39}            {t_D:^39}           {{ κ }}"
            lines.append(header_bot)
            lines.append("-" * 105)

    # --- PART 4: APPLIED LOAD SUMMARY ---
    if load_vals is not None:
        lines.append("\n" + "=" * 105)
        lines.append("SECTION 4: LOAD SUMMARY")
        lines.append("=" * 105)
        lines.append(f"  Applied N: {{{load_vals[0]}, {load_vals[1]}, {load_vals[2]}}} {u.force}")
        lines.append(f"  Applied M: {{{load_vals[3]}, {load_vals[4]}, {load_vals[5]}}} {u.moment}")

    return "\n".join(lines)


# ==============================================================================
# 3. INTERFACE
# ==============================================================================

class App:
    def __init__(self):
        self.units = UnitSystem()
        self.solver = LaminateSolver(self.units)
        self.materials = {}
        self.lib_file = "mat_lib_v20.json"
        self.load_lib()

    def clear(self):
        os.system('cls' if os.name == 'nt' else 'clear')

    def load_lib(self):
        if os.path.exists(self.lib_file):
            try:
                with open(self.lib_file, 'r') as f:
                    data = json.load(f)
                    for k, v in data.items():
                        m = v.get('matrix')
                        self.materials[k] = Material(v['name'], v['mode'], v.get('data'), m)
            except:
                pass

    def save_lib(self):
        data = {k: v.to_dict() for k, v in self.materials.items()}
        with open(self.lib_file, 'w') as f: json.dump(data, f, indent=4)

    def get_float(self, prompt, allow_default=False):
        while True:
            val = input(prompt)
            if not val and allow_default: return 0.0
            try:
                return float(val)
            except ValueError:
                print("Invalid number.")

    def setup_units(self):
        self.clear()
        print("--- UNIT SYSTEM SETUP ---")
        print("1. SI [MPa, mm, N/mm]")
        print("2. SI [GPa, mm, kN/mm] (High Stiffness)")
        print("3. SI [Pa, m, N/m]")
        print("4. US [psi, in, lb/in]")
        print("5. US [Msi, in, kips/in]")
        c = input("Select: ")
        self.units.set_units(c)

    def define_material(self):
        print("\n--- MATERIAL WIZARD ---")
        name = input("Material Name: ").upper()
        print("A) Engineering Constants (E1, E2, ν12, G12)")
        print("B) Plane Stress Q Matrix (3x3)")
        print("C) 3D Stiffness C Matrix (6x6 -> Reduced)")
        mode = input("Choice: ").upper()

        if mode == 'A':
            p = {}
            p['E1'] = self.get_float("  E1: ")
            p['E2'] = self.get_float("  E2: ")
            print("  [REMINDER: For Transversely Isotropic: E3=E2, G12=G13]")
            p['nu12'] = self.get_float("  ν12: ")
            p['G12'] = self.get_float("  G12: ")
            print("  --- Optional 3D Props (Press Enter to skip) ---")
            p['G13'] = self.get_float("  G13 [Default=G12]: ", True) or p['G12']
            p['nu23'] = self.get_float("  ν23 [Default=0]: ", True)
            self.materials[name] = Material(name, 'props_2d', props=p)
        elif mode == 'B':
            print("Enter 3x3 Q Matrix:")
            m = []
            for i in range(3): m.append(list(map(float, input(f"Row {i + 1}: ").split())))
            self.materials[name] = Material(name, 'direct_q', matrix=m)
        elif mode == 'C':
            print("Enter 6x6 C Matrix:")
            m = []
            for i in range(6):
                t = input(f"Row {i + 1}: ")
                if not t: break
                m.append(list(map(float, t.split())))
            self.materials[name] = Material(name, 'stiffness_c', matrix=m)
        self.save_lib()

    def build_stack(self):
        while True:
            self.clear()
            print(f"--- STACK BUILDER [{self.units.length}] ---")
            print(f"Total Layers: {len(self.solver.layers)}")
            print(f"Notation: {self.solver.get_notation()}")
            if self.solver.layers:
                print(f"{'ID':<4} {'Mat':<10} {'DispAng':<8} {'Thick':<8}")
                for i, l in enumerate(self.solver.layers):
                    print(f"{i + 1:<4} {l['mat'].name:<10} {l['disp_angle']:<8} {l['thick']:<8}")

            print("\n1. Add Layer")
            print("2. Remove Last")
            print("3. Make Symmetric")
            print("4. Clear All")
            print("5. Done")

            ans = input("Choice: ")
            if ans == '1':
                if not self.materials:
                    input("Define material first!");
                    return
                print("Materials:", list(self.materials.keys()))
                m = input("Name: ").upper()
                if m in self.materials:
                    disp_ang = self.get_float("Physical Angle (e.g. 45): ")
                    thk = self.get_float(f"Thickness ({self.units.length}): ")
                    pre = False
                    if self.materials[m].mode == 'direct_q':
                        print("  [?] PRE-ROTATED CHECK:")
                        print("           - Say NO if this is raw material data (0-degree properties.")
                        print("           - Say YES if this matrix is ALREADY calculated for a specific angle (e.g. Q_45).")
                        chk = input("Is this Q-matrix PRE-ROTATED? (y/n): ").lower()
                        if chk == 'y': pre = True
                    self.solver.add_layer(self.materials[m], disp_ang, thk, is_pre_rotated=pre)
            elif ans == '2':
                self.solver.remove_last_layer()
            elif ans == '3':
                self.solver.make_symmetric()
            elif ans == '4':
                self.solver.layers = []
            elif ans == '5':
                break

    def run_optimization_workflow(self):
        if not self.solver.layers: return
        print("\n[OPTIMIZATION] Initializing Monte Carlo Permutation...")
        new_layers, r_old, r_new = self.solver.optimize_stack_iterative()

        def temp_notation(lst):
            return "[" + "/".join(
                [f"{l['disp_angle']:.0f}" if l['disp_angle'].is_integer() else f"{l['disp_angle']:.1f}" for l in
                 lst]) + "]"

        print("\n" + "=" * 40)
        print("      OPTIMIZATION RESULTS")
        print("=" * 40)
        print(f"Original D16/D11 Ratio: {r_old:.4f}")
        print(f"Optimized D16/D11 Ratio: {r_new:.4f}")
        print(f"\nOld Stack: {self.solver.get_notation()}")
        print(f"New Stack: {temp_notation(new_layers)}")
        if input("\nApply this configuration? (y/n): ").lower() == 'y':
            self.solver.layers = new_layers
            print("Stack updated.")

    def run(self):
        self.setup_units()
        while True:
            self.clear()
            print(f"ADVANCED LAMINATE CALCULATOR| Units: {self.units.stress}")
            print(f"Stack: {self.solver.get_notation()}")
            print("-" * 50)
            print("1. Material Wizard")
            print("2. Stack Builder")
            print("3. Run Analysis (Report)")
            print("4. Run Deformation (Apply Loads)")
            print("5. Optimize Stack")
            print("6. Toggle Z-Axis")
            print("7. Exit")

            c = input("Select: ")
            if c == '1':
                self.define_material()
            elif c == '2':
                self.build_stack()
            elif c == '3':
                if not self.solver.layers: continue
                # Opt: Calculate once here
                A, B, D, h = self.solver.calculate_stiffness()
                rep = generate_report(self.solver, A, B, D, h)
                print(rep)
                with open("Report_Stiff.txt", "w", encoding='utf-8') as f:
                    f.write(rep)
                input("\nReport Saved. Press Enter...")
            elif c == '4':
                if not self.solver.layers: continue
                print("\n--- DEFINE LOADS ---")
                print(f"Enter Force per width ({self.units.force}) and Moment ({self.units.moment})")
                Nx = self.get_float("Nx:  ", True)
                Ny = self.get_float("Ny:  ", True)
                Nxy = self.get_float("Nxy: ", True)
                Mx = self.get_float("Mx:  ", True)
                My = self.get_float("My:  ", True)
                Mxy = self.get_float("Mxy: ", True)
                loads = np.array([Nx, Ny, Nxy, Mx, My, Mxy])

                # Opt: Calculate once here
                A, B, D, h = self.solver.calculate_stiffness()
                deformations = self.solver.solve_deformations(A, B, D, loads)

                if deformations is None:
                    print("Error: Stiffness Matrix Singular (Unstable)")
                else:
                    rep = generate_report(self.solver, A, B, D, h, deformations, loads)
                    print(rep)
                    with open("Report_Deform.txt", "w", encoding='utf-8') as f:
                        f.write(rep)
                input("\nDone. Press Enter...")
            elif c == '5':
                self.run_optimization_workflow()
            elif c == '6':
                self.solver.z_sign *= -1
                print("Z-Axis.")
            elif c == '7':
                break


if __name__ == "__main__":
    App().run()