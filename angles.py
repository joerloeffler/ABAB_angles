import json
import MDAnalysis as mda
import numpy as np
import pandas as pd

def load_config(config_file):
    with open(config_file, 'r') as f:
        return json.load(f)

def parse_residue_ranges(ranges):
    residues = []
    for r in ranges:
        if "-" in r:
            start, end = map(int, r.split("-"))
            residues.extend(range(start, end + 1))
        else:
            residues.append(int(r))
    return residues

def calculate_com(universe, residue_ranges):
    residues = parse_residue_ranges(residue_ranges)
    selection = " or ".join([f"resid {r}" for r in residues])
    atoms = universe.select_atoms(f"protein and ({selection})")
    return atoms.center_of_mass()

def calculate_angle(com1, com2, com3):
    vec1 = com1 - com3  # Vector from antigen to antibody 1
    vec2 = com2 - com3  # Vector from antigen to antibody 2
    cosine_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
    return np.degrees(np.arccos(np.clip(cosine_angle, -1.0, 1.0)))

def write_pdb(coms, output_pdb):
    with open(output_pdb, 'w') as f:
        for i, (name, com) in enumerate(coms.items(), start=1):
            f.write(f"HETATM{i:5d}  CA  {name[:3].upper()} A   1    {com[0]:8.3f}{com[1]:8.3f}{com[2]:8.3f}  1.00  0.00           C\n")

def main(config_file):
    config = load_config(config_file)
    u = mda.Universe(config["pdb_file"])

    com_antibody1 = calculate_com(u, config["antibody1_residues"])
    com_antibody2 = calculate_com(u, config["antibody2_residues"])
    com_antigen = calculate_com(u, config["antigen_residues"])

    angle = calculate_angle(com_antibody1, com_antibody2, com_antigen)

    coms = {"ANT1": com_antibody1, "ANT2": com_antibody2, "ANTG": com_antigen}
    write_pdb(coms, "com_visualization.pdb")

    df = pd.DataFrame([{"Angle": angle}])
    df.to_csv("angle_output.csv", index=False)

    print(f"Angle between antibodies: {angle:.2f} degrees")

if __name__ == "__main__":
    main("config.json")
