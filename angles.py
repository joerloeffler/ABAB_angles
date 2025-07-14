#!/usr/bin/env python3
"""
angles_distances.py
Compute COMs for two antibodies + antigen, output the angle and distances,
and (optionally) write ChimeraX or PyMOL visualisation scripts.

Usage:
    python calculate_angle_com.py config.json [--chimera] [--pymol]
"""
import os
import json
import argparse
import MDAnalysis as mda
import numpy as np
import pandas as pd


# ──────────────────────────── helpers ────────────────────────────
def load_config(path):
    with open(path, "r") as fh:
        return json.load(fh)


def parse_residue_ranges(ranges):
    """Expand list like ['1-5', '8'] → [1,2,3,4,5,8]."""
    out = []
    for r in ranges:
        if "-" in r:
            a, b = map(int, r.split("-"))
            out.extend(range(a, b + 1))
        else:
            out.append(int(r))
    return out


def build_selection(block):
    """
    Build an MDAnalysis selection string.
      - full chains: {"segid": "E"}
      - segid+residues: {"segid": "A", "residues": ["1-122"]}
    """
    clauses = []
    if "segments" not in block:
        raise ValueError("Missing 'segments' in antibody/antigen block.")
    for seg in block["segments"]:
        sid = seg["segid"]
        if "residues" in seg:
            res_clause = " or ".join(f"resid {r}" for r in parse_residue_ranges(seg["residues"]))
            clauses.append(f"(segid {sid} and ({res_clause}))")
        else:
            clauses.append(f"(segid {sid})")
    return " or ".join(clauses)


def com(universe, selection):
    atoms = universe.select_atoms(f"({selection}) and name CA")  # backbone CA only
    if not len(atoms):
        raise ValueError(f"No atoms selected for: {selection}")
    return atoms.center_of_mass()


def angle(a, b, c):
    v1, v2 = a - c, b - c
    cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))


def distance(a, b):
    return np.linalg.norm(a - b)


# ──────────────────────────── writers ────────────────────────────
def write_pdb(com_dict, out_pdb):
    """Write COMs as three pseudo-atoms; resnames A1, A2, AG (3 chars)."""
    with open(out_pdb, "w") as fh:
        for i, coord in enumerate(com_dict.values(), start=1):
            resname = ("A1", "A2", "AG")[i - 1]
            fh.write(
                f"HETATM{i:5d}  CA  {resname:<3s} A{i:4d}    "
                f"{coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}  1.00  0.00           C\n"
            )


def write_chimerax(angle_deg, pdb_file, com_pdb, out_cxc):
    with open(out_cxc, "w") as fh:
        fh.write(f"""\
open {pdb_file}
open {com_pdb}
set bgColor white

# Show cartoon for structure
show #1 cartoon
hide #1 atoms

# Color and style
style #2 sphere
color #2/A:1@CA black
color #2/A:2@CA black
color #2/A:3@CA black

# Show distances and angle
distance #2/A:1@CA #2/A:3@CA
distance #2/A:2@CA #2/A:3@CA
distance #2/A:1@CA #2/A:2@CA

distance style dashes 1
distance style color black

hide #3.1 models

angle #2/A:1@CA #2/A:3@CA #2/A:2@CA
lighting soft
""")


def write_pymol_script(angle_deg, pdb_file, com_pdb, out_pml):
    with open(out_pml, "w") as fh:
        fh.write(f"""\
load {pdb_file}, prot
load {com_pdb}, coms

hide everything, prot
show cartoon, prot

hide everything, coms
show spheres, coms
color black, coms

# Draw distances
distance dist1, coms///1/CA, coms///3/CA
distance dist2, coms///2/CA, coms///3/CA
distance dist3, coms///1/CA, coms///2/CA

# Show angle
angle ang1, coms///1/CA, coms///3/CA, coms///2/CA

set dash_color, black
set dash_width, 2
set sphere_scale, 0.6, coms

bg_color white
""")

# ──────────────────────────── main ────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config_file", help="JSON config")
    ap.add_argument("--chimera", action="store_true", help="write ChimeraX script")
    ap.add_argument("--pymol",   action="store_true", help="write PyMOL script")
    args = ap.parse_args()

    cfg = load_config(args.config_file)
    pdb_path = cfg["pdb_file"]
    tag      = os.path.splitext(os.path.basename(pdb_path))[0]

    out_pdb  = f"{tag}_coms.pdb"
    out_csv  = f"{tag}_angle.csv"
    out_cxc  = f"{tag}_visualize.cxc"
    out_pml  = f"{tag}_visualize.pml"

    u = mda.Universe(pdb_path)

    sel1 = build_selection(cfg["antibody1"])
    sel2 = build_selection(cfg["antibody2"])
    sel3 = build_selection(cfg["antigen"])

    c1, c2, c3 = com(u, sel1), com(u, sel2), com(u, sel3)
    ang = angle(c1, c2, c3)

    # Distances
    d1 = distance(c1, c3)  # AB1–AG
    d2 = distance(c2, c3)  # AB2–AG
    d12 = distance(c1, c2) # AB1–AB2

    # Write outputs
    write_pdb({"ANT1": c1, "ANT2": c2, "ANTG": c3}, out_pdb)
    pd.DataFrame([{
        "Angle_deg": round(ang, 2),
        "Distance_AB1_AG": round(d1, 2),
        "Distance_AB2_AG": round(d2, 2),
        "Distance_AB1_AB2": round(d12, 2)
    }]).to_csv(out_csv, index=False)

    print(f"Angle AB1–AB2 about AG: {ang:.2f}°")
    print(f"Distances: AB1–AG = {d1:.2f} Å, AB2–AG = {d2:.2f} Å, AB1–AB2 = {d12:.2f} Å")

    if args.chimera:
        write_chimerax(ang, pdb_path, out_pdb, out_cxc)
        print(f"ChimeraX script → {out_cxc}")

    if args.pymol:
        write_pymol_script(ang, pdb_path, out_pdb, out_pml)
        print(f"PyMOL script → {out_pml}")


if __name__ == "__main__":
    main()
