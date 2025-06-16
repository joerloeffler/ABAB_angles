#!/usr/bin/env python3
"""
calculate_angle_com.py
Compute COMs for two antibodies + antigen, output the angle,
and (optionally) write ChimeraX visualisation script.

Usage:
    python calculate_angle_com.py config.json [--chimera] 
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


# ──────────────────────────── writers ────────────────────────────
def write_pdb(com_dict, out_pdb):
    """Write COMs as three pseudo-atoms; resnames A1, A2, AG (3 chars)."""
    with open(out_pdb, "w") as fh:
        for i, coord in enumerate(com_dict.values(), start=1):
            resname = ("A1", "A2", "AG")[i - 1]  # ANT1/ANT2/ANTG labels in scripts
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

# Color
style #2 sphere

color #2/A:1@CA black
color #2/A:2@CA black
color #2/A:3@CA black

# Show distances + angle
distance #2/A:1@CA #2/A:3@CA
distance #2/A:2@CA #2/A:3@CA

distance style dashes 1
distance style color black

hide #3.1 models

angle #2/A:1@CA #2/A:3@CA #2/A:2@CA
lighting soft


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

    u = mda.Universe(pdb_path)

    sel1 = build_selection(cfg["antibody1"])
    sel2 = build_selection(cfg["antibody2"])
    sel3 = build_selection(cfg["antigen"])

    c1, c2, c3 = com(u, sel1), com(u, sel2), com(u, sel3)
    ang = angle(c1, c2, c3)

    # write outputs
    write_pdb({"ANT1": c1, "ANT2": c2, "ANTG": c3}, out_pdb)
    #pd.DataFrame([{"Angle_deg": ang}]).to_csv(out_csv, index=False)
    pd.DataFrame([{"Angle_deg": round(ang, 2)}]).to_csv(out_csv, index=False)
    print(f"Angle ANT1–ANT2 about ANTG: {ang:.2f}°")

    if args.chimera:
        write_chimerax(ang, pdb_path, out_pdb, out_cxc)
        print(f"ChimeraX script → {out_cxc}")



if __name__ == "__main__":
    main()
