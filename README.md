# angles.py

A Python script to compute the angle and distances between two antibodies and an antigen based on center-of-mass (COM) calculations. Optionally, it generates visualization scripts for **ChimeraX** or **PyMOL**.

---

##  Usage

```bash
python angles.py config.json [--chimera] [--pymol]
```
config.json: JSON file specifying the input PDB and selection definitions.

--chimera: Optional. Generate a ChimeraX script (.cxc) for visualization.

--pymol: Optional. Generate a PyMOL script (.pml) for visualization.

```
| File                 | Description                            |
| -------------------- | -------------------------------------- |
| `<pdbtag>_angle.csv` | Computed angle and distances (Å)       |
| `<pdbtag>_coms.pdb`  | Pseudo-atoms at computed COM positions |
| -------------------- | -------------------------------------- |
```
## Example jsopn file 
```json
{
  "pdb_file": "structure.pdb",
  "antibody1": {
    "segments": [{ "segid": "A", "residues": ["1-122"] }]
  },
  "antibody2": {
    "segments": [{ "segid": "B", "residues": ["1-122"] }]
  },
  "antigen": {
    "segments": [{ "segid": "C" }]
  }
}
```
