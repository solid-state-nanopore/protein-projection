# Protein Projected Area Analysis

This script estimates how the two-dimensional projected area of a protein varies with its orientation. It reads Cα coordinates from a PDB structure, rotates them about the x and y axes, and estimates the area covered by circular disks centered on the projected Cα positions. The result is a matrix of areas for all sampled orientations.

## Requirements

- Python 3
- NumPy

Install the dependency with:

```bash
python -m pip install numpy
```

## Input

Provide a PDB file containing standard `ATOM` records with Cα atoms. The script reads atom names from columns 13–16 and Cartesian coordinates from columns 31–54 of each `ATOM` record. Other atoms, including records beginning with `HETATM`, are ignored. Coordinates are assumed to be in ångströms (Å). A file without readable Cα records does not produce an output matrix.

## Usage

If the uploaded script is named `protein projection(1).py`, run:

```bash
python "protein projection(1).py" protein.pdb -o area_matrix.dat
```

For a repository-friendly filename, rename the script to `protein_projection.py` and run:

```bash
python protein_projection.py protein.pdb -o area_matrix.dat
```

Optional arguments:

```bash
python protein_projection.py protein.pdb --step 0.5 --radius 1.7 -o area_matrix.dat
```

| Argument | Meaning | Default |
| --- | --- | --- |
| `input_pdb` | Path to the input PDB structure | Required |
| `-o`, `--output` | Output matrix path | `area_matrix.dat` |
| `--step` | Spacing of the area estimation grid, in Å | `1.0` |
| `--radius` | Radius assigned to each projected Cα disk, in Å | `1.7` |

## Calculation and output

The script calculates the mean Cα position and rotates the coordinates about the x axis and then the y axis around that center. Each angle ranges from 0° to 179.5° in 0.5° increments. The rotated Cα positions are projected onto the xy plane. Their disks are rasterized on a square grid, and the number of covered grid points multiplied by `step²` gives the estimated projected area in Å².

The output is a tab-delimited **360 × 360** matrix. Rows correspond to x-axis rotation angles and columns to y-axis rotation angles; index 0 is 0°, index 1 is 0.5°, and index 359 is 179.5°. The first line is a comment-style header listing the y-axis angles. The matrix contains areas formatted to four decimal places. For example:

```python
import numpy as np

areas = np.loadtxt("area_matrix.dat")  # shape: (360, 360)
area_at_x30_y45 = areas[60, 90]         # x = 30°, y = 45°
```

The script also prints the center coordinates and each angle pair's area while running. Its full scan includes 129,600 angle pairs, so processing can take time and print substantial terminal output, particularly for large structures or small grid steps.

## Interpretation

These values are **model-based estimates** of the union of fixed-radius disks centered on projected Cα atoms. They are not solvent-accessible surface areas or atomic-resolution silhouettes: side-chain atoms, van der Waals radii, and atoms present only in `HETATM` records are not included. Use the same `--step` and `--radius` when comparing structures.
