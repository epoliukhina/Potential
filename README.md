# cryoPOT

Interactive Pluto.jl notebook for analyzing 3D particle coordinate datasets (e.g., cryo-ET / cryo-TEM picks) and extracting:
- radial distribution function **g(r)**
- potential of mean force **W(r) = -ln g(r)**
- structure factor **S(q)** (Debye histogram method)
- Kirkwood–Buff integral **G₂₂** with optional surface-to-volume (S/V) extrapolation fit

Designed to be usable by non-programmers (interactive forms) and still transparent for advanced users (toggle code).

---

## Features

- **Load coordinates** from `.txt`, `.tsv`, or `.csv` (3 columns: x y z)
- **Rotation/alignment** of coordinates using `Potential.rotate_particles` and visual fit diagnostics
- **Cropping** (box cut in XY and Z) with 2D/3D previews
- **Box summary**: particle count, box size, number density, concentration (mg/mL, from MW)
- **RDF + PMF** computation and plotting; saves `*_RDF.txt` and `*_Gr.pdf`
- **Structure factor S(q)** computation and plotting; optional SAXS window; saves `*_Sq.txt`, `*_Sq.pdf`
- **KBI G₂₂** computation across box sizes/steps; saves `*_KBI.txt`, `*_KBI.pdf`
- **Auto-fit window** for G₂₂ vs S/V extrapolation; saves `*_KBI_SVfit.pdf`

Outputs are written to a per-dataset folder inside the chosen base output directory.

---

## Prerequisites

**Julia** (recommended: Julia ≥ 1.9)

## Installation / Setup

1. **Clone or copy the repository** 
```
git clone https://github.com/epoliukhina/Potential.jl.git
```

2. Start Julia and install Pluto (once):
```julia
import Pkg
Pkg.add("Pluto")
```

3. Launch Pluto:
```julia
using Pluto
Pluto.run()
```

4. Open the `cryoPOT_pluto_notebook.jl` notebook in Pluto. Follow the subsequent steps to load the file and perform analysis. 
