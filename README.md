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

**Julia** (recommended: Julia >= 1.9). On macOS with Homebrew:
```bash
brew install julia
```

## Installation / Setup

1. **Clone or copy the repository** 
```
git clone https://github.com/epoliukhina/Potential.jl.git
```

2. Enter the repository and install the package dependencies:
```bash
cd Potential.jl
julia --project=. -e 'import Pkg; isempty(Pkg.Registry.reachable_registries()) && Pkg.Registry.add("General"); Pkg.resolve(); Pkg.instantiate()'
```

3. Install Pluto (once, in your default Julia environment):
```julia
import Pkg
Pkg.add("Pluto")
```

4. Launch Pluto:
```julia
using Pluto
Pluto.run()
```

5. Open `example/cryoPOT_pluto_notebook.jl` in Pluto. The notebook activates this repository environment automatically. Follow the notebook steps to load a file and run analysis.

If you have the BSA example dataset, put it at `data/BSA_20mgmL_Position_1_2.txt`; the notebook will use that file by default. If that file is absent, it falls back to the bundled small sample in `example/data/`.
