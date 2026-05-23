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

**Julia 1.12** is recommended — that is the version the package is tested against on CI. Earlier versions back to 1.6.7 will still install and run, but `using Potential` prints a warning every time it loads on an untested Julia.

Install Julia from https://julialang.org/downloads/ (the official `juliaup` installer works on Linux, macOS, and Windows), or on macOS with Homebrew:
```bash
brew install julia
```

After installing, confirm Julia is on your `PATH`:
```bash
julia --version
```

## Installation / Setup

All commands below are run in a terminal unless the prompt is shown as `julia>`, which means you are inside the Julia REPL.

### 1. Clone the repository

```bash
git clone https://github.com/epoliukhina/Potential.git
cd Potential
```

### 2. Install the package's dependencies

From inside the `Potential` directory:
```bash
julia --project=. -e 'import Pkg; isempty(Pkg.Registry.reachable_registries()) && Pkg.Registry.add("General"); Pkg.resolve(); Pkg.instantiate()'
```
This downloads and precompiles everything the package needs (~150 packages, a few minutes on first run).

### 3. Install Pluto into your default Julia environment (one-time)

Pluto must live in your **default** environment, not in this project's environment. Start a fresh Julia session *without* `--project=.`:
```bash
julia
```
Then, at the `julia>` prompt:
```julia
import Pkg
Pkg.add("Pluto")
```
You only need to do this once per machine.

### 4. Launch Pluto

In the same Julia session (or any new `julia` session started without `--project=.`):
```julia
using Pluto
Pluto.run()
```
Pluto opens a tab in your browser automatically.

### 5. Open the notebook

In the Pluto tab, open `example/cryoPOT_pluto_notebook.jl` from this repository. The notebook activates this repository's environment automatically. Follow the notebook steps to load a file and run analysis.

If you have the BSA example dataset, put it at `data/BSA_20mgmL_Position_1_2.txt`; the notebook will use that file by default. If that file is absent, it falls back to the bundled small sample in `example/data/`.
