<div align="center">

# BLADE

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/license/gpl-3-0)
[![Python](https://img.shields.io/badge/python-3.10%2B-brightgreen.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-0.2.0-orange.svg)](pyproject.toml)

**Boride Learning and Design Engine — automated CALPHAD thermodynamic database generation for multicomponent materials systems.**

[Report a Bug](https://github.com/ichasekatz/BLADE/issues/new?labels=bug) · [Request a Feature](https://github.com/ichasekatz/BLADE/issues/new?labels=enhancement)

</div>

---

<div align="center">
<img width="700" alt="BLADE_Framework" src="https://github.com/user-attachments/assets/8a64271d-0427-4e2e-b964-df0af7ff18d0" />
</div>

---

## Overview

BLADE automates the full CALPHAD workflow from start to finish: given an element pool and a crystal structure prototype, it enumerates every valid N-component system, generates SQS supercells, relaxes them with an ML interatomic potential, fits CALPHAD parameters, and produces `.tdb` files — without manual intervention at any step.

Inspired by [MaterialsFramework](https://github.com/dogusariturk/MaterialsFramework)'s modular design, BLADE builds on top of it as a computational backend. **BLADE requires the [ichasekatz fork of MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework)** — install it before installing BLADE.

- **Composition generation** — enumerate binary, ternary, or N-component systems from primary and secondary element pools with configurable count bounds
- **SQS generation** — drive ATAT `mcsqs` in parallel across all compositions and supercell sizes
- **TDB fitting** — relax structures with any MLIP supported by MaterialsFramework, compute formation energies and thermal properties, fit CALPHAD parameters, write `.tdb` files
- **Visualization** — combine CONTCAR structure images, stitch phase-diagram PNGs, render relaxation movies
- **Volume analysis** — scan POSCAR files and extract lattice parameters and per-atom volumes into a DataFrame

---

## What BLADE Does Differently

Existing tools handle pieces of this workflow in isolation — a calculator relaxes one structure, a fitting tool processes one system, a visualization tool reads one file. BLADE connects all of these into a single automated pipeline driven by a composition list.

| Capability | BLADE |
|---|---|
| Structure relaxation | Batched across every composition and SQS level, instead of a single structure at a time |
| Property calculation | Automated per-SQS-level pipeline, instead of per-structure manual invocation |
| Composition enumeration | `BladeCompositions` traverses the full N-component search space automatically |
| SQS generation | `BladeSQS` runs parallel `mcsqs` across all compositions simultaneously |
| CALPHAD database output | `BladeTDBGen` processes all systems and phases in one call, instead of one system per run |
| Fixed sublattice sites | `sqsgen_levels2` holds any species fixed on any sublattice, enabling ceramics and intermetallics |
| Phase diagrams | Built-in pseudo-binary plots via `pycalphad`, instead of a separate workflow |
| Visualization | `BladeVisualizer` combines structures, diagrams, and relaxation movies in one place |

The fixed-sublattice support is particularly significant: BLADE is not limited to simple alloys. Any crystal structure where some sites have fixed occupancy (e.g., boron in diborides, oxygen in oxides) can be described through a lattice coordinate string and a phase prototype dictionary, making BLADE applicable across a wide range of ceramic and intermetallic systems.

---

## Modules

### Tools

| Class | Description |
|---|---|
| `BladeCompositions` | Enumerate compositions from primary and secondary element pools |
| `BladeSQS` | Generate special quasirandom structures using ATAT `mcsqs` |
| `BladeTDBGen` | Relax structures and fit CALPHAD thermodynamic databases |
| `BladeCutoff` | Compute radial cutoff distances from SQS supercell geometry |

### Analysis

| Class | Description |
|---|---|
| `BladeVisualizer` | Combine structure images, phase diagrams, and relaxation movies |
| `BLADEVolume` | Extract volumes and lattice parameters from POSCAR files |

All classes are available via lazy top-level import — heavy dependencies are not loaded until first use:

```python
from blade import BladeCompositions, BladeSQS, BladeTDBGen
from blade import BladeVisualizer, BLADEVolume
```

---

## Getting Started

### Prerequisites

- [ichasekatz/MaterialsFramework](https://github.com/ichasekatz/MaterialsFramework) — required fork, install before BLADE
- [ATAT](https://www.brown.edu/Departments/Engineering/Labs/avdw/atat/) (`mcsqs` must be on `$PATH`)

### Installation

Install dependencies with [pixi](https://pixi.sh) and pip:

```bash
pixi add git python==3.12 cmake pip psutil pandarallel orb-models pytorch tcsh pymatgen sqsgenerator
pip install matplotlib pillow imageio pycalphad torch ase pynanoflann
pip install tensorflow
pip install tensorpotential
```

Install the required MaterialsFramework fork:

```bash
git clone https://github.com/ichasekatz/MaterialsFramework.git
cd MaterialsFramework
pip install -e .
cd ..
```

Install BLADE:

```bash
git clone https://github.com/ichasekatz/BLADE.git
cd BLADE
pip install -e .
```

---

## Quick Start

### 1. Generate compositions

```python
from blade.tools.blade_compositions import BladeCompositions

composer = BladeCompositions(
    primary_elements=["Hf", "Cr", "Ta", "Ti", "Mo"],
    secondary_elements=[],
    system_size=3,
    primary_min=3, primary_max=3,
    secondary_min=0, secondary_max=0,
    allow_lower_order=False,
)

composition_list = composer.generate_compositions()
print(f"{len(composition_list)} ternary systems: {composition_list[:5]}")
```

### 2. Generate SQS structures

```python
from blade.tools.blade_sqsgen import BladeSQS

sqs_gen = BladeSQS(
    phases_dict=phases["HEDB1"],
    sqsgen_levels=sqsgen_levels,
    level=5,
    len_comp=2,
    skip_existing_sqs=True,
)
sqs_gen.sqs_gen(phase=phase, paths=paths, iter=1_000_000, params=mcsqs_params)
```

### 3. Fit thermodynamic databases

```python
from blade.tools.blade_tdb_gen import BladeTDBGen

gen = BladeTDBGen(
    phases=phase_list,
    liquid=False,
    paths=paths,
    composition_list=composition_list,
    level=5,
    sqsgen_levels2=sqsgen_levels2,
    skip_existing=False,
)
gen.fit()  # explicit call — no side effects on construction
```

### 4. Visualize results

```python
from blade.analysis.blade_visual import BladeVisualizer

viz = BladeVisualizer()
viz.contcar(contcars, save="CrHf/Combined_CONTCARs.png")
viz.phase_diagram(pngs, save="Combined_Phase_Diagrams.png")
```

---

## Example Scripts

| Script | Description |
|---|---|
| [`01_compositions.py`](examples/01_compositions.py) | Enumerate binary/ternary/N-component compositions |
| [`02_sqs_generation.py`](examples/02_sqs_generation.py) | Generate SQS structures with ATAT `mcsqs` |
| [`03_tdb_generation.py`](examples/03_tdb_generation.py) | Relax structures and fit TDB databases |
| [`04_visualization.py`](examples/04_visualization.py) | Combine structures, phase diagrams, and movies |
| [`05_volume_analysis.py`](examples/05_volume_analysis.py) | Extract per-atom volumes from POSCAR files |
| [`tdb_gen.py`](examples/tdb_gen.py) | Full end-to-end workflow (HPC driver script) |

---

## Project Structure

```
BLADE/
├── src/
│   └── blade/
│       ├── __init__.py              # lazy loader — all public classes
│       ├── tools/
│       │   ├── __init__.py
│       │   ├── blade_compositions.py
│       │   ├── blade_cutoff.py
│       │   ├── blade_sqsgen.py
│       │   └── blade_tdb_gen.py
│       └── analysis/
│           ├── __init__.py
│           ├── blade_visual.py
│           └── blade_volume.py
├── examples/
│   ├── 01_compositions.py
│   ├── 02_sqs_generation.py
│   ├── 03_tdb_generation.py
│   ├── 04_visualization.py
│   ├── 05_volume_analysis.py
│   └── tdb_gen.py
└── pyproject.toml
```

---

## License

Distributed under the GPL-3.0-or-later License.

## Contact

Chase Katz — [ichasekatz@tamu.edu](mailto:ichasekatz@tamu.edu)
