"""Full BLADE workflow: compositions → SQS → TDB → phase diagrams → visualization.

This is the main driver script. It mirrors the original BLADE/BLADE/examples/tdb_gen.py
but uses the refactored BLADE API:

  - BladeTDBGen stores config on construction; call gen.fit() to run
  - BladeVisualizer (not BLADEVisualizer)
  - All classes imported from blade.*

Run on HPC:
    pkill -f "BLADE/examples/tdb_gen.py"
    rm tdb_gen.log
    nohup python -u BLADE/examples/tdb_gen.py > tdb_gen.log 2>&1 &
    tail -f tdb_gen.log

Mount HPC over SSHFS (optional):
    diskutil unmount force ~/HPC 2>/dev/null
    rm -rf ~/HPC && mkdir -p ~/HPC
    sshfs ichasekatz@10.125.153.116:/home/ichasekatz ~/HPC \\
      -o reconnect,ServerAliveInterval=15,ServerAliveCountMax=3
"""

import os

os.environ["OMP_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"

import torch

torch.set_num_threads(8)
torch.set_num_interop_threads(2)

from pathlib import Path

import matplotlib.pyplot as plt
from pycalphad import Database, binplot
from pycalphad import variables as v

from blade.tools.blade_compositions import BladeCompositions
from blade.tools.blade_sqsgen import BladeSQS
from blade.tools.blade_tdb_gen import BladeTDBGen
from blade.analysis.blade_visual import BladeVisualizer

# ------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------
path0 = Path("/Users/chasekatz/Desktop/School/Research")
path1 = path0 / "BLADE"
path2 = path0 / "PhaseForge" / "PhaseForge" / "atat" / "data" / "sqsdb"
paths = [path0, path1, path2]

# ------------------------------------------------------------------
# Run flags
# ------------------------------------------------------------------
level = 5
sqs_iter = 1_000_000
run_sqs = True
skip_existing_sqs = True

run_tdb = True
skip_existing_tdb = False

run_movie = False

# ------------------------------------------------------------------
# Elements and composition constraints
# ------------------------------------------------------------------
transition_metals = ["Zr", "Hf", "Ta", "Cr", "Ti", "V", "Nb", "Mo", "W"]
rare_earths: list[str] = []
# transition_metals = ["Ti", "Cr", "W"]
transition_metals = ["Hf", "Cr"]
# transition_metals = ["Ni", "Re"]

system_size = 2
tm_element_range = [2, 2]
re_element_range = [0, 0]
allow_lower_order = False

# ------------------------------------------------------------------
# Phase prototypes
# ------------------------------------------------------------------
phases: dict[str, dict] = {
    "HEDB1": {
        "a": 1,
        "b": 1,
        "c": 1,
        "alpha": 90,
        "beta": 90,
        "gamma": 120,
        "vectors": "1 0 0\n0 1 0\n0 0 1\n",
        "coords": (
            "0.000000 0.000000 0.000000 a\n"
            "0.333333 0.666667 0.500000 B\n"
            "0.666667 0.333333 0.500000 B\n"
        ),
    },
}

phase_list = [
    # supercell_size controls n_atoms = unit_cell_sites × product(supercell_size)
    # (4,3,2) → 72 atoms   (4,4,3) → 144   (6,6,2) → 216   (8,6,2) → 288
    {"generator_name": "HEDB", "lattice": "HEDB1", "supercell_size": (4, 3, 2)},
]

liquid = False

# ------------------------------------------------------------------
# mcsqs run parameters
# ------------------------------------------------------------------
mcsqs_params = {
    "use_time": True,
    "time": 30,          # seconds per sqsdb directory
    "2": 5,              # pair cutoff  = 5th-nearest shell
    "3": 4,              # triplet cutoff
    "4": 3,              # quadruplet cutoff
    "wr": 20,
    "wn": 0.75,
    "wd": 1,
    "parallel_runs": 10,
    "stop_grace": 30,
}

# ------------------------------------------------------------------
# Fixed-sublattice species (B sublattice in HEDB1)
# ------------------------------------------------------------------
sqsgen_levels2 = [
    {"element": "B", "compositions": "1.0", "letter": "b", "count": "2"},
    # {"element": "O", "compositions": "1.0", "letter": "c", "count": "1"},
]

# ------------------------------------------------------------------
# SQS composition levels
# ------------------------------------------------------------------
sqsgen_levels = [
    {"level": 0, "compositions": [[1.0, 0.0]],                            "letter": ["a"]},
    {"level": 1, "compositions": [[0.5, 0.5]],                             "letter": ["a"]},
    {"level": 2, "compositions": [[0.75, 0.25]],                           "letter": ["a"]},
    {"level": 3, "compositions": [[0.33333, 0.33333, 0.33333]],            "letter": ["a"]},
    {"level": 4, "compositions": [[0.5, 0.25, 0.25]],                      "letter": ["a"]},
    {"level": 5, "compositions": [[0.875, 0.125], [0.625, 0.375]],         "letter": ["a"]},
    {"level": 6, "compositions": [[0.75, 0.125, 0.125]],                   "letter": ["a"]},
]

# ------------------------------------------------------------------
# 1. Generate compositions
# ------------------------------------------------------------------
composer = BladeCompositions(
    primary_elements=transition_metals,
    secondary_elements=rare_earths,
    system_size=system_size,
    primary_min=tm_element_range[0],
    primary_max=tm_element_range[1],
    secondary_min=re_element_range[0],
    secondary_max=re_element_range[1],
    allow_lower_order=allow_lower_order,
)

composition_list = composer.generate_compositions()
unique_len_comps = composer.get_systems()

print(f"Compositions ({len(composition_list)} total): {composition_list}")
print(f"System sizes: {unique_len_comps}")

# ------------------------------------------------------------------
# 2. Generate SQS structures
# ------------------------------------------------------------------
if run_sqs:
    for specific_phase in phase_list:
        for len_comp in unique_len_comps:
            sqs_gen = BladeSQS(
                phases_dict=phases[specific_phase["lattice"]],
                sqsgen_levels=sqsgen_levels,
                level=level,
                len_comp=len_comp,
                skip_existing_sqs=skip_existing_sqs,
            )
            params = mcsqs_params | {"super_cell_size": specific_phase["supercell_size"]}
            sqs_gen.sqs_gen(phase=specific_phase, paths=paths, iter=sqs_iter, params=params)
            # sqs_gen.rename_files(specific_phase, paths, sqsgen_levels2)

# ------------------------------------------------------------------
# 3. Fit TDB databases
# ------------------------------------------------------------------
if run_tdb:
    gen = BladeTDBGen(
        phases=phase_list,
        liquid=liquid,
        paths=paths,
        composition_list=composition_list,
        level=level,
        sqsgen_levels2=sqsgen_levels2,
        skip_existing=skip_existing_tdb,
    )
    gen.fit()  # explicit call — no side effects on construction

# ------------------------------------------------------------------
# 4. Plot binary phase diagrams (uncomment to enable)
# ------------------------------------------------------------------
PHASE_DIAGRAM_SYSTEM_SIZE = 2

def plot_phase_diagram(tdb, elements, num_elements, output_path):
    """
    Plot a pseudo-binary metal-sublattice diagram for an HEDB1_2 diboride phase.

    For (M1,M2)B2, boron is fixed at X(B)=2/3 and the varied metal
    has a total-system mole fraction from 0 to 1/3.
    """
    if num_elements != 2:
        print(f"Skipping phase diagram: expected binary metal system, got {elements}")
        return

    metal_1 = elements[0].upper()
    metal_2 = elements[1].upper()

    # B must be included because HEDB1_2 contains a fixed boron sublattice.
    components = [metal_1, metal_2, "B"]

    phases = ['HEDB1_2']

    print(f"Plotting system: ({metal_1},{metal_2})B2")
    print(f"Active components: {components}")
    print(f"Phases: {phases}")

    conditions = {
        v.N: 1,
        v.P: 101325,
        v.T: (300, 4500, 50),

        # Fixed boron fraction for (Cr,Hf)B2
        v.X("B"): 2 / 3,

        # Hf total mole fraction varies from CrB2 to HfB2
        v.X(metal_2): (1e-6, (1 / 3) - 1e-6, 0.005),
    }

    fig, ax = plt.subplots(figsize=(10, 7))

    binplot(
        tdb,
        components,
        phases,
        conditions,
        plot_kwargs={"ax": ax},
    )

    # Convert total-system mole fraction to metal-sublattice fraction:
    # y = metal_2 / (metal_1 + metal_2) = 3 * X(metal_2)
    ticks = [0.00, 1 / 12, 1 / 6, 1 / 4, 1 / 3]
    labels = ["0.00", "0.25", "0.50", "0.75", "1.00"]

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=14)
    ax.set_xlabel(
        rf"Metal-site fraction $x_{{{metal_2}}}$ in $({metal_1}_{{1-x}}{metal_2}_x)B_2$",
        fontsize=16,
    )
    ax.set_ylabel("Temperature (K)", fontsize=16)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_title(rf"$({metal_1},{metal_2})B_2$ Pseudo-Binary Phase Diagram", fontsize=18)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved phase diagram to: {output_path}")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from pycalphad import equilibrium, variables as v


def plot_phase_diagram(tdb, elements, num_elements, output_path):
    """
    Plot Gibbs energy of the single HEDB1_2 phase across the
    (M1,M2)B2 pseudo-binary section.

    This is an energy diagnostic, not a phase diagram, because the
    current TDB contains only one active phase.
    """
    if num_elements != 2:
        print(f"Skipping Gibbs-energy plot: expected two metals, got {elements}")
        return

    metal_1 = elements[0].upper()
    metal_2 = elements[1].upper()

    components = [metal_1, metal_2, "B"]
    phases = ["HEDB1_2"]

    # Metal-site fraction of metal_2 in (metal_1_(1-x) metal_2_x)B2
    x_metal_site = np.linspace(1e-5, 1.0 - 1e-5, 201)

    # In MB2, metals are 1/3 of all atoms.
    x_total_metal_2 = x_metal_site / 3.0

    temperatures = [300, 1000, 2000, 3000, 4000]

    fig, ax = plt.subplots(figsize=(10, 7))

    for temperature in temperatures:
        conditions = {
            v.N: 1,
            v.P: 101325,
            v.T: temperature,
            v.X("B"): 2 / 3,
            v.X(metal_2): x_total_metal_2,
        }

        result = equilibrium(
            tdb,
            components,
            phases,
            conditions,
            output="GM",
        )

        gm = np.squeeze(result.GM.values)

        ax.plot(
            x_metal_site,
            gm,
            linewidth=2,
            label=f"{temperature} K",
        )

    ax.set_xlabel(
        rf"Metal-site fraction $x_{{{metal_2}}}$ in "
        rf"$({metal_1}_{{1-x}}{metal_2}_x)B_2$",
        fontsize=15,
    )
    ax.set_ylabel("Molar Gibbs Energy (J/mol-atom)", fontsize=15)
    ax.set_title(
        rf"Gibbs Energy of $({metal_1},{metal_2})B_2$ HEDB1_2 Phase",
        fontsize=17,
    )
    ax.tick_params(axis="both", labelsize=13)
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.25)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved Gibbs-energy plot to: {output_path}")

# Build filtered composition list (strip fixed-sublattice elements like B)
remove_elements = {item["element"] for item in sqsgen_levels2}
filt_comp_list = [
    [el for el in comp if el not in remove_elements]
    for comp in composition_list
]

for comp, comp_filt in zip(composition_list, filt_comp_list):
    comp_name = "".join(comp_filt)
    comp_dir = path1 / comp_name
    if not comp_dir.exists():
        continue
    os.chdir(comp_dir)
    elements = [el.upper() for el in comp]
    for tdb_file in comp_dir.glob("*.tdb"):
        tdb = Database(str(tdb_file))
        plot_phase_diagram(tdb, elements, len(comp), comp_dir / f"{comp_name}_Phase_Diagram.png")

os.chdir(path0)

# ------------------------------------------------------------------
# 5. Visualize
# ------------------------------------------------------------------
viz = BladeVisualizer()

# Combine phase diagram PNGs
pngs = []
for comp_filt in filt_comp_list:
    comp_dir = path1 / "".join(comp_filt)
    pngs.extend(comp_dir.glob("*_Phase_Diagram.png"))

if pngs:
    out_pd = path1 / "Combined_Phase_Diagrams.png"
    viz.phase_diagram(pngs, save=out_pd)
    print(f"Combined phase diagrams → {out_pd}")

# Combine CONTCAR structures for each phase in every composition
for comp_filt in filt_comp_list:
    comp_name = "".join(comp_filt)
    comp_dir = path1 / comp_name
    if not comp_dir.exists():
        continue
    for phase_dir in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
        contcars = sorted(phase_dir.glob("sqs_lev=*/CONTCAR"))
        if not contcars:
            continue
        out = comp_dir / f"Combined_CONTCARs_{comp_name}_{phase_dir.name}.png"
        viz.contcar(contcars, save=out)
        print(f"Saved combined CONTCARs → {out}")

# Relaxation movies (requires trajectory files from TDB fitting)
if run_movie:
    viz.make_combined_relaxation_movie(
        composition_list=filt_comp_list,
        path1=path1,
        traj_name="relaxation_live.xyz",
        fps=10,
    )
