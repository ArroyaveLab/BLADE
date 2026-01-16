#%conda activate materialsframework-main
#%conda activate materialsframework-main

from materialsframework.calculators import GraceCalculator as Calculator
from materialsframework.tools.sqs2tdb import Sqs2tdb
import matplotlib.pyplot as plt
from pycalphad import Database, binplot, ternplot
from pycalphad import variables as v, binplot
import pycalphad.variables as v
import torch
import os
import itertools
import subprocess
import fractions
import math
from pathlib import Path
import threading
import shutil
import re
from fractions import Fraction

N_THREADS = 8
os.environ["OMP_NUM_THREADS"] = str(N_THREADS)
os.environ["MKL_NUM_THREADS"] = str(N_THREADS)
os.environ["OPENBLAS_NUM_THREADS"] = str(N_THREADS)
torch.set_num_threads(N_THREADS)

# Specify settings

# Specify elements and system size (Total # elements)
transition_metals = ["Ti", "Zr", "Hf", "V", "Nb", "Ta", "Cr", "Mo", "W"]
#transition_metals = ["Zr", "Hf"]
rare_earths = ["Sc", "Y", "La", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]
system_size = 2
tm_element_range = [2, 2]
re_element_range = [0, 0]
allow_lower_order = True
allow_lower_order = False

phase = 'HEDB1'
path = '/Users/chasekatz/Desktop/School/Research'
path1 = os.path.join(path, 'PhaseForge/PhaseForge/atat/data/sqsdb/')
path2 = os.path.join(path, 'MF_tests/p8/')
level = 6
time = 30

a = 1
b = 1
c = 1.63299
alpha = 90
beta = 90
gamma = 120

rndstr1 = f"""
{a} {b} {c} {alpha} {beta} {gamma}
{a} 0 0
{a/2} {a*math.sqrt(3)/2} 0
0 0 {c}
"""

rndstr2 = """
0.000000 0.000000 0.000000  a
0.333333 0.666667 0.000000  a
0.666667 0.333333 0.000000  a
0.500000 0.000000 0.000000  a
0.000000 0.500000 0.000000  a
0.500000 0.500000 0.000000  a
0.250000 0.750000 0.000000  a
0.166667 0.833333 0.5       B
0.833333 0.166667 0.5       B
0.250000 0.250000 0.5       B
0.750000 0.750000 0.5       B
0.600000 0.300000 0.5       B
0.300000 0.600000 0.5       B
0.000000 0.000000 1.000000  a
0.333333 0.666667 1.000000  a
0.666667 0.333333 1.000000  a
0.500000 0.000000 1.000000  a
0.000000 0.500000 1.000000  a
0.500000 0.500000 1.000000  a
0.250000 0.750000 1.000000  a
"""

sqsgen_levels = [
"""level=0         a=1""", 
"""level=1         a=0.5,0.5""",
"""level=2         a=0.75,0.25""",
"""level=3         a=0.33333,0.33333,0.33333""",
"""level=4         a=0.5,0.25,0.25""",
"""level=5         a=0.875,0.125\nlevel=5         a=0.625,0.375""",
"""level=6         a=0.75,0.125,0.125"""]

sqsgen = ""
for i in range((level+1)):
    sqsgen += sqsgen_levels[i] + "\n"

rndstr = rndstr1.strip() + "\n" + rndstr2.strip()

# Generate compositions

# Generate all possible transition metal compositions
tm_elements = list(range(tm_element_range[0], tm_element_range[1] + 1))
re_elements = list(range(re_element_range[0], re_element_range[1] + 1))

tm_combos = []
re_combos = []
combined_comps = []
compositions = []

# Generate all possible combinations of transition metals
for i in tm_elements:
    # If only transition metals are 0, add all rare earth combinations
    if i == 0:
        for j in re_elements:
            # If both are 0, only add Boron
            if j == 0:
                compositions += [[""]]
            else:
                combined_comps += [list(c) for c in itertools.combinations(rare_earths, j)]
    if i != 0:
        tm_combos += [list(c) for c in itertools.combinations(transition_metals, i)]

# Generate all possible combinations of rare earths
for j in re_elements:
    # If rare earths are 0, add all transition metal combinations
    if j == 0:
        combined_comps += tm_combos
    else:
        re_combos += [list(c) for c in itertools.combinations(rare_earths, j)]

# Combine transition metal and rare earth combinations
for tm_comp in tm_combos:
    for re_comp in re_combos:
        combined_comps.append(list(tm_comp) + re_comp)

# Remove compositions that exceed system size
combined_comps = [c for c in combined_comps if len(c) <= system_size]

# Add Boron to each combination, remove lower order compositions if needed, and sort
for i in combined_comps:
    compositions += [sorted(i)]
if allow_lower_order == False:
    compositions = [c for c in compositions if len(c) == system_size]
compositions = sorted(compositions)
print(compositions)
print(f'Total compositions: {len(compositions)}')

# Ensure no duplicate compositions
unique_comps = {tuple(sorted(c)) for c in compositions}
print("Unique:", len(unique_comps))

# Ensure no duplicate compositions
len_comps = []
for i in compositions:
    len_comps += [len(i)]
print(len_comps)
if len(len_comps) >= 2:
    unique_len_comps = set(len_comps)
else:
    unique_len_comps = {len_comps[0]}
print("Unique Systems:", unique_len_comps)



def make_stopsqs(sqsdir, time, comp_vals):
    if not completed.is_set():
        Path(sqsdir, "stopsqs").touch()
        print(f"touch stopsqs in {comp_vals} after {time} seconds")


def exact_supercell_size(fractions, base_sites):
    # Convert floats to exact rational fractions
    fracs = [Fraction(str(f)).limit_denominator() for f in fractions]

    # Find LCM of denominators of composition fractions
    denoms = [f.denominator for f in fracs]
    denom_lcm = math.lcm(*denoms)

    # Supercell must be multiple of both base_sites and denom_lcm
    supercell_size = math.lcm(base_sites, denom_lcm)

    # Exact integer atom counts
    counts = [int(supercell_size * f) for f in fracs]

    return supercell_size, counts

# Generate sqs
for length in unique_len_comps:
    dir_name = os.path.join(path1, (phase+"_"+str(length)))
    os.makedirs(dir_name, exist_ok=True)
    file_path = os.path.join(dir_name, "rndstr.skel")
    with open(file_path, "w") as f:
        f.write(rndstr)
    print(f"File created at: {file_path}")

    dir_name = os.path.join(path1, (phase+"_"+str(length)))
    os.makedirs(dir_name, exist_ok=True)
    file_path = os.path.join(dir_name, "sqsgen.in")
    with open(file_path, "w") as f:
        f.write(sqsgen)
    print(f"File created at: {file_path}")

    cmd = ['sqs2tdb', '-mk']
    # Run the command in the target directory
    result = subprocess.run(cmd, cwd=dir_name, capture_output=True, text=True)
    # Output results (optional)
    print(result.stdout)
    if result.stderr:
        print("Error:", result.stderr)

    parent_dir = Path(os.path.join(path1, (phase+"_"+str(length))))
    for sqsdir in parent_dir.glob('sqsdb_lev=*/'):
        folder_name = sqsdir.name
        folder_lev = folder_name.split('=')[1]
        folder_lev = folder_lev.split('_')[0]
        comp_str = folder_name.split('=')[-1]   # e.g. "0.75,0.125,0.125"
        comp_vals = [float(x) for x in comp_str.split(',')]
        # Skip pure-species (end-member) folders
        if len(comp_vals) == 1:    # Only one species, no mixing
            print(comp_vals)
            print(f"Skipping pure species directory: {sqsdir}")
            continue
        base_sites = 14
        N_atoms, counts = exact_supercell_size(comp_vals, base_sites)
        if N_atoms>=100:
            N_atoms = base_sites*len(comp_vals)
        N_atoms += 6
        N_atoms = ((N_atoms + 19) // 20) * 20

        print(f"Running corrdump in {comp_vals}")
        try:
            subprocess.run(
                [
                    "corrdump",
                    "-l=rndstr.in", "-ro", "-noe", "-nop", "-clus",
                    "-2=1,2,3", "-3=1"
                ],
                cwd=sqsdir, check=True
            )

            print(f"Running mcsqs with {N_atoms} atoms in {comp_vals}")

            # Start the timer BEFORE you launch mcsqs
            completed = threading.Event()
            timer = threading.Timer(time, make_stopsqs(sqsdir, time, comp_vals))
            timer.start()

            try:
                subprocess.run(
                    ["mcsqs", f"-n={N_atoms}"],
                    cwd=sqsdir,
                    check=True
                )
                completed.set()  # Mark process as completed, so timer won't fire
            except subprocess.CalledProcessError:
                print(f"corrdump or mcsqs failed in {sqsdir}, continuing.")
            finally:
                timer.cancel()
        except subprocess.CalledProcessError:
            print(f"corrdump or mcsqs failed in {sqsdir}, continuing.")
            continue


# Create phase diagrams

# Loop through each composition
for comp in compositions:
    # Create and move to output folder
    file = f'{path2}{"".join(comp)}'
    os.makedirs(file, exist_ok=True)

    main_dir = os.getcwd()
    os.chdir(file)

    # Specify inputs
    comp_leng = len(comp)
    #phases = [(phase+"_"+str(length)), 'LIQUID']
    phases = [(phase+"_"+str(length))]

    # Initialize the wrapper
    #sqs = Sqs2tdb(fmax=0.001, verbose=True, calculator=Calculator(device="cpu"), md_timestep=100, md_temperature=2000)
    sqs = Sqs2tdb(fmax=0.001, verbose=True, calculator=Calculator(device="cpu"))

    # Fit SQS for Ni-Re system on FCC_A1 and HCP_A3 lattices
    sqs.fit(
        species=comp,
        lattices=phases,
        level=level,
        phonon=True#,
        #T_max=4000,
    )

    # Editing names
    elements = [el.upper() for el in comp]
    file_names = ["_".join(p) for p in itertools.permutations(elements)]
    print(file_names)

    # Load database and choose the phases that will be considered
    for file in file_names:
        if os.path.isfile(f'{file}.tdb'):
            tdb = Database(f'{file}.tdb')
    
            fig = plt.figure(figsize=(9,7))
            axes = fig.gca()

            # Compute the phase diagram and plot it
            print(elements)
            print(phases)
            print(tdb)
            binplot(
                tdb, elements, phases,
                {v.X(elements[0]): (0, 1, 0.02), v.T: (300, 4000, 10), v.P: 101325, v.N: 1},
                plot_kwargs={'ax': axes}
            )
            plt.tight_layout()
            plt.savefig(f'{file}_Phase_Diagram.png', dpi=300)

    # Change back to main folder
    os.chdir(main_dir)
