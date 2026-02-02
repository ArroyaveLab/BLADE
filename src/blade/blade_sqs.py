import math
import subprocess
import threading
from fractions import Fraction
from pathlib import Path


def supercell_size(self, fractions, min_a_sites=16, max_den=96):
    # safer rationalization
    fracs = [Fraction(float(f)).limit_denominator(max_den) for f in fractions]
    denom_lcm = math.lcm(*[f.denominator for f in fracs])

    # Count sites in prototype
    a_sites = 0
    b_sites = 0
    sites = 0
    for line in self.unit_cell.splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue
        sites += 1
        if parts[-1] == "a":
            a_sites += 1
        else:
            b_sites += 1

    if a_sites == 0:
        raise ValueError("No 'a' sites found in unit_cell")

    # Smallest n_a that can represent fractions exactly
    n_a0 = math.lcm(a_sites, denom_lcm)

    # Enforce a minimum number of metal sites (quality)
    if n_a0 < min_a_sites:
        mult = math.ceil(min_a_sites / n_a0)
        n_a = n_a0 * mult
    else:
        n_a = n_a0

    # Number of prototype repeats
    k = n_a // a_sites

    # Total sites follow prototype ratio exactly
    n_total = k * (a_sites + b_sites)

    # Integer counts per species on a-sublattice
    counts = [int(n_a * f) for f in fracs]
    diff = n_a - sum(counts)
    counts[0] += diff

    return n_total, counts


class BladeSQS:
    def __init__(self, phases_dict, sqsgen_levels, level):
        self.a = phases_dict["a"]
        self.b = phases_dict["b"]
        self.c = phases_dict["c"]
        self.alpha = phases_dict["alpha"]
        self.beta = phases_dict["beta"]
        self.gamma = phases_dict["gamma"]
        self.unit_cell = phases_dict["coords"]
        self.sqsgen_levels = sqsgen_levels
        self.level = level

    def sqs_struct(self):
        rndstr1 = f"""
        {self.a} {self.b} {self.c} {self.alpha} {self.beta} {self.gamma}
        1 0 0
        0 1 0
        0 0 1
        """
        sqsgen = ""
        for i in range(self.level + 1):
            sqsgen += self.sqsgen_levels[i] + "\n"

        #unit_cell_multiplied = repeat_unit_cell(self.unit_cell, nx=2, ny=2, nz=2)
        rndstr = rndstr1.strip() + "\n" + self.unit_cell.strip()
        print(rndstr)
        self.sqsgen_text = sqsgen
        self.rndstr = rndstr

        return sqsgen, rndstr

    def supercell_size(self, fractions, min_a_sites=16, max_den=96):
        # safer rationalization
        fracs = [Fraction(float(f)).limit_denominator(max_den) for f in fractions]
        denom_lcm = math.lcm(*[f.denominator for f in fracs])

        # Count sites in prototype
        a_sites = 0
        b_sites = 0
        sites = 0
        for line in self.unit_cell.splitlines():
            parts = line.split()
            if len(parts) != 4:
                continue
            sites += 1
            if parts[-1] == "a":
                a_sites += 1
            else:
                b_sites += 1

        if a_sites == 0:
            raise ValueError("No 'a' sites found in unit_cell")

        # Smallest n_a that can represent fractions exactly
        n_a0 = math.lcm(a_sites, denom_lcm)

        # Enforce a minimum number of metal sites (quality)
        if n_a0 < min_a_sites:
            mult = math.ceil(min_a_sites / n_a0)
            n_a = n_a0 * mult
        else:
            n_a = n_a0

        # Number of prototype repeats
        k = n_a // a_sites

        # Total sites follow prototype ratio exactly
        n_total = k * (a_sites + b_sites)

        # Integer counts per species on a-sublattice
        counts = [int(n_a * f) for f in fracs]
        diff = n_a - sum(counts)
        counts[0] += diff

        return n_total, counts

    def sqs_gen(self, unique_len_comps, phase, path1, time):
        # Generate sqs
        for length in unique_len_comps:
            dir_name = Path(path1) / (phase + "_" + str(length))
            dir_name.mkdir(parents=True, exist_ok=True)
            file_path = dir_name / "rndstr.skel"
            sqsgen, rndstr = self.sqs_struct()
            with file_path.open("w") as f:
                f.write(rndstr)
            print(f"File created at: {file_path}")

            dir_name = Path(path1) / (phase + "_" + str(length))
            dir_name.mkdir(parents=True, exist_ok=True)
            file_path = dir_name / "sqsgen.in"
            with file_path.open("w") as f:
                f.write(sqsgen)
            print(f"File created at: {file_path}")

            cmd = ["sqs2tdb", "-mk"]
            # Run the command in the target directory
            result = subprocess.run(cmd, cwd=dir_name, capture_output=True, text=True, check=False)
            # Output results (optional)
            print(result.stdout)
            if result.stderr:
                print("Error:", result.stderr)

            parent_dir = Path(path1) / (phase + "_" + str(length))
            for sqsdir in parent_dir.glob("sqsdb_lev=*/"):
                folder_name = sqsdir.name
                folder_lev = folder_name.split("=")[1]
                folder_lev = folder_lev.split("_")[0]
                comp_str = folder_name.split("=")[-1]  # e.g. "0.75,0.125,0.125"
                fractions = [float(x) for x in comp_str.split(",")]
                # Skip pure-species (end-member) folders
                if len(fractions) == 1:  # Only one species, no mixing
                    print(fractions)
                    print(f"Skipping pure species directory: {sqsdir}")
                    continue
                n_atoms, counts = self.supercell_size(fractions)

                print(f"Running corrdump in {fractions}")
                try:
                    subprocess.run(
                        ["corrdump", "-l=rndstr.in", "-ro", "-noe", "-nop", "-clus", "-2=1,2,3", "-3=1"],
                        cwd=sqsdir,
                        check=True,
                    )

                    print(f"Running mcsqs with {n_atoms} atoms in {fractions}")
                    completed = threading.Event()

                    def make_stopsqs(sqsdir, time, fractions):
                        if not completed.is_set():
                            Path(sqsdir, "stopsqs").touch()
                            print(f"touch stopsqs in {fractions} after {time} seconds")

                    # Start the timer BEFORE you launch mcsqs
                    timer = threading.Timer(time, make_stopsqs, args=(sqsdir, time, fractions))
                    timer.start()

                    try:
                        subprocess.run(["mcsqs", f"-n={n_atoms}"], cwd=sqsdir, check=True)
                        completed.set()  # Mark process as completed, so timer won't fire
                    except subprocess.CalledProcessError:
                        print(f"corrdump or mcsqs failed in {sqsdir}, continuing.")
                    finally:
                        timer.cancel()
                except subprocess.CalledProcessError:
                    print(f"corrdump or mcsqs failed in {sqsdir}, continuing.")
                    continue
