"""Example: structure and phase-diagram visualization with BladeVisualizer.

This script shows how to use BladeVisualizer for three common tasks:

1. Combining CONTCAR structures from multiple SQS levels into one image.
2. Stitching phase-diagram PNGs side by side.
3. Building a relaxation movie from trajectory files (requires ffmpeg for MP4).

Prerequisites:
    - matplotlib, ase, Pillow, and imageio must be installed.
    - CONTCAR / PNG files must exist at the paths used below.
    - ffmpeg is optional (needed for MP4 output only).

    python examples/04_visualization.py
"""

from pathlib import Path

from blade.analysis.blade_visual import BladeVisualizer

viz = BladeVisualizer()

# ------------------------------------------------------------------
# Paths — adjust to your environment
# ------------------------------------------------------------------
path1 = Path("/Users/chasekatz/Desktop/School/Research/BLADE")

# ------------------------------------------------------------------
# Example 1: combine CONTCAR structures for one composition + phase
# ------------------------------------------------------------------
print("=" * 60)
print("Example 1: Combined CONTCAR structures")
print("=" * 60)

comp_name = "CrHfTa"
phase_name = "HEDB1_3"
comp_dir = path1 / comp_name
phase_dir = comp_dir / phase_name

contcars = sorted(phase_dir.glob("sqs_lev=*/CONTCAR"))
if contcars:
    out = comp_dir / f"Combined_CONTCARs_{comp_name}_{phase_name}.png"
    viz.contcar(contcars, save=out)
    print(f"Saved {len(contcars)} structures -> {out}")
else:
    print(f"No CONTCARs found in {phase_dir} — skipping.")

# ------------------------------------------------------------------
# Example 2: combine phase diagram PNGs
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 2: Combined phase diagrams")
print("=" * 60)

pngs = sorted(path1.glob("**/*_Phase_Diagram.png"))
if pngs:
    out_pd = path1 / "Combined_Phase_Diagrams.png"
    viz.phase_diagram(pngs, save=out_pd)
    print(f"Combined {len(pngs)} phase diagrams -> {out_pd}")
else:
    print("No *_Phase_Diagram.png files found — skipping.")

# ------------------------------------------------------------------
# Example 3: relaxation movie (GIF + optional MP4)
# ------------------------------------------------------------------
print()
print("=" * 60)
print("Example 3: Relaxation movie")
print("=" * 60)

composition_list = [["Cr", "Hf", "Ta"]]

viz.make_combined_relaxation_movie(
    composition_list=composition_list,
    path1=path1,
    traj_name="relaxation_live.xyz",
    fps=10,
)
