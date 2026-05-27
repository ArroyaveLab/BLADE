"""Structure and phase-diagram visualization utilities.

This module provides :class:`BladeVisualizer`, which offers lightweight
helpers for inspecting BLADE workflow outputs:

- Side-by-side rendering of CONTCAR structures.
- Horizontal concatenation of phase-diagram images.
- Frame-by-frame relaxation movies exported as GIF and (optionally) MP4.

Example::

    from blade.analysis.blade_visual import BladeVisualizer
    from pathlib import Path

    viz = BladeVisualizer()

    # Combine several CONTCAR files into one image
    contcars = list(Path("MyComp/Phase1_3").glob("sqs_lev=*/CONTCAR"))
    viz.contcar(contcars, save="combined_structures.png")

    # Combine phase diagram PNGs side by side
    pngs = list(Path(".").glob("*_Phase_Diagram.png"))
    viz.phase_diagram(pngs, save="Combined_Phase_Diagrams.png")
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.io import read
from ase.visualize.plot import plot_atoms
from PIL import Image

if TYPE_CHECKING:
    pass

__author__ = "Chase Katz"


class BladeVisualizer:
    """Visualization utilities for BLADE structures and phase diagrams.

    Provides convenience methods for rendering structures and combining
    images. Intended for quick inspection and reporting rather than
    publication-quality figures.
    """

    def __init__(self) -> None:
        """Initialize BladeVisualizer."""

    def contcar(
        self,
        contcars: list[str | Path],
        save: str | Path | None = None,
    ) -> None:
        """Visualize one or more CONTCAR structures side by side.

        Each CONTCAR is read with inline element symbols (see
        :meth:`read_contcar_inline_symbols`) and rendered as a 2-D
        projection using ASE's plotting utilities.

        Args:
            contcars (list[str | Path]): Paths to CONTCAR files.
            save (str | Path | None, optional): Output path for the figure.
                If ``None``, the figure is shown interactively.
                Defaults to ``None``.
        """
        fig, axes = plt.subplots(1, len(contcars), figsize=(4 * len(contcars), 4))
        if len(contcars) == 1:
            axes = [axes]

        for ax, p in zip(axes, contcars):
            atoms = self.read_contcar_inline_symbols(p)
            plot_atoms(atoms, ax, rotation="65x,45y,0z", show_unit_cell=True)
            ax.set_axis_off()

        if save is not None:
            fig.savefig(save, dpi=200, bbox_inches="tight")
            plt.close(fig)
        else:
            plt.show()

    def phase_diagram(
        self,
        images: list[str | Path],
        save: str | Path,
    ) -> None:
        """Concatenate phase-diagram images horizontally into one file.

        All images are aligned at the top edge and padded with white
        background where heights differ.

        Args:
            images (list[str | Path]): Paths to image files to combine.
            save (str | Path): Output path for the combined image.
        """
        imgs = [Image.open(p) for p in images]
        widths, heights = zip(*(img.size for img in imgs))
        combined = Image.new("RGB", (sum(widths), max(heights)), "white")
        x_offset = 0
        for img in imgs:
            combined.paste(img, (x_offset, 0))
            x_offset += img.width
        combined.save(save)

    def read_contcar_inline_symbols(self, contcar_path: str | Path) -> Atoms:
        """Read a CONTCAR file with inline element symbols per atom line.

        Parses a VASP CONTCAR-style file where each atomic coordinate line
        ends with an explicit element symbol.  Supports both Direct and
        Cartesian coordinate formats.

        Args:
            contcar_path (str | Path): Path to the CONTCAR file.

        Returns:
            ase.Atoms: Atoms object with positions, cell, and periodic
            boundary conditions set.
        """
        with open(contcar_path) as f:
            lines = [ln.strip() for ln in f if ln.strip()]

        scale = float(lines[1])
        cell = (
            np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)], dtype=float)
            * scale
        )

        i = 5
        toks = lines[i].split()

        def _all_int(ts: list[str]) -> bool:
            try:
                [int(x) for x in ts]
                return True
            except ValueError:
                return False

        if _all_int(toks):
            counts = [int(x) for x in toks]
            i += 1
        else:
            i += 1
            counts = [int(x) for x in lines[i].split()]
            i += 1

        if lines[i].lower().startswith("s"):
            i += 1

        direct = lines[i].lower().startswith("d")
        i += 1

        n_atoms = sum(counts)
        symbols: list[str] = []
        positions: list[list[float]] = []
        for ln in lines[i : i + n_atoms]:
            parts = ln.split()
            positions.append([float(parts[0]), float(parts[1]), float(parts[2])])
            symbols.append(parts[-1])

        pos_arr = np.array(positions, dtype=float)
        if direct:
            pos_arr = pos_arr @ cell

        return Atoms(symbols=symbols, positions=pos_arr, cell=cell, pbc=True)

    def write_inline_contcar_from_ase(self, atoms: Atoms, out_path: str | Path) -> None:
        """Write an ASE Atoms object as a CONTCAR with inline element symbols.

        Each coordinate line is appended with the element symbol of that
        atom, which is required by :meth:`read_contcar_inline_symbols`.

        Args:
            atoms (ase.Atoms): Atoms object to serialize.
            out_path (str | Path): Output file path.
        """
        out_path = Path(out_path)
        cell = atoms.get_cell()
        frac_coords = atoms.get_scaled_positions()
        symbols = atoms.get_chemical_symbols()

        species_order: list[str] = []
        for s in symbols:
            if s not in species_order:
                species_order.append(s)
        counts = [symbols.count(s) for s in species_order]

        with open(out_path, "w") as f:
            f.write("Generated relaxation frame\n")
            f.write("1.0\n")
            for row in cell:
                f.write(" ".join(f"{x:.16f}" for x in row) + "\n")
            f.write(" ".join(species_order) + "\n")
            f.write(" ".join(str(c) for c in counts) + "\n")
            f.write("Direct\n")
            for frac, symbol in zip(frac_coords, symbols):
                f.write(f"{frac[0]:.16f} {frac[1]:.16f} {frac[2]:.16f} {symbol}\n")

    def make_combined_relaxation_movie(
        self,
        composition_list: list[list[str]],
        path1: str | Path,
        traj_name: str = "relaxation_live.xyz",
        fps: int = 5,
    ) -> None:
        """Generate per-composition relaxation movies as GIF (and MP4 if ffmpeg is available).

        For each composition and each phase sub-directory, reads trajectory
        files from ``sqs_lev=*`` folders, renders every relaxation frame as
        a side-by-side CONTCAR image, and assembles a GIF.  If ``ffmpeg``
        is on the system PATH, the GIF is also converted to MP4.

        Args:
            composition_list (list[list[str]]): List of chemical systems,
                each given as a list of element symbols.
            path1 (str | Path): BLADE staging directory containing
                per-composition sub-directories.
            traj_name (str, optional): Trajectory filename to search for
                inside each ``sqs_lev=*`` folder.
                Defaults to ``"relaxation_live.xyz"``.
            fps (int, optional): Frames per second for the output movie.
                Defaults to ``5``.
        """
        path1 = Path(path1)

        for comp in composition_list:
            comp_name = "".join(comp)
            comp_dir = path1 / comp_name
            if not comp_dir.exists():
                continue

            for phase_dir in sorted(p for p in comp_dir.iterdir() if p.is_dir()):
                traj_files = sorted(phase_dir.glob(f"sqs_lev=*/{traj_name}"))
                if not traj_files:
                    print(f"No trajectories found in {phase_dir}")
                    continue

                print(f"\nFound {len(traj_files)} trajectories in {phase_dir}")

                all_trajs: list[list] = []
                labels: list[str] = []
                for traj_file in traj_files:
                    try:
                        frames = read(traj_file, index=":")
                        if frames:
                            all_trajs.append(frames)
                            labels.append(traj_file.parent.name)
                    except Exception as exc:
                        print(f"Skipping {traj_file}: {exc}")

                if not all_trajs:
                    continue

                n_frames = max(len(frames) for frames in all_trajs)
                movie_dir = comp_dir / f"movie_frames_{comp_name}_{phase_dir.name}"
                movie_dir.mkdir(parents=True, exist_ok=True)

                frame_paths: list[Path] = []
                for frame_idx in range(n_frames):
                    step_dir = movie_dir / f"step_{frame_idx:04d}"
                    step_dir.mkdir(parents=True, exist_ok=True)

                    step_contcars: list[Path] = []
                    for traj_i, frames in enumerate(all_trajs):
                        out_contcar = step_dir / f"{labels[traj_i]}_step_{frame_idx:04d}.out"
                        self.write_inline_contcar_from_ase(
                            frames[min(frame_idx, len(frames) - 1)],
                            out_contcar,
                        )
                        step_contcars.append(out_contcar)

                    out_png = movie_dir / f"frame_{frame_idx:04d}.png"
                    self.contcar(step_contcars, save=out_png)
                    frame_paths.append(out_png)
                    print(f"Saved frame {frame_idx + 1}/{n_frames}: {out_png}")

                out_gif = comp_dir / f"Combined_relaxation_{comp_name}_{phase_dir.name}.gif"
                imageio.mimsave(
                    out_gif,
                    [imageio.imread(p) for p in frame_paths],
                    fps=fps,
                )
                print(f"Saved GIF -> {out_gif}")

                if shutil.which("ffmpeg") is not None:
                    out_mp4 = comp_dir / f"Combined_relaxation_{comp_name}_{phase_dir.name}.mp4"
                    subprocess.run(
                        [
                            "ffmpeg", "-y",
                            "-i", str(out_gif),
                            "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
                            "-movflags", "faststart",
                            "-pix_fmt", "yuv420p",
                            str(out_mp4),
                        ],
                        check=True,
                    )
                    print(f"Saved MP4 -> {out_mp4}")
                else:
                    print("ffmpeg not found. Install with: brew install ffmpeg or pixi add ffmpeg")
