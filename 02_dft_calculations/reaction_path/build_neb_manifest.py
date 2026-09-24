#!/usr/bin/env python3
"""Build the canonical CI-NEB pathway manifest from the archived XYZ profiles.

The trajectory filenames define the pathway identifiers:

    W_MM_B_NN_<initial-state>_<pathway>_movie_contcar.xyz

Each XYZ comment line contains the image energy relative to the initial image.
The script combines those relative energies with the absolute initial-state
energy retained in NEB_wb.csv, derives every other energetic quantity, and
classifies the transferred H atom and its relaxed final W/B anchor directly
from the endpoint geometries.

Run from this directory with:

    python build_neb_manifest.py
"""

from __future__ import annotations

import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean, median, stdev


HERE = Path(__file__).resolve().parent
TRAJECTORY_DIR = HERE / "neb_trajectories_xyz"
MANIFEST = HERE / "NEB_wb.csv"
SUMMARY = HERE / "NEB_sampling_by_cluster.csv"
TABLE_ROWS_TEX = HERE / "NEB_table_rows.tex"
SAMPLING_ROWS_TEX = HERE / "NEB_sampling_rows.tex"

FILENAME_RE = re.compile(
    r"^(W_(\d+)_B_(\d+))_(\d+)_(\d+)_movie_contcar\.xyz$"
)
ENERGY_RE = re.compile(r"\bE:\s*([-+0-9.eE]+)")


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def subtract(a: tuple[float, float, float], b: tuple[float, float, float]):
    return tuple(x - y for x, y in zip(a, b))


def norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(sum(x * x for x in a))


def read_xyz(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    frames = []
    energies = []
    cursor = 0
    while cursor < len(lines):
        n_atoms = int(lines[cursor].strip())
        match = ENERGY_RE.search(lines[cursor + 1])
        if match is None:
            raise ValueError(f"Missing relative energy in {path.name}, frame {len(frames)}")
        energies.append(float(match.group(1)))
        atoms = []
        for line in lines[cursor + 2 : cursor + 2 + n_atoms]:
            fields = line.split()
            atoms.append((fields[0], tuple(map(float, fields[1:4]))))
        frames.append(atoms)
        cursor += n_atoms + 2
    if cursor != len(lines):
        raise ValueError(f"Malformed XYZ record in {path.name}")
    if len(frames) != 10:
        raise ValueError(f"Expected 10 images in {path.name}; found {len(frames)}")
    if any([atom[0] for atom in frame] != [atom[0] for atom in frames[0]] for frame in frames):
        raise ValueError(f"Atom order changes between frames in {path.name}")
    return frames, energies


def read_initial_energies(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    energies: dict[tuple[str, int], list[float]] = defaultdict(list)
    for row in rows:
        energies[(row["Cluster"], int(row["Structure_ID"]))].append(
            float(row["E_IS (eV)"])
        )
    result = {}
    for key, values in energies.items():
        if max(values) - min(values) > 0.002:
            raise ValueError(f"Inconsistent E_IS values for {key}: {values}")
        result[key] = mean(values)
    return result


def endpoint_metadata(frames):
    initial = frames[0]
    final = frames[-1]
    symbols = [atom[0] for atom in initial]
    hosts = [i for i, symbol in enumerate(symbols) if symbol in {"W", "B"}]
    hydrogens = [i for i, symbol in enumerate(symbols) if symbol == "H"]
    nitrogens = [i for i, symbol in enumerate(symbols) if symbol == "N"]
    if len(nitrogens) != 1 or len(hydrogens) != 3:
        raise ValueError(f"Expected one N and three H atoms; found N={len(nitrogens)}, H={len(hydrogens)}")

    host_displacements = [
        subtract(final[i][1], initial[i][1]) for i in hosts
    ]
    drift = tuple(median(values) for values in zip(*host_displacements))
    corrected_h_displacements = {
        i: norm(subtract(subtract(final[i][1], initial[i][1]), drift))
        for i in hydrogens
    }
    transferred_h = max(hydrogens, key=corrected_h_displacements.get)

    n_atom = nitrogens[0]
    binding_w = min(
        (i for i in hosts if symbols[i] == "W"),
        key=lambda i: distance(initial[n_atom][1], initial[i][1]),
    )
    final_anchor = min(
        hosts,
        key=lambda i: distance(final[transferred_h][1], final[i][1]),
    )

    h_local = hydrogens.index(transferred_h) + 1
    anchor_element = symbols[final_anchor]
    anchor_local = sum(1 for i in range(final_anchor + 1) if symbols[i] == anchor_element)
    binding_w_local = sum(1 for i in range(binding_w + 1) if symbols[i] == "W")

    candidates = [i for i in hosts if i != binding_w]
    ranked = sorted(
        candidates,
        key=lambda i: distance(initial[transferred_h][1], initial[i][1]),
    )
    initial_rank = 0 if final_anchor == binding_w else ranked.index(final_anchor) + 1
    h_to_anchor_initial = distance(initial[transferred_h][1], initial[final_anchor][1])
    h_to_anchor_final = distance(final[transferred_h][1], final[final_anchor][1])
    binding_w_to_anchor_initial = distance(initial[binding_w][1], initial[final_anchor][1])
    closest_h = min(
        hydrogens,
        key=lambda i: distance(initial[i][1], initial[final_anchor][1]),
    )

    return {
        "Transferred_H": f"H{h_local}",
        "Initial_NH3_W_anchor": f"W{binding_w_local}",
        "Final_H_anchor": f"{anchor_element}{anchor_local}",
        "Final_anchor_element": anchor_element,
        "Final_anchor_rank_initial_excluding_binding_W": initial_rank,
        "Transferred_H_is_closest_H_to_final_anchor_initially": (
            "yes" if closest_h == transferred_h else "no"
        ),
        "d_H_anchor_initial (A)": h_to_anchor_initial,
        "d_H_anchor_final (A)": h_to_anchor_final,
        "d_bindingW_anchor_initial (A)": binding_w_to_anchor_initial,
    }


def build_rows(initial_energies):
    rows = []
    for path in sorted(TRAJECTORY_DIR.glob("*.xyz")):
        match = FILENAME_RE.fullmatch(path.name)
        if match is None:
            raise ValueError(f"Unexpected trajectory filename: {path.name}")
        cluster, w_count, b_count, structure_id, pathway_id = match.groups()
        structure_id = int(structure_id)
        pathway_id = int(pathway_id)
        frames, relative_energies = read_xyz(path)
        e_is = initial_energies[(cluster, structure_id)]
        delta_e = relative_energies[-1]
        e_act_fwd = max(relative_energies)
        ts_image = relative_energies.index(e_act_fwd)
        e_act_rev = e_act_fwd - delta_e
        metadata = endpoint_metadata(frames)
        rows.append(
            {
                "Cluster": cluster,
                "W": int(w_count),
                "B": int(b_count),
                "Structure_ID": structure_id,
                "Pathway_ID": pathway_id,
                "N_images": len(frames),
                "TS_image": ts_image,
                "E_IS (eV)": e_is,
                "E_FS (eV)": e_is + delta_e,
                "E_TS (eV)": e_is + e_act_fwd,
                "E_act_fwd (eV)": e_act_fwd,
                "E_act_rev (eV)": e_act_rev,
                "Delta_E (eV)": delta_e,
                **metadata,
                "Trajectory_file": path.name,
            }
        )
    rows.sort(key=lambda r: (r["W"], r["B"], r["Structure_ID"], r["Pathway_ID"]))
    return rows


def write_manifest(rows):
    fieldnames = list(rows[0])
    with MANIFEST.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            formatted = {
                key: (f"{value:.6f}" if isinstance(value, float) else value)
                for key, value in row.items()
            }
            writer.writerow(formatted)


def write_summary(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["W"], row["B"], row["Cluster"])].append(row)
    fieldnames = [
        "Cluster",
        "N_initial_states",
        "N_pathways",
        "N_H_on_B",
        "N_H_on_W",
        "E_act_fwd_mean (eV)",
        "E_act_fwd_SD (eV)",
        "E_act_fwd_min (eV)",
        "E_act_fwd_max (eV)",
        "E_act_rev_mean (eV)",
        "Delta_E_min (eV)",
        "Delta_E_max (eV)",
    ]
    with SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for (_, _, cluster), group in sorted(grouped.items()):
            fwd = [r["E_act_fwd (eV)"] for r in group]
            rev = [r["E_act_rev (eV)"] for r in group]
            delta = [r["Delta_E (eV)"] for r in group]
            writer.writerow(
                {
                    "Cluster": cluster,
                    "N_initial_states": len({r["Structure_ID"] for r in group}),
                    "N_pathways": len(group),
                    "N_H_on_B": sum(r["Final_anchor_element"] == "B" for r in group),
                    "N_H_on_W": sum(r["Final_anchor_element"] == "W" for r in group),
                    "E_act_fwd_mean (eV)": f"{mean(fwd):.6f}",
                    "E_act_fwd_SD (eV)": f"{stdev(fwd):.6f}" if len(fwd) > 1 else "",
                    "E_act_fwd_min (eV)": f"{min(fwd):.6f}",
                    "E_act_fwd_max (eV)": f"{max(fwd):.6f}",
                    "E_act_rev_mean (eV)": f"{mean(rev):.6f}",
                    "Delta_E_min (eV)": f"{min(delta):.6f}",
                    "Delta_E_max (eV)": f"{max(delta):.6f}",
                }
            )


def latex_cluster(w_count: int, b_count: int) -> str:
    return f"\\wbcluster{{{w_count}}}{{{b_count}}}"


def write_latex_rows(rows):
    table_lines = []
    for row in rows:
        table_lines.append(
            " & ".join(
                [
                    latex_cluster(row["W"], row["B"]),
                    str(row["Structure_ID"]),
                    str(row["Pathway_ID"]),
                    row["Transferred_H"],
                    row["Final_H_anchor"],
                    f"{row['E_IS (eV)']:.3f}",
                    f"{row['E_FS (eV)']:.3f}",
                    f"{row['E_TS (eV)']:.3f}",
                    f"{row['E_act_fwd (eV)']:.3f}",
                    f"{row['E_act_rev (eV)']:.3f}",
                    f"{row['Delta_E (eV)']:.3f}",
                ]
            )
            + r" \\"
        )
    TABLE_ROWS_TEX.write_text("\n".join(table_lines) + "\n", encoding="utf-8")

    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["W"], row["B"])].append(row)
    summary_lines = []
    for (w_count, b_count), group in sorted(grouped.items()):
        fwd = [r["E_act_fwd (eV)"] for r in group]
        rev = [r["E_act_rev (eV)"] for r in group]
        delta = [r["Delta_E (eV)"] for r in group]
        sd = f"{stdev(fwd):.3f}" if len(fwd) > 1 else "--"
        n_b = sum(r["Final_anchor_element"] == "B" for r in group)
        n_w = sum(r["Final_anchor_element"] == "W" for r in group)
        summary_lines.append(
            " & ".join(
                [
                    latex_cluster(w_count, b_count),
                    str(len({r["Structure_ID"] for r in group})),
                    str(len(group)),
                    f"{n_b}/{n_w}",
                    f"{mean(fwd):.3f}",
                    sd,
                    f"{min(fwd):.3f}",
                    f"{max(fwd):.3f}",
                    f"{mean(rev):.3f}",
                    f"${min(delta):+.3f}$ to ${max(delta):+.3f}$",
                ]
            )
            + r" \\"
        )
    SAMPLING_ROWS_TEX.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def validate(rows):
    keys = [(r["Cluster"], r["Structure_ID"], r["Pathway_ID"]) for r in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate cluster/initial-state/pathway identifiers")
    for row in rows:
        if abs((row["E_TS (eV)"] - row["E_IS (eV)"]) - row["E_act_fwd (eV)"]) > 1e-9:
            raise ValueError(f"Forward-barrier inconsistency in {row['Trajectory_file']}")
        if abs((row["E_TS (eV)"] - row["E_FS (eV)"]) - row["E_act_rev (eV)"]) > 1e-9:
            raise ValueError(f"Reverse-barrier inconsistency in {row['Trajectory_file']}")
        if abs((row["E_FS (eV)"] - row["E_IS (eV)"]) - row["Delta_E (eV)"]) > 1e-9:
            raise ValueError(f"Reaction-energy inconsistency in {row['Trajectory_file']}")


def main():
    initial_energies = read_initial_energies(MANIFEST)
    rows = build_rows(initial_energies)
    validate(rows)
    write_manifest(rows)
    write_summary(rows)
    write_latex_rows(rows)
    print(
        f"Wrote {len(rows)} internally consistent pathways over "
        f"{len({r['Cluster'] for r in rows})} compositions."
    )


if __name__ == "__main__":
    main()
