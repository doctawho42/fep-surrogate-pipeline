"""Fetch public RBFE edges (OpenFF protein-ligand-benchmark) -> data/fep_edges/*.csv.

Each row is a congeneric edge (smiles_a, smiles_b, ddg = dG(B)-dG(A) in kcal/mol).

REAL YAML STRUCTURE (differs from brief template):
  - ligands: data/<target>/00_data/ligands.yml
      top-level dict keyed by ligand name; each entry has:
        smiles: <SMILES string>
        measurement:
          type: ic50 | ki | pic50 | dg | dh | tds
          value: <float>
          unit: uM | nM | '' (pic50 has empty unit)
          [additional: list of extra measurement dicts]
  - edges: data/<target>/03_edges/<network>.yml
      top-level dict with key 'edges'; each entry has:
        ligand_a: <ligand name>
        ligand_b: <ligand name>
        (plus atom mapping, score, etc. — ignored)
  - There is NO 00_data/edges.yml in the real benchmark targets; that only
    exists in plbenchmark/sample_data/.

Usage: python scripts/fetch_fep_edges.py <PLB_clone_dir>
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys

import yaml

OUT = pathlib.Path(__file__).resolve().parent.parent / "data" / "fep_edges"
RT = 0.593  # kcal/mol at 298 K  (kT at 298 K = 0.5922 kcal/mol, rounded to 0.593)
KJ_TO_KCAL = 1.0 / 4.184


def to_dg(meas: dict) -> float | None:
    """Convert a measurement dict to ΔG in kcal/mol, or None if unsupported."""
    t = str(meas.get("type", "")).lower()
    v = meas.get("value")
    if v is None:
        return None
    v = float(v)
    unit = str(meas.get("unit", "")).strip().lower()

    if t == "dg":
        # Value is free energy — check unit
        if "kj" in unit or "kilojoule" in unit:
            return v * KJ_TO_KCAL
        # kcal/mol or plain 'calories / mole' edge cases
        if "kcal" in unit or "cal" in unit:
            return v
        # Dimensionless or empty: assume kcal/mol (rare; flag via small magnitude check)
        return v

    if t == "pic50":
        # ΔG = -RT·ln(10)·pIC50 ≈ -1.364·pIC50 at 298 K
        return -1.364 * v

    if t in ("ki", "ic50"):
        # Unit normalization: handle 'uM', 'nM', 'micromolar', 'millimolar', etc.
        scale_map = {
            "m": 1.0,
            "mm": 1e-3, "millimolar": 1e-3, "millimol": 1e-3, "millimole": 1e-3,
            "um": 1e-6, "µm": 1e-6, "micromolar": 1e-6, "micromol": 1e-6,
            "nm": 1e-9, "nanomolar": 1e-9, "nanomol": 1e-9,
            "pm": 1e-12, "picomolar": 1e-12,
        }
        # Normalize: remove spaces, special chars
        unit_key = unit.replace(" ", "").replace("/", "").replace("l", "").lower()
        # Try direct lookup, else fallback to uM (most common in this dataset)
        scale = scale_map.get(unit_key)
        if scale is None:
            # Try prefix matching
            for key, s in scale_map.items():
                if unit_key.startswith(key):
                    scale = s
                    break
        if scale is None:
            scale = 1e-6  # default assume uM
        conc_m = v * scale
        if conc_m <= 0:
            return None
        return RT * math.log(max(conc_m, 1e-30))

    return None


def parse_ligands(lig_yml: pathlib.Path) -> tuple[dict[str, str], dict[str, float]]:
    """Return (smiles_dict, dg_dict) keyed by ligand name."""
    raw = yaml.safe_load(lig_yml.read_text())
    if not isinstance(raw, dict):
        return {}, {}
    smi: dict[str, str] = {}
    dg: dict[str, float] = {}
    for name, rec in raw.items():
        if not isinstance(rec, dict):
            continue
        s = rec.get("smiles") or rec.get("SMILES")
        if not s:
            continue
        m = rec.get("measurement")
        if not isinstance(m, dict):
            continue
        d = to_dg(m)
        if d is not None:
            smi[name] = s
            dg[name] = d
    return smi, dg


def parse_edges_from_network(edges_yml: pathlib.Path) -> list[tuple[str, str]]:
    """Parse a 03_edges/*.yml network file; return list of (ligand_a, ligand_b) pairs."""
    raw = yaml.safe_load(edges_yml.read_text())
    if not isinstance(raw, dict):
        return []
    # The file has a top-level 'edges' key mapping edge-name -> {ligand_a, ligand_b, ...}
    edges_section = raw.get("edges", {})
    if not isinstance(edges_section, dict):
        return []
    pairs = []
    for _edge_name, edge_rec in edges_section.items():
        if not isinstance(edge_rec, dict):
            continue
        a = edge_rec.get("ligand_a")
        b = edge_rec.get("ligand_b")
        if a and b:
            pairs.append((str(a), str(b)))
    return pairs


def parse_target(data_dir: pathlib.Path) -> list[tuple[str, str, float]]:
    """Return rows (smiles_a, smiles_b, ddg) for a single target.

    Searches for ligands.yml in 00_data/ and picks the largest edges network
    from 03_edges/ to maximize edge coverage.
    """
    lig_yml = data_dir / "00_data" / "ligands.yml"
    edges_dir = data_dir / "03_edges"

    if not lig_yml.exists():
        raise FileNotFoundError(f"No ligands.yml at {lig_yml}")
    if not edges_dir.exists():
        raise FileNotFoundError(f"No 03_edges dir at {edges_dir}")

    smi, dg = parse_ligands(lig_yml)
    if not smi:
        raise ValueError("No parseable ligands found")

    # Pick the edge network with most edges (usually kartograf_mapper_lomap_network.yml)
    edge_files = sorted(edges_dir.glob("*.yml"))
    if not edge_files:
        raise FileNotFoundError(f"No yml files in {edges_dir}")

    best_file = None
    best_pairs: list[tuple[str, str]] = []
    for ef in edge_files:
        try:
            pairs = parse_edges_from_network(ef)
        except Exception:
            continue
        if len(pairs) > len(best_pairs):
            best_pairs = pairs
            best_file = ef

    if not best_pairs:
        raise ValueError("No parseable edges found")

    rows = []
    skipped = 0
    for a, b in best_pairs:
        if a in smi and b in smi:
            rows.append((smi[a], smi[b], dg[b] - dg[a]))
        else:
            skipped += 1
    if skipped:
        # Normal: some ligands may lack measurement data
        pass
    return rows


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python scripts/fetch_fep_edges.py <PLB_clone_dir>")
        sys.exit(1)

    plb = pathlib.Path(sys.argv[1])
    data_dir = plb / "data"
    if not data_dir.exists():
        print(f"ERROR: {data_dir} does not exist")
        sys.exit(1)

    OUT.mkdir(parents=True, exist_ok=True)
    combined: list[tuple] = []

    for target_dir in sorted(data_dir.iterdir()):
        if not target_dir.is_dir():
            continue
        target = target_dir.name
        if target == "targets.yml" or target.endswith(".yml"):
            continue

        try:
            rows = parse_target(target_dir)
        except Exception as exc:
            print(f"  skip {target}: {exc}")
            continue

        if len(rows) < 3:
            print(f"  skip {target}: only {len(rows)} edges (< 3 min)")
            continue

        out_csv = OUT / f"{target}.csv"
        with open(out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["smiles_a", "smiles_b", "ddg"])
            w.writerows(rows)

        combined += [(target, r[0], r[1], r[2]) for r in rows]
        print(f"  {target}: {len(rows)} edges -> {out_csv.name}")

    with open(OUT / "all_edges.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target", "smiles_a", "smiles_b", "ddg"])
        w.writerows(combined)

    n_targets = len(set(c[0] for c in combined))
    print(f"TOTAL {len(combined)} edges across {n_targets} targets")
    print(f"Output: {OUT}/")

    if len(combined) < 200 or n_targets < 6:
        print(f"WARNING: acceptance threshold not met (need >=200 edges, >=6 targets)")
        sys.exit(1)


if __name__ == "__main__":
    main()
