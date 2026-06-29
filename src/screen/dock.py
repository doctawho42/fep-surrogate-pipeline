"""Structure-based docking orchestration (smina) — the target-finding arm's scorer.

We ORCHESTRATE existing tools, we do not rebuild them (docs/target_finding_plan.md §4):
download a holo PDB, auto-detect its co-crystal ligand as the pocket reference, prep a
rigid receptor (Open Babel), and score a ligand with smina (AutoDock Vina fork). The
honest methodology (fold-disjoint split, recovery, calibration) lives in `recovery.py`.
"""
from __future__ import annotations

import pathlib
import subprocess
import urllib.request

SMINA = "/Users/nikitapolomosnov/anaconda3/envs/synfp_dock/bin/smina"
OBABEL = "obabel"

# residue names that are NOT the co-crystal ligand (water / ions / common buffers/cryo)
_NON_LIGAND = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "PI", "GOL", "EDO", "PEG", "PG4", "ACT", "DMS",
    "MES", "TRS", "FMT", "CIT", "MPD", "BME", "IPA", "NO3", "NH4", "CO3", "EPE", "TLA",
    "NA", "CL", "MG", "ZN", "CA", "K", "MN", "FE", "CD", "NI", "CU", "HG", "BR", "IOD",
}


def download_pdb(pdb_id: str, dest_dir: str) -> str:
    pdb_id = pdb_id.upper()
    path = pathlib.Path(dest_dir) / f"{pdb_id}.pdb"
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        urllib.request.urlretrieve(url, path)
    return str(path)


def extract_pocket(pdb_file: str, out_dir: str, lig_resname: str | None = None) -> tuple[str, str]:
    """Split into (receptor.pdb, ref_ligand.pdb). Receptor = all ATOM of the first
    protein chain; ref ligand = the HETATM residue named ``lig_resname`` if given, else
    the largest non-buffer HETATM (auto-detected). Use ``lig_resname`` when a cofactor
    (FMN/heme) is larger than the inhibitor of interest."""
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    receptor_lines, het = [], {}
    first_chain = None
    for ln in pathlib.Path(pdb_file).read_text().splitlines():
        rec = ln[:6].strip()
        if rec == "ATOM":
            ch = ln[21]
            if first_chain is None:
                first_chain = ch
            if ch == first_chain:
                receptor_lines.append(ln)
        elif rec == "HETATM":
            resn = ln[17:20].strip()
            if lig_resname is None and resn in _NON_LIGAND:
                continue
            key = (resn, ln[21], ln[22:26].strip())
            het.setdefault(key, []).append(ln)
    if lig_resname is not None:
        cand = {k: v for k, v in het.items() if k[0] == lig_resname}
        if not cand:
            raise ValueError(f"ligand {lig_resname} not in {pdb_file}")
        lig_key = max(cand, key=lambda k: len(cand[k]))
    elif not het:
        raise ValueError(f"no co-crystal ligand found in {pdb_file}")
    else:
        lig_key = max(het, key=lambda k: len(het[k]))
    rec_path = out / "receptor.pdb"
    lig_path = out / "ref_ligand.pdb"
    rec_path.write_text("\n".join(receptor_lines) + "\n")
    lig_path.write_text("\n".join(het[lig_key]) + "\n")
    return str(rec_path), str(lig_path)


def prep_receptor(receptor_pdb: str) -> str:
    """Rigid receptor PDBQT via Open Babel (add H at pH 7.4)."""
    out = str(pathlib.Path(receptor_pdb).with_suffix(".pdbqt"))
    subprocess.run([OBABEL, receptor_pdb, "-O", out, "-xr", "-p", "7.4"],
                   check=True, capture_output=True)
    return out


def prepare_target(pdb_id: str, work_dir: str, lig_resname: str | None = None) -> dict:
    """Download + pocket + receptor prep. ``lig_resname`` boxes on a specific inhibitor
    (use when a cofactor/lipid is the largest HETATM). Returns paths for docking."""
    d = pathlib.Path(work_dir) / pdb_id.upper()
    pdb = download_pdb(pdb_id, str(d))
    rec_pdb, ref_lig = extract_pocket(pdb, str(d), lig_resname=lig_resname)
    rec_pdbqt = prep_receptor(rec_pdb)
    return {"pdb_id": pdb_id.upper(), "receptor": rec_pdbqt, "ref_ligand": ref_lig}


def dock(receptor_pdbqt: str, ligand_sdf: str, ref_ligand_pdb: str,
         exhaustiveness: int = 8, box_pad: float = 6.0, seed: int = 1,
         timeout_s: int = 60) -> float:
    """Best smina affinity (kcal/mol) of the ligand in the autoboxed pocket; +inf on
    failure. Lower = better."""
    try:
        r = subprocess.run(
            [SMINA, "--receptor", receptor_pdbqt, "--ligand", ligand_sdf,
             "--autobox_ligand", ref_ligand_pdb, "--autobox_add", str(box_pad),
             "--exhaustiveness", str(exhaustiveness), "--num_modes", "1",
             "--seed", str(seed), "--cpu", "1", "-o", "/dev/null"],
            check=True, capture_output=True, text=True, timeout=timeout_s)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return float("inf")
    best = float("inf")
    for ln in r.stdout.splitlines():
        parts = ln.split()
        if len(parts) >= 2 and parts[0] == "1":
            try:
                best = float(parts[1])
            except ValueError:
                pass
    return best
