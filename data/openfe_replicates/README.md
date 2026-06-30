# OpenFE IndustryBenchmarks2024 — per-replicate RBFE edge data

`combined_pymbar4_edge_data.csv`: per transformation (ligand_A, ligand_B), the solvent and
complex legs each run as **3 independent replicates** (`*_repeat_{0,1,2}_DG`,
`*_dDG` = per-replicate MBAR/pymbar4 uncertainty, `*_smallest_overlap`), plus partner/system.

Binding ΔΔG per replicate = complex_DG − solvent_DG; reported single-replicate se =
sqrt(complex_dDG² + solvent_dDG²) (the MBAR uncertainty, ≡ sandwich B/I² to leading order).

Source: https://github.com/OpenFreeEnergy/IndustryBenchmarks2024
(industry_benchmarks/analysis/processed_results/combined_pymbar4_edge_data.csv).
Public precomputed FEP (no MD run here). Used by figs/make_figA_replicates.py.
