"""Target-finding arm (the cage merge) — structure-based reverse screening.

SECOND arm of the architecture (the RBFE-calibration arm is `src/bar/`). This arm is
under validation: the load-bearing gate is whether retrospective, fold-disjoint,
structure-based target-ID beats the ligand-only (Fig F) and docking-only baselines.
See `docs/target_finding_plan.md`. Built lazily — only what the gate needs.
"""
