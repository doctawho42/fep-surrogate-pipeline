"""Audit every figure the manuscript includes against the ACS artwork specification.

Four requirements, checked on the rendered PDF rather than on the plotting code, because each
of them has been breached at some point by something the code did not say:

* width -- single-column artwork at most 240 pt, double-column between 300 and 504 pt;
* stroked lines at least 0.5 pt, which marker and bar edges quietly fall below;
* set type at least 4.5 pt, which mathtext sub- and superscripts reach at 0.7 of their base;
* Helvetica or Arial, which a mathtext font set overrides for everything inside ``$...$``.

The font check is the reason this reads the PDF. Naming Arial in ``font.sans-serif`` does not
make a figure Arial: a missing glyph falls back silently, and a local ``rcParams`` override in
one script put DejaVu Sans ahead of it for the contents graphic alone.

Run: ``make acscheck``. Exits non-zero on any breach, listing the figure and the value.
"""
from __future__ import annotations

import pathlib
import re
import subprocess
import sys

from pypdf import PdfReader

ROOT = pathlib.Path(__file__).resolve().parent.parent
SINGLE_MAX = 240.0
DOUBLE_MIN, DOUBLE_MAX = 300.0, 504.0
MIN_RULE_PT = 0.5
MIN_TYPE_PT = 4.5
FOREIGN = ("Cm", "DejaVu", "STIX", "LastResort")


def included_figures() -> list[str]:
    """Every graphic the manuscript actually includes, main text, supplement and wrappers."""
    tex = "".join(
        (ROOT / "docs" / name).read_text()
        for name in ("paper_body.tex", "paper_si.tex", "paper_jctc.tex")
    )
    return sorted({m.group(1) for m in re.finditer(r"includegraphics\[[^\]]*\]\{([^}]+)\}", tex)})


def audit(path: pathlib.Path) -> dict:
    page = PdfReader(str(path)).pages[0]
    stream = page.get_contents().get_data().decode("latin-1", "replace")
    # a zero width is the device's thinnest line and is set separately, so it is not a breach
    rules = [float(m.group(1)) for m in re.finditer(r"(?<![\d.])([\d.]+)\s+w[\s\n]", stream)]
    rules = [w for w in rules if w > 0]
    type_pt = [float(m.group(1)) for m in re.finditer(r"/F\d+\s+([\d.]+)\s+Tf", stream)]
    listing = subprocess.run(["pdffonts", str(path)], capture_output=True, text=True).stdout
    fonts = {line.split()[0].split("+")[-1] for line in listing.split("\n")[2:] if line.strip()}
    return {
        "width": float(page.mediabox.width),
        "rule": min(rules) if rules else float("inf"),
        "type": min(type_pt) if type_pt else float("inf"),
        "foreign": sorted(f for f in fonts if f.startswith(FOREIGN)),
    }


def main() -> int:
    breaches: list[str] = []
    for name in included_figures():
        path = ROOT / "figs" / name
        if not path.exists():
            breaches.append(f"{name}: included by the manuscript but not on disk")
            continue
        a = audit(path)
        w = a["width"]
        if not (w <= SINGLE_MAX or DOUBLE_MIN <= w <= DOUBLE_MAX):
            breaches.append(
                f"{name}: {w:.1f} pt wide, outside <={SINGLE_MAX:.0f} and "
                f"{DOUBLE_MIN:.0f}-{DOUBLE_MAX:.0f}")
        if a["rule"] < MIN_RULE_PT:
            breaches.append(f"{name}: thinnest line {a['rule']:.2f} pt, floor {MIN_RULE_PT}")
        if a["type"] < MIN_TYPE_PT:
            breaches.append(f"{name}: smallest type {a['type']:.2f} pt, floor {MIN_TYPE_PT}")
        if a["foreign"]:
            breaches.append(f"{name}: not Helvetica or Arial: {', '.join(a['foreign'])}")
        print(f"{name:44s} {w:6.1f} pt  line {a['rule']:.2f}  type {a['type']:.2f}  "
              f"{'ok' if not a['foreign'] else ','.join(a['foreign'])}")

    print()
    if breaches:
        print(f"{len(breaches)} breach(es) of the ACS artwork specification:")
        for b in breaches:
            print(f"  {b}")
        return 1
    print(f"{len(included_figures())} figures, all within the ACS artwork specification.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
