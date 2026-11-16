#!/usr/bin/env python3
"""POST-4 §9 item 607 — census of Lagrange-P1 functionspace sites in examples/.

Stdlib only (no dolfinx). Walks ``examples/**/*.py`` with ``ast``, finds every
``functionspace(<mesh>, ("Lagrange"|"CG"|"P", 1[, shape]))`` call, and prints one
row per site: file:line, the space variable, the ``fem.Function`` built on it in
the same function body (if any) and the first ``.interpolate(<src>)`` fed into
that function (if any).

Asserted:
  (i)  each of the ten file paths in the 2026-08-12 known-issues table still
       exists (or is reported removed), with its current P1 site line(s);
  (ii) the AST row count equals an independent line-regex count printed in the
       same window (negative control: a spelling the parser misses, or a
       multi-line construction the regex misses, shows as a mismatch).

Run from the repo root:  python3 scripts/probes/post4_p1_site_census.py
"""

from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"
FAMILIES = {"Lagrange", "CG", "P"}

AUGUST_TABLE = [
    "examples/mri/01_coil_phantom_fields.py",
    "examples/magnetostatics/01_straight_wire.py",
    "examples/magnetostatics/02_circular_loop.py",
    "examples/magnetostatics/04_helmholtz_analytic_comparison.py",
    "examples/magnetostatics/05_gauge_cross_check.py",
    "examples/time_harmonic/01_lossy_plane_wave.py",
    "examples/time_harmonic/02_pec_cavity_resonances.py",
    "examples/time_harmonic/03_dielectric_sphere_in_uniform_field.py",
    "examples/time_harmonic/04_evanescent_waveguide_decay.py",
    "examples/time_harmonic/05_resonance_guard_sweep.py",
]

# Independent control: one source line that both calls functionspace( and
# names a P1 family tuple. Docstring/comment mentions without a call are not
# constructions and are excluded by requiring "functionspace(" on the line.
LINE_RE = re.compile(r'functionspace\(.*\(\s*"(Lagrange|CG|P)"\s*,\s*1\s*[,)]')
# The brief's literal grep, printed for the record (it also matches prose).
BRIEF_RE = re.compile(r'"Lagrange", 1|\("CG", 1|\("P", 1')


def _call_name(node: ast.AST) -> str | None:
    if isinstance(node, ast.Attribute):
        return node.attr
    if isinstance(node, ast.Name):
        return node.id
    return None


def _is_p1_space(call: ast.Call) -> bool:
    if _call_name(call.func) not in ("functionspace", "FunctionSpace"):
        return False
    if len(call.args) < 2 or not isinstance(call.args[1], ast.Tuple):
        return False
    elts = call.args[1].elts
    if len(elts) < 2:
        return False
    fam, deg = elts[0], elts[1]
    return (
        isinstance(fam, ast.Constant)
        and fam.value in FAMILIES
        and isinstance(deg, ast.Constant)
        and deg.value == 1
        and not isinstance(deg.value, bool)
    )


def _scopes(tree: ast.Module):
    """Yield (scope_node, direct-and-nested nodes excluding inner defs)."""
    yield tree
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            yield node


def _innermost_scope(tree: ast.Module, lineno: int) -> ast.AST:
    best = tree
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            end = getattr(node, "end_lineno", node.lineno)
            if node.lineno <= lineno <= end:
                if best is tree or node.lineno >= best.lineno:
                    best = node
    return best


def _assigned_name(tree: ast.Module, call: ast.Call) -> str | None:
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and node.value is call:
            tgt = node.targets[0]
            if isinstance(tgt, ast.Name):
                return tgt.id
            return ast.unparse(tgt)
    return None


def _functions_on(scope: ast.AST, space: str, after: int) -> list[tuple[str, int]]:
    out = []
    for node in ast.walk(scope):
        if (
            isinstance(node, ast.Assign)
            and isinstance(node.value, ast.Call)
            and _call_name(node.value.func) == "Function"
            and node.value.args
            and isinstance(node.value.args[0], ast.Name)
            and node.value.args[0].id == space
            and node.lineno >= after
            and isinstance(node.targets[0], ast.Name)
        ):
            out.append((node.targets[0].id, node.lineno))
    return sorted(out, key=lambda t: t[1])


def _first_interpolate(scope: ast.AST, fname: str, after: int) -> tuple[str, int] | None:
    hits = []
    for node in ast.walk(scope):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "interpolate"
            and isinstance(node.func.value, ast.Name)
            and node.func.value.id == fname
            and node.lineno >= after
        ):
            src = ast.unparse(node.args[0]) if node.args else "?"
            hits.append((src, node.lineno))
    return min(hits, key=lambda t: t[1]) if hits else None


def census() -> list[dict]:
    rows = []
    for path in sorted(EXAMPLES.rglob("*.py")):
        rel = path.relative_to(ROOT).as_posix()
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=rel)
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and _is_p1_space(node):
                space = _assigned_name(tree, node)
                scope = _innermost_scope(tree, node.lineno)
                funcs = _functions_on(scope, space, node.lineno) if space else []
                interp = []
                for fname, fline in funcs:
                    hit = _first_interpolate(scope, fname, fline)
                    interp.append((fname, fline, hit))
                rows.append(
                    {
                        "file": rel,
                        "line": node.lineno,
                        "space": space,
                        "scope": getattr(scope, "name", "<module>"),
                        "interp": interp,
                    }
                )
    rows.sort(key=lambda r: (r["file"], r["line"]))
    return rows


def regex_counts() -> tuple[int, int, list[str]]:
    n_ctrl = n_brief = 0
    brief_only = []
    for path in sorted(EXAMPLES.rglob("*.py")):
        rel = path.relative_to(ROOT).as_posix()
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            c = bool(LINE_RE.search(line))
            b = bool(BRIEF_RE.search(line))
            n_ctrl += c
            n_brief += b
            if b and not c:
                brief_only.append(f"{rel}:{i}: {line.strip()}")
    return n_ctrl, n_brief, brief_only


def main() -> int:
    rows = census()
    print(f"POST-4 P1-site census — {len(list(EXAMPLES.rglob('*.py')))} example .py files")
    print("file:line | space var | scope | Function <- interpolate(src)")
    for r in rows:
        if r["interp"]:
            parts = []
            for fname, fline, hit in r["interp"]:
                if hit:
                    parts.append(f"{fname}(L{fline}) <- {hit[0]} (L{hit[1]})")
                else:
                    parts.append(f"{fname}(L{fline}) <- (no interpolate)")
            desc = "; ".join(parts)
        else:
            desc = "(no fem.Function on this space in scope)"
        print(f"ROW {r['file']}:{r['line']} | {r['space']} | {r['scope']} | {desc}")

    n_ctrl, n_brief, brief_only = regex_counts()
    print()
    print(f"AST rows                         : {len(rows)}")
    print(f"line-regex control count         : {n_ctrl}")
    print(f"brief's literal grep line count  : {n_brief}")
    for line in brief_only:
        print(f"  brief-grep-only (not a construction): {line}")

    print()
    print("August 2026-08-12 table, current state:")
    ok_i = True
    for rel in AUGUST_TABLE:
        p = ROOT / rel
        if not p.exists():
            print(f"  REMOVED/MOVED {rel}")
            continue
        lines = [r["line"] for r in rows if r["file"] == rel]
        print(f"  present {rel}: P1 site line(s) {lines}")
        if not lines:
            ok_i = False
            print(f"    !! present but no P1 site found")

    ok_ii = len(rows) == n_ctrl
    print()
    print(f"(i)  August paths accounted for : {'PASS' if ok_i else 'FAIL'}")
    print(f"(ii) AST rows == regex control  : {'PASS' if ok_ii else 'FAIL'} ({len(rows)} vs {n_ctrl})")
    return 0 if (ok_i and ok_ii) else 1


if __name__ == "__main__":
    sys.exit(main())
