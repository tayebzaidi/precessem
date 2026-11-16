#!/usr/bin/env python3
"""`OPS-78` census reader: band names against grep, the classification, the control, the reads.

Stdlib only; run from the repo root.  No test, band or log is written.

(1) Each module under tests/, src/, examples/ is parsed with ``ast``.  A name resolves to its
definition through the module's own imports, because one name can carry different values in
different files (``IDENTITY_BAND`` is 1e-12 in ``test_birdcage_leg_offset.py`` and 1e-9 in
``test_geo38_f_clinical_build.py``).  Every line ``grep -rnw`` returns for a name <= 1e-11 lands in
one kind (definition, import, prose, out of scope, site, unresolved); a difference exits 2.

(2) Every site gets a class under item 767's rule.  ``RULES`` is the reading of each call: S (two
factorisations in one window), R (one factorisation, one field read twice or several drives through
one held factor), N (a mesh, area or schedule identity, or a closed form with no solve).  A module
with no ``RULES`` entry and no ``SOLVE_TOKENS`` hit is N in bulk.  The tokens are a filter, not a
proof: the modules they missed carry ``RULES`` entries.  A site nothing classifies prints UNREAD.

(3) The brief's three controls, read from ``RULES``, the source and the two logs; a split exits 3.
The source check is not a count of MUMPS factor calls.

(4) Largest recorded read of each separate-solve gate, through a ``.gz``-safe route; a read above
its band exits 1.  A pattern with two groups is a one-sided form, read as left / right against 1.
"""

from __future__ import annotations

import ast
import gzip
import re
import subprocess
import sys
import warnings
from pathlib import Path

warnings.simplefilter("ignore", SyntaxWarning)  # invalid escapes in docstrings of scanned modules

ROOTS = ("tests", "src", "examples")
LIMIT = 1.0e-11
LOGS = Path("docs/testing/logs")
STUDY = Path("src/fem_em_solver/study.py")
FLOAT = r"([0-9.]+e[+-][0-9]+)"

M21 = r"MAT-21-step"
M17 = r"MAT-17-step3|EX-99"
M20 = r"MAT-20-step2|EX-110"
# (label, builder, band, read patterns, log-name pattern).  Bands as the gates print them.
READS = (
    ("mri:37 / WF-15 s10 (S), Study vs direct", "different", 1e-12,
     (r"^(?!.*\[WF-16\]).*\|\|S_study - S_direct\|\|/\|\|S_direct\|\|\s+" + FLOAT,),
     r"EX-167|WF-15-s10|OPS-77-mri37"),
    ("WF-17 s25 (S), block vs direct", "different", 1e-12,
     (r"\|\|S_block - S_direct\|\|/\|\|S_direct\|\| " + FLOAT,), r"WF-17-s25"),
    ("WF-17 s25 (P=Y), YAML vs block", "same", 1e-12,
     (r"\|\|S_yaml - S_block\|\|/\|\|S_block\|\| " + FLOAT,), r"WF-17-s25"),
    ("mri:38 / WF-15 s11 (S), sweep vs point", "same", 1e-12,
     (r"sweep vs point rel " + FLOAT,), r"EX-170|WF-15-s11"),
    ("PORT-8 anchor S', flipped sheets vs DSD", "same", 1e-12,
     (r"‖S' − DSD‖_F/‖S‖_F = " + FLOAT,), r"PORT-8"),
    ("PORT-8 anchor Z', flipped sheets vs DZD", "same", 1e-12,
     (r"‖Z' − DZD‖_F/‖Z‖_F = " + FLOAT,), r"PORT-8"),
    ("PORT-21 step 5, option vs hand-rolled", "same", 1e-12, (r"max\|dS\|/max\|S\| = " + FLOAT,),
     r"PORT-2[12]"),
    ("MAT-17 step 3 / mat:15 (Delta)", "same", 1e-12,
     (r"MAT17s3 \(Delta\).* rel L2 " + FLOAT, r"\[anchor\] \(Delta\) Duhamel superposition " + FLOAT), M17),
    ("MAT-17 step 3 / mat:15 (D')", "same", 1e-13,
     (r"MAT17s3 \(D'\).* rel L2 " + FLOAT, r"\[anchor\] \(D'\) chi=1 vs step 1 " + FLOAT), M17),
    ("mat:15 (Hist), history runs vs the gate's", "same", 1e-13,
     (r"\(Hist\) history duty run vs gate's " + FLOAT, r"chi=1 vs T_star\[0\] " + FLOAT,
      r"(?<=[ ])n=\d+: " + FLOAT), M17),
    ("MAT-20 / mat:18 (H0)", "same", 1e-12,
     (r"MAT20fem \(H0\).* = " + FLOAT, r"\[anchor\] \(H0\) " + FLOAT), M20),
    ("MAT-20 / mat:18 (FEM)", "same", 1e-13,
     (r"MAT20fem \(FEM\).* = " + FLOAT, r"\[anchor\] \(FEM\) W = .* windows " + FLOAT), M20),
    ("MAT-20 / mat:18 (FEM-2d)", "same", 1e-13,
     (r"MAT20fem \(FEM-2d\) 2-D call: " + FLOAT, r"\[anchor\] \(FEM-2d\) " + FLOAT), M20),
    ("mat:18 (Cm), per-vertex field vs series", "same", 1e-13,
     (r"\[anchor\] \(Cm\) .* series " + FLOAT,), r"EX-110"),
    ("MAT-21 steps 2, 4 (H0)", "same", 1e-12, (r"MAT21-[24] \(H0\) max_n.* = " + FLOAT,), M21),
    ("MAT-21 steps 2, 4 (Delta)", "same", 1e-12, (r"MAT21-[24] \(Δ\).* = " + FLOAT,), M21),
    ("MAT-21 steps 2, 4, 11 (D')", "same", 1e-13,
     (r"MAT21-[24] \(D'\).* = " + FLOAT, r"MAT21-11 \(D'-quad\).* rel L2 " + FLOAT), M21),
    ("MAT-21 steps 2, 4 (FEM'), residual subtracted", "same", 1e-13,
     (r"MAT21-[24] \(FEM'\) max.* = " + FLOAT,), M21),
    ("MAT-21 steps 2, 4 (FEM'-2d), residual subtracted", "same", 1e-13,
     (r"MAT21-[24] \(FEM'\) max.*\(FEM'-2d\) " + FLOAT,), M21),
    ("MAT-21 steps 2, 4 (FEM-c), left / right", "same", 1.0,
     (r"MAT21-[24] \(FEM-c\) max \|a - b\| = " + FLOAT + r" K <= .* = " + FLOAT + " K",), M21),
    ("MAT-21 step 13 (C100-h)", "same", 1e-13, (r"\(C100-h\).*\|ratio - 1\| " + FLOAT,), M21),
    ("WF-15 s13 (Id), caps Study vs direct", "different", 1e-12,
     (r"\[WF-15 s13\] .*\(Id\) \|\|S_study - S_direct\|\|/\|\|S_direct\|\| " + FLOAT,), r"WF-15-s13"),
    ("WF-17 s28 (Id), caps block vs direct", "different", 1e-12,
     (r"\[WF-17 s28\] \(Id\) \|\|S_block - S_direct\|\|/\|\|S_direct\|\| " + FLOAT,), r"WF-17-s28"),
    ("WF-16 s18 (I), fields_at sweep vs private route", "same", 1e-12,
     (r"\[WF-16 s18\] \(I\) \|B1\+\| P\d: valid \d+/\d+, max rel " + FLOAT,
      r"\[WF-16 s18\] \(I\) SAR10g .* rel " + FLOAT), r"WF-16-s18"),
)

SOLVE_TOKENS = (".solve(", "solve_sweep", "run_n_port_sparameter_sweep", "build_four_port_sweep(", "_sweep(",
                "solve_time_harmonic", "solve_transient", "_solve_", "solver.solve", ".solve_ports",
                "run_two_port", "_run_drive")

# The (FEM) of MAT-21 steps 2 and 4 as posed is printed only since the 2026-10-13 ruling (no site).
ALL = None
LIN, DUHAMEL, CONTROL = "(Lin) Y_in, node pair held", "Duhamel / chi = 1 runs", "control, a miss"
RULES = {
    "examples/implants/05_lead_shaft_polarisation_in_gel.py": [((300, 302), "R", "(X1) map of one field"),
                                                               ((303,), "N", "(X2) closed form")],
    "examples/implants/06_insulated_lead_shaft_in_gel.py": [((344, 346), "R", "(X1) map of one field"),
                                                            ((347,), "N", "(X2) closed form")],
    "examples/materials/14_coil_driven_phantom_in_container_wall.py": [(ALL, "N", "(R) numpy round-off")],
    "examples/materials/15_duty_cycled_coil_heating_phantom.py": [((218, 262), "S", CONTROL), (ALL, "S", DUHAMEL)],
    "examples/materials/18_six_minute_window_coil_heating.py": [(ALL, "N", "name bindings")],
    "examples/materials/31_gel_box_wall_duty_cycled_transient.py": [(ALL, "N", "(S0) census")],
    "examples/materials/32_gel_box_wall_thermal_dose.py": [(ALL, "N", "(S0) census")],
    "examples/materials/36_gel_box_wait_across_schedules.py": [(ALL, "N", "(S0) census")],
    "examples/mri/17_f_human_tuned_quadrature_b1_sar.py": [(ALL, "R", LIN)],
    "examples/mri/18_f_human_four_port_quadrature.py": [(ALL, "R", LIN)],
    "examples/mri/28_f_human_four_port_sar10g_c4.py": [(ALL, "R", LIN)],
    "examples/mri/37_implant_rod_through_a_study.py": [((211, 212, 213), "S", CONTROL), (ALL, "S", "(S) rod")],
    "examples/mri/38_implant_rod_on_a_sweep.py": [((162, 164), "S", CONTROL), (ALL, "S", "(S) sweep vs point")],
    "examples/ports/03_lumped_sheet_port_widths.py": [(ALL, "R", "one solved field, two reads")],
    "tests/package/test_wf15_study_implant_rod.py": [((298, 299), "N", "sheet geometry"), (ALL, "S", "(S) rod")],
    "tests/package/test_wf15_study_rod_caps.py": [(ALL, "S", "(Id) caps, Study vs direct")],
    "tests/package/test_wf15_study_rod_sweep.py": [((138, 140), "S", CONTROL), (ALL, "S", "(S) sweep vs point")],
    "tests/package/test_wf16_study_four_port.py": [(ALL, "N", "sheet geometry")],
    "tests/package/test_wf16_sweep_fields_at.py": [(ALL, "S", "(I) fields_at sweep vs private route")],
    "tests/package/test_wf16_study_tuned_f_human.py": [((323,), "R", "(Y) one s_node"), (ALL, "N", "spec / weights")],
    "tests/package/test_wf17_implant_rod_tutorial.py": [((175, 362), "N", "page text"), ((281, 302), "S", CONTROL),
                                                        ((422, 423, 430), "S", "(Id) caps, block vs direct"),
                                                        (ALL, "S", "(S), (P=Y) rod")],
    "tests/package/test_wf17_tuned_copper_tutorial.py": [(ALL, "N", "spied centre")],
    "tests/ports/test_port_orientation_sensitivity.py": [(ALL, "S", "PORT-8 anchors, two sweeps, one mesh")],
    "tests/solver/test_energy_and_point_evaluation.py": [(ALL, "R", "(Im W) one field")],
    "tests/solver/test_single_port_excitation.py": [(ALL, "R", "Y symmetry, one held factor")],
    "tests/unit/test_ans8_v3_rows.py": [(ALL, "N", "synthetic network")],
    "tests/validation/test_ans8_v3_record_rows.py": [((260, 261), "N", "synthetic network"), (ALL, "R", LIN)],
    "tests/validation/test_geo38_f_clinical_esr_root.py": [(ALL, "N", "synthetic network")],
    "tests/validation/test_imp4_pec_cylinder_fem.py": [(ALL, "N", "closed form / quadrature")],
    "tests/validation/test_imp5_coated_cylinder_fem.py": [(ALL, "N", "closed form / quadrature")],
    "tests/validation/test_mat15_f_human_tuned_heating.py": [(ALL, "R", LIN)],
    "tests/validation/test_mat17_duty_cycle_coil_driven.py": [((188, 293), "S", CONTROL), (ALL, "S", DUHAMEL)],
    "tests/validation/test_mat18_wall_coil_driven.py": [(ALL, "N", "(R) numpy round-off")],
    "tests/validation/test_mat20_fem_history.py": [((183, 196), "S", CONTROL), (ALL, "S", DUHAMEL)],
    "tests/validation/test_mat21_box_dose_fem.py": [(ALL, "R", "(FEM-c-win) one run vs series")],
    "tests/validation/test_mat21_box_probe_dose_fem.py": [(ALL, "R", "(FEM-c-win) one run vs series")],
    "tests/validation/test_mat21_box_duty_dose_fem.py": [(ALL, "S", "(C100-h) chi = 1 vs default")],
    "tests/validation/test_mat21_box_duty_fem.py": [((296, 304), "S", CONTROL), (ALL, "S", DUHAMEL)],
    "tests/validation/test_mat21_box_duty_fem_p2.py": [((270, 278), "S", CONTROL), (ALL, "S", DUHAMEL)],
    "tests/validation/test_mat21_box_quadratic_duty_fem.py": [(ALL, "S", DUHAMEL)],
    "tests/validation/test_mat21_pre_equilibration_duty_dose_fem.py": [(ALL, "N", "(S0) census")],
    "tests/validation/test_mat21_wall_duty25_fem.py": [(ALL, "N", "(S0) census")],
    "tests/validation/test_mat21_wall_duty_dose_fem.py": [(ALL, "N", "(S0) census")],
    "tests/validation/test_mat21_wall_duty_fem.py": [(ALL, "N", "(S0) census")],
    "tests/validation/test_mat21_wall_pre_equilibration_fem.py": [(ALL, "N", "(S0) census")],
    "tests/validation/test_port22_driven_sweep_resonance.py": [(ALL, "S", "PORT-21 route identity, two sweeps")],
    "tests/validation/test_port_circuit_skrf.py": [(ALL, "N", "recorded S through algebra")],
    "tests/validation/test_port_lumped_narrowed_sheet.py": [(ALL, "R", "one solved field, two reads")],
    "tests/validation/test_port_lumped_two_torus.py": [((717,), "N", "area partition"),
                                                       (ALL, "R", "one solved field, two reads")],
    "tests/validation/test_port_package_sparameters.py": [(ALL, "R", "S vs z_to_s(Z), one sweep")],
    "tests/validation/test_th16_quarter_birdcage.py": [(ALL, "R", "reciprocity, one held factor")],
    "tests/validation/test_th20_quarter_degree2.py": [(ALL, "R", "reciprocity printed, one held factor")],
    "tests/validation/test_th20_quarter_solves.py": [(ALL, "R", "reciprocity printed, one held factor")],
    "tests/validation/test_wf10_f_human_copper_esr_root.py": [(ALL, "N", "synthetic network / planted root")],
    "tests/validation/test_wf10_f_human_four_port_sar10g.py": [(ALL, "R", LIN)],
    "tests/validation/test_wf10_f_human_refined_b1_sar.py": [(ALL, "R", LIN)],
    "tests/validation/test_wf10_f_human_tuned_b1_sar.py": [(ALL, "R", LIN)],
    "tests/validation/test_wf10_f_human_tuned_sar10g.py": [(ALL, "R", LIN)],
    "tests/validation/test_wf12_f_human_copper_c_sweep.py": [(ALL, "R", "route printed, one matrix")],
}
# The 21 lines ast cannot follow.  mat:18 binds m17s3's D_BAND / L_BAND, mri:32 binds b.P_BAND = 1e-2.
UNRESOLVED = {
    "examples/materials/18_six_minute_window_coil_heating.py": [((222, 235, 267, 269), "S", CONTROL),
                                                                (ALL, "S", DUHAMEL)],
    "examples/mri/32_femem_run_ball_probes_and_duty_quadratic_files.py": [(ALL, "out", "P_BAND 1e-2, MAT-10's")],
}


def _number(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
        return float(node.value)
    return None


def _dotted(path):
    key = ".".join(path.with_suffix("").parts)
    return key[4:] if key.startswith("src.") else key


class Module:
    """One parsed file: band definitions, import bindings, and the lines names are read on."""

    def __init__(self, path):
        self.path = path
        tree = ast.parse(path.read_text(errors="replace"))
        self.defs, self.names, self.mods = {}, {}, {}
        self.importlines, self.reads = set(), []
        for node in ast.walk(tree):
            if isinstance(node, (ast.Assign, ast.AnnAssign)):
                targets = node.targets if isinstance(node, ast.Assign) else [node.target]
                for t in targets:
                    if isinstance(t, ast.Name) and t.id.endswith("_BAND") and node.value is not None:
                        self.defs[t.id] = (_number(node.value), node.lineno)
            elif isinstance(node, ast.ImportFrom) and node.module:
                for a in node.names:
                    self.names[a.asname or a.name] = (node.module, a.name)
                    self.importlines.add((a.name, a.lineno))
                    self.mods.setdefault(a.asname or a.name, f"{node.module}.{a.name}")
            elif isinstance(node, ast.Import):
                for a in node.names:
                    self.mods[a.asname or a.name] = a.name
            elif isinstance(node, ast.Name) and node.id.endswith("_BAND"):
                if not isinstance(node.ctx, ast.Store):
                    self.reads.append((node.id, node.lineno, None))
            elif isinstance(node, ast.Attribute) and node.attr.endswith("_BAND"):
                alias = node.value.id if isinstance(node.value, ast.Name) else None
                self.reads.append((node.attr, node.lineno, alias))


def _resolve(mods, dotted, name, seen=()):
    """The value of ``name`` as read in module ``dotted``, or None if it cannot be followed."""
    if (dotted, name) in seen or dotted not in mods:
        return None
    m = mods[dotted]
    if name in m.defs:
        return m.defs[name]
    if name in m.names:
        return _resolve(mods, *m.names[name], seen + ((dotted, name),))
    return None


def _grep(name):
    cmd = ["grep", "-rnw", "--include=*.py", name, *ROOTS]
    rows = subprocess.run(cmd, capture_output=True, text=True, check=False).stdout.splitlines()
    return [(r.split(":", 2)[0], int(r.split(":", 2)[1])) for r in rows]


def _kind(mods, m, name, line):
    dotted = _dotted(m.path)
    small = lambda r: r is not None and r[0] is not None and abs(r[0]) <= LIMIT  # noqa: E731
    if name in m.defs and m.defs[name][1] == line:
        return "def" if small(m.defs[name]) else "out_of_scope"
    if (name, line) in m.importlines:
        return "import" if small(_resolve(mods, dotted, name)) else "out_of_scope"
    rows = [x for x in m.reads if x[0] == name and x[1] == line]
    if not rows:
        return "prose"  # comment or docstring text; code reads have an AST node on their line
    alias = rows[0][2]
    target = m.mods.get(alias) if alias else None
    r = _resolve(mods, dotted, name) if alias is None else (_resolve(mods, target, name) if target else None)
    if r is None:
        return "unresolved"
    return "site" if small(r) else "out_of_scope"


def census():
    """Per-name reconciliation against grep; returns (every line got a kind, sites by kind)."""
    mods = {}
    for root in ROOTS:
        for path in sorted(Path(root).rglob("*.py")):
            mods[_dotted(path)] = Module(path)
    by_path = {str(m.path): m for m in mods.values()}
    names = sorted({n for m in mods.values() for n, (v, _) in m.defs.items()
                    if v is not None and abs(v) <= LIMIT})
    print(f"RECORD {len(by_path)} modules parsed; {len(names)} names defined <= {LIMIT:g}")
    ok, sites, unresolved = True, 0, 0
    found = {}
    for name in names:
        rows = _grep(name)
        kinds = {}
        for path, line in rows:
            kind = _kind(mods, by_path[path], name, line)
            kinds[kind] = kinds.get(kind, 0) + 1
            if kind in ("site", "unresolved"):
                found.setdefault(kind, {}).setdefault(path, []).append((line, name))
        same = sum(kinds.values()) == len(rows)
        ok &= same
        sites += kinds.get("site", 0)
        unresolved += kinds.get("unresolved", 0)
        parts = " ".join(f"{k}={v}" for k, v in sorted(kinds.items()))
        print(f"NAME {name}: grep {len(rows)} lines, kinds {sum(kinds.values())} "
              f"({'equal' if same else 'DIFFER'}): {parts}")
    print(f"RESULT census names {len(names)}: {sites} use sites resolved <= {LIMIT:g}, "
          f"{unresolved} unresolved, grep totals {'equal' if ok else 'DIFFER'}")
    return ok, found


def _rule(rules, line):
    for lines, cls, gate in rules:
        if lines is None or line in lines:
            return cls, gate
    return None


def classify(found):
    """Every site and unresolved line to a class; returns (per-class counts, unread rows)."""
    counts, gates, unread = {}, {}, []
    for kind, table in (("site", RULES), ("unresolved", UNRESOLVED)):
        for path, rows in sorted(found.get(kind, {}).items()):
            has_tokens = any(t in Path(path).read_text(errors="replace") for t in SOLVE_TOKENS)
            for line, name in rows:
                hit = _rule(table[path], line) if path in table else None
                if hit is None and kind == "site" and not has_tokens and path not in table:
                    hit = ("N", "no solver token")
                if hit is None:
                    unread.append((path, line, name))
                    continue
                counts[hit[0]] = counts.get(hit[0], 0) + 1
                gates.setdefault((hit[0], hit[1]), []).append(f"{path}:{line}")
    for cls in ("S", "R", "N", "out"):
        print(f"CLASS {cls}: {counts.get(cls, 0)} sites")
    for (cls, gate), where in sorted(gates.items()):
        print(f"GATE {cls} | {gate} | {len(where)} sites in {len({w.split(':')[0] for w in where})} modules")
    print(f"RESULT classification: {sum(counts.values())} classified, {len(unread)} unread")
    for path, line, name in unread:
        print(f"UNREAD {path}:{line} {name}")
    return counts, unread


def _function_source(path, name):
    tree = ast.parse(path.read_text(errors="replace"))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return "\n".join(path.read_text(errors="replace").splitlines()[node.lineno - 1:node.end_lineno])
    sys.exit(f"{path}: no function {name}")


def _own_factor(path, name):
    """True when a function builds its own mesh and calls the port sweep on it."""
    src = _function_source(path, name)
    return "geom.build(comm)" in src and "run_n_port_sparameter_sweep(" in src


def _text(path):
    """A harness log's text, from the plain file or its ``.gz``."""
    if path.exists():
        return path.read_text(errors="replace")
    with gzip.open(f"{path}.gz", "rt", errors="replace") as fh:
        return fh.read()


def control():
    """The three controls of item 767; True when each lands as the brief requires."""
    mri38 = "examples/mri/38_implant_rod_on_a_sweep.py"
    mri37 = "examples/mri/37_implant_rod_through_a_study.py"
    torus = "tests/validation/test_port_lumped_two_torus.py"
    good = True

    cls = _rule(RULES[mri38], 144)[0]
    built = all(_own_factor(STUDY, f) for f in ("_study_solve_sweep", "_study_solve_ports"))
    pat = re.compile(r"sweep vs point rel " + FLOAT)
    zeros = [(n, float(m.group(1)))
             for n, line in enumerate(_text(LOGS / "20261102T060640Z_EX-170.log").splitlines(), 1)
             if (m := pat.search(line))]
    ok = cls == "S" and built and all(v == 0.0 for _, v in zeros)
    good &= ok
    print(f"CONTROL mri:38 -> {cls}, each side builds and factorises: {built}, builder same (Study on both "
          f"sides); EX-170 reads {zeros}: {'as required' if ok else 'SPLIT WRONGLY'}")

    cls = _rule(RULES[mri37], 177)[0]
    built = _own_factor(STUDY, "_study_solve_port_fields")
    src = Path(mri37).read_text(errors="replace")
    other = "_build_rod(" in src and "_sweep(" in src
    r2 = "20261102T121721Z_OPS-77-mri37-r2.log"
    line = _text(LOGS / r2).splitlines()[5567]
    ok = cls == "S" and built and other and "1.147e-13" in line
    good &= ok
    print(f"CONTROL mri:37 -> {cls}, Study side builds and factorises: {built}, direct side _build_rod: "
          f"{other}, builder different; {r2}:5568 {line.strip()[:70]}: {'as required' if ok else 'SPLIT WRONGLY'}")

    cls = _rule(RULES[torus], 709)[0]
    src = _function_source(Path(torus), "lumped_run")
    solves = len(re.findall(r"\.solve\(", src))
    ok = cls == "R" and solves == 1
    good &= ok
    print(f"CONTROL two_torus:709 -> {cls}, `.solve(` calls in lumped_run: {solves}, "
          f"both voltages read off that field: {'as required' if ok else 'SPLIT WRONGLY'}")
    print(f"CONTROL verdict: {'as required' if good else 'SPLIT WRONGLY'}")
    return good


def reads():
    """Largest recorded read per gate, with log and line; returns True when none is above band."""
    names = sorted({p.name.removesuffix(".gz") for p in LOGS.iterdir() if p.is_file()})
    keys = [(label, builder, band, [re.compile(x) for x in pats], re.compile(logs))
            for label, builder, band, pats, logs in READS]
    found = {k[0]: [] for k in keys}
    for name in names:
        chosen = [k for k in keys if k[4].search(name)]
        if not chosen:
            continue
        for n, line in enumerate(_text(LOGS / name).splitlines(), 1):
            for label, _, _, pats, _ in chosen:
                for pat in pats:
                    for m in pat.finditer(line):
                        g = [float(x) for x in m.groups()]
                        found[label].append((g[0] / g[1] if len(g) == 2 else g[0], name, n))
    clean = True
    rows = [(k[0], k[1], k[2], found[k[0]]) for k in keys]
    rows.sort(key=lambda r: -max([v for v, _, _ in r[3]], default=-1.0) / r[2])
    for label, builder, band, hits in rows:
        if not hits:
            print(f"READ {label}: builder {builder}, no printed read found (band {band:g})")
            continue
        top, name, n = max(hits)
        room = f"{band / top:.1f} x" if top > 0 else "floor"
        zeros = sum(v == 0.0 for v, _, _ in hits)
        print(f"READ {label}: builder {builder}, largest {top:.3e} {name}:{n}, {top / band:.4f} x band "
              f"{band:g}, headroom {room}; {len(hits)} reads, {zeros} exactly 0")
        clean &= top <= band
    return clean


def main():
    totals_ok, found = census()
    _, unread = classify(found)
    if not control():
        sys.exit(3)
    if not reads():
        sys.exit("a read stands above its band")
    if not totals_ok:
        sys.exit(2)
    if unread:
        sys.exit(4)


if __name__ == "__main__":
    main()
