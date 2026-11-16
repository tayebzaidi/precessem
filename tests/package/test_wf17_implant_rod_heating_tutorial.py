"""`WF-17` step 32 (§9 item 839, §10 U128): the implant-rod page's heating section.

``docs/guide/implant_rod_in_a_study.md`` "Heating the phantom around the rod"
carries a ``implant-rod-heating`` Python block (the rod ``Study`` and
``temperature(..., balls=..., drive="P1")``) and a ``implant-rod-heating-file``
block (``run_study_file`` on the tracked ``tests/package/wf15/rod_heating.yaml``).

Brief: docs/planning/chunks/WF-17.md, "§9 item 839" (review ruling 2026-11-08
03:00: serial on item 884, whose W1 record is the comparand).

Anchors (bands imported, never restated):
- (Y) the two blocks' tip θ̄(t*) equal at ``IDENTITY_BAND``;
- (R) each, per watt of phantom power (θ / ``source_power_w``, linear in the
  source scale), equals item 884's record printed in
  ``20261108T223537Z_WF-15-s14p-W1.log:5680–5681`` rounded to ``.9e``, within one
  unit of that print's last digit (``PRINT_LAST_DIGIT``);
- (E) ``max_balance_residual`` ≤ ``E_BAND``;
- (S2) the volume slope of one BE step at ``DT_1`` ≤ ``S2_BAND`` (item 884's
  formula, ``G P`` = the step's ``source_power_w``);
- (D) 0 dead docrefs.
Control (asserted, no mesh): the heating block with ``implant_rod`` removed is
refused before any mesh.  Measured, not as the brief predicted: the refusal is
the ``ValueError`` "drive is defined on the implant-rod four-port Study only",
because ``Study.temperature`` checks ``drive=`` before the port-Study
``TypeError`` (``study.py:1644–1645``).  (c-page): a removed or doubled tag raises.

The solve node runs only with ``FEM_EM_WF17_S32=1`` (complex build, ``-n 2``).
"""

from __future__ import annotations

import ast
import math
import os
import time

import pytest
import yaml
from mpi4py import MPI

from fem_em_solver import FourPortBirdcageGeometry, RobinCooling, Study, TransientSpec

from tests.complex_mode import IS_COMPLEX
from tests.mesh.test_birdcage_implant_rod import ROD
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.package.test_wf15_study_implant_rod import ROD_EXACT_CELL_RECORD, _phantom
from tests.package.test_wf15_study_rod_heating import _be, _cn, _yaml_mapping
from tests.package.test_wf17_first_study_tutorial import extract_tagged_block
from tests.package.test_wf17_implant_rod_tutorial import (
    PAGE,
    ROD_CALL,
    TAG,
    _exec_block,
    _page_rod,
    _scratch,
)
from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ
from tests.validation.test_imp3_rod_birdcage import TIPS
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    DT_1,
    E_BAND,
    H_GEL,
    K_GEL,
    RHO_C,
    S2_BAND,
    T_INF,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

S32_ENV = "FEM_EM_WF17_S32"
HEAT_TAG = "implant-rod-heating"
HEAT_FILE_TAG = "implant-rod-heating-file"
HEAT_BANDS_TAG = "implant-rod-heating-bands"
OFF_CALL = "FourPortBirdcageGeometry()"
DRIVE = "P1"
#: Item 884's tip θ̄(t*) per watt (arm R), printed `.12e` at
#: `20261108T223537Z_WF-15-s14p-W1.log:5680–5681` (2.840168784615e+00 /
#: 2.783507266272e+00 K), rounded here to the `.9e` print the brief names.
RECORD_9E = {"top": "2.840168785e+00", "bottom": "2.783507266e+00"}
#: One unit of the `.9e` print's last digit, relative to its decade.
PRINT_LAST_DIGIT = 1.0e-9


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S32_ENV, "").strip()
    if raw not in ("", "1"):
        raise ValueError(f"{S32_ENV}={raw!r}: only unset/empty or '1' are accepted")
    return raw == "1"


def _assigned(block: str, name: str):
    for node in ast.parse(block).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == name:
            return node.value
    raise LookupError(f"no `{name} = ...` in the block")


def _is_call(node, attr: str, owner: str | None = None) -> bool:
    if not isinstance(node, ast.Call):
        return False
    f = node.func
    if owner is None:
        return isinstance(f, ast.Name) and f.id == attr
    return (isinstance(f, ast.Attribute) and f.attr == attr
            and isinstance(f.value, ast.Name) and f.value.id == owner)


def _em_free(block: str, through_heat: bool = False) -> ast.Module:
    """The block with ``study.sar(...)`` … ``P = ...`` replaced by ``P = 1.0``.

    Everything from ``heat = ...`` on is dropped, unless ``through_heat`` keeps
    that one statement (the (c-off) control).  No EM solve, so no mesh.
    """
    body = ast.parse(block).body
    i_sar = [i for i, n in enumerate(body)
             if isinstance(n, ast.Expr) and _is_call(n.value, "sar", "study")]
    i_p = [i for i, n in enumerate(body)
           if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "P"]
    i_heat = [i for i, n in enumerate(body)
              if isinstance(n, ast.Assign) and getattr(n.targets[0], "id", "") == "heat"]
    assert len(i_sar) == 1 and len(i_p) == 1 and len(i_heat) == 1, (i_sar, i_p, i_heat)
    assert i_sar[0] < i_p[0] < i_heat[0], (i_sar, i_p, i_heat)
    stop = i_heat[0] + 1 if through_heat else i_heat[0]
    new = body[:i_sar[0]] + ast.parse("P = 1.0").body + body[i_p[0] + 1:stop]
    return ast.fix_missing_locations(ast.Module(body=new, type_ignores=[]))


def _head(block: str, through_heat: bool = False) -> dict:
    """The EM-free head (``_em_free``) exec'd: no mesh."""
    ns = {"__name__": "__main__"}
    exec(compile(_em_free(block, through_heat), f"{PAGE}#{HEAT_TAG}-head", "exec"), ns)
    return ns


def _bands_text() -> str:
    return (f"Python vs the study file, tip temperature at t*: {IDENTITY_BAND:.0e}\n"
            f"energy balance (E): {E_BAND:.0e}\n"
            f"volume slope (S2): {S2_BAND:.0e}")


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s32_heating_blocks_and_control(monkeypatch):
    text = PAGE.read_text()
    block = extract_tagged_block(text, HEAT_TAG)
    file_block = extract_tagged_block(text, HEAT_FILE_TAG)
    bands = extract_tagged_block(text, HEAT_BANDS_TAG, lang="text").strip()

    calls = []

    def _build(self, comm=MPI.COMM_WORLD):
        calls.append(self)
        raise AssertionError("a mesh was built")

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", _build)

    # The block's inputs are item 884's (its Study / cooling / CN / balls / drive).
    assert block.count(ROD_CALL) == 1
    assert _page_rod(block) == ROD
    # (In) per watt (item 908): one EM `study.sar` on P1, one `mean_sar` on the
    # phantom tag read as dissipated power, the cooling scaled by 1 / P.
    assert block.count("Study(") == 1
    assert block.count("study.sar(") == 1 and block.count("mean_sar(") == 1
    tree = ast.parse(block)
    sar_calls = [n for n in ast.walk(tree) if _is_call(n, "sar", "study")]
    assert len(sar_calls) == 1
    drives = {k.arg: k.value for k in sar_calls[0].keywords}["drives"]
    assert ast.literal_eval(drives) == [DRIVE], ast.unparse(drives)
    p_val = _assigned(block, "P")
    assert isinstance(p_val, ast.Subscript) and ast.literal_eval(p_val.slice) == "dissipated_power_w"
    assert _is_call(p_val.value, "mean_sar")
    kw = {k.arg: k.value for k in p_val.value.keywords}
    assert ast.unparse(kw["subdomain_ids"]) == "study.geometry.phantom_tag", ast.unparse(kw["subdomain_ids"])
    cool = _assigned(block, "cooling")
    assert _is_call(cool, "RobinCooling")
    assert ast.unparse({k.arg: k.value for k in cool.keywords}["source_scale"]) == "1.0 / P"
    assert FourPortBirdcageGeometry().phantom_tag == PHANTOM_CELL_TAG
    ns = _head(block)
    st, cooling, tr = ns["study"], ns["cooling"], ns["transient"]
    assert isinstance(st, Study) and st.geometry == FourPortBirdcageGeometry(implant_rod=ROD)
    assert st.material == _phantom() and st.degree == 1
    assert float(st.frequency_hz) == STEP2D_FREQUENCY_HZ
    assert isinstance(cooling, RobinCooling)
    assert (cooling.k_w_per_m_k, cooling.h_w_per_m2_k, cooling.t_inf_c, cooling.degree,
            cooling.source_scale) == (K_GEL, H_GEL, T_INF, 1, 1.0)
    cn = _cn()
    assert isinstance(tr, TransientSpec)
    assert abs(tr.dt_s / cn.dt_s - 1.0) <= IDENTITY_BAND, (tr.dt_s, cn.dt_s)
    assert (tr.n_steps, tr.rho_c_j_per_m3_k, tr.theta, tr.startup_half_steps) == (
        cn.n_steps, cn.rho_c_j_per_m3_k, cn.theta, cn.startup_half_steps)
    assert [tuple(t) for t in ns["tips"]] == list(TIPS.values())
    assert 'drive="P1"' in block and '"mass_g": 1.0' in block
    assert 'run_study_file("rod_heating.yaml")' in file_block
    assert "](../../tests/package/wf15/rod_heating.yaml)" in text
    assert bands == _bands_text(), bands
    # No YAML fence in the heating section (the file is linked, not inlined).
    section = text.split("## Heating the phantom around the rod")[1].split("## What this")[0]
    assert "```yaml" not in section
    _p0(f"\n[WF-17 s32] (In) rod {ROD}; dt {tr.dt_s!r} x {tr.n_steps}; tips {ns['tips']}; "
        f"bands {bands!r}")

    # Control (asserted, ruling 2026-11-09 03:00): the EM-free head with the rod
    # removed, run through the `heat = study.temperature(` statement, is refused
    # before any mesh.
    off = block.replace(ROD_CALL, OFF_CALL)
    assert off.count(OFF_CALL) == 1 and "implant_rod" not in off
    with pytest.raises(ValueError, match="drive is defined on the implant-rod four-port") as err:
        _head(off, through_heat=True)
    assert calls == [], calls
    _p0(f"[WF-17 s32] (c-off) refused before any mesh: {type(err.value).__name__}: {err.value}")
    # Printed, not asserted: the *full* rod-free block meshes at `study.sar(`
    # (a rod-free four-port `sar` is legal), so the build spy is reached there.
    try:
        exec(compile(off, f"{PAGE}#{HEAT_TAG}-off-full", "exec"), {"__name__": "__main__"})
        full = "no exception"
    except Exception as exc:  # noqa: BLE001 - a reading, printed
        full = f"{type(exc).__name__}: {exc}"
    _p0(f"[WF-17 s32] (c-off full, printed) rod-free full block: {full}; build spy calls "
        f"{len(calls)}")
    calls.clear()

    # (N) the page's limits stay; no external tool named.
    for token in ("No two-port drive", "No 1 g rung", "No published comparand", "No tissue"):
        assert token in text, token
    for token in ("AED", "HFSS"):
        assert token not in text, token
    # The step-25 node's `Study(` count stays on its own tag.
    assert extract_tagged_block(text, TAG).count("Study(") == 1

    # (c-page): each heating tag removed or doubled raises.
    for tag, body in ((HEAT_TAG, block), (HEAT_FILE_TAG, file_block)):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```python {tag}\n", "```python other\n"), tag)
        with pytest.raises(LookupError):
            extract_tagged_block(text + f"\n```python {tag}\n{body}```\n", tag)


# ------------------------------------------------------------- solve window (-n 2)


def _tips(balls) -> dict:
    names = list(TIPS)
    return {names[i]: float(b["theta_k"][-1]) for i, b in enumerate(balls)}


def test_wf17_s32_heating_blocks(monkeypatch):
    if not _knob():
        pytest.skip(f"{S32_ENV} unset: the heating section's solve window is environment-selected")
    import numpy as np
    import ufl
    from dolfinx import default_scalar_type, fem

    import fem_em_solver.study as study_module

    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port Study needs the complex build"
    assert comm.size == 2, comm.size
    text = PAGE.read_text()

    t0 = time.perf_counter()
    ns, out = _exec_block(monkeypatch, extract_tagged_block(text, HEAT_TAG), "heating")
    t_py = time.perf_counter() - t0
    heat, st = ns["heat"], ns["study"]

    # (S2) one BE step at DT_1 on the block's Study (fields cached: no second EM solve).
    t0 = time.perf_counter()
    be = st.temperature(ns["cooling"], transient=_be(), balls={"mass_g": 1.0, "centres": ns["tips"]},
                        drive=DRIVE, comm=comm)
    res = st.fields_result
    sub, _pid = study_module._tag_submesh(res.mesh, res.cell_tags, PHANTOM_CELL_TAG)
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(np.real(comm.allreduce(fem.assemble_scalar(fem.form(one * ufl.dx(domain=sub))),
                                         op=MPI.SUM)))
    s2 = abs(RHO_C * be.volume_mean_k[1] / (DT_1 * be.source_power_w / v_sub) - 1.0)
    t_be = time.perf_counter() - t0

    t0 = time.perf_counter()
    # `_exec_block` wipes its own scratch, so the tracked file is copied beside it
    # (W1 r1, `20261109T014948Z_WF-17-s32-W1.log:1910`, read "no such file") and
    # the block's one file name is replaced by that copy's absolute path.
    # Item 908: the copy carries `thermal.source_scale = 1 / P` with the block's
    # P (gitignored scratch only; no tracked file holds the float).
    scratch = _scratch("heating-file-src")
    if comm.rank == 0:
        m = _yaml_mapping()
        m["thermal"]["source_scale"] = 1.0 / float(ns["P"])
        (scratch / "rod_heating.yaml").write_text(yaml.safe_dump(m, sort_keys=False))
    comm.Barrier()
    file_block = extract_tagged_block(text, HEAT_FILE_TAG)
    ns_f, out_f = _exec_block(
        monkeypatch, file_block.replace('"rod_heating.yaml"', repr(str(scratch / "rod_heating.yaml")), 1),
        "heating-file")
    assert file_block.count('"rod_heating.yaml"') == 1
    t_file = time.perf_counter() - t0
    summ = ns_f["run"].summary["temperature"]
    # (Scale) precondition: the file ran at exactly the block's scale.
    scale_ok = comm.allreduce(int(summ["source_scale"] != ns["cooling"].source_scale), op=MPI.SUM) == 0
    _p0(f"[WF-17 s32] (Scale) file {summ['source_scale']!r} vs block {ns['cooling'].source_scale!r}: "
        f"{'equal' if scale_ok else 'DIFFERENT'} (precondition, ASSERTED)")
    assert scale_ok, (summ["source_scale"], ns["cooling"].source_scale)
    _p0(f"[WF-17 s32] block P {float(ns['P']):.12e} W vs item 884 1.646387247816e-06 W "
        f"(`20261108T223537Z_WF-15-s14p-W1.log:5674`); heat.source_power_w "
        f"{float(heat.source_power_w):.15e} (884 (B) read 9.999999999999829e-01)")
    _p0("[WF-17 s32] predictions (item 884, same scale, same mesh; printed, never asserted): "
        "(E) CN ~4.4e-14, BE ~1.5e-11, file CN same order as CN; (Y) ~1e-14; "
        "(R) raw 3.85e-10 / 2.72e-10; (S2) ~4.8353e-06")

    th_py = _tips(heat.balls)
    th_file = _tips(summ["balls"])
    rel_y = {n: abs(th_file[n] / th_py[n] - 1.0) for n in TIPS}
    per_w = {n: th_py[n] / float(heat.source_power_w) for n in TIPS}
    per_w_file = {n: th_file[n] / float(summ["source_power_w"]) for n in TIPS}
    unit = {n: PRINT_LAST_DIGIT * 10.0 ** math.floor(math.log10(abs(float(RECORD_9E[n])))) for n in TIPS}
    dev = {n: abs(float(f"{per_w[n]:.9e}") - float(RECORD_9E[n])) for n in TIPS}
    dev_raw = {n: abs(per_w[n] - float(RECORD_9E[n])) for n in TIPS}

    _p0(f"[WF-17 s32] -n {comm.size}: heating block {t_py:.1f} s, BE step {t_be:.1f} s, "
        f"file block {t_file:.1f} s")
    _p0(f"[WF-17 s32] (Cells) block {heat.cells}, file {ns_f['run'].summary['cells']}, record "
        f"{ROD_EXACT_CELL_RECORD}")
    _p0(f"[WF-17 s32] P_phantom block {float(heat.source_power_w):.12e} W, file "
        f"{float(summ['source_power_w']):.12e} W")
    for n in TIPS:
        _p0(f"[WF-17 s32] (Y) {n}: block theta(t*) {th_py[n]:.15e} K vs file {th_file[n]:.15e} K: "
            f"rel {rel_y[n]:.3e} (IDENTITY_BAND {IDENTITY_BAND:.0e}, ASSERTED)")
        _p0(f"[WF-17 s32] (R) {n}: per W {per_w[n]:.12e} -> {per_w[n]:.9e} vs item 884 "
            f"{RECORD_9E[n]}: |.9e - rec| {dev[n]:.1e}, raw {dev_raw[n]:.3e} (one last-digit unit "
            f"{unit[n]:.0e}, ASSERTED); file per W {per_w_file[n]:.12e}")
    _p0(f"[WF-17 s32] (E) CN {heat.max_balance_residual:.3e}, file {summ['max_balance_residual']:.3e}, "
        f"BE {be.max_balance_residual:.3e} (E_BAND {E_BAND:.0e}, ASSERTED)")
    _p0(f"[WF-17 s32] (S2) volume slope |r-1| {s2:.4e} (S2_BAND {S2_BAND:.0e}, ASSERTED; item 884 "
        f"read 4.8353e-06); V_sub {v_sub:.12e}")
    gates = {
        "cells": int(heat.cells) == ROD_EXACT_CELL_RECORD,
        "stdout": (f"cells = {ROD_EXACT_CELL_RECORD}" in out.splitlines()) if comm.rank == 0 else True,
        "E cn": heat.max_balance_residual <= E_BAND,
        "E file": summ["max_balance_residual"] <= E_BAND,
        "E be": be.max_balance_residual <= E_BAND,
        "S2": s2 <= S2_BAND,
    }
    for n in TIPS:
        gates[f"Y {n}"] = rel_y[n] <= IDENTITY_BAND
        gates[f"R {n}"] = dev_raw[n] <= unit[n]
    red = sorted(set(comm.allreduce([k for k, v in gates.items() if not v], op=MPI.SUM)))
    _p0(f"[WF-17 s32] gates red: {red or 'none'} of {len(gates)}")
    assert not red, red
