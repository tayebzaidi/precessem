"""Thermal solvers (Phase 8). `MAT-9` steps 1–2: steady and transient conduction;
`MAT-5` step 1: Picard for a temperature-dependent source."""

from .probe import ProbeBall, ProbeSeries, ball_average, probe_series
from .steady import (
    PicardConductionResult,
    SteadyConductionResult,
    solve_steady_conduction,
    solve_steady_conduction_picard,
)
from .transient import TransientConductionResult, solve_transient_conduction

__all__ = [
    "SteadyConductionResult",
    "solve_steady_conduction",
    "PicardConductionResult",
    "solve_steady_conduction_picard",
    "TransientConductionResult",
    "solve_transient_conduction",
    "ProbeBall",
    "ProbeSeries",
    "ball_average",
    "probe_series",
]
