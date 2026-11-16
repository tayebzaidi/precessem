"""`WF-15` cost probe (§9 item 240 after W2's 590 s kill): one degree-2 Study
solve on `MAT-4` step 1's fixture, phase timings flushed as they happen so a
killed run still says where the time went.  Measurement only, no assertion.
"""

import resource
import time

from mpi4py import MPI

from fem_em_solver import GelledSalinePhantomMaterial, QuasiStaticSphereField, SphereInBoxGeometry, Study

comm = MPI.COMM_WORLD
R = 0.01
t0 = time.perf_counter()


def say(msg):
    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM) / 1024.0**2
    if comm.rank == 0:
        print(f"[WF-15 probe] t = {time.perf_counter() - t0:7.1f} s  rss(sum) = {rss:6.2f} GiB  {msg}", flush=True)


geom = SphereInBoxGeometry(R, 0.02, R / 10.0, R / 5.0, sphere_tag=1)
study = Study(
    geom,
    GelledSalinePhantomMaterial(0.57, 78.0, 64.0e6),
    QuasiStaticSphereField(e0=1.0),
    frequency_hz=64.0e6,
)
say(f"degree = {study.degree}; solving")
res = study.solve(comm)
V = res.fields.e_complex.function_space
say(f"solve done: {V.dofmap.index_map.size_global} dofs, "
    f"{res.mesh.topology.index_map(3).size_global} cells")
sar = study.mean_sar()["mean_sar_w_per_kg"]
say(f"mean SAR = {sar:.10e}")
