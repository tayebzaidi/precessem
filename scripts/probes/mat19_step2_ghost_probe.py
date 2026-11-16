"""MAT-19 step 2 diagnostic (not a gate): ghost-cell counts on the MAT-15 parent mesh and tag-3 submesh.

W1 (20261005T183423Z_MAT-19-step2-W1.log:1977) read the (U3c) ghost-included Omega volume
equal to the owned one to 9.8e-15; this prints index_map(3).num_ghosts per rank and the
parent mesh's ghost mode evidence, so the control's failure has a measured cause.
"""
from mpi4py import MPI

from tests.validation.test_birdcage_sar_1g_rung import PHANTOM_RESOLUTION_1G_RUNG
from tests.validation.test_mat15_coil_driven_heating_phantom import _tag_submesh
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

comm = MPI.COMM_WORLD
sweep = build_four_port_sweep(phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG)
msh, ct = sweep["mesh"], sweep["cell_tags"]
sub, _ = _tag_submesh(msh, ct, PHANTOM_CELL_TAG)
pm, sm = msh.topology.index_map(3), sub.topology.index_map(3)
rows = comm.gather((comm.rank, pm.size_local, pm.num_ghosts, sm.size_local, sm.num_ghosts), root=0)
if comm.rank == 0:
    for r in rows:
        print(f"GHOSTPROBE rank {r[0]}: parent owned {r[1]} ghosts {r[2]}; sub owned {r[3]} ghosts {r[4]}",
              flush=True)
    print(f"GHOSTPROBE parent global {pm.size_global}, sub global {sm.size_global}", flush=True)
