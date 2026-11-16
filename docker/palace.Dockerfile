# Palace — the independent second reference for OPS-61 (operator directive
# 2026-09-20). Never the solver: nothing in src/ imports or links it.
#
# Why a CMake superbuild and not Spack (2026-09-28, interactive session):
#   - Spack's builtin `palace` package stops at 0.16.0; 0.18.x lives only in
#     Palace's in-tree spack_repo, and upstream's container recipe
#     (.github/actions/build-container) is a CI template tracking `develop`,
#     with a token-gated ghcr buildcache and a gcc-15 PPA workaround (#731).
#   - The superbuild pins, per Palace release, the MFEM / hypre / PETSc-SLEPc /
#     MUMPS / SuperLU_DIST / libCEED versions upstream tested together, in one
#     `cmake && make`.
#
# Build (needs network; ~30-90 min, capped at 12 cores on the shared box):
#   docker build -f docker/palace.Dockerfile -t fem-em-palace:0.18.1 docker
# Smoke (step 0):
#   docker run --rm fem-em-palace:0.18.1 bash -lc \
#     'cd /opt/palace/examples/spheres && palace -np 2 spheres.json'

ARG PALACE_VERSION=v0.18.1

FROM ubuntu:24.04 AS build
ARG PALACE_VERSION
ARG BUILD_JOBS=12
ENV DEBIAN_FRONTEND=noninteractive

# ubuntu 24.04: gcc 13, cmake 3.28 (Palace needs >= 3.24), OpenMPI 4.1.
# python3 is for PETSc's configure; git/curl fetch the superbuild's sources.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gfortran cmake git curl ca-certificates pkg-config \
        python3 patch unzip \
        libopenmpi-dev openmpi-bin libopenblas-dev liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 --branch "${PALACE_VERSION}" \
        https://github.com/awslabs/palace.git /src/palace

# CPU only: MAGMA is a GPU library and OpenMP stays off (MPI ranks only, as
# in our own runs). MUMPS on beside the default SuperLU_DIST, so a direct-
# solver disagreement can be split from a model disagreement. SLEPc (default)
# keeps the eigenmode driver available for a later TH-17 cross-check.
WORKDIR /src/palace/build
RUN cmake .. \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/opt/palace \
        -DPALACE_WITH_OPENMP=OFF \
        -DPALACE_WITH_MAGMA=OFF \
        -DPALACE_WITH_CUDA=OFF \
        -DPALACE_WITH_HIP=OFF \
        -DPALACE_WITH_MUMPS=ON \
        -DPALACE_WITH_SUPERLU=ON \
        -DPALACE_WITH_SLEPC=ON \
    && make -j"${BUILD_JOBS}" \
    && cp -r /src/palace/examples /opt/palace/examples


FROM ubuntu:24.04
ARG PALACE_VERSION
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        openmpi-bin libopenmpi3t64 libopenblas0 liblapack3 \
        libgfortran5 libgomp1 libstdc++6 python3 ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /opt/palace /opt/palace

# Fail the build here, not at the first run, if the runtime stage is missing
# a shared library the build stage had.
RUN set -e; for f in /opt/palace/bin/palace-*.bin; do \
        if LD_LIBRARY_PATH=/opt/palace/lib:/opt/palace/lib64 ldd "$f" | grep "not found"; then \
            echo "missing runtime library for $f" >&2; exit 1; \
        fi; \
    done

ENV PATH=/opt/palace/bin:$PATH \
    LD_LIBRARY_PATH=/opt/palace/lib:/opt/palace/lib64 \
    OMP_NUM_THREADS=1 \
    OMPI_ALLOW_RUN_AS_ROOT=1 \
    OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    PALACE_VERSION=${PALACE_VERSION}

WORKDIR /workspace
CMD ["/bin/bash"]
