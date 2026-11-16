FROM dolfinx/dolfinx:v0.11.0

# Install h5py from source to match system HDF5
# This prevents version mismatch errors
RUN pip install --no-cache-dir --no-binary=h5py h5py

# Install additional dependencies
RUN pip install --no-cache-dir \
    gmsh \
    pyvista \
    scikit-rf==2.1.0 \
    meshio \
    tqdm \
    pydantic \
    pyyaml \
    black \
    isort \
    mypy \
    flake8 \
    pytest \
    pytest-xdist \
    mkdocs \
    mkdocs-material \
    "mkdocstrings[python]" \
    jupyter \
    ipywidgets \
    matplotlib

# CadQuery, for authoring test STEP files (GEO-34 step 1c; operator 2026-10-05).
# Isolated in its own venv: cadquery pulls its own OCCT (cadquery-ocp), vtk and
# numpy, and must never touch the numpy / vtk the DolfinX builds are linked
# against. Run it as `cq-python` (PYTHONPATH cleared so the DolfinX
# dist-packages on the compose PYTHONPATH cannot leak in). The solver side reads
# the STEP files through gmsh, never imports cadquery.
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3-venv libgl1 libxrender1 libxext6 \
 && rm -rf /var/lib/apt/lists/* \
 && python3 -m venv /opt/cadquery \
 && /opt/cadquery/bin/pip install --no-cache-dir --upgrade pip \
 && /opt/cadquery/bin/pip install --no-cache-dir cadquery \
 && printf '#!/bin/sh\nexec env -u PYTHONPATH /opt/cadquery/bin/python "$@"\n' > /usr/local/bin/cq-python \
 && chmod +x /usr/local/bin/cq-python \
 && cq-python -c "import cadquery; print('cadquery', cadquery.__version__)"

# Set up working directory
WORKDIR /workspace

# Copy source code
COPY . /workspace/

# Install package in editable mode
RUN pip install -e .

# Expose Jupyter port
EXPOSE 8888

# Default command
CMD ["/bin/bash"]
