# =============================================================================
# csubst image — cdskit, IQ-TREE 3, ete4, csubst  (slheidig/csubst:01)
# =============================================================================
# Used by: BACKALIGN, PRUNE_MSA, MAKE_FOREGROUND, MAKE_CONSTRAINT, IQTREE,
#          ROOT_TREE, CSUBST  (see nextflow.config)
#
# Toolchain: IQ-TREE 3 from conda (native linux-64 AND linux-aarch64); csubst
# (1.14.x) and cdskit from their GitHub repos so behaviour matches the
# csubst_trial results run (verified: csubst gives identical output with
# iqtree3 as with iqtree2). csubst pulls in ete4/numpy/scipy.
#
# DUAL-ARCH: every dependency has native amd64 + arm64 builds, so this image is
# built for both platforms in one manifest — no emulation, no cross tricks.
#
# Build (from the repo root):
#   docker/build.sh csubst              # amd64 only, loaded locally
#   docker/build.sh --push csubst       # amd64 + arm64 manifest, pushed
# =============================================================================
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="csubst"
LABEL org.opencontainers.image.description="csubst + cdskit + IQ-TREE 3 + ete4 for the bea2 pipeline"
LABEL org.opencontainers.image.source="https://github.com/kfuku52/csubst"

# System/toolchain layer: IQ-TREE binary + Python + build tools for csubst's
# Cython/ete4 extensions (installed via pip below).
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 "iqtree>=3" pip numpy cython git \
        c-compiler cxx-compiler && \
    micromamba clean --all --yes

# Activate the base env for the following RUN steps and at container runtime.
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# csubst + cdskit from source (pins to the repos used for the analysis).
RUN pip install --no-cache-dir "git+https://github.com/kfuku52/cdskit" && \
    pip install --no-cache-dir "git+https://github.com/kfuku52/csubst"

# Smoke test: the CLIs must resolve at build time.
# 3.x ships the binary as iqtree3; accept any of the names.
RUN csubst --version && cdskit --version && \
    { iqtree3 --version || iqtree2 --version || iqtree --version; } | head -1

# Nextflow stages inputs into the task work dir and mounts it; nothing is baked in.
WORKDIR /data
CMD ["bash"]
