# =============================================================================
# b2bTools + scientific python  (slheidig/og_b2b_pca)
# =============================================================================
# Used by: B2BTOOLS, HYPHY_BYSITE, MAP_TO_MSA, COMBINE_TABLE, PLOT_OG,
#          GLOBAL_STATS, SPLIT_CUSTOM_PREDICTIONS  (see nextflow.config)
#
# Same recipe as proteinfam_pca_nf/docker/b2btools/Dockerfile, with the versions
# pinned (the image currently on Docker Hub was built with unpinned pip installs
# and resolved to b2bTools 3.0.8 / python 3.12).
#
# DUAL-ARCH: python:3.12-slim and every wheel here exist for amd64 + arm64, so
# both platforms build natively in one manifest.
#
# Build (from the repo root):
#   docker/build.sh og_b2b_pca           # amd64 only, loaded locally
#   docker/build.sh --push og_b2b_pca    # amd64 + arm64 manifest, pushed
# =============================================================================
FROM python:3.12-slim

LABEL org.opencontainers.image.title="og_b2b_pca"
LABEL org.opencontainers.image.version="3.0.8"
LABEL org.opencontainers.image.description="b2bTools + pandas + matplotlib + scikit-learn for the bea2 pipeline"
LABEL org.opencontainers.image.source="https://pypi.org/project/b2bTools/"

# gcc/build-essential: some wheels are source-only on arm64
RUN apt-get update && apt-get install -y --no-install-recommends \
        gcc \
        build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    b2bTools==3.0.8 \
    pandas==2.2.3 \
    matplotlib==3.9.4 \
    scikit-learn==1.6.1 \
    numpy==1.26.4

# fail the build if the stack does not import
RUN python -c "import b2bTools, pandas, matplotlib, sklearn, numpy; print('ok')"

# non-root user; /data must be writable by it (Nextflow mounts the task dir there)
RUN useradd --create-home --shell /bin/bash b2buser && \
    mkdir -p /data && chown b2buser:b2buser /data
USER b2buser
WORKDIR /data
CMD ["b2bTools", "--help"]
