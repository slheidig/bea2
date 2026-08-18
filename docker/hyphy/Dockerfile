# =============================================================================
# HyPhy 2.5.1 + python3 for the bea2 selection module
# =============================================================================
# Build:
#   docker build -t slheidig/hyphy:2.5.1 docker/hyphy
# HyPhy 2.5.1 predates bioconda's arm64 builds, so amd64 is pinned (runs under
# emulation on Apple Silicon; native on the clusters).
# =============================================================================
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "hyphy=2.5.1" "python>=3.10" && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# fail the build if the pinned toolchain isn't callable
RUN hyphy --version && python3 --version

WORKDIR /data
CMD ["bash"]
