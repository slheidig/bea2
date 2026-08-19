# =============================================================================
# HyPhy + python3 for the bea2 selection module   (slheidig/hyphy:2.5.101)
# =============================================================================
# Methods called by modules/hyphy.nf: FEL, FUBAR, MEME, contrast-fel, relax.
#
# 2.5.101 matches the HyPhy installed locally on the dev machine.
#
# Why not 2.5.1, as originally pinned:
#   * bioconda's hyphy 2.5.1 links libcrypto.so.1.1, and today's conda-forge
#     resolves openssl 3 -> the binary no longer starts.
#   * 2.5.1 ships SelectionAnalyses/FEL-contrast.bf, not contrast-fel.bf, so
#     `hyphy contrast-fel` does not resolve there.
#   * bioconda has no linux-aarch64 build below 2.5.60, which is what forced
#     the old amd64 platform pin.
# 2.5.101 has contrast-fel + RELAX, openssl 3, and native linux-64 AND
# linux-aarch64 builds -> genuine dual-arch image, no emulation needed.
#
# Build (from the repo root):
#   docker/build.sh hyphy           # amd64 + arm64 manifest, pushed
#   docker/build.sh --load hyphy    # single-arch local image for testing
# =============================================================================
FROM mambaorg/micromamba:1.5.8

ARG HYPHY_VERSION=2.5.101

# python is pinned rather than ">=3.10": the bin/*.py steps that run in this
# image (mask_stops, tag_foreground, parse_hyphy_persite, parse_relax) are
# stdlib-only, so any 3.x works — but an unpinned range silently drifted to
# python 3.14 and reproducibility beats novelty here.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "hyphy=${HYPHY_VERSION}" "python=3.12" && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# fail the build if the toolchain isn't callable, or if the two method files the
# pipeline needs are absent from this HyPhy build (2.5.1 had FEL-contrast.bf
# instead of contrast-fel.bf, which would only surface at pipeline runtime)
RUN hyphy --version && python3 --version && \
    prefix="$(dirname "$(dirname "$(command -v hyphy)")")" && \
    for bf in contrast-fel.bf RELAX.bf FEL.bf FUBAR.bf MEME.bf; do \
        find "$prefix" -name "$bf" -print -quit | grep -q . \
          || { echo "MISSING HyPhy batch file: $bf" >&2; exit 1; }; \
    done

WORKDIR /data
CMD ["bash"]
