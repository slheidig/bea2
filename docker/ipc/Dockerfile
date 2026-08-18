# =============================================================================
# IPC — Isoelectric Point Calculator (Kozlowski 2016, http://isoelectric.org)
# =============================================================================
# Build:
#   docker build -t slheidig/ipc:01 docker/ipc
#
# The pipeline's bin/run_ipc.py computes pI with the published IPC pKa sets
# (built-in, stdlib-only) and per-residue charge at a chosen pH, so this image
# only needs python. To use the OFFICIAL IPC script instead, download it from
# isoelectric.org (academic licence), COPY it in below, set IPC_CMD, and run
# the pipeline with --use-cli added to run_ipc.py in modules/predictors.nf.
# =============================================================================
FROM python:3.11-slim

# --- optional: official IPC CLI -------------------------------------------
# COPY ipc-2.0/ /opt/ipc/
# ENV IPC_CMD="python /opt/ipc/ipc.py"
# ---------------------------------------------------------------------------

RUN useradd --create-home ipcuser
USER ipcuser
WORKDIR /data
CMD ["python", "--version"]
