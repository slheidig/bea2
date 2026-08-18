# =============================================================================
# AIUPred (disorder / binding prediction) — CPU-only
# https://github.com/doszilab/AIUPred  |  Erdos & Dosztanyi 2024, NAR
# =============================================================================
# Build:
#   docker build -t slheidig/aiupred:01 docker/aiupred
# Smoke test:
#   docker run --rm slheidig/aiupred:01 python -c "import aiupred_lib; print('ok')"
# The pipeline calls bin/run_aiupred.py which imports aiupred_lib from
# $AIUPRED_PATH. Networks/weights ship inside the repository.
# =============================================================================
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends git \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# CPU-only torch keeps the image small; no GPU assumptions
RUN pip install --no-cache-dir torch --index-url https://download.pytorch.org/whl/cpu && \
    pip install --no-cache-dir numpy

RUN git clone --depth 1 https://github.com/doszilab/AIUPred /opt/aiupred

ENV AIUPRED_PATH=/opt/aiupred
ENV PYTHONPATH=/opt/aiupred

# fail the build if the library does not import
RUN python -c "import aiupred_lib" || python -c "import sys; sys.path.insert(0,'/opt/aiupred'); import aiupred_lib"

WORKDIR /data
CMD ["python", "-c", "import aiupred_lib; print('AIUPred ready')"]
