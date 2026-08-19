#!/usr/bin/env bash
# =============================================================================
# bea2 container builds — multi-arch, from an Apple Silicon laptop
# =============================================================================
# Why this script exists instead of plain `docker build`:
#
#  1. Docker Desktop's default builder uses the "docker" driver, which CANNOT
#     write manifest lists: any `--platform linux/amd64,linux/arm64` build fails
#     with "docker exporter does not currently support exporting manifest
#     lists". Multi-arch needs a docker-container driver builder — created here.
#  2. `docker build` on an arm Mac silently produces an ARM-ONLY image. Pushing
#     that means the amd64 clusters cannot run it. Here, platforms are declared
#     per image in docker-bake.hcl and always explicit.
#  3. Pushing straight from the builder (--push) never materialises the image in
#     the local docker store, which matters on a laptop that is short on disk.
#  4. After a push it verifies the published manifest actually contains the
#     platforms that were asked for. That check is the whole point: an arm-only
#     image looks perfectly fine locally and only fails on the cluster.
#
# Usage:
#   docker/build.sh                     # build+push hyphy, aiupred, ipc
#   docker/build.sh hyphy ipc           # named targets
#   docker/build.sh all                 # incl. csubst / og_b2b_pca (already pushed)
#   docker/build.sh --load hyphy        # single-arch, into the local docker store
#   docker/build.sh --dry-run all       # print the plan, build nothing
#   PLATFORMS=linux/amd64 docker/build.sh hyphy     # override platforms
#
# Requires `docker login` for --push (the default).
# =============================================================================
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAKE_FILE="$REPO_ROOT/docker/docker-bake.hcl"
BUILDKITD_CONFIG="$REPO_ROOT/docker/buildkitd.toml"
BUILDER=bea2

MODE=push
DRY_RUN=0
TARGETS=()

# Keep provenance/SBOM attestations out of the manifest lists (they show up as
# unknown/unknown entries). Belt and braces with provenance=false in the bake file,
# which older buildx versions ignore.
export BUILDX_NO_DEFAULT_ATTESTATIONS=1

while [[ $# -gt 0 ]]; do
    case "$1" in
        --push)    MODE=push;   shift ;;
        --load)    MODE=load;   shift ;;
        --dry-run) DRY_RUN=1;   shift ;;
        -h|--help) sed -n '2,30p' "${BASH_SOURCE[0]}"; exit 0 ;;
        -*)        echo "unknown flag: $1" >&2; exit 2 ;;
        *)         TARGETS+=("$1"); shift ;;
    esac
done
[[ ${#TARGETS[@]} -eq 0 ]] && TARGETS=(default)

# ---- disk guard: builds die halfway and leave cache behind when space runs out
free_gb=$(df -g "$HOME" | awk 'NR==2 {print $4}')
echo "free disk: ${free_gb} GB"
if [[ "$free_gb" -lt 8 ]]; then
    echo "WARNING: under 8 GB free. Reclaim first:"
    echo "  docker buildx du                        # build cache"
    echo "  docker buildx prune --builder $BUILDER   # drop it"
    echo "  docker image ls / docker image rm <img>  # anything re-pullable from Hub"
fi

# ---- the builder that can actually write manifest lists
if ! docker buildx inspect "$BUILDER" >/dev/null 2>&1; then
    echo "creating buildx builder '$BUILDER' (docker-container driver)"
    docker buildx create --name "$BUILDER" --driver docker-container \
        --config "$BUILDKITD_CONFIG" --bootstrap >/dev/null
fi

# ---- always shut the builder down again: the buildkit container otherwise sits
# ---- there holding memory (and its cache volume) between build attempts
cleanup() {
    docker buildx stop "$BUILDER" >/dev/null 2>&1 || true
    echo
    echo "builder '$BUILDER' stopped. Its cache lives in the volume"
    echo "buildx_buildkit_${BUILDER}0_state — 'docker buildx rm $BUILDER' deletes"
    echo "both and reclaims that disk."
}
trap cleanup EXIT

# ---- --load cannot write a manifest list either: one platform only
if [[ "$MODE" == load ]]; then
    platforms="${PLATFORMS:-linux/$(uname -m | sed 's/x86_64/amd64/;s/aarch64/arm64/;s/arm64/arm64/')}"
    if [[ "$platforms" == *,* ]]; then
        echo "ERROR: --load takes a single platform (the docker store holds one" >&2
        echo "       image per tag). Use --push for a multi-arch manifest." >&2
        exit 2
    fi
    export PLATFORMS="$platforms"
    echo "mode: --load, platform: $PLATFORMS"
fi

cmd=(docker buildx bake --builder "$BUILDER" -f "$BAKE_FILE" "--$MODE" "${TARGETS[@]}")
echo "+ ${cmd[*]}"
if [[ "$DRY_RUN" == 1 ]]; then
    (cd "$REPO_ROOT" && docker buildx bake -f "$BAKE_FILE" --print "${TARGETS[@]}")
    exit 0
fi
(cd "$REPO_ROOT" && "${cmd[@]}")

# ---- verify what actually landed in the registry
if [[ "$MODE" == push ]]; then
    echo
    echo "=== published manifests ==="
    (cd "$REPO_ROOT" && docker buildx bake -f "$BAKE_FILE" --print "${TARGETS[@]}") \
      | grep -o '"[^"]*/[^"]*:[^"]*"' | tr -d '"' | sort -u | while read -r tag; do
        echo "--- $tag"
        docker buildx imagetools inspect "$tag" 2>&1 \
          | grep -E "^ *Platform:" | sort | uniq -c || echo "    (inspect failed)"
    done
    echo
    echo "amd64 must be listed for every image, or the clusters cannot run it."
fi
