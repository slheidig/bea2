# =============================================================================
# bea2 container images — buildx bake definition
# =============================================================================
# One place that says, per image, WHICH PLATFORMS it is valid for. Use
# docker/build.sh (it sets up the right builder) rather than calling bake
# directly; see that script for why.
#
#   docker/build.sh                 # build+push the images that need building
#   docker/build.sh hyphy           # one target
#   docker/build.sh --load hyphy    # single-arch local image for testing
#   docker/build.sh all             # every target, incl. the ones already pushed
# =============================================================================

variable "REGISTRY" { default = "docker.io/slheidig" }

# Overridable platform list, e.g. PLATFORMS=linux/amd64 docker/build.sh hyphy
# Empty = use each target's own default.
variable "PLATFORMS" { default = "" }

variable "HYPHY_VERSION"      { default = "2.5.101" }
variable "IPC_TAG"            { default = "01" }
variable "CSUBST_TAG"         { default = "01" }
variable "OG_B2B_PCA_TAG"     { default = "latest" }

# amd64 is the platform that must always work (the clusters); arm64 is the
# convenience platform for local testing on Apple Silicon.
DUAL  = PLATFORMS != "" ? split(",", PLATFORMS) : ["linux/amd64", "linux/arm64"]
AMD64 = PLATFORMS != "" ? split(",", PLATFORMS) : ["linux/amd64"]

# Images we build ourselves and that still need pushing.
# AIUPred is NOT here: the pipeline uses the authors' official image
# ghcr.io/doszilab/aiupred:cpu (see nextflow.config), so there is nothing to build.
group "default" {
  targets = ["hyphy", "ipc"]
}

# Everything, including images already on Docker Hub — only needed for a rebuild.
group "all" {
  targets = ["hyphy", "ipc", "og_b2b_pca", "csubst"]
}

target "_common" {
  # No provenance/SBOM attestations: they add unknown/unknown entries to the
  # manifest list that buy nothing here and only confuse apptainer pulls.
  provenance = false
  sbom       = false
}

target "hyphy" {
  inherits  = ["_common"]
  context   = "docker/hyphy"
  args      = { HYPHY_VERSION = HYPHY_VERSION }
  tags      = ["${REGISTRY}/hyphy:${HYPHY_VERSION}"]
  platforms = DUAL           # 2.5.97 has native linux-64 + linux-aarch64 conda builds
}

# No aiupred target: the pipeline uses the authors' official amd64 image
# ghcr.io/doszilab/aiupred:cpu (published with github.com/doszilab/AIUPred-NF).

target "ipc" {
  inherits  = ["_common"]
  context   = "docker/ipc"
  tags      = ["${REGISTRY}/ipc:${IPC_TAG}"]
  platforms = DUAL           # stdlib-only python
}

target "og_b2b_pca" {
  inherits  = ["_common"]
  context   = "docker/og_b2b_pca"
  tags      = ["${REGISTRY}/og_b2b_pca:${OG_B2B_PCA_TAG}"]
  platforms = DUAL           # already pushed as a dual-arch manifest
}

target "csubst" {
  inherits  = ["_common"]
  context   = "docker/csubst"
  tags      = ["${REGISTRY}/csubst:${CSUBST_TAG}"]
  platforms = DUAL           # already pushed as a dual-arch manifest
}
