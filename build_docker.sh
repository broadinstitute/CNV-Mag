#!/usr/bin/env bash
set -euo pipefail

REGISTRY="us.gcr.io/tag-public"
IMAGE="cnv-mag"
TAG="${1:-}"

if [[ -z "$TAG" ]]; then
    echo "Usage: ./build_docker.sh <tag>  (e.g. ./build_docker.sh v0.5)"
    exit 1
fi

FULL_IMAGE="${REGISTRY}/${IMAGE}:${TAG}"

echo "Building ${FULL_IMAGE} ..."
docker build -t "${FULL_IMAGE}" .

echo "Pushing ${FULL_IMAGE} ..."
docker push "${FULL_IMAGE}"

echo ""
echo "Done. Update the WDL default:"
echo "  String dockerImage = \"${FULL_IMAGE}\""
