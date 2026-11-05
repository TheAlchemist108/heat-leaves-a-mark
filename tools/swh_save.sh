#!/usr/bin/env bash
# Request an independent, dated snapshot of this repository from Software Heritage.
# Usage: tools/swh_save.sh
set -euo pipefail

REPO_URL="https://github.com/TheAlchemist108/heat-leaves-a-mark"

curl -sS -X POST "https://archive.softwareheritage.org/api/1/origin/save/git/url/${REPO_URL}/"
echo
echo "Check status at: https://archive.softwareheritage.org/browse/origin/?origin_url=${REPO_URL}"
