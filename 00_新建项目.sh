#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
exec bash "$SCRIPT_DIR/00_new_project.sh" "$@"
