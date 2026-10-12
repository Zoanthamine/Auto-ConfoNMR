#!/usr/bin/env bash
# Create a fixed project snapshot; no calculation is submitted.
set -Eeuo pipefail
IFS=$'\n\t'
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
if [[ -n "${CONFONMR_PYTHON:-}" ]]; then
  PYTHON="$CONFONMR_PYTHON"
elif [[ -f "$SCRIPT_DIR/config/cluster.json" ]]; then
  # System Python only reads JSON; the configured environment runs the workflow.
  PYTHON="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["python"])' "$SCRIPT_DIR/config/cluster.json")"
else
  printf '请先按 README.md 的迁移配置部分 配置 config/cluster.json。\n' >&2
  exit 2
fi
readonly PYTHON
exec "$PYTHON" "$SCRIPT_DIR/project_console.py" prepare "$@"
