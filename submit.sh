#!/usr/bin/env bash
# Internal helper; users normally call the project's 01_submit.sh.
set -Eeuo pipefail
IFS=$'\n\t'
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly SCRIPT_DIR
if [[ $# -ne 1 ]]; then
  printf 'Usage: %s /absolute/path/job.json\n' "$0" >&2
  exit 2
fi
config="$(readlink -f -- "$1")"
if [[ "$config" == *','* || "$config" == *' '* || "$config" == *$'\n'* ]]; then
  printf 'PBS configuration path must not contain spaces, commas or newlines.\n' >&2
  exit 2
fi
if [[ ! -f "$SCRIPT_DIR/cluster_runtime.sh" ]]; then
  printf 'Use 00_new_project.sh to create a fixed project first.\n' >&2
  exit 2
fi
# shellcheck disable=SC1091
source "$SCRIPT_DIR/cluster_runtime.sh"
readonly PYTHON="${CONFONMR_PYTHON:?Missing Python in project environment}"
"$PYTHON" "$SCRIPT_DIR/nmrflow.py" check "$config"
args_text="$("$PYTHON" "$SCRIPT_DIR/cluster_config.py" scheduler-args "$config")"
mapfile -t scheduler_args <<< "$args_text"
job_id="$(qsub "${scheduler_args[@]}" -v "JOB_CONFIG=$config" "$SCRIPT_DIR/run_nmr.pbs")"
printf '%s\t%s\t%s\n' "$(date -Is)" "$job_id" "$config" >> "$SCRIPT_DIR/submissions.tsv"
printf 'Submitted %s. Monitor with qstat -f %s\n' "$job_id" "$job_id"
