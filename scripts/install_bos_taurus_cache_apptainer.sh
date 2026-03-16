#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_bos_taurus_cache_apptainer.sh [additional INSTALL.pl arguments]

Runs the Bos taurus cache/FASTA installation inside an Apptainer image.

Environment variables:
  APPTAINER_IMAGE       Path to the .sif image (default: <repo>/vep-cattle_115.2.sif)
  VEP_HOST_DATA_DIR     Host directory mounted to /data (default: <repo>/data)
  VEP_HOST_SCRIPTS_DIR  Host scripts directory mounted inside the container
  VEP_SPECIES           Species name passed into the container (default: bos_taurus)
  VEP_ASSEMBLY          Assembly passed into the container (default: ARS-UCD2.0)
  VEP_CACHE_VERSION     Cache version passed into the container (default: 115)
  APPTAINER_EXTRA_OPTS  Extra options inserted before "exec" arguments
EOF
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

APPTAINER_IMAGE="${APPTAINER_IMAGE:-${repo_root}/vep-cattle_115.2.sif}"
VEP_HOST_DATA_DIR="${VEP_HOST_DATA_DIR:-${repo_root}/data}"
VEP_HOST_SCRIPTS_DIR="${VEP_HOST_SCRIPTS_DIR:-${repo_root}/scripts}"

if [[ ! -f "${APPTAINER_IMAGE}" ]]; then
  echo "Apptainer image not found: ${APPTAINER_IMAGE}" >&2
  exit 1
fi

if [[ ! -d "${VEP_HOST_SCRIPTS_DIR}" ]]; then
  echo "Host scripts directory not found: ${VEP_HOST_SCRIPTS_DIR}" >&2
  exit 1
fi

mkdir -p "${VEP_HOST_DATA_DIR}/input" "${VEP_HOST_DATA_DIR}/output" "${VEP_HOST_DATA_DIR}/logs"

export APPTAINERENV_VEP_DATA_DIR=/data
export APPTAINERENV_VEP_SPECIES="${VEP_SPECIES:-bos_taurus}"
export APPTAINERENV_VEP_ASSEMBLY="${VEP_ASSEMBLY:-ARS-UCD2.0}"
export APPTAINERENV_VEP_CACHE_VERSION="${VEP_CACHE_VERSION:-115}"

apptainer_extra_opts=()
if [[ -n "${APPTAINER_EXTRA_OPTS:-}" ]]; then
  # shellcheck disable=SC2206
  apptainer_extra_opts=(${APPTAINER_EXTRA_OPTS})
fi

apptainer exec \
  --cleanenv \
  "${apptainer_extra_opts[@]}" \
  --bind "${VEP_HOST_DATA_DIR}:/data" \
  --bind "${VEP_HOST_SCRIPTS_DIR}:/opt/vep_cattle_scripts:ro" \
  "${APPTAINER_IMAGE}" \
  /opt/vep_cattle_scripts/install_bos_taurus_cache.sh \
  "$@"
