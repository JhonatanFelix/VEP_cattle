#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: install_bos_taurus_cache.sh [additional INSTALL.pl arguments]

Downloads the Ensembl VEP cache and FASTA for Bos taurus into /data by default.
By default this uses the VEP installer's standard download protocol selection.

Environment variables:
  VEP_DATA_DIR        Target directory for VEP cache/FASTA (default: /data)
  VEP_SPECIES         Species name for VEP (default: bos_taurus)
  VEP_ASSEMBLY        Assembly name (default: ARS-UCD2.0)
  VEP_CACHE_VERSION   Cache version matching the VEP release (default: 115)
  VEP_INSTALLER       Path to INSTALL.pl inside the container
  VEP_USE_HTTPS_PROTO Set to 1 to add --USE_HTTPS_PROTO explicitly
EOF
}

find_installer() {
  local candidate

  for candidate in \
    "${VEP_INSTALLER:-}" \
    "$(command -v INSTALL.pl 2>/dev/null || true)" \
    /opt/vep/src/ensembl-vep/INSTALL.pl \
    /opt/vep/INSTALL.pl \
    /usr/local/bin/INSTALL.pl
  do
    if [[ -n "${candidate}" && -f "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done

  return 1
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

VEP_DATA_DIR="${VEP_DATA_DIR:-/data}"
VEP_SPECIES="${VEP_SPECIES:-bos_taurus}"
VEP_ASSEMBLY="${VEP_ASSEMBLY:-ARS-UCD2.0}"
VEP_CACHE_VERSION="${VEP_CACHE_VERSION:-115}"
VEP_INSTALLER="$(find_installer || true)"
VEP_USE_HTTPS_PROTO="${VEP_USE_HTTPS_PROTO:-0}"

mkdir -p "${VEP_DATA_DIR}/input" "${VEP_DATA_DIR}/output" "${VEP_DATA_DIR}/logs"

if [[ -z "${VEP_INSTALLER}" ]]; then
  echo "VEP installer not found in the container image." >&2
  exit 1
fi

echo "Installing VEP cache and FASTA for ${VEP_SPECIES} (${VEP_ASSEMBLY}) into ${VEP_DATA_DIR}"

install_args=(
  -a cf \
  -s "${VEP_SPECIES}" \
  -y "${VEP_ASSEMBLY}" \
  -c "${VEP_DATA_DIR}" \
  --CACHE_VERSION "${VEP_CACHE_VERSION}" \
)

if [[ "${VEP_USE_HTTPS_PROTO}" == "1" ]]; then
  install_args+=(--USE_HTTPS_PROTO)
fi

perl "${VEP_INSTALLER}" \
  "${install_args[@]}" \
  "$@"
