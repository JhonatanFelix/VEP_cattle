#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_vep_bos_taurus.sh <input_vcf> <output_vcf> [additional vep arguments]

Examples:
  run_vep_bos_taurus.sh input/sample.vcf.gz output/sample.vep.vcf.gz
  run_vep_bos_taurus.sh input/sample.vcf output/sample.vep.vcf --everything --fork 8

Relative paths are resolved inside VEP_DATA_DIR (default: /data).
EOF
}

resolve_path() {
  local path="$1"

  if [[ "${path}" = /* ]]; then
    printf '%s\n' "${path}"
  else
    printf '%s/%s\n' "${VEP_DATA_DIR}" "${path}"
  fi
}

find_fasta() {
  local candidate

  while IFS= read -r candidate; do
    printf '%s\n' "${candidate}"
    return 0
  done < <(
    # In the Ensembl container /data is a symlink to /opt/vep/.vep, so we
    # must follow symlinks to discover the installed FASTA files.
    find -L "${VEP_DATA_DIR}" -type f \
      \( -iname "*.fa" -o -iname "*.fa.gz" -o -iname "*.fasta" -o -iname "*.fasta.gz" \) \
      \( -iname "*${VEP_SPECIES}*" -o -iname "*${VEP_ASSEMBLY}*" \) \
      | sort
  )

  return 1
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

if [[ $# -lt 2 ]]; then
  usage
  exit 1
fi

VEP_DATA_DIR="${VEP_DATA_DIR:-/data}"
VEP_SPECIES="${VEP_SPECIES:-bos_taurus}"
VEP_ASSEMBLY="${VEP_ASSEMBLY:-ARS-UCD2.0}"
VEP_CACHE_VERSION="${VEP_CACHE_VERSION:-115}"
VEP_FORKS="${VEP_FORKS:-4}"

input_file="$(resolve_path "$1")"
output_file="$(resolve_path "$2")"
shift 2

if [[ ! -f "${input_file}" ]]; then
  echo "Input VCF not found: ${input_file}" >&2
  exit 1
fi

if [[ ! -d "${VEP_DATA_DIR}/${VEP_SPECIES}" && ! -d "${VEP_DATA_DIR}/${VEP_SPECIES}_refseq" ]]; then
  echo "No VEP cache detected under ${VEP_DATA_DIR}. Run install_bos_taurus_cache.sh first." >&2
  exit 1
fi

mkdir -p "$(dirname "${output_file}")"

stats_file="${output_file}.stats.html"
warning_file="${output_file}.warnings.txt"

vep_args=(
  --input_file "${input_file}"
  --output_file "${output_file}"
  --species "${VEP_SPECIES}"
  --assembly "${VEP_ASSEMBLY}"
  --cache
  --offline
  --dir_cache "${VEP_DATA_DIR}"
  --format vcf
  --vcf
  --force_overwrite
  --fork "${VEP_FORKS}"
  --stats_file "${stats_file}"
  --warning_file "${warning_file}"
  --symbol
  --biotype
  --canonical
  --variant_class
  --numbers
  --total_length
  --allele_number
)

if [[ "${output_file}" == *.gz ]]; then
  vep_args+=(--compress_output bgzip)
fi

if fasta_file="$(find_fasta)"; then
  echo "Using FASTA file: ${fasta_file}"
  vep_args+=(--fasta "${fasta_file}" --hgvs)
else
  echo "No Bos taurus FASTA found under ${VEP_DATA_DIR}; continuing without --fasta/--hgvs" >&2
fi

echo "Annotating ${input_file}"
vep "${vep_args[@]}" "$@"
