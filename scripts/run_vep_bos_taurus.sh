#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: run_vep_bos_taurus.sh [--profile PROFILE] <input_vcf> <output_vcf> [additional vep arguments]

Examples:
  run_vep_bos_taurus.sh input/sample.vcf.gz output/sample.vep.vcf.gz
  run_vep_bos_taurus.sh --profile cattle input/sample.vcf.gz output/sample.vep.vcf.gz
  VEP_ENABLE_SIFT=1 VEP_ENABLE_POLYPHEN=1 run_vep_bos_taurus.sh input/sample.vcf.gz output/sample.vep.vcf.gz
  run_vep_bos_taurus.sh --profile everything input/sample.vcf output/sample.vep.vcf --fork 8

Relative paths are resolved inside VEP_DATA_DIR (default: /data).

Profiles:
  minimal     Core consequence/gene/transcript annotations
  cattle      Adds richer transcript/protein metadata for Bos taurus analysis
  everything  Delegates to VEP --everything for the broadest annotation set

Environment variables:
  VEP_PROFILE          Default profile when --profile is not given (default: cattle)
  VEP_ENABLE_SIFT      Set to 1 to add --sift b
  VEP_ENABLE_POLYPHEN  Set to 1 to add --polyphen b
EOF
}

list_profiles() {
  cat <<'EOF'
minimal
cattle
everything
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

profile_args=()

append_profile_args() {
  case "$1" in
    minimal)
      ;;
    cattle)
      profile_args+=(
        --protein
        --ccds
        --uniprot
        --domains
        --xref_refseq
        --transcript_version
        --gene_version
        --flag_pick_allele
      )
      ;;
    everything)
      profile_args+=(--everything)
      ;;
    *)
      echo "Unknown VEP profile: $1" >&2
      echo "Available profiles:" >&2
      list_profiles >&2
      exit 1
      ;;
  esac
}

if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
  exit 0
fi

if [[ "${1:-}" == "--list-profiles" ]]; then
  list_profiles
  exit 0
fi

profile="${VEP_PROFILE:-cattle}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --help|-h)
      usage
      exit 0
      ;;
    --list-profiles)
      list_profiles
      exit 0
      ;;
    --profile)
      if [[ $# -lt 2 ]]; then
        echo "Missing value after --profile" >&2
        exit 1
      fi
      profile="$2"
      shift 2
      ;;
    --)
      shift
      break
      ;;
    -*)
      echo "Unknown wrapper option before input/output arguments: $1" >&2
      echo "Pass VEP options after <output_vcf>, or use --profile for wrapper configuration." >&2
      exit 1
      ;;
    *)
      break
      ;;
  esac
done

if [[ $# -lt 2 ]]; then
  usage
  exit 1
fi

VEP_DATA_DIR="${VEP_DATA_DIR:-/data}"
VEP_SPECIES="${VEP_SPECIES:-bos_taurus}"
VEP_ASSEMBLY="${VEP_ASSEMBLY:-ARS-UCD2.0}"
VEP_CACHE_VERSION="${VEP_CACHE_VERSION:-115}"
VEP_FORKS="${VEP_FORKS:-4}"
VEP_ENABLE_SIFT="${VEP_ENABLE_SIFT:-0}"
VEP_ENABLE_POLYPHEN="${VEP_ENABLE_POLYPHEN:-0}"

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

append_profile_args "${profile}"

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
  "${profile_args[@]}"
)

if [[ "${output_file}" == *.gz ]]; then
  vep_args+=(--compress_output bgzip)
fi

if [[ "${VEP_ENABLE_SIFT}" == "1" ]]; then
  vep_args+=(--sift b)
fi

if [[ "${VEP_ENABLE_POLYPHEN}" == "1" ]]; then
  vep_args+=(--polyphen b)
fi

if fasta_file="$(find_fasta)"; then
  echo "Using FASTA file: ${fasta_file}"
  vep_args+=(--fasta "${fasta_file}" --hgvs)
else
  echo "No Bos taurus FASTA found under ${VEP_DATA_DIR}; continuing without --fasta/--hgvs" >&2
fi

echo "Annotating ${input_file}"
echo "VEP profile: ${profile}"
echo "SIFT enabled: ${VEP_ENABLE_SIFT}"
echo "PolyPhen enabled: ${VEP_ENABLE_POLYPHEN}"
vep "${vep_args[@]}" "$@"
