#!/usr/bin/env bash
#SBATCH --job-name=vep_bos_taurus
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-vep-%j.out
#SBATCH --error=slurm-vep-%j.err

set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  sbatch scripts/run_vep_bos_taurus_slurm.sh <input_vcf> <output_vcf> [additional vep arguments]

Examples:
  sbatch scripts/run_vep_bos_taurus_slurm.sh input/sample.vcf.gz output/sample.vep.vcf.gz
  sbatch --export=ALL,VEP_PROFILE=everything --cpus-per-task=8 --mem=48G scripts/run_vep_bos_taurus_slurm.sh input/sample.vcf.gz output/sample.vep.vcf.gz
  sbatch --export=ALL,VEP_ENABLE_SIFT=1,VEP_ENABLE_POLYPHEN=1 scripts/run_vep_bos_taurus_slurm.sh input/sample.vcf.gz output/sample.vep.vcf.gz

Behavior:
  - Runs the existing Apptainer wrapper from this repository
  - Uses SLURM_CPUS_PER_TASK as VEP_FORKS unless VEP_FORKS is already set
  - Forwards any extra arguments directly to VEP

Useful environment variables:
  APPTAINER_IMAGE       Path to the .sif image
  VEP_HOST_DATA_DIR     Host directory bound to /data
  VEP_HOST_SCRIPTS_DIR  Host scripts directory bound inside the container
  VEP_FORKS             Overrides the default fork count
  VEP_PROFILE           Runner profile: minimal, cattle, or everything
  VEP_ENABLE_SIFT       Set to 1 to request SIFT annotations
  VEP_ENABLE_POLYPHEN   Set to 1 to request PolyPhen annotations
  APPTAINER_EXTRA_OPTS  Extra options inserted into apptainer exec
EOF
}

pick_repo_root() {
  local candidate

  for candidate in \
    "${VEP_REPO_ROOT:-}" \
    "${SLURM_SUBMIT_DIR:-}" \
    "${SLURM_SUBMIT_DIR:-}/.." \
    "$(pwd)" \
    "$(pwd)/.." \
    "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  do
    if [[ -n "${candidate}" && -x "${candidate}/scripts/run_vep_bos_taurus_apptainer.sh" ]]; then
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

if [[ $# -lt 2 ]]; then
  usage
  exit 1
fi

repo_root="$(pick_repo_root || true)"

if [[ -z "${repo_root}" ]]; then
  echo "Could not locate the repository root." >&2
  echo "Set VEP_REPO_ROOT to the VEP_cattle directory and resubmit the job." >&2
  exit 1
fi

runner="${repo_root}/scripts/run_vep_bos_taurus_apptainer.sh"

if [[ ! -x "${runner}" ]]; then
  echo "Apptainer runner not found or not executable: ${runner}" >&2
  exit 1
fi

mkdir -p "${repo_root}/data/logs"
cd "${repo_root}"

export VEP_FORKS="${VEP_FORKS:-${SLURM_CPUS_PER_TASK:-4}}"

echo "Job ID: ${SLURM_JOB_ID:-no_slurm_id}"
echo "Job name: ${SLURM_JOB_NAME:-vep_bos_taurus}"
echo "Node list: ${SLURM_JOB_NODELIST:-unknown}"
echo "Submit dir: ${SLURM_SUBMIT_DIR:-unknown}"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK:-unknown}"
echo "VEP forks: ${VEP_FORKS}"
echo "Input: $1"
echo "Output: $2"
echo "Repo root: ${repo_root}"

if command -v srun >/dev/null 2>&1; then
  srun "${runner}" "$@"
else
  "${runner}" "$@"
fi
