# VEP_cattle

This repository builds an Apptainer image for running Ensembl VEP locally on `Bos_taurus` VCF files and provides both container-side and host-side helper scripts.

- `install_bos_taurus_cache.sh`: runs inside the container and downloads the `Bos_taurus` cache and FASTA into `/data`.
- `run_vep_bos_taurus.sh`: runs inside the container and annotates a VCF using the local cache in offline mode.
- `install_bos_taurus_cache_apptainer.sh`: host-side wrapper that launches the install script with `apptainer exec`.
- `run_vep_bos_taurus_apptainer.sh`: host-side wrapper that launches VEP annotation with `apptainer exec`.
- `run_vep_bos_taurus_slurm.sh`: Slurm batch script that runs the Apptainer wrapper with `sbatch`.

The setup is pinned to Ensembl VEP `release_115.2` and defaults to the `Bos_taurus` assembly `ARS-UCD2.0`.
The Apptainer wrappers bind the repository `scripts/` directory into the container, so helper-script fixes do not require rebuilding the `.sif` image.

## Repository layout

```text
.
├── Apptainer.def
├── README.md
└── scripts
    ├── install_bos_taurus_cache.sh
    ├── install_bos_taurus_cache_apptainer.sh
    ├── run_vep_bos_taurus_slurm.sh
    ├── run_vep_bos_taurus.sh
    └── run_vep_bos_taurus_apptainer.sh
```

## 1. Build the Apptainer image

From the root of this repository:

```bash
apptainer build vep-cattle_115.2.sif Apptainer.def
```

If your system requires rootless builds, use:

```bash
apptainer build --fakeroot vep-cattle_115.2.sif Apptainer.def
```

This produces a local `.sif` image that already contains the two container-side VEP helper scripts.

## 2. Prepare a persistent data directory

Create a local directory that will be mounted into the container at `/data`. This keeps the cache, FASTA, input VCFs, and results on your machine.

```bash
mkdir -p data/input data/output data/logs
```

Copy your `Bos_taurus` VCF file into `data/input/`.

## 3. Install the Bos taurus VEP cache and FASTA

The easiest way is to use the host-side wrapper script:

```bash
./scripts/install_bos_taurus_cache_apptainer.sh
```

This expects the image at `./vep-cattle_115.2.sif` and the mounted data directory at `./data`.
It also bind-mounts `./scripts` into the container so the latest repo-local helper script is used.

Equivalent direct Apptainer command:

```bash
apptainer exec --cleanenv \
  --bind "$(pwd)/data:/data" \
  --bind "$(pwd)/scripts:/opt/vep_cattle_scripts:ro" \
  vep-cattle_115.2.sif \
  /opt/vep_cattle_scripts/install_bos_taurus_cache.sh
```

What this does:

- Downloads the Ensembl VEP cache for `bos_taurus`
- Downloads the corresponding FASTA for `ARS-UCD2.0`
- Stores everything under your local `data/` directory

If you need to override the defaults, you can set environment variables:

```bash
VEP_ASSEMBLY=ARS-UCD2.0 \
VEP_CACHE_VERSION=115 \
./scripts/install_bos_taurus_cache_apptainer.sh
```

## 4. Run VEP on your VCF

Example with a gzipped input VCF:

```bash
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

Example with an uncompressed VCF:

```bash
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf \
  output/my_samples.vep.vcf
```

The script automatically:

- Uses the local cache in offline mode
- Detects the installed `Bos_taurus` FASTA when present
- Adds useful annotation fields such as gene symbol, biotype, canonical transcript, HGVS, and variant class
- Writes a warnings file and VEP stats HTML next to the output VCF

## 5. Pass additional VEP options

Anything after the output filename is passed directly to `vep`.

```bash
VEP_FORKS=8 \
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz \
  --everything
```

Equivalent direct Apptainer command:

```bash
APPTAINERENV_VEP_FORKS=8 \
apptainer exec --cleanenv \
  --bind "$(pwd)/data:/data" \
  --bind "$(pwd)/scripts:/opt/vep_cattle_scripts:ro" \
  vep-cattle_115.2.sif \
  /opt/vep_cattle_scripts/run_vep_bos_taurus.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz \
  --everything
```

## 6. Expected outputs

After a run, you should see files like:

```text
data/output/my_samples.vep.vcf.gz
data/output/my_samples.vep.vcf.gz.stats.html
data/output/my_samples.vep.vcf.gz.warnings.txt
```

## 7. Run on Slurm with sbatch

You can submit the annotation as a Slurm job with the included batch script:

```bash
sbatch scripts/run_vep_bos_taurus_slurm.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

To request more resources, pass standard `sbatch` options:

```bash
sbatch \
  --cpus-per-task=8 \
  --mem=48G \
  --time=48:00:00 \
  scripts/run_vep_bos_taurus_slurm.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz \
  --everything
```

The Slurm script automatically sets `VEP_FORKS` from `SLURM_CPUS_PER_TASK`, so requesting more CPUs will increase the number of VEP forks unless you explicitly export `VEP_FORKS` yourself.

If your cluster runs batch scripts from a temporary Slurm spool directory, the script will try to recover the repository root from `SLURM_SUBMIT_DIR`. You can also force it explicitly:

```bash
sbatch --export=ALL,VEP_REPO_ROOT="$(pwd)" \
  scripts/run_vep_bos_taurus_slurm.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

## 8. Useful environment variables

The host-side wrapper scripts understand:

- `APPTAINER_IMAGE`: path to the `.sif` file
- `VEP_HOST_DATA_DIR`: host directory to bind to `/data`
- `VEP_HOST_SCRIPTS_DIR`: host scripts directory to bind inside the container
- `VEP_REPO_ROOT`: repository root used by the Slurm wrapper if auto-detection is not enough
- `VEP_SPECIES`: defaults to `bos_taurus`
- `VEP_ASSEMBLY`: defaults to `ARS-UCD2.0`
- `VEP_CACHE_VERSION`: defaults to `115`
- `VEP_FORKS`: defaults to `4` for annotation
- `APPTAINER_EXTRA_OPTS`: extra flags inserted into `apptainer exec`

Example using a scratch directory and a custom image path:

```bash
APPTAINER_IMAGE=/scratch/$USER/vep-cattle_115.2.sif \
VEP_HOST_DATA_DIR=/scratch/$USER/vep_cattle_data \
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

## Notes

- Make sure your VCF assembly matches the installed `Bos_taurus` assembly. The current default in this repository is `ARS-UCD2.0`.
- If you want a different assembly or Ensembl release, update the defaults in `Apptainer.def` and rerun the cache installation.
- The first cache download can take a while and requires internet access during the container run.
- Apptainer is the continuation of Singularity. If your cluster exposes `singularity` instead of `apptainer`, the same commands usually work after replacing the executable name.
- If Apptainer prints `underlay ... required more than 50 bind mounts`, this is usually a host-side mount warning rather than the real failure. If your cluster auto-binds many filesystems, try `APPTAINER_EXTRA_OPTS="--no-mount hostfs"` with the wrapper script.
