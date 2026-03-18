# VEP_cattle

This repository builds an Apptainer image for running Ensembl VEP locally on `Bos_taurus` VCF files and provides both container-side and host-side helper scripts.

- `install_bos_taurus_cache.sh`: runs inside the container and downloads the `Bos_taurus` cache and FASTA into `/data`.
- `run_vep_bos_taurus.sh`: runs inside the container and annotates a VCF using the local cache in offline mode.
- `install_bos_taurus_cache_apptainer.sh`: host-side wrapper that launches the install script with `apptainer exec`.
- `run_vep_bos_taurus_apptainer.sh`: host-side wrapper that launches VEP annotation with `apptainer exec`.
- `run_vep_bos_taurus_slurm.sh`: Slurm batch script that runs the Apptainer wrapper with `sbatch`.
- `export_vep_table.py`: converts the annotated VCF `CSQ` field into an analysis-ready TSV with all or selected VEP fields.
- `annotate_bim_from_vep.py`: joins VEP annotations back onto a BIM file while retaining all `CSQ` columns in the output.

The setup is pinned to Ensembl VEP `release_115.2` and defaults to the `Bos_taurus` assembly `ARS-UCD2.0`.
The Apptainer wrappers bind the repository `scripts/` directory into the container, so helper-script fixes do not require rebuilding the `.sif` image.

## Repository layout

```text
.
├── Apptainer.def
├── README.md
└── scripts
    ├── annotate_bim_from_vep.py
    ├── export_vep_table.py
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

## 4. Run VEP on your VCF with the full annotation set

Recommended example with a gzipped input VCF and the broadest built-in annotation set:

```bash
VEP_PROFILE=everything \
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

Equivalent example with an uncompressed VCF:

```bash
VEP_PROFILE=everything \
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf \
  output/my_samples.vep.vcf
```

The script automatically:

- Uses the local cache in offline mode
- Detects the installed `Bos_taurus` FASTA when present
- Runs VEP with `--everything` when `VEP_PROFILE=everything`
- Keeps the annotations inside the VCF `CSQ` field so all VEP subfields remain available for export later
- Writes a warnings file and VEP stats HTML next to the output VCF

If you want the full run plus optional score-style annotations such as `SIFT` and `PolyPhen`, use:

```bash
VEP_PROFILE=everything \
VEP_ENABLE_SIFT=1 \
VEP_ENABLE_POLYPHEN=1 \
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

## 5. Annotation profiles

The container-side runner now supports three built-in profiles:

- `minimal`: consequence, impact, gene, transcript, and core variant metadata
- `cattle`: the default profile, adding richer transcript/protein metadata for Bos taurus work
- `everything`: forwards `--everything` to VEP for the broadest annotation set

If your goal is to keep all VEP annotations available for downstream analysis, use `VEP_PROFILE=everything`.
The lighter profiles are still available when you want smaller outputs or faster runs.

List the available profiles:

```bash
./scripts/run_vep_bos_taurus_apptainer.sh --list-profiles
```

Use the default `cattle` profile:

```bash
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

Use the broadest built-in profile:

```bash
VEP_PROFILE=everything \
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

Request score-style annotations such as `SIFT` and `PolyPhen`:

```bash
VEP_ENABLE_SIFT=1 \
VEP_ENABLE_POLYPHEN=1 \
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

Anything after the output filename is still passed directly to `vep`, so you can combine profiles with explicit VEP flags when needed.

## 6. Export all VEP annotations to a table

The annotated VCF keeps VEP data inside the `CSQ` field. To analyse all annotations in a TSV, use the exporter script. The most complete path is:

1. Run annotation with `VEP_PROFILE=everything`
2. Export with `--mode all`
3. Keep `--fields all` (the default) so every `CSQ` subfield is written to the table

Show the fields present in one annotated VCF:

```bash
python3 scripts/export_vep_table.py \
  --vcf data/output/my_samples.vep.vcf.gz \
  --show-fields
```

Export every matched `CSQ` row for every ALT allele, keeping all VEP subfields:

```bash
python3 scripts/export_vep_table.py \
  --vcf data/output/my_samples.vep.vcf.gz \
  --output data/output/my_samples.vep.all_csq.tsv.gz \
  --mode all
```

This `all_csq` table is the one to use when you do not want to lose any VEP information.

Export one best-ranked annotation row per allele with selected fields:

```bash
python3 scripts/export_vep_table.py \
  --vcf data/output/my_samples.vep.vcf.gz \
  --output data/output/my_samples.vep.best.tsv.gz \
  --mode best \
  --fields Consequence,IMPACT,SYMBOL,Gene,Feature,BIOTYPE,CANONICAL,HGVSc,HGVSp,SIFT,PolyPhen,Protein_position,Amino_acids,Codons,DOMAINS,PICK
```

This exporter is not limited to a hard-coded field set. If the annotated VCF contains extra `CSQ` fields, they are exported automatically when `--fields all` is used.

## 7. Expected outputs

After a run, you should see files like:

```text
data/output/my_samples.vep.vcf.gz
data/output/my_samples.vep.vcf.gz.stats.html
data/output/my_samples.vep.vcf.gz.warnings.txt
```

If you also export the annotations as a table, you will get files such as:

```text
data/output/my_samples.vep.all_csq.tsv.gz
data/output/my_samples.vep.best.tsv.gz
```

## 8. Run on Slurm with sbatch

You can submit the annotation as a Slurm job with the included batch script:

```bash
sbatch scripts/run_vep_bos_taurus_slurm.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

To request more resources, pass standard `sbatch` options:

```bash
sbatch \
  --export=ALL,VEP_PROFILE=everything \
  --cpus-per-task=8 \
  --mem=48G \
  --time=48:00:00 \
  scripts/run_vep_bos_taurus_slurm.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

The Slurm script automatically sets `VEP_FORKS` from `SLURM_CPUS_PER_TASK`, so requesting more CPUs will increase the number of VEP forks unless you explicitly export `VEP_FORKS` yourself.

If your cluster runs batch scripts from a temporary Slurm spool directory, the script will try to recover the repository root from `SLURM_SUBMIT_DIR`. You can also force it explicitly:

```bash
sbatch --export=ALL,VEP_REPO_ROOT="$(pwd)" \
  scripts/run_vep_bos_taurus_slurm.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

To request score annotations in Slurm:

```bash
sbatch --export=ALL,VEP_ENABLE_SIFT=1,VEP_ENABLE_POLYPHEN=1 \
  scripts/run_vep_bos_taurus_slurm.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

## 9. Useful environment variables

The host-side wrapper scripts understand:

- `APPTAINER_IMAGE`: path to the `.sif` file
- `VEP_HOST_DATA_DIR`: host directory to bind to `/data`
- `VEP_HOST_SCRIPTS_DIR`: host scripts directory to bind inside the container
- `VEP_REPO_ROOT`: repository root used by the Slurm wrapper if auto-detection is not enough
- `VEP_SPECIES`: defaults to `bos_taurus`
- `VEP_ASSEMBLY`: defaults to `ARS-UCD2.0`
- `VEP_CACHE_VERSION`: defaults to `115`
- `VEP_FORKS`: defaults to `4` for annotation
- `VEP_PROFILE`: defaults to `cattle`
- `VEP_ENABLE_SIFT`: set to `1` to request SIFT annotations
- `VEP_ENABLE_POLYPHEN`: set to `1` to request PolyPhen annotations
- `APPTAINER_EXTRA_OPTS`: extra flags inserted into `apptainer exec`

Example using a scratch directory and a custom image path:

```bash
APPTAINER_IMAGE=/scratch/$USER/vep-cattle_115.2.sif \
VEP_HOST_DATA_DIR=/scratch/$USER/vep_cattle_data \
VEP_ENABLE_SIFT=1 \
VEP_ENABLE_POLYPHEN=1 \
./scripts/run_vep_bos_taurus_apptainer.sh \
  input/my_samples.vcf.gz \
  output/my_samples.vep.vcf.gz
```

## Notes

- Make sure your VCF assembly matches the installed `Bos_taurus` assembly. The current default in this repository is `ARS-UCD2.0`.
- If you want a different assembly or Ensembl release, update the defaults in `Apptainer.def` and rerun the cache installation.
- The first cache download can take a while and requires internet access during the container run.
- Availability of optional prediction fields such as `SIFT` and `PolyPhen` depends on what Ensembl VEP provides for the chosen species, cache, and release. The scripts now request them when enabled, and they will appear in the exported table if VEP emits them.
- Apptainer is the continuation of Singularity. If your cluster exposes `singularity` instead of `apptainer`, the same commands usually work after replacing the executable name.
- If Apptainer prints `underlay ... required more than 50 bind mounts`, this is usually a host-side mount warning rather than the real failure. If your cluster auto-binds many filesystems, try `APPTAINER_EXTRA_OPTS="--no-mount hostfs"` with the wrapper script.
