#!/bin/bash
#SBATCH --job-name=vep_qc
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=../vep_qc_%j.out
#SBATCH --error=../vep_qc_%j.err

set -euo pipefail

module load releases/2022b
module load BCFtools/1.17-GCC-12.2.0

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$(pwd -P)}"
cd "$SUBMIT_DIR"

ANNOT_VCF="../data/output/samples.vep.vcf.gz"
ORIG_VCF="../data/input/samples.vcf.gz"
OUTDIR="../data/output/vep_qc"
LOGDIR="logs"

mkdir -p "$OUTDIR" "$LOGDIR"

log() {
    echo "[$(date '+%F %T')] $*"
}

require_file() {
    local f="$1"
    if [[ ! -f "$f" ]]; then
        echo "ERROR: file not found: $f" >&2
        exit 1
    fi
}

log "Starting VEP QC"
log "Submit directory: $SUBMIT_DIR"
log "Annotated VCF: $ANNOT_VCF"
log "Original VCF:  $ORIG_VCF"

require_file "$ANNOT_VCF"
require_file "$ORIG_VCF"

bcftools --version > "$OUTDIR/bcftools_version.txt"

# Index if needed
for f in "$ANNOT_VCF" "$ORIG_VCF"; do
    if [[ ! -f "${f}.csi" && ! -f "${f}.tbi" ]]; then
        log "Index not found for $f ; creating one"
        bcftools index -f "$f"
    fi
done

log "Saving headers"
bcftools view -h "$ORIG_VCF"  > "$OUTDIR/original_header.txt"
bcftools view -h "$ANNOT_VCF" > "$OUTDIR/annotated_header.txt"

grep '^##INFO=<ID=CSQ' "$OUTDIR/annotated_header.txt" > "$OUTDIR/csq_header.txt" || true
if [[ ! -s "$OUTDIR/csq_header.txt" ]]; then
    echo "ERROR: no CSQ header found in annotated VCF. This does not look like a standard VEP-annotated VCF." >&2
    exit 1
fi

log "Listing CSQ subfields"
bcftools +split-vep "$ANNOT_VCF" -l > "$OUTDIR/csq_fields.txt"

log "Computing basic counts"
ORIG_N=$(bcftools view -H "$ORIG_VCF"  | wc -l | awk '{print $1}')
ANNOT_N=$(bcftools view -H "$ANNOT_VCF" | wc -l | awk '{print $1}')
SAMPLES_N=$(bcftools query -l "$ORIG_VCF" | wc -l | awk '{print $1}')
SITES_WITH_CSQ=$(bcftools view -i 'INFO/CSQ!="."' -H "$ANNOT_VCF" | wc -l | awk '{print $1}')

{
    echo -e "metric\tvalue"
    echo -e "original_variants\t$ORIG_N"
    echo -e "annotated_variants\t$ANNOT_N"
    echo -e "annotated_sites_with_CSQ\t$SITES_WITH_CSQ"
    echo -e "samples_in_original_vcf\t$SAMPLES_N"
} > "$OUTDIR/variant_counts.tsv"


if ! bcftools view -h "$ANNOT_VCF" | grep -q '^##INFO=<ID=CSQ'; then
    echo "ERROR: no CSQ header found in annotated VCF. This does not look like a standard VEP-annotated VCF." >&2
    exit 1
fi

log "Counting total SNP records in the annotated VCF"
TOTAL_SNPS=$(bcftools view -H -v snps "$ANNOT_VCF" | wc -l | awk '{print $1}')
log "Shape: $TOTAL_SNPS SNP records in total"

{
    echo -e "metric\tvalue"
    echo -e "annotated_snp_records\t$TOTAL_SNPS"
} >> "$OUTDIR/variant_counts.tsv"

log "Saving a small preview of parsed VEP annotations"
{
    echo -e "CHROM\tPOS\tREF\tALT\tConsequence\tIMPACT\tSYMBOL\tGene\tFeature\tBIOTYPE"
    bcftools +split-vep "$ANNOT_VCF" \
        -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature\t%BIOTYPE\n' \
        -s worst | awk 'NR<=20'
} > "$OUTDIR/preview_worst_transcript.tsv"

log "Counting worst-transcript consequence terms"
bcftools +split-vep "$ANNOT_VCF" -f '%Consequence\n' -s worst \
| awk '
    function add_term(term){ if(term != "" && term != ".") c[term]++ }
    {
        if($0 == "" || $0 == ".") next
        n = split($0, a, /&/)
        for(i=1; i<=n; i++) add_term(a[i])
    }
    END {
        for(k in c) print k "\t" c[k]
    }
' | sort -k2,2nr > "$OUTDIR/consequence_counts_worst.tsv"

log "Counting worst-transcript IMPACT classes"
bcftools +split-vep "$ANNOT_VCF" -f '%IMPACT\n' -s worst \
| awk '
    $0 != "" && $0 != "." { c[$0]++ }
    END {
        for(k in c) print k "\t" c[k]
    }
' | sort -k2,2nr > "$OUTDIR/impact_counts_worst.tsv"

log "Computing per-site completeness using the worst transcript"
bcftools +split-vep "$ANNOT_VCF" \
    -f '%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature\t%BIOTYPE\n' \
    -s worst \
| awk '
    function present(x){ return (x != "" && x != ".") }
    BEGIN {
        FS=OFS="\t"
        cons=impact=symbol=gene=feature=biotype=all6=any_missing=0
        total=0
    }
    {
        total++
        p1=present($1); p2=present($2); p3=present($3); p4=present($4); p5=present($5); p6=present($6)
        if(p1) cons++
        if(p2) impact++
        if(p3) symbol++
        if(p4) gene++
        if(p5) feature++
        if(p6) biotype++
        if(p1 && p2 && p3 && p4 && p5 && p6) all6++
        if(!(p1 && p2 && p3 && p4 && p5 && p6)) any_missing++
    }
    END {
        print "metric","non_missing","missing","pct_non_missing"
        print "Consequence", cons, total-cons, 100*cons/total
        print "IMPACT", impact, total-impact, 100*impact/total
        print "SYMBOL", symbol, total-symbol, 100*symbol/total
        print "Gene", gene, total-gene, 100*gene/total
        print "Feature", feature, total-feature, 100*feature/total
        print "BIOTYPE", biotype, total-biotype, 100*biotype/total
        print "ALL_6_FIELDS_PRESENT", all6, total-all6, 100*all6/total
        print "ANY_MISSING_FIELD", total-all6, all6, 100*(total-all6)/total
    }
' > "$OUTDIR/completeness_worst_transcript.tsv"

log "Computing completeness across all CSQ entries (transcript/consequence level)"
bcftools +split-vep "$ANNOT_VCF" \
    -f '%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature\t%BIOTYPE\n' \
    -d \
| awk '
    function present(x){ return (x != "" && x != ".") }
    BEGIN {
        FS=OFS="\t"
        cons=impact=symbol=gene=feature=biotype=all6=0
        total=0
    }
    {
        total++
        p1=present($1); p2=present($2); p3=present($3); p4=present($4); p5=present($5); p6=present($6)
        if(p1) cons++
        if(p2) impact++
        if(p3) symbol++
        if(p4) gene++
        if(p5) feature++
        if(p6) biotype++
        if(p1 && p2 && p3 && p4 && p5 && p6) all6++
    }
    END {
        print "metric","non_missing","missing","pct_non_missing"
        print "Consequence", cons, total-cons, 100*cons/total
        print "IMPACT", impact, total-impact, 100*impact/total
        print "SYMBOL", symbol, total-symbol, 100*symbol/total
        print "Gene", gene, total-gene, 100*gene/total
        print "Feature", feature, total-feature, 100*feature/total
        print "BIOTYPE", biotype, total-biotype, 100*biotype/total
        print "ALL_6_FIELDS_PRESENT", all6, total-all6, 100*all6/total
    }
' > "$OUTDIR/completeness_all_csq_entries.tsv"

log "Saving up to 200 sites whose worst transcript is missing at least one key field"
bcftools +split-vep "$ANNOT_VCF" \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%Consequence\t%IMPACT\t%SYMBOL\t%Gene\t%Feature\t%BIOTYPE\n' \
    -s worst \
| awk '
    function present(x){ return (x != "" && x != ".") }
    BEGIN {
        FS=OFS="\t"
        print "CHROM","POS","REF","ALT","Consequence","IMPACT","SYMBOL","Gene","Feature","BIOTYPE","missing_fields"
    }
    {
        miss=""
        if(!present($5)) miss = miss "Consequence,"
        if(!present($6)) miss = miss "IMPACT,"
        if(!present($7)) miss = miss "SYMBOL,"
        if(!present($8)) miss = miss "Gene,"
        if(!present($9)) miss = miss "Feature,"
        if(!present($10)) miss = miss "BIOTYPE,"
        if(miss != "") {
            sub(/,$/, "", miss)
            print $0, miss
            n++
            if(n >= 200) exit
        }
    }
' > "$OUTDIR/sites_missing_key_fields.tsv"

log "Building a concise summary report"
{
    echo "VEP QC SUMMARY"
    echo "=============="
    echo
    echo "Annotated VCF: $ANNOT_VCF"
    echo
    echo "Shape"
    cat "$OUTDIR/variant_counts.tsv"
    echo
    echo "This run intentionally skips the pre-preview QC block."
    echo
    echo "Files generated in $OUTDIR"
    echo "- plugin_diagnostic.txt (only if plugin loading failed)"
    echo "- variant_counts.tsv"
    echo "- preview_worst_transcript.tsv"
    echo "- consequence_counts_worst.tsv"
    echo "- impact_counts_worst.tsv"
    echo "- completeness_worst_transcript.tsv"
    echo "- completeness_all_csq_entries.tsv"
    echo "- sites_missing_key_fields.tsv"
} > "$OUTDIR/summary_report.txt"

log "Done"
log "Main report: $OUTDIR/summary_report.txt"
