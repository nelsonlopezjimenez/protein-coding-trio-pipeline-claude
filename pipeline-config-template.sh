# pipeline_config.conf
# Configuration file for Trio Protein-Coding Gene Pipeline

# Reference genome (local file or URL)
REF="/path/to/reference/hg38.fa"

# Gene annotation file (GFF3 or GTF format)
# Download from: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gff3.gz
GFF="/path/to/annotation/gencode.v44.annotation.gff3"

# Trio sample IDs
FATHER_ID="NA12877"
MOTHER_ID="NA12878" 
CHILD_ID="NA12882"

# 1000 Genomes VCF URL template (use {CHR} placeholder)
VCF_URL_TEMPLATE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{CHR}.recalibrated_variants.vcf.gz"

# Output directory
OUTPUT_DIR="trio_protein_coding_analysis"

# Processing parameters
THREADS=8
MIN_GENE_LENGTH=300

# Temporary directory
TEMP_DIR="/tmp/trio_analysis"

# Chromosomes to process (space-separated)
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"

#############################################################################
# Multi-Chromosome Processing Script
# run_multi_chromosome.sh
#############################################################################

#!/bin/bash

# Multi-Chromosome Trio Analysis Runner
# Processes multiple chromosomes sequentially or in parallel

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_SCRIPT="${SCRIPT_DIR}/trio_consensus_pipeline.sh"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.conf"

# Default values
PARALLEL_MODE="false"
MAX_JOBS=4
CHROMOSOMES_TO_PROCESS=""
DOWNLOAD_MODE="true"

usage() {
    cat << EOF
Multi-Chromosome Trio Analysis Runner

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -c, --chromosomes LIST    Chromosomes to process (e.g., "chr1 chr2 chr22")
    -p, --parallel           Enable parallel processing
    -j, --jobs NUM           Maximum parallel jobs (default: 4)
    -d, --download           Enable download mode
    --config FILE            Configuration file (default: pipeline_config.conf)
    -h, --help               Show this help

EXAMPLES:
    # Process single chromosome
    $0 -c "chr22"
    
    # Process multiple chromosomes sequentially
    $0 -c "chr21 chr22"
    
    # Process chromosomes in parallel
    $0 -c "chr1 chr2 chr3" --parallel -j 3
    
    # Process all autosomes
    $0 -c "$(seq -f 'chr%.0f' 1 22 | tr '\n' ' ')"

EOF
}

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Load configuration
load_config() {
    if [[ -f "$CONFIG_FILE" ]]; then
        source "$CONFIG_FILE"
        log "Loaded configuration from $CONFIG_FILE"
    else
        log "ERROR: Configuration file not found: $CONFIG_FILE"
        exit 1
    fi
}

# Process single chromosome
process_chromosome() {
    local chr="$1"
    local vcf_url="${VCF_URL_TEMPLATE/\{CHR\}/$chr}"
    
    log "Processing chromosome $chr"
    
    local chr_output="${OUTPUT_DIR}/${chr}"
    mkdir -p "$chr_output"
    
    # Run pipeline for this chromosome
    "$PIPELINE_SCRIPT" \
        --chromosome "$chr" \
        --vcf "$vcf_url" \
        --reference "$REF" \
        --gff "$GFF" \
        --output "$chr_output" \
        --father "$FATHER_ID" \
        --mother "$MOTHER_ID" \
        --child "$CHILD_ID" \
        --download \
        --threads "$THREADS" \
        --min-length "$MIN_GENE_LENGTH" \
        --temp-dir "${TEMP_DIR}_${chr}"
    
    if [[ $? -eq 0 ]]; then
        log "Successfully completed chromosome $chr"
        return 0
    else
        log "ERROR: Failed to process chromosome $chr"
        return 1
    fi
}

# Process chromosomes in parallel
process_parallel() {
    local chromosomes=($CHROMOSOMES_TO_PROCESS)
    local pids=()
    local results=()
    
    log "Starting parallel processing of ${#chromosomes[@]} chromosomes with max $MAX_JOBS jobs"
    
    for chr in "${chromosomes[@]}"; do
        # Wait if we've reached max jobs
        while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
            for i in "${!pids[@]}"; do
                if ! kill -0 "${pids[$i]}" 2>/dev/null; then
                    wait "${pids[$i]}"
                    local exit_code=$?
                    results+=("${chromosomes[$i]}:$exit_code")
                    unset pids[$i]
                    break
                fi
            done
            # Compact arrays
            pids=("${pids[@]}")
            sleep 1
        done
        
        # Start new job
        log "Starting chromosome $chr"
        process_chromosome "$chr" &
        pids+=($!)
    done
    
    # Wait for remaining jobs
    for pid in "${pids[@]}"; do
        wait "$pid"
        local exit_code=$?
        # Find chromosome for this PID (simplified approach)
        results+=("unknown:$exit_code")
    done
    
    # Report results
    log "Parallel processing completed"
    local failed=0
    for result in "${results[@]}"; do
        local chr="${result%:*}"
        local code="${result#*:}"
        if [[ $code -ne 0 ]]; then
            log "FAILED: $chr (exit code: $code)"
            ((failed++))
        else
            log "SUCCESS: $chr"
        fi
    done
    
    log "Summary: $((${#results[@]} - failed)) succeeded, $failed failed"
    return $failed
}

# Process chromosomes sequentially  
process_sequential() {
    local chromosomes=($CHROMOSOMES_TO_PROCESS)
    local failed=0
    
    log "Starting sequential processing of ${#chromosomes[@]} chromosomes"
    
    for chr in "${chromosomes[@]}"; do
        if process_chromosome "$chr"; then
            log "SUCCESS: $chr"
        else
            log "FAILED: $chr"
            ((failed++))
        fi
    done
    
    log "Sequential processing completed"
    log "Summary: $((${#chromosomes[@]} - failed)) succeeded, $failed failed"
    return $failed
}

# Combine results from all chromosomes
combine_results() {
    log "Combining results from all chromosomes"
    
    local combined_hashes="${OUTPUT_DIR}/combined_trio_protein_coding_hashes.csv"
    local combined_summary="${OUTPUT_DIR}/combined_analysis_summary.txt"
    
    # Combine hash files
    local first_file=true
    for chr_dir in "${OUTPUT_DIR}"/chr*/; do
        if [[ -d "$chr_dir" ]]; then
            local chr=$(basename "$chr_dir")
            local hash_file="${chr_dir}/trio_${chr}_protein_coding_hashes.csv"
            
            if [[ -f "$hash_file" ]]; then
                if [[ "$first_file" == "true" ]]; then
                    # Include header from first file
                    cat "$hash_file" > "$combined_hashes"
                    first_file=false
                else
                    # Skip header for subsequent files
                    tail -n +2 "$hash_file" >> "$combined_hashes"
                fi
            fi
        fi
    done
    
    # Generate combined summary
    if [[ -f "$combined_hashes" ]]; then
        python3 << EOF > "$combined_summary"
import pandas as pd
from collections import Counter

try:
    df = pd.read_csv('$combined_hashes')
    
    print("Combined Trio Protein-Coding Gene Analysis Summary")
    print("=" * 60)
    print(f"Analysis date: $(date)")
    print(f"Trio: $FATHER_ID (father), $MOTHER_ID (mother), $CHILD_ID (child)")
    print()
    
    if len(df) == 0:
        print("No data found")
        exit()
    
    # Overall statistics
    print("Overall Statistics:")
    print(f"  Total gene sequences analyzed: {len(df):,}")
    print(f"  Unique genes: {df['gene_id'].nunique():,}")
    print(f"  Chromosomes processed: {df['chromosome'].nunique()}")
    print(f"  Unique hash values: {df['sha256_hash'].nunique():,}")
    print()
    
    # Chromosome breakdown
    chr_stats = df.groupby('chromosome').agg({
        'gene_id': 'nunique',
        'sha256_hash': 'nunique',
        'clean_length': 'mean'
    }).round(1)
    chr_stats.columns = ['Unique_Genes', 'Unique_Hashes', 'Avg_Length']
    print("Per-Chromosome Statistics:")
    print(chr_stats)
    print()
    
    # Genetic diversity analysis
    hash_diversity = df.groupby('gene_id')['sha256_hash'].nunique()
    print("Genetic Diversity Across Trio:")
    print(f"  Genes with identical sequences (all 6 haplotypes): {sum(hash_diversity == 1):,}")
    print(f"  Genes with 2-3 variants: {sum((hash_diversity >= 2) & (hash_diversity <= 3)):,}")
    print(f"  Genes with 4-6 variants: {sum(hash_diversity >= 4):,}")
    print(f"  Most diverse gene: {hash_diversity.idxmax()} ({hash_diversity.max()} variants)")
    print()
    
    # Sample comparison
    sample_stats = df.groupby('sample').agg({
        'gene_id': 'nunique',
        'clean_length': 'mean'
    }).round(1)
    print("Sample Statistics:")
    print(sample_stats)
    print()
    
    # Length distribution
    print("Sequence Length Statistics:")
    print(f"  Mean length: {df['clean_length'].mean():.1f} bp")
    print(f"  Median length: {df['clean_length'].median():.1f} bp") 
    print(f"  Range: {df['clean_length'].min()}-{df['clean_length'].max()} bp")
    print(f"  Std deviation: {df['clean_length'].std():.1f} bp")
    
except Exception as e:
    print(f"Error generating summary: {e}")
EOF
        
        log "Combined results saved:"
        log "  Hashes: $combined_hashes"
        log "  Summary: $combined_summary"
        
        local total_records=$(tail -n +2 "$combined_hashes" | wc -l)
        log "Total records: $total_records"
    else
        log "WARNING: No hash files found to combine"
    fi
}

# Main function for multi-chromosome runner
main() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -c|--chromosomes)
                CHROMOSOMES_TO_PROCESS="$2"
                shift 2
                ;;
            -p|--parallel)
                PARALLEL_MODE="true"
                shift
                ;;
            -j|--jobs)
                MAX_JOBS="$2"
                shift 2
                ;;
            -d|--download)
                DOWNLOAD_MODE="true"
                shift
                ;;
            --config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                log "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done
    
    # Load configuration
    load_config
    
    # Use chromosomes from config if not specified
    if [[ -z "$CHROMOSOMES_TO_PROCESS" ]]; then
        CHROMOSOMES_TO_PROCESS="$CHROMOSOMES"
    fi
    
    if [[ -z "$CHROMOSOMES_TO_PROCESS" ]]; then
        log "ERROR: No chromosomes specified"
        exit 1
    fi
    
    # Check if pipeline script exists
    if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
        log "ERROR: Pipeline script not found: $PIPELINE_SCRIPT"
        exit 1
    fi
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    log "Starting multi-chromosome analysis"
    log "Chromosomes: $CHROMOSOMES_TO_PROCESS"
    log "Mode: $([ "$PARALLEL_MODE" == "true" ] && echo "parallel ($MAX_JOBS jobs)" || echo "sequential")"
    
    # Process chromosomes
    if [[ "$PARALLEL_MODE" == "true" ]]; then
        process_parallel
    else
        process_sequential
    fi
    
    local exit_code=$?
    
    # Combine results
    combine_results
    
    log "Multi-chromosome analysis completed with exit code: $exit_code"
    return $exit_code
}

#############################################################################
# Quick Start Script for 1000 Genomes
# quick_start_1000g.sh
#############################################################################

#!/bin/bash

# Quick start script for 1000 Genomes trio analysis
# Downloads reference and annotation files, then processes chromosome 22

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${SCRIPT_DIR}/1000g_quickstart"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

setup_1000g_analysis() {
    log "Setting up 1000 Genomes trio analysis"
    
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"
    
    # Download reference genome (chromosome 22 only for quick start)
    log "Downloading reference genome (chr22)..."
    if [[ ! -f "chr22.fa" ]]; then
        wget -O chr22.fa.gz "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_chr22_reference.fa.gz"
        gunzip chr22.fa.gz
        samtools faidx chr22.fa
    fi
    
    # Download gene annotations
    log "Downloading gene annotations..."
    if [[ ! -f "gencode.v44.annotation.gff3" ]]; then
        wget -O gencode.v44.annotation.gff3.gz "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gff3.gz"
        gunzip gencode.v44.annotation.gff3.gz
    fi
    
    # Create config file
    log "Creating configuration file..."
    cat > pipeline_config.conf << EOF
# Quick start configuration for 1000 Genomes trio analysis
REF="$WORK_DIR/chr22.fa"
GFF="$WORK_DIR/gencode.v44.annotation.gff3"

# CEU trio (Ashkenazi Jewish)
FATHER_ID="NA12877"
MOTHER_ID="NA12878"
CHILD_ID="NA12882"

# 1000G VCF URL
VCF_URL_TEMPLATE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{CHR}.recalibrated_variants.vcf.gz"

OUTPUT_DIR="$WORK_DIR/results"
THREADS=4
MIN_GENE_LENGTH=300
TEMP_DIR="$WORK_DIR/temp"
CHROMOSOMES="chr22"
EOF
    
    log "Setup complete! Ready to run analysis."
    log "Run: ${SCRIPT_DIR}/trio_consensus_pipeline.sh --config $WORK_DIR/pipeline_config.conf -c chr22 -d"
}

# Run setup
setup_1000g_analysis

#############################################################################
# Example Usage Commands
#############################################################################

# Save these as executable scripts:

# 1. Single chromosome analysis:
# ./trio_consensus_pipeline.sh -c chr22 -d \
#   -v "ftp://ftp.1000genomes.ebi.ac.uk/.../chr22.vcf.gz" \
#   -r hg38_chr22.fa -g gencode.v44.gff3 \
#   -f NA12877 -m NA12878 -k NA12882

# 2. Multiple chromosomes sequentially:
# ./run_multi_chromosome.sh -c "chr21 chr22"

# 3. Parallel processing:
# ./run_multi_chromosome.sh -c "chr1 chr2 chr3" --parallel -j 3

# 4. All autosomes:
# ./run_multi_chromosome.sh -c "$(seq -f 'chr%.0f' 1 22 | tr '\n' ' ')" --parallel

# 5. Quick start with chromosome 22:
# ./quick_start_1000g.sh
# ./trio_consensus_pipeline.sh --config 1000g_quickstart/pipeline_config.conf -c chr22 -d