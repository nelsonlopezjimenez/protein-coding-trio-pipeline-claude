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
    $0 -c "\$(seq -f 'chr%.0f' 1 22 | tr '\n' ' ')"

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
    local chr_for_pid=()
    local results=()
    
    log "Starting parallel processing of ${#chromosomes[@]} chromosomes with max $MAX_JOBS jobs"
    
    for chr in "${chromosomes[@]}"; do
        # Wait if we've reached max jobs
        while [[ ${#pids[@]} -ge $MAX_JOBS ]]; do
            for i in "${!pids[@]}"; do
                if ! kill -0 "${pids[$i]}" 2>/dev/null; then
                    wait "${pids[$i]}"
                    local exit_code=$?
                    results+=("${chr_for_pid[$i]}:$exit_code")
                    unset pids[$i]
                    unset chr_for_pid[$i]
                fi
            done
            # Compact arrays
            local new_pids=()
            local new_chr_for_pid=()
            for i in "${!pids[@]}"; do
                new_pids+=("${pids[$i]}")
                new_chr_for_pid+=("${chr_for_pid[$i]}")
            done
            pids=("${new_pids[@]}")
            chr_for_pid=("${new_chr_for_pid[@]}")
            sleep 2
        done
        
        # Start new job
        log "Starting chromosome $chr"
        process_chromosome "$chr" &
        local pid=$!
        pids+=($pid)
        chr_for_pid+=("$chr")
    done
    
    # Wait for remaining jobs
    for i in "${!pids[@]}"; do
        wait "${pids[$i]}"
        local exit_code=$?
        results+=("${chr_for_pid[$i]}:$exit_code")
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
                log "Combined results from $chr"
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

# Check if pipeline script exists
check_pipeline() {
    if [[ ! -f "$PIPELINE_SCRIPT" ]]; then
        log "ERROR: Pipeline script not found: $PIPELINE_SCRIPT"
        log "Make sure trio_consensus_pipeline.sh is in the same directory"
        exit 1
    fi
    
    if [[ ! -x "$PIPELINE_SCRIPT" ]]; then
        log "ERROR: Pipeline script is not executable: $PIPELINE_SCRIPT"
        log "Run: chmod +x $PIPELINE_SCRIPT"
        exit 1
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
    
    # Check pipeline script
    check_pipeline
    
    # Load configuration
    load_config
    
    # Use chromosomes from config if not specified
    if [[ -z "$CHROMOSOMES_TO_PROCESS" ]]; then
        CHROMOSOMES_TO_PROCESS="$CHROMOSOMES"
    fi
    
    if [[ -z "$CHROMOSOMES_TO_PROCESS" ]]; then
        log "ERROR: No chromosomes specified"
        log "Use -c option or set CHROMOSOMES in config file"
        exit 1
    fi
    
    # Create output directory
    mkdir -p "$OUTPUT_DIR"
    
    log "Starting multi-chromosome analysis"
    log "Chromosomes: $CHROMOSOMES_TO_PROCESS"
    log "Mode: $([ "$PARALLEL_MODE" == "true" ] && echo "parallel ($MAX_JOBS jobs)" || echo "sequential")"
    log "Config: $CONFIG_FILE"
    
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

# Run main function if script is executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi