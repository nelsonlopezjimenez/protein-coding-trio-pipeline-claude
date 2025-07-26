#!/bin/bash

# Chromosome-Specific Protein-Coding Trio Pipeline
# Processes one chromosome at a time focusing only on protein-coding genes
# Optimized for large datasets like 1000 Genomes

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.conf"

# Default values
REF=""
INPUT_VCF=""
GFF=""
OUTPUT_DIR="trio_protein_coding_output"
CHROMOSOME=""
FATHER_ID=""
MOTHER_ID=""
CHILD_ID=""
DOWNLOAD_MODE="false"
THREADS=4
MIN_GENE_LENGTH=300
TEMP_DIR="/tmp/trio_pipeline_$$"

# Function to print usage
usage() {
    cat << EOF
Chromosome-Specific Protein-Coding Trio Pipeline

USAGE:
    $0 [OPTIONS]

OPTIONS:
    -c, --chromosome CHROM     Chromosome to process (e.g., chr22, 22)
    -v, --vcf FILE            Input VCF file (can be URL for download mode)
    -r, --reference FILE      Reference genome FASTA file
    -g, --gff FILE            Gene annotation file (GFF3/GTF)
    -o, --output DIR          Output directory (default: trio_protein_coding_output)
    -f, --father ID           Father sample ID
    -m, --mother ID           Mother sample ID
    -k, --child ID            Child sample ID
    -d, --download            Enable download mode for remote files
    -t, --threads NUM         Number of threads (default: 4)
    -l, --min-length NUM      Minimum gene length (default: 300)
    --temp-dir DIR            Temporary directory (default: /tmp)
    --config FILE             Configuration file
    -h, --help                Show this help

EXAMPLES:
    # Process chromosome 22 from local files
    $0 -c chr22 -v chr22.vcf.gz -r ref.fa -g genes.gff3 \\
       -f NA12877 -m NA12878 -k NA12882

    # Download and process 1000G chromosome 22
    $0 -c chr22 -d \\
       -v "ftp://ftp.1000genomes.ebi.ac.uk/.../chr22.vcf.gz" \\
       -r hg38.fa -g gencode.gff3 \\
       -f NA12877 -m NA12878 -k NA12882

    # Process multiple chromosomes with config file
    $0 --config trio_config.conf -c chr1
    $0 --config trio_config.conf -c chr2

CONFIGURATION FILE FORMAT:
    REF=/path/to/reference.fa
    GFF=/path/to/annotation.gff3
    FATHER_ID=NA12877
    MOTHER_ID=NA12878
    CHILD_ID=NA12882
    OUTPUT_DIR=my_output
    THREADS=8
EOF
}

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

# Function to check dependencies
check_dependencies() {
    local deps=("bcftools" "bedtools" "awk" "python3")
    
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            log "ERROR: Required dependency '$dep' not found"
            exit 1
        fi
    done
    
    # Check bcftools plugins
    if ! bcftools plugin -l | grep -q "trio-phase"; then
        log "WARNING: bcftools trio-phase plugin not available"
        log "Trio phasing will be skipped"
    fi
}

# Function to load configuration
load_config() {
    if [[ -f "$CONFIG_FILE" ]]; then
        log "Loading configuration from $CONFIG_FILE"
        source "$CONFIG_FILE"
    fi
}

# Function to validate inputs
validate_inputs() {
    local errors=0
    
    if [[ -z "$CHROMOSOME" ]]; then
        log "ERROR: Chromosome must be specified (-c/--chromosome)"
        ((errors++))
    fi
    
    if [[ -z "$REF" ]]; then
        log "ERROR: Reference genome must be specified (-r/--reference)"
        ((errors++))
    fi
    
    if [[ -z "$INPUT_VCF" ]]; then
        log "ERROR: Input VCF must be specified (-v/--vcf)"
        ((errors++))
    fi
    
    if [[ -z "$GFF" ]]; then
        log "ERROR: Gene annotation must be specified (-g/--gff)"
        ((errors++))
    fi
    
    if [[ -z "$FATHER_ID" || -z "$MOTHER_ID" || -z "$CHILD_ID" ]]; then
        log "ERROR: All trio sample IDs must be specified"
        ((errors++))
    fi
    
    if [[ "$DOWNLOAD_MODE" == "false" ]]; then
        if [[ ! -f "$REF" ]]; then
            log "ERROR: Reference file not found: $REF"
            ((errors++))
        fi
        
        if [[ ! -f "$INPUT_VCF" ]]; then
            log "ERROR: VCF file not found: $INPUT_VCF"
            ((errors++))
        fi
        
        if [[ ! -f "$GFF" ]]; then
            log "ERROR: GFF file not found: $GFF"
            ((errors++))
        fi
    fi
    
    if [[ $errors -gt 0 ]]; then
        log "Found $errors error(s). Exiting."
        exit 1
    fi
}

# Function to download file if needed
download_if_needed() {
    local url="$1"
    local output="$2"
    
    if [[ "$url" =~ ^(http|ftp):// ]]; then
        log "Downloading $url to $output"
        if command -v wget &> /dev/null; then
            wget -O "$output" "$url"
        elif command -v curl &> /dev/null; then
            curl -o "$output" "$url"
        else
            log "ERROR: Neither wget nor curl available for downloading"
            exit 1
        fi
    else
        # Local file - create symlink or copy
        if [[ ! -f "$output" ]]; then
            ln -s "$(realpath "$url")" "$output" 2>/dev/null || cp "$url" "$output"
        fi
    fi
}

# Function to extract protein-coding genes for chromosome
extract_protein_coding_genes() {
    local gff_file="$1"
    local chromosome="$2"
    local output_bed="$3"
    
    log "Extracting protein-coding genes for $chromosome"
    
    # Normalize chromosome name (handle both "chr22" and "22" formats)
    local chr_patterns="$chromosome"
    if [[ "$chromosome" =~ ^chr ]]; then
        chr_patterns="$chromosome|${chromosome#chr}"
    else
        chr_patterns="$chromosome|chr$chromosome"
    fi
    
    # Extract protein-coding genes from GFF
    awk -v chr_pat="^($chr_patterns)$" -v min_len="$MIN_GENE_LENGTH" '
    BEGIN { OFS="\t" }
    $1 ~ chr_pat && $3 == "gene" && $9 ~ /gene_type[= ]"?protein_coding"?/ {
        # Extract gene ID
        match($9, /gene_id[= ]"?([^";]+)"?/, id)
        gene_id = id[1]
        if (gene_id == "") {
            match($9, /ID=([^;]+)/, id)
            gene_id = id[1]
        }
        if (gene_id == "") gene_id = "gene_" NR
        
        # Check minimum length
        gene_length = $5 - $4 + 1
        if (gene_length >= min_len) {
            print $1, $4-1, $5, gene_id, gene_length, $7
        }
    }' "$gff_file" > "$output_bed"
    
    local gene_count=$(wc -l < "$output_bed")
    log "Found $gene_count protein-coding genes for $chromosome (min length: $MIN_GENE_LENGTH bp)"
    
    if [[ $gene_count -eq 0 ]]; then
        log "WARNING: No protein-coding genes found for $chromosome"
        log "Check chromosome naming in GFF file (chr22 vs 22)"
    fi
}

# Function to process VCF for chromosome
process_chromosome_vcf() {
    local input_vcf="$1"
    local chromosome="$2"
    local output_vcf="$3"
    
    log "Processing VCF for $chromosome"
    
    # Check if samples exist in VCF
    if ! bcftools query -l "$input_vcf" | grep -q "$FATHER_ID"; then
        log "ERROR: Father ID '$FATHER_ID' not found in VCF"
        exit 1
    fi
    
    if ! bcftools query -l "$input_vcf" | grep -q "$MOTHER_ID"; then
        log "ERROR: Mother ID '$MOTHER_ID' not found in VCF"
        exit 1
    fi
    
    if ! bcftools query -l "$input_vcf" | grep -q "$CHILD_ID"; then
        log "ERROR: Child ID '$CHILD_ID' not found in VCF"
        exit 1
    fi
    
    # Extract trio samples for the specific chromosome
    bcftools view -r "$chromosome" -s "${FATHER_ID},${MOTHER_ID},${CHILD_ID}" \
        "$input_vcf" --threads "$THREADS" -O z -o "$output_vcf"
    
    bcftools index --threads "$THREADS" "$output_vcf"
    
    local variant_count=$(bcftools view -H "$output_vcf" | wc -l)
    log "Extracted $variant_count variants for trio on $chromosome"
}

# Function to phase trio variants
phase_trio_variants() {
    local input_vcf="$1"
    local output_vcf="$2"
    local ped_file="$3"
    
    log "Phasing trio variants"
    
    # Create pedigree file
    cat > "$ped_file" << EOF
FAM001	$FATHER_ID	0	0	1	1
FAM001	$MOTHER_ID	0	0	2	1
FAM001	$CHILD_ID	$FATHER_ID	$MOTHER_ID	2	1
EOF
    
    # Phase variants using trio information if plugin available
    if bcftools plugin -l | grep -q "trio-phase"; then
        bcftools +trio-phase "$input_vcf" -p "$ped_file" \
            --threads "$THREADS" -O z -o "$output_vcf"
    else
        log "WARNING: Trio-phase plugin not available, using statistical phasing"
        bcftools phase "$input_vcf" --threads "$THREADS" -O z -o "$output_vcf"
    fi
    
    bcftools index --threads "$THREADS" "$output_vcf"
}

# Function to generate consensus sequences
generate_consensus_sequences() {
    local vcf_file="$1"
    local ref_file="$2"
    local output_dir="$3"
    
    log "Generating consensus sequences for trio"
    
    # Create individual VCFs
    for sample in "$FATHER_ID" "$MOTHER_ID" "$CHILD_ID"; do
        local sample_vcf="${output_dir}/${sample}_${CHROMOSOME}.vcf.gz"
        bcftools view -s "$sample" "$vcf_file" --threads "$THREADS" -O z -o "$sample_vcf"
        bcftools index --threads "$THREADS" "$sample_vcf"
        
        # Generate consensus for each haplotype
        log "Generating consensus for $sample"
        bcftools consensus -f "$ref_file" -H 1 "$sample_vcf" > "${output_dir}/${sample}_${CHROMOSOME}_hap1.fa"
        bcftools consensus -f "$ref_file" -H 2 "$sample_vcf" > "${output_dir}/${sample}_${CHROMOSOME}_hap2.fa"
    done
}

# Function to extract gene sequences
extract_gene_sequences() {
    local genes_bed="$1"
    local output_dir="$2"
    
    log "Extracting protein-coding gene sequences"
    
    for sample in "$FATHER_ID" "$MOTHER_ID" "$CHILD_ID"; do
        for hap in hap1 hap2; do
            local consensus_fa="${output_dir}/${sample}_${CHROMOSOME}_${hap}.fa"
            local genes_fa="${output_dir}/${sample}_${CHROMOSOME}_${hap}_genes.fa"
            
            if [[ -f "$consensus_fa" && -s "$genes_bed" ]]; then
                bedtools getfasta -fi "$consensus_fa" -bed "$genes_bed" -name -s > "$genes_fa"
                local gene_count=$(grep -c "^>" "$genes_fa" || echo "0")
                log "Extracted $gene_count gene sequences for ${sample}_${hap}"
            else
                log "WARNING: Missing consensus file or empty gene list for ${sample}_${hap}"
                touch "$genes_fa"  # Create empty file
            fi
        done
    done
}

# Function to generate hashes
generate_gene_hashes() {
    local output_dir="$1"
    local final_output="$2"
    
    log "Generating SHA-256 hashes for gene sequences"
    
    # Create CSV header
    echo "chromosome,sample,haplotype,gene_id,gene_length,gene_strand,raw_length,clean_length,sha256_hash" > "$final_output"
    
    for sample in "$FATHER_ID" "$MOTHER_ID" "$CHILD_ID"; do
        for hap in hap1 hap2; do
            local genes_fa="${output_dir}/${sample}_${CHROMOSOME}_${hap}_genes.fa"
            
            if [[ -f "$genes_fa" && -s "$genes_fa" ]]; then
                log "Hashing ${sample}_${hap} genes"
                
                # Process FASTA and generate hashes
                awk -v chr="$CHROMOSOME" -v sample="$sample" -v hap="$hap" '
                /^>/ {
                    if (gene_id && sequence) {
                        # Parse gene info from header
                        split(gene_info, parts, ":")
                        gene_length = parts[3] - parts[2] + 1
                        gene_strand = parts[4] ? parts[4] : "+"
                        
                        # Clean sequence
                        raw_len = length(sequence)
                        gsub(/[-N]/, "", sequence)
                        sequence = toupper(sequence)
                        clean_len = length(sequence)
                        
                        if (clean_len > 0) {
                            printf "%s,%s,%s,%s,%d,%s,%d,%d,%s\n", 
                                chr, sample, hap, gene_id, gene_length, gene_strand, 
                                raw_len, clean_len, sequence
                        }
                    }
                    
                    # Parse new header: >gene_id::chr:start-end or >gene_id
                    gene_line = substr($0, 2)
                    if (match(gene_line, /^([^:]+)::/)) {
                        gene_id = substr(gene_line, 1, RLENGTH-2)
                        gene_info = substr(gene_line, RLENGTH+1)
                    } else {
                        gene_id = gene_line
                        gene_info = ""
                    }
                    sequence = ""
                    next
                }
                {sequence = sequence $0}
                END {
                    if (gene_id && sequence) {
                        split(gene_info, parts, ":")
                        gene_length = parts[3] - parts[2] + 1
                        gene_strand = parts[4] ? parts[4] : "+"
                        
                        raw_len = length(sequence)
                        gsub(/[-N]/, "", sequence)
                        sequence = toupper(sequence)
                        clean_len = length(sequence)
                        
                        if (clean_len > 0) {
                            printf "%s,%s,%s,%s,%d,%s,%d,%d,%s\n", 
                                chr, sample, hap, gene_id, gene_length, gene_strand, 
                                raw_len, clean_len, sequence
                        }
                    }
                }' "$genes_fa" | \
                while IFS=',' read -r chr sample hap gene_id gene_length gene_strand raw_len clean_len sequence; do
                    if [[ -n "$sequence" ]]; then
                        hash=$(echo -n "$sequence" | sha256sum | cut -d' ' -f1)
                        echo "$chr,$sample,$hap,$gene_id,$gene_length,$gene_strand,$raw_len,$clean_len,$hash"
                    fi
                done >> "$final_output"
                
            else
                log "WARNING: No gene sequences found for ${sample}_${hap}"
            fi
        done
    done
}

# Function to generate summary
generate_summary() {
    local hash_file="$1"
    local summary_file="$2"
    
    log "Generating summary statistics"
    
    python3 << EOF > "$summary_file"
import pandas as pd
import sys
from collections import Counter

try:
    df = pd.read_csv('$hash_file')
    
    print(f"Chromosome-Specific Protein-Coding Gene Analysis Summary")
    print(f"=" * 60)
    print(f"Chromosome: $CHROMOSOME")
    print(f"Processing date: $(date)")
    print(f"")
    
    if len(df) == 0:
        print("No gene sequences processed")
        sys.exit(0)
    
    print(f"Sample Summary:")
    print(f"  Father: $FATHER_ID")
    print(f"  Mother: $MOTHER_ID") 
    print(f"  Child: $CHILD_ID")
    print(f"")
    
    print(f"Processing Statistics:")
    print(f"  Total gene sequences: {len(df)}")
    print(f"  Unique genes: {df['gene_id'].nunique()}")
    print(f"  Unique hash values: {df['sha256_hash'].nunique()}")
    print(f"  Average clean sequence length: {df['clean_length'].mean():.1f} bp")
    print(f"  Sequence length range: {df['clean_length'].min()}-{df['clean_length'].max()} bp")
    print(f"")
    
    # Hash diversity per gene
    hash_diversity = df.groupby('gene_id')['sha256_hash'].nunique()
    print(f"Genetic Diversity:")
    print(f"  Genes with 1 unique sequence: {sum(hash_diversity == 1)}")
    print(f"  Genes with 2-3 unique sequences: {sum((hash_diversity >= 2) & (hash_diversity <= 3))}")
    print(f"  Genes with 4+ unique sequences: {sum(hash_diversity >= 4)}")
    print(f"  Most diverse gene: {hash_diversity.idxmax()} ({hash_diversity.max()} variants)")
    print(f"")
    
    # Sample-specific stats
    sample_stats = df.groupby(['sample', 'haplotype']).agg({
        'gene_id': 'count',
        'clean_length': 'mean'
    }).round(1)
    print(f"Sample-Haplotype Statistics:")
    print(sample_stats)
    
except Exception as e:
    print(f"Error in analysis: {e}")
    print(f"Basic summary: {$(wc -l < "$hash_file")} total records")
EOF
}

# Function to cleanup
cleanup() {
    if [[ -d "$TEMP_DIR" ]]; then
        log "Cleaning up temporary files"
        rm -rf "$TEMP_DIR"
    fi
}

# Main function
main() {
    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            -c|--chromosome)
                CHROMOSOME="$2"
                shift 2
                ;;
            -v|--vcf)
                INPUT_VCF="$2"
                shift 2
                ;;
            -r|--reference)
                REF="$2"
                shift 2
                ;;
            -g|--gff)
                GFF="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -f|--father)
                FATHER_ID="$2"
                shift 2
                ;;
            -m|--mother)
                MOTHER_ID="$2"
                shift 2
                ;;
            -k|--child)
                CHILD_ID="$2"
                shift 2
                ;;
            -d|--download)
                DOWNLOAD_MODE="true"
                shift
                ;;
            -t|--threads)
                THREADS="$2"
                shift 2
                ;;
            -l|--min-length)
                MIN_GENE_LENGTH="$2"
                shift 2
                ;;
            --temp-dir)
                TEMP_DIR="$2"
                shift 2
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
    
    # Load configuration if specified
    load_config
    
    # Check dependencies
    check_dependencies
    
    # Validate inputs
    validate_inputs
    
    # Create directories
    mkdir -p "$OUTPUT_DIR"
    mkdir -p "$TEMP_DIR"
    
    # Set up cleanup trap
    trap cleanup EXIT
    
    log "Starting chromosome-specific protein-coding gene analysis"
    log "Chromosome: $CHROMOSOME"
    log "Output directory: $OUTPUT_DIR"
    
    # Download files if needed
    local working_vcf="$TEMP_DIR/input.vcf.gz"
    local working_ref="$TEMP_DIR/reference.fa"
    local working_gff="$TEMP_DIR/annotation.gff3"
    
    if [[ "$DOWNLOAD_MODE" == "true" ]]; then
        download_if_needed "$INPUT_VCF" "$working_vcf"
        download_if_needed "$REF" "$working_ref"
        download_if_needed "$GFF" "$working_gff"
    else
        working_vcf="$INPUT_VCF"
        working_ref="$REF"
        working_gff="$GFF"
    fi
    
    # Extract protein-coding genes for chromosome
    local genes_bed="$TEMP_DIR/protein_coding_genes_${CHROMOSOME}.bed"
    extract_protein_coding_genes "$working_gff" "$CHROMOSOME" "$genes_bed"
    
    # Process VCF for chromosome and trio
    local trio_vcf="$TEMP_DIR/trio_${CHROMOSOME}.vcf.gz"
    process_chromosome_vcf "$working_vcf" "$CHROMOSOME" "$trio_vcf"
    
    # Phase trio variants
    local phased_vcf="$TEMP_DIR/trio_${CHROMOSOME}_phased.vcf.gz"
    local ped_file="$TEMP_DIR/trio.ped"
    phase_trio_variants "$trio_vcf" "$phased_vcf" "$ped_file"
    
    # Generate consensus sequences
    generate_consensus_sequences "$phased_vcf" "$working_ref" "$OUTPUT_DIR"
    
    # Extract gene sequences
    extract_gene_sequences "$genes_bed" "$OUTPUT_DIR"
    
    # Generate hashes
    local hash_output="$OUTPUT_DIR/trio_${CHROMOSOME}_protein_coding_hashes.csv"
    generate_gene_hashes "$OUTPUT_DIR" "$hash_output"
    
    # Generate summary
    local summary_output="$OUTPUT_DIR/trio_${CHROMOSOME}_summary.txt"
    generate_summary "$hash_output" "$summary_output"
    
    # Copy gene bed file to output for reference
    cp "$genes_bed" "$OUTPUT_DIR/protein_coding_genes_${CHROMOSOME}.bed"
    
    log "Pipeline completed successfully!"
    log "Results:"
    log "  Hashes: $hash_output"
    log "  Summary: $summary_output"
    log "  Gene coordinates: $OUTPUT_DIR/protein_coding_genes_${CHROMOSOME}.bed"
    
    # Print final statistics
    if [[ -f "$hash_output" ]]; then
        local total_hashes=$(tail -n +2 "$hash_output" | wc -l)
        log "Total gene sequence hashes generated: $total_hashes"
    fi
}

# Run main function
main "$@"