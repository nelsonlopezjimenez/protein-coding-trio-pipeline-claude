#!/bin/bash

# Quick Start Script for 1000 Genomes Trio Analysis
# Downloads reference and annotation files, then sets up configuration for chromosome 22

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK_DIR="${SCRIPT_DIR}/1000g_quickstart"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >&2
}

check_dependencies() {
    local deps=("wget" "samtools" "bcftools" "bedtools")
    local missing=()
    
    for dep in "${deps[@]}"; do
        if ! command -v "$dep" &> /dev/null; then
            missing+=("$dep")
        fi
    done
    
    if [[ ${#missing[@]} -gt 0 ]]; then
        log "ERROR: Missing dependencies: ${missing[*]}"
        log "Please install the missing tools and try again"
        exit 1
    fi
    
    log "All dependencies found"
}

setup_1000g_analysis() {
    log "Setting up 1000 Genomes trio analysis workspace"
    
    mkdir -p "$WORK_DIR"
    cd "$WORK_DIR"
    
    # Use existing GRCh38 reference file
    log "Setting up reference genome..."
    
    # Prompt user for their GRCh38 file location
    if [[ -z "${GRCH38_REF:-}" ]]; then
        echo
        echo "Please provide the path to your GRCh38 reference file:"
        echo "Examples:"
        echo "  /path/to/GRCh38.fa"
        echo "  /path/to/hg38.fa.gz"
        echo "  /path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        echo
        read -p "GRCh38 reference file path: " GRCH38_REF
    fi
    
    if [[ ! -f "$GRCH38_REF" ]]; then
        log "ERROR: Reference file not found: $GRCH38_REF"
        exit 1
    fi
    
    log "Using reference file: $GRCH38_REF"
    
    # Handle compressed vs uncompressed reference
    if [[ "$GRCH38_REF" == *.gz ]]; then
        log "  Detected compressed reference, creating symlink..."
        ln -sf "$GRCH38_REF" "hg38.fa.gz"
        
        # Check if index exists
        if [[ -f "${GRCH38_REF}.gzi" && -f "${GRCH38_REF}.fai" ]]; then
            log "  Using existing index files"
            ln -sf "${GRCH38_REF}.gzi" "hg38.fa.gz.gzi"
            ln -sf "${GRCH38_REF}.fai" "hg38.fa.gz.fai"
        else
            log "  Creating index for compressed reference..."
            samtools faidx "hg38.fa.gz"
        fi
        REF_FILE="$WORK_DIR/hg38.fa.gz"
    else
        log "  Detected uncompressed reference, creating symlink..."
        ln -sf "$GRCH38_REF" "hg38.fa"
        
        # Check if index exists
        if [[ -f "${GRCH38_REF}.fai" ]]; then
            log "  Using existing index file"
            ln -sf "${GRCH38_REF}.fai" "hg38.fa.fai"
        else
            log "  Creating index..."
            samtools faidx "hg38.fa"
        fi
        REF_FILE="$WORK_DIR/hg38.fa"
    fi
    
    # Verify reference has chr22
    log "  Verifying chromosome 22 is present..."
    if samtools faidx "$REF_FILE" -l | grep -q -E "^(chr22|22)$"; then
        log "  ✓ Chromosome 22 found in reference"
    else
        log "  WARNING: chr22/22 not found in reference. Available chromosomes:"
        samtools faidx "$REF_FILE" -l | head -10
        log "  You may need to adjust chromosome naming (chr22 vs 22)"
    fi
    
    # Try multiple sources for VCF data
    log "Downloading 1000 Genomes VCF data..."
    if [[ ! -f "chr22_test.vcf.gz" ]]; then
        log "  Trying HTTPS source..."
        if wget -O chr22_test.vcf.gz "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.recalibrated_variants.vcf.gz" 2>/dev/null; then
            log "  Successfully downloaded from HTTPS"
        else
            log "  HTTPS failed, trying Phase 3 data..."
            if wget -O chr22_test.vcf.gz "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" 2>/dev/null; then
                log "  Successfully downloaded Phase 3 data"
            else
                log "  Download failed, will use remote URL in pipeline"
                rm -f chr22_test.vcf.gz
            fi
        fi
        
        if [[ -f "chr22_test.vcf.gz" ]]; then
            bcftools index chr22_test.vcf.gz
            local variant_count=$(bcftools view -H chr22_test.vcf.gz | wc -l)
            log "  VCF ready: $variant_count variants"
        fi
    else
        log "  VCF file already exists"
    fi
        log "  Fetching GENCODE v44 annotations..."
        wget -O gencode.v44.annotation.gff3.gz "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gff3.gz"
        
        log "  Extracting annotations..."
        gunzip gencode.v44.annotation.gff3.gz
        
        # Count protein-coding genes on chr22
        local gene_count=$(awk '$1=="chr22" && $3=="gene" && $9~/gene_type[= ]"?protein_coding"?/' gencode.v44.annotation.gff3 | wc -l)
        log "  Annotations ready: $gene_count protein-coding genes on chr22"
    else
        log "  Annotations already exist"
    fi
    
    # Create config file
    log "Creating configuration file..."
    cat > pipeline_config.conf << EOF
# Quick start configuration for 1000 Genomes trio analysis
# Generated on $(date)

# Local reference files
REF="$REF_FILE"
GFF="$WORK_DIR/gencode.v44.annotation.gff3"

# CEU trio (Ashkenazi Jewish family)
# Father: NA12877, Mother: NA12878, Child: NA12882
FATHER_ID="NA12877"
MOTHER_ID="NA12878"
CHILD_ID="NA12882"

# 1000 Genomes VCF URL template (with fallbacks)
VCF_URL_TEMPLATE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_{CHR}.recalibrated_variants.vcf.gz"

# Fallback URL for Phase 3 data if high-coverage fails
VCF_FALLBACK_TEMPLATE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.{CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

# Output and processing settings
OUTPUT_DIR="$WORK_DIR/results"
THREADS=4
MIN_GENE_LENGTH=300
TEMP_DIR="$WORK_DIR/temp"

# Chromosomes for testing (start with chr22)
CHROMOSOMES="chr22"
EOF
    
    # Create run script for easy execution
    log "Creating convenient run script..."
    cat > run_chr22_analysis.sh << 'EOF'
#!/bin/bash

# Convenient script to run chromosome 22 analysis
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE="${SCRIPT_DIR}/../trio_consensus_pipeline.sh"

if [[ ! -f "$PIPELINE" ]]; then
    echo "ERROR: Pipeline script not found at $PIPELINE"
    echo "Make sure trio_consensus_pipeline.sh is in the parent directory"
    exit 1
fi

echo "Running chromosome 22 trio analysis..."
echo "This will download ~47MB VCF file and process protein-coding genes"
echo

"$PIPELINE" \
    --config "$SCRIPT_DIR/pipeline_config.conf" \
    --chromosome chr22 \
    --download \
    --threads 4

echo
echo "Analysis complete! Check results in: $SCRIPT_DIR/results/chr22/"
EOF
    
    chmod +x run_chr22_analysis.sh
    
    # Create summary of what was set up
    log "Creating setup summary..."
    cat > setup_summary.txt << EOF
1000 Genomes Trio Analysis - Quick Start Setup
==============================================

Setup completed on: $(date)
Working directory: $WORK_DIR

Files created:
- hg38.fa (or hg38.fa.gz)       : Full GRCh38 reference genome (symlink)
- hg38.fa.fai (or hg38.fa.gz.fai) : Reference index
- gencode.v44.annotation.gff3   : Gene annotations
- pipeline_config.conf          : Configuration file
- run_chr22_analysis.sh         : Convenient run script
- setup_summary.txt             : This file

Reference source: $GRCH38_REF

Trio configuration:
- Father:  NA12877 (CEU population)
- Mother:  NA12878 (CEU population) 
- Child:   NA12882 (CEU population)

Next steps:
1. Run the analysis:
   ./run_chr22_analysis.sh

2. Or run manually:
   ../trio_consensus_pipeline.sh --config pipeline_config.conf -c chr22 -d

3. Check results:
   ls results/chr22/

Expected output:
- ~500 protein-coding genes on chromosome 22
- 6 haplotype sequences per gene (3 samples × 2 haplotypes)
- SHA-256 hashes for genetic diversity analysis
- Summary statistics

File sizes:
- Reference: Full GRCh38 genome (~3.1GB)
- Annotations: ~$(du -h gencode.v44.annotation.gff3 2>/dev/null | cut -f1 || echo "~40MB")
- VCF download: ~47MB (downloaded automatically)
EOF
    
    log "Setup complete!"
    log ""
    log "Summary:"
    log "  Working directory: $WORK_DIR"
    log "  Reference genome: chr22.fa"
    log "  Gene annotations: gencode.v44.annotation.gff3"
    log "  Configuration: pipeline_config.conf"
    log ""
    log "To run the analysis:"
    log "  cd $WORK_DIR"
    log "  ./run_chr22_analysis.sh"
    log ""
    log "Or manually:"
    log "  ${SCRIPT_DIR}/trio_consensus_pipeline.sh --config $WORK_DIR/pipeline_config.conf -c chr22 -d"
    log ""
    log "Setup summary saved to: $WORK_DIR/setup_summary.txt"
}

# Function to test the setup
test_setup() {
    log "Testing setup..."
    
    if [[ ! -f "$WORK_DIR/chr22.fa" ]]; then
        log "ERROR: Reference file not found"
        return 1
    fi
    
    if [[ ! -f "$WORK_DIR/gencode.v44.annotation.gff3" ]]; then
        log "ERROR: Annotation file not found"
        return 1
    fi
    
    if [[ ! -f "$WORK_DIR/pipeline_config.conf" ]]; then
        log "ERROR: Configuration file not found"
        return 1
    fi
    
    # Test that we can extract genes
    local gene_count=$(awk '$1=="chr22" && $3=="gene" && $9~/gene_type[= ]"?protein_coding"?/' "$WORK_DIR/gencode.v44.annotation.gff3" | wc -l)
    
    if [[ $gene_count -eq 0 ]]; then
        log "WARNING: No protein-coding genes found on chr22 in annotation file"
        log "This might indicate a format issue"
        return 1
    fi
    
    log "Setup test passed!"
    log "Found $gene_count protein-coding genes on chr22"
    return 0
}

# Main function
main() {
    log "Starting 1000 Genomes quick start setup"
    
    # Check dependencies
    check_dependencies
    
    # Set up analysis environment
    setup_1000g_analysis
    
    # Test the setup
    if test_setup; then
        log "Setup completed successfully!"
        log "Ready to run trio analysis on chromosome 22"
        
        echo
        echo "================================================================"
        echo "Setup complete! To run the analysis:"
        echo "  cd $WORK_DIR"
        echo "  ./run_chr22_analysis.sh"
        echo "================================================================"
        echo
    else
        log "ERROR: Setup test failed"
        exit 1
    fi
}

# Show usage if requested
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    cat << EOF
Quick Start Setup for 1000 Genomes Trio Analysis

This script downloads and sets up everything needed to run a trio analysis
on chromosome 22 using 1000 Genomes Project data.

USAGE:
    $0

WHAT IT DOES:
1. Downloads GRCh38 chromosome 22 reference genome
2. Downloads GENCODE v44 gene annotations  
3. Creates configuration file for CEU trio (NA12877, NA12878, NA12882)
4. Sets up convenient run script

REQUIREMENTS:
- Internet connection for downloads (~100MB total)
- Tools: wget, samtools, bcftools, bedtools
- ~200MB free disk space

OUTPUT:
- Working directory: 1000g_quickstart/
- Ready-to-run configuration for chromosome 22 analysis

EOF
    exit 0
fi

# Run main function
main "$@"