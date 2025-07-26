#!/bin/bash

# Trio Genome Consensus and Gene Hashing Pipeline
# Generates full genome consensus sequences from trio VCF and hashes gene sequences

# Configuration
REF="reference.fa"
VCF="trio_phased.vcf.gz"
GFF="annotation.gff3"
OUTPUT_DIR="trio_consensus_output"

# Create output directory
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "=== Trio Genome Consensus Pipeline Started ==="
echo "Reference: $REF"
echo "VCF: $VCF"
echo "Annotation: $GFF"
echo "Output directory: $OUTPUT_DIR"
echo

# Step 1: Generate consensus sequences for each sample and haplotype
echo "Step 1: Generating consensus sequences..."
for sample in father mother child; do
    echo "Processing $sample..."
    
    # Extract sample-specific VCF
    bcftools view -s $sample ../$VCF -O z -o ${sample}.vcf.gz
    bcftools index ${sample}.vcf.gz
    
    # Generate consensus for each haplotype
    echo "  Generating haplotype 1..."
    bcftools consensus -f ../$REF -H 1 ${sample}.vcf.gz > ${sample}_hap1.fa
    
    echo "  Generating haplotype 2..."
    bcftools consensus -f ../$REF -H 2 ${sample}.vcf.gz > ${sample}_hap2.fa
    
    echo "  $sample consensus sequences complete"
done
echo "Consensus generation complete!"
echo

# Step 2: Extract gene coordinates from annotation
echo "Step 2: Extracting gene coordinates..."
awk '$3=="gene"' ../$GFF | \
awk '{
    # Extract gene ID from attributes (assuming standard GFF format)
    match($9, /ID=([^;]+)/, id)
    gene_id = id[1]
    if (gene_id == "") {
        match($9, /gene_id[= ]"?([^";]+)"?/, id)
        gene_id = id[1]
    }
    if (gene_id == "") gene_id = "gene_" NR
    
    print $1 "\t" $4-1 "\t" $5 "\t" gene_id
}' > genes.bed

gene_count=$(wc -l < genes.bed)
echo "Extracted $gene_count genes"
echo

# Step 3: Extract gene sequences from each haplotype
echo "Step 3: Extracting gene sequences..."
for sample in father mother child; do
    for hap in hap1 hap2; do
        echo "  Extracting genes from ${sample}_${hap}..."
        bedtools getfasta -fi ${sample}_${hap}.fa -bed genes.bed -name > ${sample}_${hap}_genes.fa
    done
done
echo "Gene extraction complete!"
echo

# Step 4: Generate SHA-256 hashes for all gene sequences
echo "Step 4: Generating SHA-256 hashes..."

# Create hash function
create_gene_hashes() {
    local fasta_file=$1
    local sample=$2
    local haplotype=$3
    
    awk '
    /^>/ {
        if (gene_id && sequence) {
            # Remove gaps and Ns, convert to uppercase
            gsub(/[-N]/, "", sequence)
            sequence = toupper(sequence)
            printf "%s,%s,%s,%d,", "'$sample'", "'$haplotype'", gene_id, length(sequence)
            # Print sequence for hashing
            print sequence
        }
        gene_id = substr($0, 2)  # Remove >
        sequence = ""
        next
    }
    {sequence = sequence $0}
    END {
        if (gene_id && sequence) {
            gsub(/[-N]/, "", sequence)
            sequence = toupper(sequence)
            printf "%s,%s,%s,%d,", "'$sample'", "'$haplotype'", gene_id, length(sequence)
            print sequence
        }
    }' $fasta_file | \
    while IFS=',' read sample hap gene_id length sequence; do
        if [ ! -z "$sequence" ]; then
            hash=$(echo -n "$sequence" | sha256sum | cut -d' ' -f1)
            echo "$sample,$hap,$gene_id,$length,$hash"
        fi
    done
}

# Generate header
echo "sample,haplotype,gene_id,sequence_length,sha256_hash" > trio_gene_hashes.csv

# Process all samples and haplotypes
total_hashes=0
for sample in father mother child; do
    for hap in hap1 hap2; do
        echo "  Hashing ${sample}_${hap} genes..."
        gene_file="${sample}_${hap}_genes.fa"
        if [ -f "$gene_file" ]; then
            hashes=$(create_gene_hashes $gene_file $sample $hap | tee -a trio_gene_hashes.csv | wc -l)
            total_hashes=$((total_hashes + hashes))
            echo "    Generated $hashes hashes"
        else
            echo "    Warning: $gene_file not found"
        fi
    done
done

echo "Hashing complete! Total hashes generated: $total_hashes"
echo

# Step 5: Generate summary statistics
echo "Step 5: Generating summary statistics..."

# Count unique hashes per gene across all samples/haplotypes
echo "Generating hash diversity analysis..."
python3 -c "
import pandas as pd
import sys

try:
    df = pd.read_csv('trio_gene_hashes.csv')
    
    # Summary by sample
    print('=== Summary by Sample ===')
    summary = df.groupby(['sample', 'haplotype']).agg({
        'gene_id': 'count',
        'sequence_length': ['mean', 'min', 'max'],
        'sha256_hash': 'nunique'
    }).round(2)
    print(summary)
    print()
    
    # Hash diversity per gene
    print('=== Hash Diversity per Gene (top 10) ===')
    hash_diversity = df.groupby('gene_id')['sha256_hash'].nunique().sort_values(ascending=False).head(10)
    print(hash_diversity)
    print()
    
    # Identical sequences across samples
    print('=== Genes with Identical Sequences Across All Haplotypes ===')
    identical = df.groupby('gene_id')['sha256_hash'].nunique()
    identical_genes = identical[identical == 1].index.tolist()
    print(f'Number of genes with identical sequences: {len(identical_genes)}')
    if len(identical_genes) <= 10:
        print('Genes:', ', '.join(identical_genes))
    else:
        print('First 10 genes:', ', '.join(identical_genes[:10]))
    
    # Save summary
    with open('trio_summary.txt', 'w') as f:
        f.write('Trio Gene Hashing Summary\\n')
        f.write('=' * 30 + '\\n')
        f.write(f'Total genes processed: {df[\"gene_id\"].nunique()}\\n')
        f.write(f'Total hashes generated: {len(df)}\\n')
        f.write(f'Unique hash values: {df[\"sha256_hash\"].nunique()}\\n')
        f.write(f'Genes with identical sequences: {len(identical_genes)}\\n')
        f.write(f'Average sequence length: {df[\"sequence_length\"].mean():.2f}\\n')
        
except ImportError:
    print('pandas not available, skipping detailed analysis')
    print('Basic summary:')
    with open('trio_gene_hashes.csv', 'r') as f:
        lines = sum(1 for line in f) - 1  # subtract header
        print(f'Total hash records: {lines}')
except Exception as e:
    print(f'Error in analysis: {e}')
"

# Step 6: Cleanup intermediate files (optional)
read -p "Remove intermediate files? (y/N): " cleanup
if [[ $cleanup =~ ^[Yy]$ ]]; then
    echo "Cleaning up intermediate files..."
    rm -f *.vcf.gz *.vcf.gz.csi
    rm -f *_hap1.fa *_hap2.fa
    rm -f *_genes.fa
    rm -f genes.bed
    echo "Cleanup complete"
fi

echo
echo "=== Pipeline Complete ==="
echo "Results saved in: $OUTPUT_DIR/"
echo "Main output file: trio_gene_hashes.csv"
echo "Summary: trio_summary.txt"
echo
echo "Hash file contains $(wc -l < trio_gene_hashes.csv) total records"
echo "Pipeline execution finished: $(date)"