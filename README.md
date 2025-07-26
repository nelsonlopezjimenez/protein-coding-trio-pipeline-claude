## Setup Instructions
```sh
# 1. Save all three scripts and make executable
chmod +x trio_consensus_pipeline.sh
chmod +x run_multi_chromosome.sh  
chmod +x quick_start_1000g.sh

# 2. Run quick setup
./quick_start_1000g.sh

# 3. Test with single chromosome
cd 1000g_quickstart
./run_chr22_analysis.sh

# 4. Or run multiple chromosomes
./run_multi_chromosome.sh -c "chr21 chr22"

# 5. Or parallel processing
./run_multi_chromosome.sh -c "chr1 chr2 chr3" --parallel -j 3
```

## Usage
1. Single chromosome test:
```bash
./run_multi_chromosome.sh -c "chr22"
```
1. Multiple chromosomes sequentially:
```bash
./run_multi_chromosome.sh -c "chr21 chr22 chrX"
```
1. Parallel processing:
```bash
./run_multi_chromosome.sh -c "chr1 chr2 chr3 chr4" --parallel -j 4
```
1. All autosomes (for production runs):
```bash
./run_multi_chromosome.sh -c "$(seq -f 'chr%.0f' 1 22 | tr '\n' ' ')" --parallel -j 6
```
Now you have all the scripts separated and ready to use! The multi-chromosome runner will automatically combine results from all chromosomes into summary files.

## What the Quick Start Creates:
```
1000g_quickstart/
├── chr22.fa                          # Reference genome (chr22)
├── chr22.fa.fai                      # Reference index
├── gencode.v44.annotation.gff3       # Gene annotations
├── pipeline_config.conf              # Configuration file
├── run_chr22_analysis.sh             # Convenient run script
├── setup_summary.txt                 # What was set up
└── results/                          # Results will go here
    └── chr22/
        ├── trio_chr22_protein_coding_hashes.csv
        ├── trio_chr22_summary.txt
        └── [other output files]
```

You're right! The quick start script was embedded within the second artifact. Let me create it as a separate, standalone file:Perfect! Now you have the `quick_start_1000g.sh` script as a standalone file. Here's how to use it:

## 🚀 **Setup Instructions**

**1. Save the Scripts:**
```bash
# Save all three scripts:
# - trio_consensus_pipeline.sh (from first artifact)
# - run_multi_chromosome.sh (from second artifact, if needed)
# - quick_start_1000g.sh (from this artifact)

# Make them executable
chmod +x trio_consensus_pipeline.sh
chmod +x quick_start_1000g.sh
```

**2. Run the Quick Start Setup:**
```bash
# This will download reference files and create config
./quick_start_1000g.sh
```

**3. Run the Analysis:**
```bash
# Go to the setup directory
cd 1000g_quickstart

# Run the convenient script
./run_chr22_analysis.sh

# OR run manually
../trio_consensus_pipeline.sh --config pipeline_config.conf -c chr22 -d
```

## 📁 **What the Quick Start Creates:**

```
1000g_quickstart/
├── chr22.fa                          # Reference genome (chr22)
├── chr22.fa.fai                      # Reference index
├── gencode.v44.annotation.gff3       # Gene annotations
├── pipeline_config.conf              # Configuration file
├── run_chr22_analysis.sh             # Convenient run script
├── setup_summary.txt                 # What was set up
└── results/                          # Results will go here
    └── chr22/
        ├── trio_chr22_protein_coding_hashes.csv
        ├── trio_chr22_summary.txt
        └── [other output files]
```

## 💡 **What Happens During Analysis:**

1. **Downloads** the 1000G chr22 VCF (~47MB)
2. **Extracts** trio samples (NA12877, NA12878, NA12882)  
3. **Finds** ~500 protein-coding genes on chr22
4. **Generates** consensus sequences for each haplotype
5. **Creates** SHA-256 hashes for genetic analysis
6. **Produces** summary statistics

The quick start script makes it much easier to get going with a test run before scaling up!

## 📁 Complete File List

trio_consensus_pipeline.sh - Main pipeline (from first artifact)
run_multi_chromosome.sh - Multi-chromosome runner (this artifact)
quick_start_1000g.sh - Setup script (from previous artifact)

## 📁 Complete File Structure You'll Have
```
your_project/
├── trio_consensus_pipeline.sh    # Main pipeline script
├── hash_genes.py                 # Python hashing utility
├── reference.fa                  # Your reference genome
├── trio_phased.vcf.gz            # Your phased VCF
├── annotation.gff3               # Gene annotations
└── trio_consensus_output/        # Results directory
    ├── trio_gene_hashes.csv      # Final hash results
    ├── trio_summary.txt          # Summary statistics
    └── [intermediate files]
```

## 📊 Expected Output Format

csvsample,haplotype,gene_id,sequence_length,sha256_hash
father,hap1,GENE001,1500,a1b2c3d4e5f6...
father,hap2,GENE001,1503,b2c3d4e5f6a1...
mother,hap1,GENE001,1500,c3d4e5f6a1b2...
mother,hap2,GENE001,1497,d4e5f6a1b2c3...
child,hap1,GENE001,1500,a1b2c3d4e5f6...
child,hap2,GENE001,1503,b2c3d4e5f6a1...