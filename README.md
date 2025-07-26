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