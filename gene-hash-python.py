#!/usr/bin/env python3
"""
Gene Sequence SHA-256 Hasher for Trio Analysis

This script generates SHA-256 hashes for gene sequences extracted from
consensus genome sequences. It handles cleaning sequences (removing gaps,
Ns) and provides detailed output with sequence statistics.

Usage:
    python3 hash_genes.py input.fasta output.csv sample_name haplotype_name

Example:
    python3 hash_genes.py father_hap1_genes.fa father_hap1_hashes.csv father hap1
"""

import hashlib
import sys
import csv
import argparse
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def clean_sequence(sequence):
    """
    Clean DNA sequence by removing gaps, Ns, and converting to uppercase
    
    Args:
        sequence (str): Raw DNA sequence
        
    Returns:
        str: Cleaned sequence
    """
    # Remove gaps (-), Ns, and other non-ATCG characters, convert to uppercase
    cleaned = ''.join(c.upper() for c in sequence if c.upper() in 'ATCG')
    return cleaned

def parse_fasta_simple(fasta_file):
    """
    Simple FASTA parser that yields (header, sequence) tuples
    
    Args:
        fasta_file (str): Path to FASTA file
        
    Yields:
        tuple: (gene_id, sequence)
    """
    gene_id = None
    sequence = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Output previous sequence if exists
                if gene_id and sequence:
                    yield gene_id, ''.join(sequence)
                
                # Start new sequence
                gene_id = line[1:]  # Remove '>'
                sequence = []
            else:
                sequence.append(line)
        
        # Output last sequence
        if gene_id and sequence:
            yield gene_id, ''.join(sequence)

def generate_sha256_hash(sequence):
    """
    Generate SHA-256 hash of DNA sequence
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        str: SHA-256 hash in hexadecimal
    """
    return hashlib.sha256(sequence.encode('utf-8')).hexdigest()

def hash_gene_sequences(fasta_file, output_file, sample_name, haplotype, 
                       min_length=0, verbose=False):
    """
    Generate SHA-256 hashes for all gene sequences in a FASTA file
    
    Args:
        fasta_file (str): Input FASTA file path
        output_file (str): Output CSV file path  
        sample_name (str): Sample identifier
        haplotype (str): Haplotype identifier
        min_length (int): Minimum sequence length to include
        verbose (bool): Print detailed progress
        
    Returns:
        list: List of dictionaries with hash results
    """
    results = []
    sequences_processed = 0
    sequences_skipped = 0
    
    logger.info(f"Processing FASTA file: {fasta_file}")
    logger.info(f"Sample: {sample_name}, Haplotype: {haplotype}")
    
    try:
        for gene_id, raw_sequence in parse_fasta_simple(fasta_file):
            sequences_processed += 1
            
            # Clean sequence
            clean_seq = clean_sequence(raw_sequence)
            
            # Skip if too short
            if len(clean_seq) < min_length:
                sequences_skipped += 1
                if verbose:
                    logger.warning(f"Skipping {gene_id}: length {len(clean_seq)} < {min_length}")
                continue
            
            # Generate hash
            sha256_hash = generate_sha256_hash(clean_seq)
            
            # Store result
            result = {
                'sample': sample_name,
                'haplotype': haplotype,
                'gene_id': gene_id,
                'raw_length': len(raw_sequence),
                'clean_length': len(clean_seq),
                'sha256_hash': sha256_hash
            }
            results.append(result)
            
            if verbose and sequences_processed % 100 == 0:
                logger.info(f"Processed {sequences_processed} sequences...")
    
    except FileNotFoundError:
        logger.error(f"FASTA file not found: {fasta_file}")
        return []
    except Exception as e:
        logger.error(f"Error processing FASTA file: {e}")
        return []
    
    # Write results to CSV
    try:
        with open(output_file, 'w', newline='') as csvfile:
            fieldnames = ['sample', 'haplotype', 'gene_id', 'raw_length', 
                         'clean_length', 'sha256_hash']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        
        logger.info(f"Results written to: {output_file}")
        logger.info(f"Total sequences processed: {sequences_processed}")
        logger.info(f"Sequences included: {len(results)}")
        logger.info(f"Sequences skipped: {sequences_skipped}")
        
    except Exception as e:
        logger.error(f"Error writing output file: {e}")
        return []
    
    return results

def validate_inputs(fasta_file, output_file, sample_name, haplotype):
    """
    Validate input parameters
    
    Args:
        fasta_file (str): Input FASTA file path
        output_file (str): Output CSV file path
        sample_name (str): Sample identifier  
        haplotype (str): Haplotype identifier
        
    Returns:
        bool: True if inputs are valid
    """
    # Check if FASTA file exists
    if not Path(fasta_file).exists():
        logger.error(f"FASTA file does not exist: {fasta_file}")
        return False
    
    # Check if FASTA file is readable
    try:
        with open(fasta_file, 'r') as f:
            first_line = f.readline()
            if not first_line.startswith('>'):
                logger.error(f"File does not appear to be a valid FASTA: {fasta_file}")
                return False
    except Exception as e:
        logger.error(f"Cannot read FASTA file: {e}")
        return False
    
    # Check output directory is writable
    output_path = Path(output_file)
    try:
        output_path.parent.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Cannot create output directory: {e}")
        return False
    
    # Validate sample and haplotype names
    if not sample_name or not sample_name.strip():
        logger.error("Sample name cannot be empty")
        return False
    
    if not haplotype or not haplotype.strip():
        logger.error("Haplotype name cannot be empty") 
        return False
    
    return True

def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(
        description='Generate SHA-256 hashes for gene sequences from FASTA file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('fasta_file', help='Input FASTA file with gene sequences')
    parser.add_argument('output_file', help='Output CSV file for hash results')
    parser.add_argument('sample_name', help='Sample identifier (e.g., father, mother, child)')
    parser.add_argument('haplotype', help='Haplotype identifier (e.g., hap1, hap2)')
    
    parser.add_argument('--min-length', type=int, default=0,
                       help='Minimum sequence length to include (default: 0)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not validate_inputs(args.fasta_file, args.output_file, 
                          args.sample_name, args.haplotype):
        sys.exit(1)
    
    # Set logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Process sequences
    results = hash_gene_sequences(
        args.fasta_file,
        args.output_file, 
        args.sample_name,
        args.haplotype,
        args.min_length,
        args.verbose
    )
    
    if results:
        logger.info("Gene hashing completed successfully!")
        
        # Print summary statistics
        lengths = [r['clean_length'] for r in results]
        if lengths:
            logger.info(f"Sequence length stats - Min: {min(lengths)}, "
                       f"Max: {max(lengths)}, Mean: {sum(lengths)/len(lengths):.1f}")
    else:
        logger.error("Gene hashing failed!")
        sys.exit(1)

if __name__ == "__main__":
    main()