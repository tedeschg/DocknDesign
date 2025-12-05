#!/usr/bin/env python3
"""
Calculate mutation fraction from LigandMPNN output.
Compares original and mutated protein sequences.
"""

import argparse
import sys
import os

def parse_fasta(filename):
    """Parse FASTA file and return sequences with their headers."""
    if not os.path.exists(filename):
        print(f"Error: File not found: {filename}")
        sys.exit(1)
    
    sequences = []
    headers = []
    current_seq = []
    current_header = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq and current_header:
                    sequences.append(''.join(current_seq))
                    headers.append(current_header)
                current_header = line[1:]  # Remove '>'
                current_seq = []
            elif line:  # Skip empty lines
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_seq and current_header:
            sequences.append(''.join(current_seq))
            headers.append(current_header)
    
    return headers, sequences

def calculate_mutations(seq1, seq2):
    """
    Calculate mutation statistics between two sequences.
    
    Args:
        seq1: Original sequence
        seq2: Mutated sequence
    
    Returns:
        Dictionary with mutation statistics
    """
    if len(seq1) != len(seq2):
        print(f"Error: Sequences have different lengths: {len(seq1)} vs {len(seq2)}")
        sys.exit(1)
    
    # Count mutations and total residues (excluding 'X' positions)
    mutations = 0
    total_residues = 0
    
    for i, (aa1, aa2) in enumerate(zip(seq1, seq2)):
        # Skip positions with X (ligand/masked residues)
        if aa1 == 'X' or aa2 == 'X':
            continue
        
        total_residues += 1
        if aa1 != aa2:
            mutations += 1
    
    if total_residues == 0:
        print("Error: No valid residues found (all positions are 'X')")
        sys.exit(1)
    
    # Calculate mutation fraction: mutations / total residues
    mut_fraction = mutations / total_residues
    
    return {
        'total_residues': total_residues,
        'mutations': mutations,
        'mut_fraction': mut_fraction,
        'sequence_length': len(seq1)
    }

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Calculate mutation fraction from LigandMPNN output',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
  python lmpnn_score.py -f output.fa
  python lmpnn_score.py --fasta output.fa --verbose
  
The script compares the first two sequences in the FASTA file.
Positions marked with 'X' (ligand/masked residues) are excluded from analysis.
Mutation fraction = mutations / total_residues
        """
    )
    parser.add_argument('-f', '--fasta', required=True,
                        help='Input FASTA file with original and mutated sequences')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed mutation list')
    
    args = parser.parse_args()
    
    # Parse sequences
    headers, sequences = parse_fasta(args.fasta)
    
    if len(sequences) < 2:
        print(f"Error: Expected at least 2 sequences in the FASTA file, found {len(sequences)}")
        sys.exit(1)
    
    original_seq = sequences[0]
    mutated_seq = sequences[1]
    
    # Calculate mutations
    stats = calculate_mutations(original_seq, mutated_seq)
    
    # Print results
    print("=" * 60)
    print("LigandMPNN Mutation Analysis")
    print("=" * 60)
    print(f"Input file: {args.fasta}")
    print(f"Sequence 1: {headers[0][:50]}...")
    print(f"Sequence 2: {headers[1][:50]}...")
    print(f"Sequence length: {stats['sequence_length']} residues")
    print(f"\nMutation Statistics:")
    print(f"  Total residues (non-X):    {stats['total_residues']}")
    print(f"  Actual mutations:          {stats['mutations']}")
    print(f"  Mutation fraction (m/n):   {stats['mut_fraction']:.4f}")
    print(f"  ({stats['mut_fraction']*100:.2f}%)")
    print("=" * 60)
    
    # Show mutations if verbose
    if args.verbose:
        print("\nMutations found:")
        mutation_count = 0
        for i, (aa1, aa2) in enumerate(zip(original_seq, mutated_seq), 1):
            if aa1 != 'X' and aa2 != 'X' and aa1 != aa2:
                print(f"  Position {i}: {aa1} → {aa2}")
                mutation_count += 1
        if mutation_count == 0:
            print("  No mutations found")

if __name__ == "__main__":
    main()
