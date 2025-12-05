from typing import List, Dict, Tuple
import argparse
import sys

# ===============================
#   1. Hydrophobicity & Occupancy tables
# ===============================

HYDRO_KD: Dict[str, float] = {
    "A": 1.8,  "R": -4.5, "N": -3.5, "D": -3.5,
    "C": 2.5,  "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5,  "L": 3.8,  "K": -3.9,
    "M": 1.9,  "F": 2.8,  "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

SC_HEAVY_ATOMS: Dict[str, int] = {
    "G": 0, "A": 1,
    "S": 2, "C": 2,
    "P": 3, "V": 3, "T": 3,
    "L": 4, "I": 4, "M": 4,
    "D": 4, "N": 4,
    "E": 5, "Q": 5, "K": 5,
    "F": 7, "R": 7,
    "H": 6,
    "Y": 8,
    "W": 10,
}

MAX_SC = 10.0      # max side-chain length (Trp)
ALPHA  = 0.5       # weight for dist_score in final_score


# ===============================
#   2. Parse PDB and build residue map
# ===============================

def parse_pdb_residues(pdb_file: str, target_chains: List[str] = None) -> Dict[Tuple[str, int], int]:
    """
    Parse PDB file and create a mapping from (chain, resid) to sequential position (1-based).
    
    Args:
        pdb_file: Path to PDB file
        target_chains: List of chains to include (None = all chains)
    
    Returns:
        Dict mapping (chain, resid) -> sequence_position (1-based)
    """
    residues = []
    seen = set()
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                chain = line[21].strip()
                if not chain:
                    chain = 'A'
                
                # Skip if we're filtering by chains
                if target_chains and chain not in target_chains:
                    continue
                
                resid = int(line[22:26].strip())
                
                res_key = (chain, resid)
                if res_key not in seen:
                    residues.append(res_key)
                    seen.add(res_key)
    
    # Build mapping: (chain, resid) -> sequential position
    residue_map = {}
    for idx, res_key in enumerate(residues, start=1):
        residue_map[res_key] = idx
    
    return residue_map


# ===============================
#   3. Read pocket residues file
# ===============================

def read_pocket_residues(pocket_file: str) -> List[Tuple[str, int]]:
    """
    Read pocket_residues.txt file in format "A123 B45 C67..."
    
    Returns:
        List of (chain, resid) tuples
    """
    with open(pocket_file, 'r') as f:
        content = f.read().strip()
    
    residues = []
    for token in content.split():
        # Parse format like "A123" -> chain='A', resid=123
        if len(token) < 2:
            continue
        
        chain = token[0]
        try:
            resid = int(token[1:])
            residues.append((chain, resid))
        except ValueError:
            print(f"Warning: Could not parse residue '{token}'")
            continue
    
    return residues


# ===============================
#   4. Map pocket residues to sequence positions
# ===============================

def map_pocket_to_sequence(
    pocket_residues: List[Tuple[str, int]],
    residue_map: Dict[Tuple[str, int], int]
) -> List[int]:
    """
    Map pocket residues (chain, resid) to sequence positions using PDB mapping.
    
    Args:
        pocket_residues: List of (chain, resid) from pocket file
        residue_map: Mapping from (chain, resid) to sequence position
    
    Returns:
        List of sequence positions (1-based)
    """
    positions = []
    not_found = []
    
    for chain, resid in pocket_residues:
        res_key = (chain, resid)
        if res_key in residue_map:
            positions.append(residue_map[res_key])
        else:
            not_found.append(f"{chain}{resid}")
    
    if not_found:
        print(f"Warning: {len(not_found)} pocket residues not found in PDB:")
        print(f"  {', '.join(not_found[:10])}" + (" ..." if len(not_found) > 10 else ""))
    
    return sorted(positions)


# ===============================
#   5. Read complex.fa (WT + DESIGN)
# ===============================

def read_complex_fasta(path: str) -> Tuple[str, str]:
    """
    Read complex.fa and return (wt_seq, design_seq).

    Assumes:
      - first sequence  = WT
      - second sequence = designed
    """
    sequences = []
    current = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                    current = []
            else:
                current.append(line)

    if current:
        sequences.append("".join(current))

    if len(sequences) < 2:
        raise ValueError("complex.fa must contain at least two sequences (WT and design).")

    wt = sequences[0].upper()
    des = sequences[1].upper()
    return wt, des


# ===============================
#   6. Score on a single sequence
# ===============================

def compute_sequence_score(sequence: str, pocket_positions: List[int]) -> float:
    """
    Compute score = H_pocket + Occ_pocket for a given sequence.

    pocket_positions must be 1-based indices on the sequence.
    """
    H_values = []
    Occ_values = []

    N = len(sequence)

    for pos in pocket_positions:
        idx = pos - 1  # convert to 0-based
        if idx < 0 or idx >= N:
            # position outside sequence, skip
            continue

        aa = sequence[idx]
        if aa not in HYDRO_KD or aa not in SC_HEAVY_ATOMS:
            # unknown / non-standard residue, skip
            continue

        H_values.append(HYDRO_KD[aa])
        Occ_values.append(SC_HEAVY_ATOMS[aa] / MAX_SC)

    if not H_values:
        raise ValueError("No valid pocket residues found for this sequence.")

    H_pocket = sum(H_values) / len(H_values)
    Occ_pocket = sum(Occ_values) / len(Occ_values)

    score = H_pocket + Occ_pocket
    return score


# ===============================
#   7. Count mutations in the pocket
# ===============================

def count_pocket_mutations(
    wt_seq: str,
    des_seq: str,
    pocket_positions: List[int],
) -> int:
    """
    Count how many residues differ between WT and DESIGN,
    considering only the positions in pocket_positions (1-based).
    """
    n_mut = 0
    N_wt = len(wt_seq)
    N_des = len(des_seq)

    for pos in pocket_positions:
        idx = pos - 1
        if idx < 0 or idx >= N_wt or idx >= N_des:
            # outside one of the sequences
            continue

        aa_wt = wt_seq[idx]
        aa_des = des_seq[idx]

        if aa_wt != aa_des:
            n_mut += 1

    return n_mut


# ===============================
#   8. WT vs DESIGN + distance + final score
# ===============================

def compute_scores_and_final(
    wt_seq: str,
    des_seq: str,
    pocket_positions: List[int],
) -> dict:
    """
    Compute:
      - score_wt       = H_pocket + Occ_pocket on WT
      - score_design   = H_pocket + Occ_pocket on DESIGN
      - dist_score     = abs(score_design - score_wt)
      - n_mut_pocket   = number of mutations in the pocket
      - final_score    = n_mut_pocket + ALPHA * dist_score

    Returns a dict with all values.
    """
    score_wt = compute_sequence_score(wt_seq, pocket_positions)
    score_design = compute_sequence_score(des_seq, pocket_positions)
    dist_score = abs(score_design - score_wt)

    n_mut_pocket = count_pocket_mutations(wt_seq, des_seq, pocket_positions)

    final_score = n_mut_pocket + ALPHA * (100*dist_score)

    return {
        "score_wt": score_wt,
        "score_design": score_design,
        "dist_score": dist_score,
        "n_mut_pocket": n_mut_pocket,
        "final_score": final_score,
    }


# ===============================
#   9. Main
# ===============================

def main():
    parser = argparse.ArgumentParser(
        description='Compute pocket-based scores from fpocket output and sequences.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -f complex.fa -p protein.pdb -r pocket_residues.txt
  %(prog)s -f complex.fa -p protein.pdb -r pocket_residues.txt --chains A B
        """
    )
    
    parser.add_argument('-f', '--fasta',
                        required=True,
                        help='FASTA file with WT and design sequences (complex.fa)')
    
    parser.add_argument('-p', '--pdb',
                        required=True,
                        help='PDB file used for fpocket analysis')
    
    parser.add_argument('-r', '--residues',
                        default='pocket_residues.txt',
                        help='Pocket residues file from fpocket (default: pocket_residues.txt)')
    
    parser.add_argument('--chains',
                        nargs='*',
                        default=None,
                        help='Chains to consider (default: all chains)')
    
    args = parser.parse_args()
    
    try:
        # 1. Read WT and designed sequences
        print("Reading sequences from:", args.fasta)
        wt_seq, des_seq = read_complex_fasta(args.fasta)
        print(f"WT sequence length: {len(wt_seq)}")
        print(f"Design sequence length: {len(des_seq)}")
        
        # 2. Parse PDB and build residue mapping
        print("\nParsing PDB file:", args.pdb)
        residue_map = parse_pdb_residues(args.pdb, args.chains)
        print(f"Total residues in PDB: {len(residue_map)}")
        
        # 3. Read pocket residues from fpocket output
        print("\nReading pocket residues from:", args.residues)
        pocket_residues = read_pocket_residues(args.residues)
        print(f"Pocket residues found: {len(pocket_residues)}")
        
        # 4. Map pocket residues to sequence positions
        print("\nMapping pocket residues to sequence positions...")
        pocket_positions = map_pocket_to_sequence(pocket_residues, residue_map)
        print(f"Mapped positions: {len(pocket_positions)}")
        print(f"Positions: {pocket_positions}")
        
        # 5. Compute scores
        print("\n" + "="*60)
        print("COMPUTING SCORES")
        print("="*60)
        
        result = compute_scores_and_final(wt_seq, des_seq, pocket_positions)
        
        print(f"\nWT score         : {result['score_wt']:.4f}")
        print(f"DESIGN score     : {result['score_design']:.4f}")
        print(f"dist_score       : {result['dist_score']:.4f}")
        print(f"n_mut_pocket     : {result['n_mut_pocket']}")
        print(f"\nFINAL SCORE      : {result['final_score']:.4f}")
        print(f"(final_score = n_mut_pocket + {ALPHA} * dist_score)")
        
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()