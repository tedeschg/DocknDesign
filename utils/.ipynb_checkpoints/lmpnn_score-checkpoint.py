from typing import List, Dict, Tuple
import argparse
import sys

# ===============================
#   BLOSUM62 Matrix (compressed)
# ===============================

BLOSUM62 = {
    'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 0, 'W': -3, 'Y': -2, 'V': 0},
    'R': {'A': -1, 'R': 5, 'N': 0, 'D': -2, 'C': -3, 'Q': 1, 'E': 0, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 2, 'M': -1, 'F': -3, 'P': -2, 'S': -1, 'T': -1, 'W': -3, 'Y': -2, 'V': -3},
    'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -3, 'L': -3, 'K': 0, 'M': -2, 'F': -3, 'P': -2, 'S': 1, 'T': 0, 'W': -4, 'Y': -2, 'V': -3},
    'D': {'A': -2, 'R': -2, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': -1, 'H': -1, 'I': -3, 'L': -4, 'K': -1, 'M': -3, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -4, 'Y': -3, 'V': -3},
    'C': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -4, 'G': -3, 'H': -3, 'I': -1, 'L': -1, 'K': -3, 'M': -1, 'F': -2, 'P': -3, 'S': -1, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
    'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -2, 'H': 0, 'I': -3, 'L': -2, 'K': 1, 'M': 0, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -2, 'Y': -1, 'V': -2},
    'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -4, 'Q': 2, 'E': 5, 'G': -2, 'H': 0, 'I': -3, 'L': -3, 'K': 1, 'M': -2, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'G': {'A': 0, 'R': -2, 'N': 0, 'D': -1, 'C': -3, 'Q': -2, 'E': -2, 'G': 6, 'H': -2, 'I': -4, 'L': -4, 'K': -2, 'M': -3, 'F': -3, 'P': -2, 'S': 0, 'T': -2, 'W': -2, 'Y': -3, 'V': -3},
    'H': {'A': -2, 'R': 0, 'N': 1, 'D': -1, 'C': -3, 'Q': 0, 'E': 0, 'G': -2, 'H': 8, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -1, 'P': -2, 'S': -1, 'T': -2, 'W': -2, 'Y': 2, 'V': -3},
    'I': {'A': -1, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -3, 'E': -3, 'G': -4, 'H': -3, 'I': 4, 'L': 2, 'K': -3, 'M': 1, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -3, 'Y': -1, 'V': 3},
    'L': {'A': -1, 'R': -2, 'N': -3, 'D': -4, 'C': -1, 'Q': -2, 'E': -3, 'G': -4, 'H': -3, 'I': 2, 'L': 4, 'K': -2, 'M': 2, 'F': 0, 'P': -3, 'S': -2, 'T': -1, 'W': -2, 'Y': -1, 'V': 1},
    'K': {'A': -1, 'R': 2, 'N': 0, 'D': -1, 'C': -3, 'Q': 1, 'E': 1, 'G': -2, 'H': -1, 'I': -3, 'L': -2, 'K': 5, 'M': -1, 'F': -3, 'P': -1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -2},
    'M': {'A': -1, 'R': -1, 'N': -2, 'D': -3, 'C': -1, 'Q': 0, 'E': -2, 'G': -3, 'H': -2, 'I': 1, 'L': 2, 'K': -1, 'M': 5, 'F': 0, 'P': -2, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1},
    'F': {'A': -2, 'R': -3, 'N': -3, 'D': -3, 'C': -2, 'Q': -3, 'E': -3, 'G': -3, 'H': -1, 'I': 0, 'L': 0, 'K': -3, 'M': 0, 'F': 6, 'P': -4, 'S': -2, 'T': -2, 'W': 1, 'Y': 3, 'V': -1},
    'P': {'A': -1, 'R': -2, 'N': -2, 'D': -1, 'C': -3, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -3, 'L': -3, 'K': -1, 'M': -2, 'F': -4, 'P': 7, 'S': -1, 'T': -1, 'W': -4, 'Y': -3, 'V': -2},
    'S': {'A': 1, 'R': -1, 'N': 1, 'D': 0, 'C': -1, 'Q': 0, 'E': 0, 'G': 0, 'H': -1, 'I': -2, 'L': -2, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -3, 'Y': -2, 'V': -2},
    'T': {'A': 0, 'R': -1, 'N': 0, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0},
    'W': {'A': -3, 'R': -3, 'N': -4, 'D': -4, 'C': -2, 'Q': -2, 'E': -3, 'G': -2, 'H': -2, 'I': -3, 'L': -2, 'K': -3, 'M': -1, 'F': 1, 'P': -4, 'S': -3, 'T': -2, 'W': 11, 'Y': 2, 'V': -3},
    'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -3, 'C': -2, 'Q': -1, 'E': -2, 'G': -3, 'H': 2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 3, 'P': -3, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1},
    'V': {'A': 0, 'R': -3, 'N': -3, 'D': -3, 'C': -1, 'Q': -2, 'E': -2, 'G': -3, 'H': -3, 'I': 3, 'L': 1, 'K': -2, 'M': 1, 'F': -1, 'P': -2, 'S': -2, 'T': 0, 'W': -3, 'Y': -1, 'V': 4}
}

def blosum_score(aa1: str, aa2: str) -> int:
    """Get BLOSUM62 score for amino acid pair."""
    if aa1 not in BLOSUM62 or aa2 not in BLOSUM62[aa1]:
        return 0
    return BLOSUM62[aa1][aa2]


# ===============================
#   Parse PDB and extract coordinates
# ===============================

def parse_pdb_residue_map(
    pdb_file: str,
    target_chains: List[str] = None
) -> Dict[Tuple[str, int], int]:
    """Map (chain, resid) to sequence position (1-based)."""
    residues: List[Tuple[str, int]] = []
    seen = set()

    with open(pdb_file, "r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            chain = line[21].strip() or "A"
            if target_chains and chain not in target_chains:
                continue
            resid = int(line[22:26].strip())
            res_key = (chain, resid)
            if res_key not in seen:
                residues.append(res_key)
                seen.add(res_key)

    return {res: idx for idx, res in enumerate(residues, 1)}


# ===============================
#   Read pocket residues file
# ===============================

def read_pocket_residues(pocket_file: str) -> List[Tuple[str, int]]:
    """Read pocket_residues.txt: 'A123 B45 ...' -> [(chain, resid), ...]"""
    with open(pocket_file, "r") as f:
        content = f.read().strip()

    residues: List[Tuple[str, int]] = []
    for token in content.split():
        if len(token) < 2:
            continue
        chain = token[0]
        try:
            resid = int(token[1:])
            residues.append((chain, resid))
        except ValueError:
            print(f"Warning: Could not parse '{token}'")
    return residues


# ===============================
#   Map pocket to sequence positions
# ===============================

def map_pocket_to_sequence(
    pocket_residues: List[Tuple[str, int]],
    residue_map: Dict[Tuple[str, int], int]
) -> Tuple[List[int], List[Tuple[str, int]]]:
    """Map pocket residues to sequence positions."""
    positions: List[int] = []
    mapped_residues: List[Tuple[str, int]] = []
    not_found: List[str] = []

    for chain, resid in pocket_residues:
        res_key = (chain, resid)
        if res_key in residue_map:
            positions.append(residue_map[res_key])
            mapped_residues.append(res_key)
        else:
            not_found.append(f"{chain}{resid}")

    if not_found:
        print(f"Warning: {len(not_found)} pocket residues not in PDB:")
        print(f"  {', '.join(not_found[:10])}" + (" ..." if len(not_found) > 10 else ""))

    return positions, mapped_residues


# ===============================
#   Read FASTA
# ===============================

def read_complex_fasta(path: str) -> Tuple[str, str]:
    """Read FASTA with WT (first) and DESIGN (second) sequences."""
    sequences: List[str] = []
    current: List[str] = []

    with open(path, "r") as f:
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
        raise ValueError("FASTA must contain WT and design sequences.")

    return sequences[0].upper(), sequences[1].upper()


# ===============================
#   Compute pocket scores
# ===============================

def compute_pocket_scores(
    wt_seq: str,
    des_seq: str,
    pocket_positions: List[int],
    mapped_residues: List[Tuple[str, int]],
) -> Tuple[int, float, float, List[Dict]]:
    """
    Compute pocket mutation metrics.
    
    Returns:
        - n_mut: number of mutations
        - blosum_sum: sum of BLOSUM62 scores
        - final_score: n_mut - blosum_sum (higher = worse)
        - mutations: detailed list
    """
    n_mut = 0
    blosum_sum = 0.0
    mutations: List[Dict] = []

    for pos, res_key in zip(pocket_positions, mapped_residues):
        idx = pos - 1
        if idx < 0 or idx >= len(wt_seq) or idx >= len(des_seq):
            continue

        aa_wt = wt_seq[idx]
        aa_des = des_seq[idx]

        if aa_wt != aa_des:
            score = blosum_score(aa_wt, aa_des)
            n_mut += 1
            blosum_sum += score
            mutations.append({
                "position": pos,
                "residue_key": f"{res_key[0]}{res_key[1]}",
                "wt": aa_wt,
                "design": aa_des,
                "blosum": score
            })

    final_score = n_mut - blosum_sum
    return n_mut, blosum_sum, final_score, mutations


# ===============================
#   Main
# ===============================

def main():
    parser = argparse.ArgumentParser(
        description="Compute pocket mutation scores with BLOSUM62",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -f complex.fa -p protein.pdb -r pocket_residues.txt
  %(prog)s -f complex.fa -p protein.pdb -r pocket_residues.txt --show-mutations

Final Score:
  FINAL_SCORE = mutation_count - BLOSUM62_sum
  
  Components:
  - More mutations → higher score (worse)
  - Dissimilar mutations (negative BLOSUM) → higher score (worse)
  - Similar mutations (positive BLOSUM) → lower score (better)
  
  Higher score = worse overall (more and/or dissimilar mutations)
        """
    )

    parser.add_argument("-f", "--fasta", required=True,
                        help="FASTA with WT and design sequences")
    parser.add_argument("-p", "--pdb", required=True,
                        help="PDB file for fpocket analysis")
    parser.add_argument("-r", "--residues", default="pocket_residues.txt",
                        help="Pocket residues file (default: pocket_residues.txt)")
    parser.add_argument("--chains", nargs="*", default=None,
                        help="Chains to consider (default: all)")
    parser.add_argument("--show-mutations", action="store_true",
                        help="Show detailed mutation information")

    args = parser.parse_args()

    try:
        # Read sequences
        print("Reading sequences:", args.fasta)
        wt_seq, des_seq = read_complex_fasta(args.fasta)
        print(f"WT length    : {len(wt_seq)}")
        print(f"Design length: {len(des_seq)}")

        # Parse PDB
        print(f"\nParsing PDB: {args.pdb}")
        residue_map = parse_pdb_residue_map(args.pdb, args.chains)
        print(f"PDB residues : {len(residue_map)}")

        # Read pocket
        print(f"\nReading pocket: {args.residues}")
        pocket_residues = read_pocket_residues(args.residues)
        print(f"Pocket residues: {len(pocket_residues)}")

        # Map pocket to sequence
        print("\nMapping pocket to sequence...")
        pocket_positions, mapped_residues = map_pocket_to_sequence(pocket_residues, residue_map)
        print(f"Mapped positions: {len(pocket_positions)}")

        # Compute scores
        print("\n" + "=" * 60)
        print("POCKET MUTATION SCORES")
        print("=" * 60)

        n_mut, blosum_sum, final_score, mutations = compute_pocket_scores(
            wt_seq, des_seq, pocket_positions, mapped_residues
        )

        print(f"\nMutation count : {n_mut}")
        print(f"BLOSUM62 sum   : {blosum_sum:.1f}")
        print(f"\n{'='*60}")
        print(f"FINAL SCORE    : {final_score:.1f}")
        print(f"{'='*60}")

        # Show details if requested
        if args.show_mutations and mutations:
            print("\n" + "=" * 60)
            print("MUTATION DETAILS")
            print("=" * 60)
            print(f"\n{'Pos':<6} {'Residue':<10} {'Change':<8} {'BLOSUM':<8}")
            print("-" * 40)
            for m in mutations:
                print(f"{m['position']:<6} {m['residue_key']:<10} "
                      f"{m['wt']}->{m['design']:<6} {m['blosum']:<8}")

    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()