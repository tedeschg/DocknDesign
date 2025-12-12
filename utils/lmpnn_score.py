from typing import List, Dict, Tuple
import argparse
import sys
import numpy as np

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

# Density penalty parameters
DENSITY_THRESHOLD = 8.0  # Angstroms - mutations within this distance are considered "clustered"
DENSITY_WEIGHT = 1.0     # Weight factor for density penalty


# ===============================
#   2. Parse PDB and extract coordinates
# ===============================

def parse_pdb_with_coords(pdb_file: str, target_chains: List[str] = None) -> Tuple[Dict[Tuple[str, int], int], Dict[Tuple[str, int], np.ndarray]]:
    """
    Parse PDB file and create:
    1. Mapping from (chain, resid) to sequential position (1-based)
    2. Mapping from (chain, resid) to CA coordinates (x, y, z)
    
    Args:
        pdb_file: Path to PDB file
        target_chains: List of chains to include (None = all chains)
    
    Returns:
        Tuple of:
        - Dict mapping (chain, resid) -> sequence_position (1-based)
        - Dict mapping (chain, resid) -> np.array([x, y, z])
    """
    residues = []
    coords = {}
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
                atom_name = line[12:16].strip()
                
                res_key = (chain, resid)
                
                # Store CA coordinates
                if atom_name == 'CA':
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords[res_key] = np.array([x, y, z])
                
                if res_key not in seen:
                    residues.append(res_key)
                    seen.add(res_key)
    
    # Build mapping: (chain, resid) -> sequential position
    residue_map = {}
    for idx, res_key in enumerate(residues, start=1):
        residue_map[res_key] = idx
    
    return residue_map, coords


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
#   5. Map pocket residues to sequence positions
# ===============================

def map_pocket_to_sequence(
    pocket_residues: List[Tuple[str, int]],
    residue_map: Dict[Tuple[str, int], int]
) -> Tuple[List[int], List[Tuple[str, int]]]:
    """
    Map pocket residues (chain, resid) to sequence positions using PDB mapping.
    
    Args:
        pocket_residues: List of (chain, resid) from pocket file
        residue_map: Mapping from (chain, resid) to sequence position
    
    Returns:
        Tuple of:
        - List of sequence positions (1-based)
        - List of (chain, resid) that were successfully mapped
    """
    positions = []
    mapped_residues = []
    not_found = []
    
    for chain, resid in pocket_residues:
        res_key = (chain, resid)
        if res_key in residue_map:
            positions.append(residue_map[res_key])
            mapped_residues.append(res_key)
        else:
            not_found.append(f"{chain}{resid}")
    
    if not_found:
        print(f"Warning: {len(not_found)} pocket residues not found in PDB:")
        print(f"  {', '.join(not_found[:10])}" + (" ..." if len(not_found) > 10 else ""))
    
    return positions, mapped_residues


# ===============================
#   6. Read complex.fa (WT + DESIGN)
# ===============================

def read_complex_fasta(path: str) -> Tuple[str, str]:
    """
    Read complex.fa and return (wt_seq, design_seq).
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
#   7. Score on a single sequence
# ===============================

def compute_sequence_score(sequence: str, pocket_positions: List[int]) -> float:
    """
    Compute score = H_pocket + Occ_pocket for a given sequence.
    """
    H_values = []
    Occ_values = []

    N = len(sequence)

    for pos in pocket_positions:
        idx = pos - 1  # convert to 0-based
        if idx < 0 or idx >= N:
            continue

        aa = sequence[idx]
        if aa not in HYDRO_KD or aa not in SC_HEAVY_ATOMS:
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
#   8. Compute mutation density metrics
# ===============================

def compute_mutation_density(
    mutation_coords: List[np.ndarray],
    threshold: float = DENSITY_THRESHOLD
) -> Tuple[float, float, float]:
    """
    Compute density metrics for mutations.
    
    Args:
        mutation_coords: List of 3D coordinates for mutated residues
        threshold: Distance threshold for considering mutations "clustered"
    
    Returns:
        Tuple of:
        - avg_pairwise_dist: Average pairwise distance between mutations
        - cluster_score: Fraction of mutation pairs within threshold (0=dispersed, 1=clustered)
        - density_weight_factor: Weight factor (higher = more clustered = higher weight)
    """
    n_mut = len(mutation_coords)
    
    if n_mut <= 1:
        # Single or no mutations: full weight
        return 0.0, 1.0, 1.0
    
    # Compute all pairwise distances
    pairwise_distances = []
    within_threshold = 0
    total_pairs = 0
    
    for i in range(n_mut):
        for j in range(i + 1, n_mut):
            dist = np.linalg.norm(mutation_coords[i] - mutation_coords[j])
            pairwise_distances.append(dist)
            total_pairs += 1
            
            if dist <= threshold:
                within_threshold += 1
    
    avg_pairwise_dist = np.mean(pairwise_distances)
    
    # Cluster score: fraction of pairs within threshold
    # 1.0 = all mutations clustered, 0.0 = all mutations dispersed
    cluster_score = within_threshold / total_pairs if total_pairs > 0 else 0.0
    
    # Density weight factor: same as cluster score
    # High cluster_score (clustered) -> high weight
    # Low cluster_score (dispersed) -> low weight
    density_weight_factor = cluster_score
    
    return avg_pairwise_dist, cluster_score, density_weight_factor


# ===============================
#   9. Count mutations with density weighting
# ===============================

def count_pocket_mutations_with_density(
    wt_seq: str,
    des_seq: str,
    pocket_positions: List[int],
    mapped_residues: List[Tuple[str, int]],
    coords: Dict[Tuple[str, int], np.ndarray],
    threshold: float = DENSITY_THRESHOLD,
    density_weight: float = DENSITY_WEIGHT
) -> Tuple[int, float, Dict]:
    """
    Count mutations in the pocket with density-based weighting.
    More clustered/dense mutations -> higher weight
    More dispersed mutations -> lower weight
    
    Args:
        wt_seq: Wild-type sequence
        des_seq: Design sequence
        pocket_positions: List of sequence positions (1-based)
        mapped_residues: List of (chain, resid) corresponding to positions
        coords: Dict mapping (chain, resid) to coordinates
        threshold: Distance threshold for clustering (Angstroms)
        density_weight: Multiplier for density effect
    
    Returns:
        Tuple of:
        - n_mut_simple: Simple count of mutations
        - n_mut_weighted: Weighted count based on mutation density
        - mutation_details: Dict with detailed information
    """
    n_mut_simple = 0
    N_wt = len(wt_seq)
    N_des = len(des_seq)
    
    mutations = []
    mutation_coords = []
    
    # First pass: identify all mutations and their coordinates
    for pos, res_key in zip(pocket_positions, mapped_residues):
        idx = pos - 1
        if idx < 0 or idx >= N_wt or idx >= N_des:
            continue

        aa_wt = wt_seq[idx]
        aa_des = des_seq[idx]

        if aa_wt != aa_des:
            n_mut_simple += 1
            
            if res_key in coords:
                coord = coords[res_key]
                mutation_coords.append(coord)
                
                mutations.append({
                    'position': pos,
                    'residue_key': f"{res_key[0]}{res_key[1]}",
                    'wt': aa_wt,
                    'design': aa_des,
                    'coords': coord
                })
    
    # Compute density metrics
    if mutation_coords:
        avg_dist, cluster_score, density_weight_factor = compute_mutation_density(
            mutation_coords, threshold
        )
    else:
        avg_dist = 0.0
        cluster_score = 0.0
        density_weight_factor = 1.0
    
    # Apply density weight to get final weighted mutation count
    # n_mut_weighted = n_mut_simple * (density_weight_factor ^ density_weight)
    # This means:
    #   - Clustered mutations (high density_weight_factor): weight close to n_mut_simple
    #   - Dispersed mutations (low density_weight_factor): weight reduced
    # Using power allows for tuning the strength of the effect
    weight_multiplier = density_weight_factor ** density_weight
    n_mut_weighted = n_mut_simple * weight_multiplier
    
    # Add pairwise distances to mutation details
    for i, mut in enumerate(mutations):
        mut['avg_distance_to_others'] = 0.0
        if len(mutation_coords) > 1:
            distances_to_others = []
            for j, other_coord in enumerate(mutation_coords):
                if i != j:
                    dist = np.linalg.norm(mutation_coords[i] - other_coord)
                    distances_to_others.append(dist)
            if distances_to_others:
                mut['avg_distance_to_others'] = np.mean(distances_to_others)
    
    mutation_details = {
        'mutations': mutations,
        'n_mutations': n_mut_simple,
        'avg_pairwise_distance': avg_dist,
        'cluster_score': cluster_score,
        'density_weight_factor': density_weight_factor,
        'weight_multiplier': weight_multiplier
    }
    
    return n_mut_simple, n_mut_weighted, mutation_details


# ===============================
#   10. WT vs DESIGN + density penalty
# ===============================

def compute_scores_and_final(
    wt_seq: str,
    des_seq: str,
    pocket_positions: List[int],
    mapped_residues: List[Tuple[str, int]],
    coords: Dict[Tuple[str, int], np.ndarray],
    threshold: float = DENSITY_THRESHOLD,
    density_weight: float = DENSITY_WEIGHT,
) -> dict:
    """
    Compute scores with density-based weighting for mutations.
    """
    score_wt = compute_sequence_score(wt_seq, pocket_positions)
    score_design = compute_sequence_score(des_seq, pocket_positions)
    dist_score = abs(score_design - score_wt)

    n_mut_simple, n_mut_weighted, mut_details = count_pocket_mutations_with_density(
        wt_seq, des_seq, pocket_positions, mapped_residues, coords, threshold, density_weight
    )

    # Final score uses weighted mutation count
    final_score = n_mut_weighted + ALPHA * (100 * dist_score)

    return {
        "score_wt": score_wt,
        "score_design": score_design,
        "dist_score": dist_score,
        "n_mut_simple": n_mut_simple,
        "n_mut_weighted": n_mut_weighted,
        "mutation_details": mut_details,
        "final_score": final_score,
    }


# ===============================
#   11. Main
# ===============================

def main():
    parser = argparse.ArgumentParser(
        description='Compute pocket-based scores with mutation density penalty from fpocket output and sequences.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -f complex.fa -p protein.pdb -r pocket_residues.txt
  %(prog)s -f complex.fa -p protein.pdb -r pocket_residues.txt --chains A B
  %(prog)s -f complex.fa -p protein.pdb -r pocket_residues.txt --threshold 10.0 --density-weight 1.5
  
Density weighting concept:
  - Clustered mutations (close together) = HIGHER weight in final score
  - Dispersed mutations (far apart) = LOWER weight in final score
  - Use --threshold to define "close" (default: 8.0 Angstroms)
  - Use --density-weight to control effect strength (default: 1.0, higher = stronger effect)
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
    
    parser.add_argument('--threshold',
                        type=float,
                        default=DENSITY_THRESHOLD,
                        help=f'Distance threshold for mutation clustering in Angstroms (default: {DENSITY_THRESHOLD})')
    
    parser.add_argument('--density-weight',
                        type=float,
                        default=DENSITY_WEIGHT,
                        help=f'Weight exponent for density effect (default: {DENSITY_WEIGHT}, higher = stronger effect: clustered mutations get more weight, dispersed get less)')
    
    parser.add_argument('--show-mutations',
                        action='store_true',
                        help='Show detailed information about each mutation')
    
    args = parser.parse_args()
    
    try:
        # 1. Read WT and designed sequences
        print("Reading sequences from:", args.fasta)
        wt_seq, des_seq = read_complex_fasta(args.fasta)
        print(f"WT sequence length: {len(wt_seq)}")
        print(f"Design sequence length: {len(des_seq)}")
        
        # 2. Parse PDB and build residue mapping + coordinates
        print("\nParsing PDB file:", args.pdb)
        residue_map, coords = parse_pdb_with_coords(args.pdb, args.chains)
        print(f"Total residues in PDB: {len(residue_map)}")
        print(f"CA coordinates found: {len(coords)}")
        
        # 3. Read pocket residues from fpocket output
        print("\nReading pocket residues from:", args.residues)
        pocket_residues = read_pocket_residues(args.residues)
        print(f"Pocket residues found: {len(pocket_residues)}")
        
        
        # 5. Map pocket residues to sequence positions
        print("\nMapping pocket residues to sequence positions...")
        pocket_positions, mapped_residues = map_pocket_to_sequence(pocket_residues, residue_map)
        print(f"Mapped positions: {len(pocket_positions)}")
    
        # 6. Compute scores
        print("\n" + "="*60)
        print("COMPUTING SCORES WITH MUTATION DENSITY WEIGHTING")
        print("="*60)
        
        result = compute_scores_and_final(
            wt_seq, des_seq, pocket_positions, mapped_residues, coords,
            args.threshold, args.density_weight
        )
        
        print(f"\nWT score              : {result['score_wt']:.4f}")
        print(f"DESIGN score          : {result['score_design']:.4f}")
        print(f"dist_score            : {result['dist_score']:.4f}")
        
        print(f"\nMutations (simple)    : {result['n_mut_simple']}")
        print(f"Mutations (weighted)  : {result['n_mut_weighted']:.4f}")
        
        mut_details = result['mutation_details']
        print(f"\nDensity Analysis:")
        print(f"  Avg pairwise distance : {mut_details['avg_pairwise_distance']:.2f} Å")
        print(f"  Cluster score         : {mut_details['cluster_score']:.4f} (1.0 = clustered, 0.0 = dispersed)")
        print(f"  Density weight factor : {mut_details['density_weight_factor']:.4f} (higher = more clustered)")
        print(f"  Weight multiplier     : {mut_details['weight_multiplier']:.4f}")
        print(f"  (threshold = {args.threshold} Å, density_weight = {args.density_weight})")
        print(f"\n  → Clustered mutations get HIGHER weight")
        print(f"  → Dispersed mutations get LOWER weight")
        
        print(f"\nFINAL SCORE           : {result['final_score']:.4f}")
        print(f"  (n_mut_weighted + {ALPHA} * 100 * dist_score)")
        
        # Show mutation details if requested
        if args.show_mutations and mut_details['mutations']:
            print("\n" + "="*60)
            print("MUTATION DETAILS (sorted by avg distance to other mutations)")
            print("="*60)
            mutations = sorted(mut_details['mutations'], 
                             key=lambda x: x['avg_distance_to_others'], reverse=True)
            
            print(f"\n{'Pos':<6} {'Residue':<10} {'WT->Des':<8} {'Avg Dist(Å)':<12}")
            print("-" * 50)
            for mut in mutations:
                print(f"{mut['position']:<6} {mut['residue_key']:<10} "
                      f"{mut['wt']}->{mut['design']:<6} {mut['avg_distance_to_others']:<12.2f}")
        
    except Exception as e:
        print(f"\nError: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()