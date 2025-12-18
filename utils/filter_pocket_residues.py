 #!/usr/bin/env python3
"""
Filter pocket residues by distance from ligand.
Useful for reducing noise by keeping only residues within a certain distance.
"""

import argparse
import numpy as np
from typing import List, Tuple, Optional


def extract_ligand_center(pdb_file: str, ligand_chain: str, ligand_resid: int) -> Optional[np.ndarray]:
    """Extract center coordinates of ligand."""
    coords = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            
            chain = line[21].strip() or "A"
            resid = int(line[22:26].strip())
            
            if chain == ligand_chain and resid == ligand_resid:
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append(np.array([x, y, z]))
                except ValueError:
                    continue
    
    if not coords:
        return None
    
    return np.mean(coords, axis=0)


def extract_residue_center(pdb_file: str, chain: str, resid: int) -> Optional[np.ndarray]:
    """Extract center coordinates of a residue."""
    coords = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            
            res_chain = line[21].strip() or "A"
            res_resid = int(line[22:26].strip())
            
            if res_chain == chain and res_resid == resid:
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append(np.array([x, y, z]))
                except ValueError:
                    continue
    
    if not coords:
        return None
    
    return np.mean(coords, axis=0)


def read_pocket_residues(pocket_file: str) -> List[Tuple[str, int]]:
    """Read pocket residues from file (format: A123 B45 ...)."""
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
            print(f"Warning: Could not parse '{token}'")
            continue
    
    return residues


def filter_residues_by_distance(
    pdb_file: str,
    pocket_residues: List[Tuple[str, int]],
    ligand_center: np.ndarray,
    max_distance: float
) -> List[Tuple[str, int, float]]:
    """
    Filter pocket residues by distance from ligand.
    
    Returns:
        List of (chain, resid, distance) for residues within max_distance
    """
    filtered = []
    
    for chain, resid in pocket_residues:
        res_center = extract_residue_center(pdb_file, chain, resid)
        
        if res_center is None:
            print(f"Warning: Could not find residue {chain}{resid}")
            continue
        
        distance = np.linalg.norm(res_center - ligand_center)
        
        if distance <= max_distance:
            filtered.append((chain, resid, distance))
    
    return filtered


def main():
    parser = argparse.ArgumentParser(
        description="Filter pocket residues by distance from ligand",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Keep only residues within 5Å of ligand
  %(prog)s -p protein.pdb -r pocket_residues.txt \\
      --ligand-chain A --ligand-resid 500 \\
      --distance 5.0 -o pocket_filtered.txt
  
  # Keep residues within 8Å and show statistics
  %(prog)s -p protein.pdb -r pocket_residues.txt \\
      --ligand-chain A --ligand-resid 500 \\
      --distance 8.0 --show-stats
        """
    )
    
    parser.add_argument("-p", "--pdb", required=True,
                       help="PDB file")
    parser.add_argument("-r", "--residues", required=True,
                       help="Pocket residues file (from fpocket)")
    parser.add_argument("--ligand-chain", required=True,
                       help="Ligand chain ID")
    parser.add_argument("--ligand-resid", type=int, required=True,
                       help="Ligand residue ID")
    parser.add_argument("--distance", type=float, default=6.0,
                       help="Maximum distance from ligand (Angstroms, default: 6.0)")
    parser.add_argument("-o", "--output", default="pocket_filtered.txt",
                       help="Output file (default: pocket_filtered.txt)")
    parser.add_argument("--show-stats", action="store_true",
                       help="Show distance statistics")
    parser.add_argument("--verbose", action="store_true",
                       help="Show detailed filtering information")
    
    args = parser.parse_args()
    
    print(f"Reading pocket residues from: {args.residues}")
    pocket_residues = read_pocket_residues(args.residues)
    print(f"Found {len(pocket_residues)} pocket residues")
    
    print(f"\nExtracting ligand center from: {args.pdb}")
    print(f"Ligand: chain {args.ligand_chain}, residue {args.ligand_resid}")
    ligand_center = extract_ligand_center(args.pdb, args.ligand_chain, args.ligand_resid)
    
    if ligand_center is None:
        print(f"Error: Ligand not found at {args.ligand_chain}{args.ligand_resid}")
        return 1
    
    print(f"Ligand center: [{ligand_center[0]:.2f}, {ligand_center[1]:.2f}, {ligand_center[2]:.2f}]")
    
    print(f"\nFiltering residues within {args.distance} Å...")
    filtered = filter_residues_by_distance(
        args.pdb, pocket_residues, ligand_center, args.distance
    )
    
    print(f"Residues within {args.distance} Å: {len(filtered)}/{len(pocket_residues)}")
    
    if args.show_stats and filtered:
        distances = [d for _, _, d in filtered]
        print(f"\nDistance statistics:")
        print(f"  Min distance : {min(distances):.2f} Å")
        print(f"  Max distance : {max(distances):.2f} Å")
        print(f"  Mean distance: {np.mean(distances):.2f} Å")
        print(f"  Median distance: {np.median(distances):.2f} Å")
        
        # Distance distribution
        bins = [0, 3.5, 6.0, 10.0, args.distance + 1]
        labels = ["0-3.5Å (direct contact)", "3.5-6Å (2nd shell)", 
                 f"6-10Å (peripheral)", f"10-{args.distance}Å (distant)"]
        
        print(f"\nDistance distribution:")
        for i in range(len(bins)-1):
            count = sum(1 for d in distances if bins[i] <= d < bins[i+1])
            if count > 0:
                print(f"  {labels[i]}: {count} residues")
    
    if args.verbose and filtered:
        print(f"\nFiltered residues:")
        print(f"{'Residue':<10} {'Distance (Å)':<12}")
        print("-" * 25)
        for chain, resid, dist in sorted(filtered, key=lambda x: x[2]):
            print(f"{chain}{resid:<9} {dist:<12.2f}")
    
    # Write output
    with open(args.output, 'w') as f:
        residues_str = " ".join([f"{chain}{resid}" for chain, resid, _ in filtered])
        f.write(residues_str)
    
    print(f"\nFiltered residues written to: {args.output}")
    print(f"Use this file with pocket_mutation_score_enhanced.py:")
    print(f"  python pocket_mutation_score_enhanced.py -r {args.output} ...")
    
    return 0


if __name__ == "__main__":
    exit(main())