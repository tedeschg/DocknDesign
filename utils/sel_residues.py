#!/usr/bin/env python3
"""
Script to run fpocket, identify binding pocket residues and output them to a text file.
Also generates a grid.txt file with pocket center coordinates for docking.
Requires: fpocket installed and available in PATH
"""

import argparse
import subprocess
import os
import sys
import re

import os
import sys
import subprocess
import shutil

def run_fpocket(pdb_file):
    """
    Run fpocket on the input PDB file and ensure that the *_out directory
    ends up nel working directory corrente (os.getcwd()).
    """
    try:
        pdb_file_abs = os.path.abspath(pdb_file)
        pdb_dir = os.path.dirname(pdb_file_abs)
        base_name = os.path.splitext(os.path.basename(pdb_file_abs))[0]

        cwd = os.getcwd()
        print(f"Running fpocket on {pdb_file_abs}...")
        print(f"Current working directory: {cwd}")

        # Lancia fpocket
        result = subprocess.run(
            ['fpocket', '-f', pdb_file_abs],
            capture_output=True,
            text=True,
            check=True
        )

        # Possibili posizioni della cartella di output
        out_name = f"{base_name}_out"
        expected_in_cwd = os.path.join(cwd, out_name)
        expected_near_pdb = os.path.join(pdb_dir, out_name)

        # 1) Se fpocket ha già creato la cartella dove siamo, usa quella
        if os.path.exists(expected_in_cwd):
            output_dir = expected_in_cwd

        # 2) Altrimenti se l'ha creata accanto al PDB, spostala nel cwd
        elif os.path.exists(expected_near_pdb):
            print(f"Found output in {expected_near_pdb}, moving to {expected_in_cwd}")
            # Se per qualche motivo esiste già, fermiamoci per sicurezza
            if os.path.exists(expected_in_cwd):
                print(f"Error: target output directory already exists: {expected_in_cwd}")
                sys.exit(1)
            shutil.move(expected_near_pdb, expected_in_cwd)
            output_dir = expected_in_cwd

        else:
            print("Error: Expected output directory not found in any of:")
            print(f"  - {expected_in_cwd}")
            print(f"  - {expected_near_pdb}")
            print("stdout from fpocket:")
            print(result.stdout)
            print("stderr from fpocket:")
            print(result.stderr)
            sys.exit(1)

        print(f"fpocket completed. Output directory: {output_dir}")
        return output_dir

    except subprocess.CalledProcessError as e:
        print(f"Error running fpocket: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        print("Error: fpocket not found. Please ensure fpocket is installed and in your PATH")
        sys.exit(1)


def get_pockets_ranked_by_druggability(info_file):
    """
    Parse fpocket info file and rank all pockets by druggability score.
    
    Args:
        info_file: Path to the fpocket info.txt file
        
    Returns:
        List of tuples (rank, original_pocket_number, druggability_score) sorted by score (highest first)
    """
    if not os.path.exists(info_file):
        print(f"Error: Info file not found: {info_file}")
        sys.exit(1)
    
    with open(info_file, 'r') as f:
        content = f.read()
    
    # Find all pockets and their druggability scores
    # Pattern: Pocket X : ... Druggability Score : Y.YY
    pocket_pattern = r"Pocket\s+(\d+)\s*:.*?Druggability Score\s*:\s*([\d.]+)"
    matches = re.findall(pocket_pattern, content, re.DOTALL | re.IGNORECASE)
    
    if not matches:
        print("Warning: No druggability scores found. Using pocket 1 as default.")
        return [(1, 1, None)]
    
    # Sort by druggability score (highest first)
    sorted_pockets = sorted(matches, key=lambda x: float(x[1]), reverse=True)
    
    # Create ranked list with (rank, original_pocket_num, score)
    ranked_pockets = [
        (rank + 1, int(pocket), float(score)) 
        for rank, (pocket, score) in enumerate(sorted_pockets)
    ]
    
    print(f"\nFound {len(ranked_pockets)} pockets ranked by druggability:")
    print(f"{'Rank':<6} {'Pocket':<8} {'Druggability Score'}")
    print("-" * 40)
    for rank, pocket, score in ranked_pockets:
        print(f"{rank:<6} {pocket:<8} {score:.3f}")
    
    return ranked_pockets

def parse_fpocket_info(info_file, pocket_number=1):
    """
    Parse fpocket info file to extract residues from a specific pocket.
    
    Args:
        info_file: Path to the fpocket info.txt file
        pocket_number: Pocket number to extract (default: 1)
        
    Returns:
        List of tuples (chain, resid, resname)
    """
    if not os.path.exists(info_file):
        print(f"Error: Info file not found: {info_file}")
        sys.exit(1)
    
    with open(info_file, 'r') as f:
        content = f.read()
    
    # Find the pocket section
    pocket_pattern = rf"Pocket\s+{pocket_number}\s*:(.*?)(?=Pocket\s+\d+\s*:|$)"
    pocket_match = re.search(pocket_pattern, content, re.DOTALL)
    
    if not pocket_match:
        print(f"Error: Pocket {pocket_number} not found in info file")
        sys.exit(1)
    
    pocket_section = pocket_match.group(1)
    
    # Extract residues - looking for patterns like "A 123 LEU" or similar
    residues = []
    
    # Try to find residue list section
    if "Residues" in pocket_section or "residue" in pocket_section.lower():
        lines = pocket_section.split('\n')
        for line in lines:
            # Look for lines with chain, residue number, and residue name
            # Common formats: "A 123 LEU", "A:123:LEU", etc.
            parts = line.strip().split()
            
            # Try to identify chain, resid, resname patterns
            if len(parts) >= 3:
                # Check if second element is a number (residue ID)
                try:
                    resid = int(parts[1])
                    chain = parts[0]
                    resname = parts[2] if len(parts) > 2 else ""
                    residues.append((chain, resid, resname))
                except ValueError:
                    continue
    
    return residues

def parse_pocket_pdb(pocket_pdb):
    """
    Parse pocket PDB file to extract residues and atom coordinates.
    
    Args:
        pocket_pdb: Path to the pocket PDB file
        
    Returns:
        Tuple of (residues_list, coordinates_list)
        - residues_list: List of tuples (chain, resid, resname)
        - coordinates_list: List of tuples (x, y, z) for all atoms
    """
    if not os.path.exists(pocket_pdb):
        print(f"Error: Pocket PDB file not found: {pocket_pdb}")
        sys.exit(1)
    
    residues = set()
    coordinates = []
    
    with open(pocket_pdb, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain = line[21].strip()
                if not chain:
                    chain = 'A'  # Default chain if empty
                resid = int(line[22:26].strip())
                resname = line[17:20].strip()
                residues.add((chain, resid, resname))
                
                # Extract coordinates
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                coordinates.append((x, y, z))
    
    return sorted(list(residues)), coordinates

def calculate_pocket_center(coordinates):
    """
    Calculate the geometric center (center of mass) of pocket atoms.
    
    Args:
        coordinates: List of tuples (x, y, z)
        
    Returns:
        Tuple (center_x, center_y, center_z)
    """
    if not coordinates:
        return (0.0, 0.0, 0.0)
    
    n = len(coordinates)
    center_x = sum(coord[0] for coord in coordinates) / n
    center_y = sum(coord[1] for coord in coordinates) / n
    center_z = sum(coord[2] for coord in coordinates) / n
    
    return (center_x, center_y, center_z)

def write_grid_file(center, output_file='grid.txt', box_size=20.0):
    """
    Write grid configuration file for docking.
    
    Args:
        center: Tuple (x, y, z) coordinates of pocket center
        output_file: Output file path
        box_size: Size of the grid box in each dimension (default: 20.0 Å)
    """
    center_x, center_y, center_z = center
    
    with open(output_file, 'w') as f:
        f.write(f"center_x {center_x:.3f}\n")
        f.write(f"center_y {center_y:.3f}\n")
        f.write(f"center_z {center_z:.3f}\n")
        f.write(f"size_x {box_size:.1f}\n")
        f.write(f"size_y {box_size:.1f}\n")
        f.write(f"size_z {box_size:.1f}\n")
    
    print(f"\nGrid file written to: {output_file}")
    print(f"Center: ({center_x:.3f}, {center_y:.3f}, {center_z:.3f})")
    print(f"Box size: {box_size:.1f} Å in each dimension")

def write_residues_output(residues, output_file='pocket_residues.txt', format_type='simple'):
    """
    Write residues to output file.
    
    Args:
        residues: List of tuples (chain, resid, resname)
        output_file: Output file path
        format_type: 'simple' for "A123 B45" format, 'detailed' for more info
    """
    with open(output_file, 'w') as f:
        if format_type == 'simple':
            # Write as single line: A123 B45 C67
            residue_list = [f'{chain}{resid}' for chain, resid, _ in residues]
            f.write(' '.join(residue_list))
        else:
            # Write detailed format: one per line
            for chain, resid, resname in residues:
                f.write(f'{chain} {resid} {resname}\n')
    
    print(f"\nOutput written to: {output_file}")
    print(f"Total pocket residues: {len(residues)}")

def main():
    parser = argparse.ArgumentParser(
        description='Run fpocket and extract binding pocket residues with grid coordinates.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i protein.pdb
  %(prog)s -i protein.pdb -p 2 -o pocket2_residues.txt
  %(prog)s -i protein.pdb -p 3 --format detailed --box-size 25.0
  
Pockets are automatically ranked by druggability score (highest first).
-p 1 selects the best pocket, -p 2 the second best, etc.

Outputs:
  - pocket_residues.txt: List of pocket residues
  - grid.txt: Grid box configuration for docking (center and size)
        """
    )
    
    parser.add_argument('-i', '--input', 
                        required=True,
                        help='Input PDB file (protein structure)')
    
    parser.add_argument('-p', '--pocket',
                        type=int,
                        default=1,
                        help='Pocket rank by druggability (1=best, 2=second best, etc.) (default: 1)')
    
    parser.add_argument('-o', '--output',
                        default='pocket_residues.txt',
                        help='Output text file for pocket residues (default: pocket_residues.txt)')
    
    parser.add_argument('-g', '--grid',
                        default='grid.txt',
                        help='Output grid file for docking (default: grid.txt)')
    
    parser.add_argument('--box-size',
                        type=float,
                        default=20.0,
                        help='Grid box size in Angstroms (default: 20.0)')
    
    parser.add_argument('--format',
                        choices=['simple', 'detailed'],
                        default='simple',
                        help='Output format: simple (A123 B45) or detailed (one per line with resname)')
    
    parser.add_argument('--skip-fpocket',
                        action='store_true',
                        help='Skip running fpocket (use existing output directory)')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Run fpocket
    if not args.skip_fpocket:
        output_dir = run_fpocket(args.input)
    else:
        base_name = os.path.splitext(os.path.basename(args.input))[0]
        output_dir = f"{base_name}_out"
        print(f"Skipping fpocket run, using existing directory: {output_dir}")
    
    # Get pockets ranked by druggability
    base_name = os.path.splitext(os.path.basename(args.input))[0]
    info_file = os.path.join(output_dir, f"{base_name}_info.txt")
    
    ranked_pockets = get_pockets_ranked_by_druggability(info_file)
    
    # Select pocket by rank
    if args.pocket < 1 or args.pocket > len(ranked_pockets):
        print(f"\nError: Rank {args.pocket} is out of range. Available ranks: 1-{len(ranked_pockets)}")
        sys.exit(1)
    
    rank, pocket_number, drug_score = ranked_pockets[args.pocket - 1]
    print(f"\nSelected: Rank {rank}, Pocket {pocket_number} (Druggability: {drug_score:.3f})")
    
    # Try to parse pocket PDB file first (more reliable)
    pocket_pdb = os.path.join(output_dir, 'pockets', f'pocket{pocket_number}_atm.pdb')
    
    if os.path.exists(pocket_pdb):
        print(f"\nParsing pocket PDB file: {pocket_pdb}")
        residues, coordinates = parse_pocket_pdb(pocket_pdb)
        
        # Calculate pocket center
        center = calculate_pocket_center(coordinates)
        
        # Write grid file
        write_grid_file(center, args.grid, args.box_size)
    else:
        # Fallback to info file (no coordinates available)
        print(f"\nPocket PDB not found, trying info file...")
        print("Warning: Cannot generate grid.txt without pocket PDB coordinates")
        residues = parse_fpocket_info(info_file, pocket_number)
        center = None
    
    if not residues:
        print(f"Warning: No residues found for pocket {pocket_number}")
        sys.exit(1)
    
    # Write residues output
    write_residues_output(residues, args.output, args.format)
    
    # Print summary
    print(f"\nSummary:")
    print(f"Pocket number: {pocket_number}")
    print(f"Residues found: {len(residues)}")
    if len(residues) > 0:
        chains = set(chain for chain, _, _ in residues)
        print(f"Chains involved: {', '.join(sorted(chains))}")

if __name__ == '__main__':
    main()