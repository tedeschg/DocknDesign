#!/usr/bin/env python3
"""
Assemble a protein–ligand complex using the best docking pose from GNINA output.

- Protein: PDB
- Ligand: SDF (multi-pose from GNINA) - automatically extracts the best (first) pose
- Can also handle single-molecule SDF/MOL/PDB files

Requires:
    pip install MDAnalysis rdkit
"""

import argparse
import os
import sys
import tempfile

import MDAnalysis as mda

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None
    AllChem = None


def extract_best_pose_from_sdf(sdf_file, pose_index=0):
    """
    Extract a specific pose from an SDF file (default: first/best pose).
    
    Args:
        sdf_file: Path to SDF file with multiple poses
        pose_index: Index of pose to extract (0 = best/first)
    
    Returns:
        Path to temporary PDB file containing the selected pose
    """
    if Chem is None or AllChem is None:
        print("ERROR: RDKit is required to handle SDF ligands.")
        print("Install it with: pip install rdkit")
        sys.exit(1)

    if not os.path.exists(sdf_file):
        print(f"ERROR: SDF file not found: {sdf_file}")
        sys.exit(1)

    # Read all molecules from SDF
    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mols = [m for m in suppl if m is not None]
    
    if not mols:
        print(f"ERROR: Could not read any molecules from SDF: {sdf_file}")
        sys.exit(1)
    
    print(f"Found {len(mols)} pose(s) in SDF file")
    
    if pose_index >= len(mols):
        print(f"ERROR: Requested pose index {pose_index} but only {len(mols)} poses available")
        sys.exit(1)
    
    mol = mols[pose_index]
    print(f"Extracting pose {pose_index + 1}/{len(mols)} (index {pose_index})")
    
    # Print properties if available (GNINA scores)
    if mol.HasProp('_Name'):
        print(f"  Ligand name: {mol.GetProp('_Name')}")
    if mol.HasProp('minimizedAffinity'):
        print(f"  Affinity: {mol.GetProp('minimizedAffinity')} kcal/mol")
    if mol.HasProp('CNNscore'):
        print(f"  CNN score: {mol.GetProp('CNNscore')}")
    if mol.HasProp('CNNaffinity'):
        print(f"  CNN affinity: {mol.GetProp('CNNaffinity')}")
    
    # Ensure we have 3D coordinates
    if mol.GetNumConformers() == 0:
        print("WARNING: No conformers found; generating 3D coordinates...")
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

    # Write to temporary PDB
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmp_pdb = tmp.name
    tmp.close()

    Chem.MolToPDBFile(mol, tmp_pdb)
    print(f"Temporary ligand PDB written to: {tmp_pdb}")

    return tmp_pdb


def sdf_to_pdb_temp(sdf_file, pose_index=0):
    """
    Convert an SDF/MOL file to a temporary PDB file using RDKit.
    For multi-pose SDF files, extracts the specified pose (default: first/best).
    
    Args:
        sdf_file: Path to SDF file
        pose_index: Which pose to extract (0 = best/first)
    
    Returns:
        Path to temporary PDB file
    """
    return extract_best_pose_from_sdf(sdf_file, pose_index)


def assemble_complex(protein_file, ligand_file, output_file, pose_index=0):
    """
    Assemble protein-ligand complex.
    
    Args:
        protein_file: Path to protein PDB
        ligand_file: Path to ligand file (PDB, SDF, MOL, MOL2)
        output_file: Path for output complex PDB
        pose_index: For multi-pose SDF files, which pose to use (0 = best)
    """
    if not os.path.exists(protein_file):
        print(f"ERROR: Protein file not found: {protein_file}")
        sys.exit(1)
    if not os.path.exists(ligand_file):
        print(f"ERROR: Ligand file not found: {ligand_file}")
        sys.exit(1)

    lig_ext = os.path.splitext(ligand_file)[1].lower()
    tmp_pdb_to_delete = None

    # --- Load protein ---
    print(f"\nLoading protein from: {protein_file}")
    u_prot = mda.Universe(protein_file)

    # --- Load ligand ---
    if lig_ext in [".pdb", ".ent", ".mol2"]:
        print(f"Ligand detected as PDB-like format: {ligand_file}")
        u_lig = mda.Universe(ligand_file)
    elif lig_ext in [".sdf", ".mol"]:
        print(f"Ligand detected as SDF/MOL: {ligand_file}")
        tmp_pdb = sdf_to_pdb_temp(ligand_file, pose_index)
        tmp_pdb_to_delete = tmp_pdb
        u_lig = mda.Universe(tmp_pdb)
    else:
        print(f"ERROR: Unsupported ligand format: {lig_ext}")
        print("Supported ligand formats: .pdb, .ent, .mol2, .sdf, .mol")
        sys.exit(1)

    print(f"\nProtein atoms: {u_prot.atoms.n_atoms}")
    print(f"Ligand atoms:  {u_lig.atoms.n_atoms}")

    # --- Merge ---
    print("\nMerging protein and ligand into one Universe...")
    u_complex = mda.Merge(u_prot.atoms, u_lig.atoms)

    # Set ligand residue name to LIG for clarity
    ligand_start_idx = u_prot.atoms.n_atoms
    u_complex.atoms[ligand_start_idx:].residues.resnames = "LIG"
    
    print(f"Complex atoms total: {u_complex.atoms.n_atoms}")
    print(f"Writing complex to: {output_file}")
    u_complex.atoms.write(output_file)

    # Cleanup temporary file if used
    if tmp_pdb_to_delete is not None:
        try:
            os.remove(tmp_pdb_to_delete)
            print(f"Temporary file removed: {tmp_pdb_to_delete}")
        except OSError:
            print(f"Warning: could not remove temporary file: {tmp_pdb_to_delete}")

    print("\n✓ Complex assembly completed successfully!")


def main():
    parser = argparse.ArgumentParser(
        description="Assemble a protein–ligand complex using the best docking pose.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract best pose from GNINA output
  %(prog)s -p protein.pdb -l docked_poses.sdf -o complex.pdb

  # Extract 2nd best pose (index 1)
  %(prog)s -p protein.pdb -l docked_poses.sdf -o complex.pdb --pose 1

  # Use a single-pose ligand file
  %(prog)s -p protein.pdb -l ligand.pdb -o complex.pdb

Note: For multi-pose SDF files (e.g., from GNINA), pose 0 is typically the best pose.
        """
    )
    parser.add_argument(
        "-p", "--protein",
        required=True,
        help="Protein PDB file"
    )
    parser.add_argument(
        "-l", "--ligand",
        required=True,
        help="Ligand file (.pdb, .ent, .mol2, .sdf, .mol). For SDF, extracts best pose by default."
    )
    parser.add_argument(
        "-o", "--output",
        default="complex.pdb",
        help="Output PDB file (default: complex.pdb)"
    )
    parser.add_argument(
        "--pose",
        type=int,
        default=0,
        help="For multi-pose SDF files, which pose to extract (0 = best/first, 1 = second, etc.). Default: 0"
    )

    args = parser.parse_args()
    
    print("=" * 60)
    print("Protein-Ligand Complex Assembly")
    print("=" * 60)
    
    assemble_complex(args.protein, args.ligand, args.output, args.pose)


if __name__ == "__main__":
    main()