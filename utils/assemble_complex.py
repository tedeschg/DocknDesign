#!/usr/bin/env python3
"""
Assemble a protein–ligand complex using MDAnalysis.

- Protein: PDB
- Ligand: PDB or SDF/MOL (SDF/MOL via RDKit → temporary PDB)

Requires:
    pip install MDAnalysis rdkit-pypi
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


def sdf_to_pdb_temp(sdf_file):
    """
    Convert an SDF/MOL file to a temporary PDB file using RDKit.

    Returns:
        path to temporary PDB file
    """
    if Chem is None or AllChem is None:
        print("ERROR: RDKit is required to handle SDF/MOL ligands.")
        print("Install it with: pip install rdkit-pypi")
        sys.exit(1)

    if not os.path.exists(sdf_file):
        print(f"ERROR: SDF/MOL file not found: {sdf_file}")
        sys.exit(1)

    suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
    mol = None
    for m in suppl:
        if m is not None:
            mol = m
            break

    if mol is None:
        print(f"ERROR: Could not read any molecule from SDF: {sdf_file}")
        sys.exit(1)

    # Ensure we have 3D coordinates
    if mol.GetNumConformers() == 0:
        print("No conformers found in SDF; generating a 3D conformer with RDKit...")
        AllChem.EmbedMolecule(mol)
        AllChem.UFFOptimizeMolecule(mol)

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmp_pdb = tmp.name
    tmp.close()

    Chem.MolToPDBFile(mol, tmp_pdb)
    print(f"Temporary ligand PDB written to: {tmp_pdb}")

    return tmp_pdb


def assemble_complex(protein_file, ligand_file, output_file):
    if not os.path.exists(protein_file):
        print(f"ERROR: Protein file not found: {protein_file}")
        sys.exit(1)
    if not os.path.exists(ligand_file):
        print(f"ERROR: Ligand file not found: {ligand_file}")
        sys.exit(1)

    lig_ext = os.path.splitext(ligand_file)[1].lower()
    tmp_pdb_to_delete = None

    # --- Load protein ---
    print(f"Loading protein from: {protein_file}")
    u_prot = mda.Universe(protein_file)

    # --- Load ligand ---
    if lig_ext in [".pdb", ".ent", ".mol2"]:
        print(f"Ligand detected as PDB-like format: {ligand_file}")
        u_lig = mda.Universe(ligand_file)
    elif lig_ext in [".sdf", ".mol"]:
        print(f"Ligand detected as SDF/MOL: {ligand_file}")
        tmp_pdb = sdf_to_pdb_temp(ligand_file)
        tmp_pdb_to_delete = tmp_pdb
        u_lig = mda.Universe(tmp_pdb)
    else:
        print(f"ERROR: Unsupported ligand format: {lig_ext}")
        print("Supported ligand formats: .pdb, .ent, .mol2, .sdf, .mol")
        sys.exit(1)

    print(f"Protein atoms: {u_prot.atoms.n_atoms}")
    print(f"Ligand atoms:  {u_lig.atoms.n_atoms}")

    # --- Merge ---
    print("Merging protein and ligand into one Universe...")
    u_complex = mda.Merge(u_prot.atoms, u_lig.atoms)

    # You can optionally set a name / residue / segment for the ligand atoms here
    # e.g. u_complex.atoms[u_prot.atoms.n_atoms:].resnames = "LIG"

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

    print("Done.")


def main():
    parser = argparse.ArgumentParser(
        description="Assemble a protein–ligand complex using MDAnalysis."
    )
    parser.add_argument("-p", "--protein", required=True, help="Protein PDB file")
    parser.add_argument("-l", "--ligand", required=True,
                        help="Ligand file (.pdb, .ent, .mol2, .sdf, .mol)")
    parser.add_argument("-o", "--output", default="complex.pdb",
                        help="Output PDB file (default: complex.pdb)")

    args = parser.parse_args()
    assemble_complex(args.protein, args.ligand, args.output)


if __name__ == "__main__":
    main()
