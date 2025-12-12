#!/usr/bin/env python3
"""
Script to convert SMILES from a file into a YAML configuration file
"""

import yaml
import argparse
from pathlib import Path


def read_smiles_file(smiles_file):
    """
    Read SMILES from a file. Extracts only the SMILES string (first column).
    
    Args:
        smiles_file: Path to the SMILES file
        
    Returns:
        List of SMILES strings
    """
    smiles_list = []
    
    with open(smiles_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Extract only the first column (SMILES)
            if '\t' in line or ' ' in line:
                # Split by tab or space and take first element
                smiles = line.split()[0]
                smiles_list.append(smiles)
            else:
                # Simple format: one SMILES per line
                smiles_list.append(line)
    
    return smiles_list


def create_yaml_config(smiles_list, complex_name, protein_path, output_file):
    """
    Create YAML configuration file with SMILES data
    
    Args:
        smiles_list: List of SMILES strings
        complex_name: Name for the complex
        protein_path: Path to the protein PDB file
        output_file: Output YAML file path
    """
    config = {
        'complex_name': complex_name,
        'ligand_smiles': smiles_list,
        'protein_path': protein_path,
        'protein_sequence': ""
    }
    
    with open(output_file, 'w') as f:
        yaml.dump(config, f, default_flow_style=False, sort_keys=False, allow_unicode=True)
    
    print(f"✓ Created YAML file: {output_file}")
    print(f"✓ Total SMILES: {len(smiles_list)}")


def main():
    parser = argparse.ArgumentParser(
        description='Convert SMILES file to YAML configuration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python script.py smile.smi -o config.yaml
  
  # Specify complex name and protein path
  python script.py smile.smi -o config.yaml -n "my_complex" -p "/path/to/protein.pdb"
  
  # Custom all parameters
  python script.py smile.smi -o output.yaml -n "test_multi" -p "/home/user/protein.pdb"
        """
    )
    
    parser.add_argument('smiles_file', type=str,
                       help='Input SMILES file (.smi)')
    parser.add_argument('-o', '--output', type=str, default='config.yaml',
                       help='Output YAML file (default: config.yaml)')
    parser.add_argument('-n', '--name', type=str, default='test_multi',
                       help='Complex name (default: test_multi)')
    parser.add_argument('-p', '--protein', type=str, 
                       default='/home/tedeschg/prj/dockndesign/pdbs/pdk1/pdk1.pdb',
                       help='Path to protein PDB file')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.smiles_file).exists():
        print(f"Error: SMILES file not found: {args.smiles_file}")
        return 1
    
    # Read SMILES
    print(f"Reading SMILES from: {args.smiles_file}")
    smiles_list = read_smiles_file(args.smiles_file)
    
    if not smiles_list:
        print("Error: No SMILES found in input file")
        return 1
    
    # Create YAML
    create_yaml_config(smiles_list, args.name, args.protein, args.output)
    
    return 0


if __name__ == '__main__':
    exit(main())