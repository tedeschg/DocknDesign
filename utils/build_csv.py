import pandas as pd
import os
import yaml
import argparse

def load_yaml(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"YAML file not found: {path}")
    with open(path, "r") as f:
        return yaml.safe_load(f)

def create_csv(config, output_csv):
    # Required fields
    ligands = config.get("ligand_smiles")
    protein = config.get("protein_path")

    if ligands is None:
        raise ValueError("YAML is missing required field: ligand_smiles")
    if protein is None:
        raise ValueError("YAML is missing required field: protein_path")

    # Ensure ligand_smiles is a list
    if isinstance(ligands, str):
        # Split by comma only if user still gives a single string
        ligands = [x.strip() for x in ligands.split(",")]
    elif not isinstance(ligands, list):
        raise ValueError("ligand_smiles must be a string or a list")

    # Optional fields
    complex_name = config.get("complex_name", "complex")
    protein_sequence = config.get("protein_sequence", "")

    # Build rows
    rows = []
    for idx, smi in enumerate(ligands, start=1):
        # Se ci sono più ligandi, aggiungi il numero sequenziale
        if len(ligands) > 1:
            name = f"{complex_name}_{idx}"
        else:
            name = complex_name
            
        rows.append({
            "complex_name": name,
            "protein_path": protein,
            "ligand_description": smi,
            "protein_sequence": protein_sequence
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)

    print(f"CSV created with {len(df)} ligands: {output_csv}")

def main():
    parser = argparse.ArgumentParser(description="Create CSV for DiffDock from YAML input")

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input YAML file with ligand and protein info"
    )
    parser.add_argument(
        "-o", "--output",
        default="dataset.csv",
        help="Output CSV filename (default: dataset.csv)"
    )

    args = parser.parse_args()

    # Load config
    config = load_yaml(args.input)

    # Create CSV
    create_csv(config, args.output)

if __name__ == "__main__":
    main()