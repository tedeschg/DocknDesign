#!/usr/bin/env python

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def check_dependencies(gnina_path=None):
    """Check if required external tools are available."""
    missing = []
    
    # Check Meeko scripts
    meeko_receptor = shutil.which("mk_prepare_receptor.py")
    meeko_ligand = shutil.which("mk_prepare_ligand.py")
    
    if not meeko_receptor:
        missing.append("mk_prepare_receptor.py (Meeko)")
    if not meeko_ligand:
        missing.append("mk_prepare_ligand.py (Meeko)")
    
    # Check GNINA
    if gnina_path:
        gnina_exe = Path(gnina_path)
        if not gnina_exe.exists():
            missing.append(f"gnina (specified path: {gnina_path})")
        elif not os.access(gnina_exe, os.X_OK):
            raise RuntimeError(f"GNINA at {gnina_path} is not executable")
        else:
            logger.info(f"Using GNINA from: {gnina_path}")
    else:
        if not shutil.which("gnina"):
            missing.append("gnina")
    
    if missing:
        raise RuntimeError(
            f"Required tools not found: {', '.join(missing)}\n"
            f"Please install Meeko: pip install meeko\n"
            f"And ensure GNINA is available."
        )
    
    logger.info("All required dependencies found (Meeko + GNINA).")


def read_grid_file(grid_path):
    """Read grid parameters (center_x, center_y, center_z, size_x, size_y, size_z) from a txt file."""
    if grid_path is None:
        return None
        
    params = {}
    try:
        with open(grid_path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) != 2:
                    raise ValueError(f"Invalid line in grid file: {line}")
                key, value = parts
                params[key] = float(value)
    except IOError as e:
        raise RuntimeError(f"Error reading grid file: {e}")

    required = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
    for k in required:
        if k not in params:
            raise ValueError(f"Missing parameter in grid file: {k}")

    return params


def read_smiles_file(smiles_path):
    """Return list of (name, smiles) tuples."""
    ligands = []
    try:
        with open(smiles_path, "r") as f:
            for idx, line in enumerate(f, start=1):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                parts = line.split()
                if len(parts) == 1:
                    smi = parts[0]
                    name = f"lig_{idx}"
                elif len(parts) >= 2:
                    smi = parts[0]
                    name = parts[1]
                else:
                    logger.warning(f"Skipping invalid line {idx}: {line}")
                    continue
                
                # Basic SMILES validation
                if len(smi) < 2:
                    logger.warning(f"Suspicious SMILES at line {idx}: '{smi}' - skipping")
                    continue
                
                ligands.append((name, smi))
    except IOError as e:
        raise RuntimeError(f"Error reading SMILES file: {e}")
    
    return ligands


def smiles_to_3d_sdf(smiles, name, out_sdf_path, random_seed=0xF00D):
    """Generate a 3D SDF from a SMILES using RDKit."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    
    # Try multiple times if embedding fails
    max_attempts = 10
    for attempt in range(max_attempts):
        result = AllChem.EmbedMolecule(mol, params)
        if result == 0:
            break
        if attempt < max_attempts - 1:
            params.randomSeed = random_seed + attempt + 1
            logger.warning(f"3D embedding attempt {attempt+1} failed for {name}, retrying...")
    
    if result != 0:
        raise ValueError(f"3D embedding failed for {name} after {max_attempts} attempts")

    # Optimize geometry more thoroughly
    result = AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    if result != 0:
        logger.warning(f"UFF optimization did not fully converge for {name}")

    # Write SDF
    writer = Chem.SDWriter(str(out_sdf_path))
    writer.write(mol)
    writer.close()
    
    # Verify SDF was created
    if not Path(out_sdf_path).exists():
        raise RuntimeError(f"SDF file was not created: {out_sdf_path}")
    
    sdf_size = Path(out_sdf_path).stat().st_size
    if sdf_size == 0:
        raise RuntimeError(f"SDF file is empty: {out_sdf_path}")
    
    logger.debug(f"Generated 3D structure for {name}: {mol.GetNumAtoms()} atoms, {sdf_size} bytes")


def sdf_to_pdbqt_with_meeko(in_sdf, out_pdbqt):
    """Convert SDF to PDBQT with Meeko."""
    meeko_script = shutil.which("mk_prepare_ligand.py")
    if meeko_script is None:
        raise RuntimeError("mk_prepare_ligand.py not found. Install meeko: pip install meeko")
    
    cmd = [
        meeko_script,
        "-i",
        str(in_sdf),
        "-o",
        str(out_pdbqt),
    ]
    
    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        if result.stdout:
            logger.debug(f"Meeko ligand prep stdout: {result.stdout}")
        if result.stderr:
            logger.debug(f"Meeko ligand prep stderr: {result.stderr}")
        
        # Verify output
        if not Path(out_pdbqt).exists():
            raise RuntimeError(f"PDBQT file was not created: {out_pdbqt}")
        
        if Path(out_pdbqt).stat().st_size == 0:
            raise RuntimeError(f"PDBQT file is empty: {out_pdbqt}")
            
        logger.debug(f"Created PDBQT with Meeko: {out_pdbqt} ({Path(out_pdbqt).stat().st_size} bytes)")
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Meeko ligand preparation failed.\n"
        error_msg += f"Command: {' '.join(cmd)}\n"
        error_msg += f"Return code: {e.returncode}\n"
        if e.stdout:
            error_msg += f"Stdout: {e.stdout}\n"
        if e.stderr:
            error_msg += f"Stderr: {e.stderr}\n"
        logger.error(error_msg)
        raise RuntimeError(f"SDF to PDBQT conversion with Meeko failed for {in_sdf}")


def clean_receptor_pdb(input_pdb, output_pdb):
    """
    Basic receptor cleanup:
    - keep ATOM and HETATM records
    - remove water residues (HOH)
    Everything else is copied as-is (TER, END, etc.).
    """
    try:
        with open(input_pdb, "r") as fin, open(output_pdb, "w") as fout:
            for line in fin:
                if line.startswith(("ATOM", "HETATM")):
                    resname = line[17:20].strip()
                    if resname == "HOH":
                        continue  # skip waters
                fout.write(line)
    except IOError as e:
        raise RuntimeError(f"Error cleaning receptor PDB: {e}")


def prepare_receptor_pdb_to_pdbqt(receptor_pdb, out_dir, ph=7.4):
    """
    Prepare receptor PDBQT using Meeko.
    Returns path to receptor PDBQT.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cleaned_pdb = out_dir / "receptor_nowat.pdb"
    clean_receptor_pdb(receptor_pdb, cleaned_pdb)

    meeko_script = shutil.which("mk_prepare_receptor.py")
    if meeko_script is None:
        raise RuntimeError(
            "mk_prepare_receptor.py not found.\n"
            "Please install Meeko: pip install meeko"
        )
    
    out_basename = out_dir / "receptor_prepared"
    cmd = [
        meeko_script,
        "-i",
        str(cleaned_pdb),
        "-o",
        str(out_basename),
        "-p",
    ]
    logger.info("Preparing receptor with Meeko (mk_prepare_receptor.py)...")
    try:
        result = subprocess.run(
            cmd, 
            check=True, 
            capture_output=True, 
            text=True
        )
        logger.debug(f"Meeko receptor prep output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Meeko receptor preparation failed: {e.stderr}")
        raise RuntimeError("Receptor preparation with Meeko failed")
    
    receptor_pdbqt = out_basename.with_suffix(".pdbqt")

    if not receptor_pdbqt.exists():
        raise RuntimeError(f"Receptor PDBQT not created: {receptor_pdbqt}")

    logger.info(f"Receptor PDBQT created: {receptor_pdbqt}")
    return receptor_pdbqt


def run_gnina(receptor_pdbqt, ligand_pdbqt, grid_params, out_dir,
              base_name, exhaustiveness=8, num_modes=10, cnn_scoring="rescore",
              gnina_path=None, seed=None, autobox_ligand=None, autobox_add=None,
              save_best_only=False):
    """Run GNINA docking for one ligand."""
    out_dir = Path(out_dir)
    out_sdf = out_dir / f"{base_name}_docked.sdf"
    log_file = out_dir / f"{base_name}_docking.log"

    # Use custom GNINA path if provided
    gnina_exe = str(gnina_path) if gnina_path else "gnina"

    cmd = [
        gnina_exe,
        "--receptor",
        str(receptor_pdbqt),
        "--ligand",
        str(ligand_pdbqt),
    ]
    
    # Use either autobox or manual grid
    if autobox_ligand:
        cmd.extend(["--autobox_ligand", str(autobox_ligand)])
        if autobox_add:
            cmd.extend(["--autobox_add", str(autobox_add)])
        logger.debug(f"Using autobox with reference: {autobox_ligand}")
    elif grid_params:
        cmd.extend([
            "--center_x", str(grid_params["center_x"]),
            "--center_y", str(grid_params["center_y"]),
            "--center_z", str(grid_params["center_z"]),
            "--size_x", str(grid_params["size_x"]),
            "--size_y", str(grid_params["size_y"]),
            "--size_z", str(grid_params["size_z"]),
        ])
        logger.debug(f"Using manual grid: center=({grid_params['center_x']}, {grid_params['center_y']}, {grid_params['center_z']})")
    else:
        raise ValueError("Either grid_params or autobox_ligand must be provided")
    
    cmd.extend([
        "--out", str(out_sdf),
        "--log", str(log_file),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--cnn_scoring", cnn_scoring,
    ])
    
    # Add seed if specified
    if seed is not None:
        cmd.extend(["--seed", str(seed)])

    logger.info(f"Docking {base_name}...")
    logger.debug(f"GNINA command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            logger.debug(f"GNINA stdout: {result.stdout}")
        if result.stderr:
            logger.debug(f"GNINA stderr: {result.stderr}")
    except subprocess.CalledProcessError as e:
        logger.error(f"GNINA docking failed for {base_name}")
        logger.error(f"Stderr: {e.stderr}")
        logger.error(f"Stdout: {e.stdout}")
        raise RuntimeError(f"Docking failed for {base_name}")
    
    # Verify output was created
    if not out_sdf.exists():
        raise RuntimeError(f"GNINA output file not created: {out_sdf}")
    
    # Filter to keep only best pose if requested
    if save_best_only:
        try:
            suppl = Chem.SDMolSupplier(str(out_sdf))
            mols = [m for m in suppl if m is not None]
            
            if mols:
                # First molecule has the best CNN score
                best_mol = mols[0]
                
                # Overwrite with only the best pose
                writer = Chem.SDWriter(str(out_sdf))
                writer.write(best_mol)
                writer.close()
                
                logger.info(f"Kept only best pose for {base_name} (out of {len(mols)} generated)")
            else:
                logger.warning(f"No valid molecules found in {out_sdf}")
        except Exception as e:
            logger.warning(f"Could not filter to best pose for {base_name}: {e}")
    
    logger.debug(f"Docking completed: {out_sdf} ({out_sdf.stat().st_size} bytes)")

    return out_sdf, log_file


def process_single_ligand(args_tuple):
    """Process a single ligand (for parallel execution)."""
    (name, smi, receptor_pdbqt, grid_params, out_dir, 
     exhaustiveness, num_modes, cnn_scoring, random_seed, gnina_path, gnina_seed,
     autobox_ligand, autobox_add, save_best_only) = args_tuple
    
    try:
        logger.info(f"Processing ligand: {name}")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            sdf_path = tmpdir / f"{name}.sdf"
            pdbqt_path = tmpdir / f"{name}.pdbqt"

            # 1) SMILES -> 3D SDF
            smiles_to_3d_sdf(smi, name, sdf_path, random_seed=random_seed)

            # 2) SDF -> PDBQT with Meeko
            sdf_to_pdbqt_with_meeko(sdf_path, pdbqt_path)

            # 3) Docking with GNINA
            out_sdf, log_file = run_gnina(
                receptor_pdbqt,
                pdbqt_path,
                grid_params,
                out_dir,
                base_name=name,
                exhaustiveness=exhaustiveness,
                num_modes=num_modes,
                cnn_scoring=cnn_scoring,
                gnina_path=gnina_path,
                seed=gnina_seed,
                autobox_ligand=autobox_ligand,
                autobox_add=autobox_add,
                save_best_only=save_best_only,
            )

            logger.info(f"✓ Completed {name}: {out_sdf}")
            return (name, True, str(out_sdf), str(log_file))

    except Exception as e:
        logger.error(f"✗ Failed {name}: {str(e)}")
        return (name, False, None, str(e))


def main():
    parser = argparse.ArgumentParser(
        description="Automated GNINA pipeline using Meeko: receptor prep + SMILES -> 3D -> PDBQT -> docking.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  %(prog)s --receptor_pdb protein.pdb --smiles ligands.smi --grid grid.txt --outdir results
  %(prog)s --receptor_pdb protein.pdb --smiles ligands.smi --grid grid.txt --best_only --num_modes 20

Requirements:
  - RDKit (for SMILES to 3D)
  - Meeko (pip install meeko)
  - GNINA (molecular docking)

Grid file format (space-separated):
  center_x 10.5
  center_y 20.3
  center_z 15.0
  size_x 20.0
  size_y 20.0
  size_z 20.0

SMILES file format (one per line):
  CCO ethanol
  c1ccccc1 benzene
  CCCC  (name will be auto-generated)
        """
    )
    parser.add_argument(
        "--receptor_pdb",
        required=True,
        help="Receptor file in PDB format (e.g., protein.pdb)",
    )
    parser.add_argument(
        "--smiles",
        required=True,
        help="SMILES file (one per line, optional name)",
    )
    parser.add_argument(
        "--grid",
        default=None,
        help="Grid txt file (center_x, center_y, center_z, size_x, size_y, size_z). Optional if using --autobox_ligand",
    )
    parser.add_argument(
        "--autobox_ligand",
        default=None,
        help="Reference ligand file (PDB/SDF/MOL2) for automatic box generation",
    )
    parser.add_argument(
        "--autobox_add",
        type=float,
        default=4.0,
        help="Buffer to add around autobox ligand in Angstroms (default: 4.0)",
    )
    parser.add_argument(
        "--outdir",
        default="gnina_out",
        help="Output directory (default: gnina_out)",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=8,
        help="GNINA exhaustiveness (default: 8)",
    )
    parser.add_argument(
        "--num_modes",
        type=int,
        default=10,
        help="Number of poses per ligand (default: 10)",
    )
    parser.add_argument(
        "--cnn_scoring",
        choices=["rescore", "all", "none"],
        default="rescore",
        help="CNN scoring mode (default: rescore)",
    )
    parser.add_argument(
        "--best_only",
        action="store_true",
        help="Save only the best CNN-scored pose per ligand (after computing all num_modes)",
    )
    parser.add_argument(
        "--ph",
        type=float,
        default=7.4,
        help="pH for receptor preparation (default: 7.4)",
    )
    parser.add_argument(
        "--gnina_path",
        default=None,
        help="Full path to GNINA executable (default: search in PATH)",
    )
    parser.add_argument(
        "--random_seed",
        type=int,
        default=0xF00D,
        help="Random seed for 3D conformer generation (default: 0xF00D)",
    )
    parser.add_argument(
        "--gnina_seed",
        type=int,
        default=None,
        help="Random seed for GNINA docking (default: None, uses GNINA's default)",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=1,
        help="Number of parallel processes (default: 1, sequential)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )
    
    args = parser.parse_args()

    # Setup logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    logger.info("=" * 60)
    logger.info("GNINA Automated Docking Pipeline (Meeko-based)")
    logger.info("=" * 60)

    # Check dependencies first
    try:
        check_dependencies(gnina_path=args.gnina_path)
    except RuntimeError as e:
        logger.error(str(e))
        sys.exit(1)

    receptor_pdb = Path(args.receptor_pdb).resolve()
    smiles_file = Path(args.smiles).resolve()
    grid_file = Path(args.grid).resolve() if args.grid else None
    autobox_ligand = Path(args.autobox_ligand).resolve() if args.autobox_ligand else None
    out_dir = Path(args.outdir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Basic checks
    if not receptor_pdb.exists():
        logger.error(f"Receptor PDB not found: {receptor_pdb}")
        sys.exit(1)
    if not smiles_file.exists():
        logger.error(f"SMILES file not found: {smiles_file}")
        sys.exit(1)
    if grid_file and not grid_file.exists():
        logger.error(f"Grid file not found: {grid_file}")
        sys.exit(1)
    if autobox_ligand and not autobox_ligand.exists():
        logger.error(f"Autobox ligand file not found: {autobox_ligand}")
        sys.exit(1)
    if not grid_file and not autobox_ligand:
        logger.error("Either --grid or --autobox_ligand must be provided")
        sys.exit(1)

    # Prepare receptor
    try:
        receptor_pdbqt = prepare_receptor_pdb_to_pdbqt(
            receptor_pdb, out_dir, ph=args.ph
        )
    except Exception as e:
        logger.error(f"Receptor preparation failed: {e}")
        sys.exit(1)

    # Read grid and ligands
    try:
        grid_params = read_grid_file(grid_file) if grid_file else None
        ligands = read_smiles_file(smiles_file)
    except Exception as e:
        logger.error(f"Error reading input files: {e}")
        sys.exit(1)

    logger.info(f"Found {len(ligands)} ligands in SMILES file")
    logger.info(f"Output directory: {out_dir}")
    if grid_params:
        logger.info(f"Using manual grid: center=({grid_params['center_x']:.2f}, {grid_params['center_y']:.2f}, {grid_params['center_z']:.2f})")
    if autobox_ligand:
        logger.info(f"Using autobox with reference: {autobox_ligand} (buffer: {args.autobox_add} Å)")
    logger.info(f"Parallelization: {args.parallel} processes")
    if args.best_only:
        logger.info(f"Mode: Save only best pose (will generate {args.num_modes} poses, keep best)")
    else:
        logger.info(f"Mode: Save all {args.num_modes} poses per ligand")
    logger.info("-" * 60)

    # Prepare arguments for parallel processing
    task_args = [
        (name, smi, receptor_pdbqt, grid_params, out_dir,
         args.exhaustiveness, args.num_modes, args.cnn_scoring, args.random_seed, 
         args.gnina_path, args.gnina_seed, autobox_ligand, args.autobox_add, args.best_only)
        for name, smi in ligands
    ]

    results = []
    if args.parallel > 1:
        # Parallel execution
        with ProcessPoolExecutor(max_workers=args.parallel) as executor:
            futures = {executor.submit(process_single_ligand, task): task[0] 
                      for task in task_args}
            
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
    else:
        # Sequential execution
        for task in task_args:
            result = process_single_ligand(task)
            results.append(result)

    # Summary
    logger.info("=" * 60)
    logger.info("DOCKING SUMMARY")
    logger.info("=" * 60)
    
    successful = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]
    
    logger.info(f"Total ligands: {len(results)}")
    logger.info(f"Successful: {len(successful)}")
    logger.info(f"Failed: {len(failed)}")
    
    if successful:
        logger.info("\nSuccessful ligands:")
        for name, _, output, log in successful:
            logger.info(f"  ✓ {name}")
            logger.info(f"    Output: {output}")
            logger.info(f"    Log: {log}")
    
    if failed:
        logger.warning("\nFailed ligands:")
        for name, _, _, error in failed:
            logger.warning(f"  ✗ {name}: {error}")
    
    logger.info("\n✓ Pipeline completed!")


if __name__ == "__main__":
    main()