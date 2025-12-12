1) sel_residues.py #Automaticamente seleziona il pocket con il druggability-score migliore

  python /home/tedeschg/prj/dockndesign/utils/sel_residues.py -i protein.pdb
  python /home/tedeschg/prj/dockndesign/utils/sel_residues.py -i protein.pdb -p 2 -o pocket2_residues.txt
  python /home/tedeschg/prj/dockndesign/utils/sel_residues.py -i protein.pdb -p 3 --format detailed

Note: -p 3 -i protein.pdb #Si specifica -p (ig: 3) per selezionare il terzo pocket con il druggability-score piu alto

2) build_yaml.py

python /home/tedeschg/prj/dockndesign/utils/build_yaml.py smile.smi -o config.yaml
python /home/tedeschg/prj/dockndesign/utils/script.py smile.smi -o config.yaml -n "test_multi" -p "/path/to/protein.pdb"

3) build_csv.py #for DiffDock

  python /home/tedeschg/prj/dockndesign/utils/build_csv.py -i xxx.yaml -o xxx.csv

"""""""""""""""""""""""""""""""""
./run_diffdoc k.sh

OR

python ../../utils/docking_gnina.py --receptor_pdb ../../pdbs/pdk1/pdk1.pdb --smiles ligand.smi --autobox_ligand ../../pdbs/pdk1/pyrazoloquinazoline.pdb --gnina_path /home/tedeschg/miniforge3/envs/easydock/bin/gnina --best_only --verbose

python ../../utils/docking_gnina.py --receptor_pdb ../../pdbs/pdk1/pdk1.pdb --smiles ligand.smi --grid grid.txt --gnina_path /home/tedeschg/miniforge3/envs/easydock/bin/gnina --verbose
"""""""""""""""""""""""""""""""""

3) assamble_complex.py #Per assemblare la proteina.pdb con il ligando.sdf dopo il docking con boltz

  python /home/tedeschg/prj/dockndesign/utils/python assemble_complex.py -p protein.pdb -l ligand.sdf -o complex.pdb

4) lmpnn_score.py

  python /home/tedeschg/prj/dockndesign/utils/lmpnn_score.py -f complex.fa -p protein.pdb -r pocket_residues.txt 
  python /home/tedeschg/prj/dockndesign/utils/lmpnn_score.py -f complex.fa -p protein.pdb -r pocket_residues.txt --chains A B
--best_only --verbose
