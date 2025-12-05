1) sel_residues.py #Automaticamente seleziona il pocket con il druggability-score migliore

  python /home/tedeschg/prj/dockndesign/utils/sel_residues.py -i protein.pdb
  python /home/tedeschg/prj/dockndesign/utils/sel_residues.py -i protein.pdb -p 2 -o pocket2_residues.txt
  python /home/tedeschg/prj/dockndesign/utils/sel_residues.py -i protein.pdb -p 3 --format detailed
  python /home/tedeschg/prj/dockndesign/utils/sel_residues.py -i protein.pdb --fasta seq.fasta --fasta-chain A 

Note: -p 3 -i protein.pdb #Si specifica -p (ig: 3) per selezionare il terzo pocket con il druggability-score piu alto

2) build_csv.py #for DiffDock

  python /home/tedeschg/prj/dockndesign/utils/build_csv.py -i input.yaml -o hiv.csv

3) assamble_complex.py #Per assemblare la proteina.pdb con il ligando.sdf dopo il docking con boltz

  python /home/tedeschg/prj/dockndesign/utils/python assemble_complex.py -p protein.pdb -l ligand.sdf -o complex.pdb

4) lmpnn_score.py

  python /home/tedeschg/prj/dockndesign/utils/lmpnn_score.py -f complex.fa -p protein.pdb -r pocket_residues.txt
  python /home/tedeschg/prj/dockndesign/utils/lmpnn_score.py -f complex.fa -p protein.pdb -r pocket_residues.txt --chains A B

#################################################################################################
  The scoring procedure quantifies how different a designed protein sequence is from 
the wild-type (WT) sequence, focusing exclusively on the residues that constitute the 
binding pocket identified by fpocket. For each sequence (WT and design), the method 
computes a pocket score based on two averaged physicochemical properties: hydrophobicity 
and side-chain occupancy.

For every pocket position, the corresponding amino acid is retrieved from the sequence. 
Its hydrophobicity value is taken from a predefined Kyte--Doolittle scale, denoted 
as $KD(\mathrm{AA}_i)$. Side-chain occupancy is instead estimated from the number of 
heavy atoms in each amino-acid side chain, normalized by the maximum value (10 atoms 
for tryptophan). For residue $i$ in the pocket, the occupancy term is therefore:

\[
\mathrm{occ}_i = \frac{\#\text{heavy atoms of } \mathrm{AA}_i}{10}.
\]

The hydrophobicity and occupancy contributions are averaged over all valid 
pocket residues. If the pocket contains $N$ residues, the respective averages are:

\[
H_{\text{pocket}} = \frac{1}{N} \sum_{i=1}^{N} KD(\mathrm{AA}_i),
\qquad
\mathrm{Occ}_{\text{pocket}} = \frac{1}{N} \sum_{i=1}^{N} \mathrm{occ}_i.
\]

The pocket score for a given sequence is then defined as:

\[
\mathrm{score} = H_{\text{pocket}} + \mathrm{Occ}_{\text{pocket}}.
\]

After computing this score for both the WT and the designed sequence, the absolute 
difference between the two is used to quantify how much the physicochemical profile 
of the pocket has changed:

\[
\mathrm{dist\_score} = \big| \mathrm{score}_{\text{design}} - 
\mathrm{score}_{\text{WT}} \big|.
\]

In addition, the number of mutations occurring specifically at pocket positions is 
computed by comparing the WT and designed sequences residue by residue. This yields 
the quantity:

\[
n_{\text{mut,pocket}} = \#\{ i \in \text{pocket} \; | \; 
\mathrm{AA}^{\text{WT}}_i \neq \mathrm{AA}^{\text{design}}_i \}.
\]

The final score combines the pocket mutation count with the physicochemical deviation. 
The formula implemented in the script is:

\[
\mathrm{final\_score} = 
n_{\text{mut,pocket}} + \alpha \cdot \left( 100 \times \mathrm{dist\_score} \right),
\qquad \text{with } \alpha = 0.5.
\]

Thus, the final score increases both with the number of mutations in the binding 
pocket and with the magnitude of the hydrophobicity/occupancy change. Lower values 
indicate designs whose pocket properties remain closer to those of the wild type.
##########################################################################################