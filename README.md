# AFanalysis
Quick analysis for AlphaFold models. 

Creates an json file from the provided pkl file containing plddt, pae, max_pae, ptm and iptm scores. Draws a PAE graph showing iptm and ptm values. 

Creates two figures from the provided PDB file, one is colored by chains and other is colored by AF coloring. 


Usage: 
"""
	python AFanalysis.py --input_pkl <pkl file>
	python AFanalysis.py --input_pdb <pdb file>
	or 
	python AFanalysis.py --input_pkl <pkl file> --input_pdb <pdb file>
"""
	
Example: 
	python AFanalysis.py --input_pkl result_model_5_multimer_v2_pred_4.pkl
	python AFanalysis.py --input_pdb rank0.pdb
	python AFanalysis.py --input_pkl result_model_5_multimer_v2_pred_4.pkl --input_pdb rank0.pdb


Recommended to use this script by creating a conda environment.
	conda create -n AFanalysis python=3.7
	conda activate AFanalysis
	conda install -c schrodinger pymol -y
	conda install matplotlib



