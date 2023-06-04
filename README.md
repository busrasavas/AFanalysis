# AFanalysis

Here, you can use the provided python script for the analysis of AlphaFold models, including both local predictions that produce .pkl files and ColabFold predictions that produce .json files storing quality evaluation scores. 

With this script, you can only analyze one AlphaFold model at a time.

	If only pkl/json file is provided, this script outputs:
	   - a json file containing plddt, pae, max_pae, ranking confidence, ptm and iptm scores (%s.json)
	   - draws a PAE graph showing iptm and ptm values (%s_PAE.jpeg)
	   - draws a pLDDT graph showing the plddt score for each residue (%s_Plddt.jpeg)
	   
	If only pdb file is provided, this script outputs:
	- a PyMOL-made figure showing the structure colored by chains (%s_cbc.png)
	- a PyMOL-made figure showing the structure colored by AlphaFold coloring (%s_cbaf.png)
	
	If both pkl/json and pdb files are provided, this script outputs:
	- All mentioned above and updated versions of PAE and pLDDT graphs showing the residue numbers of each monomers as lines. 

Usage: 
```
	python AFanalysis.py --data_file <pkl file or json file>
	python AFanalysis.py --pdb_file <pdb file>
	or 
	python AFanalysis.py --data_file <pkl file or json file> --pdb_file <pdb file>
```
	
Example: 
```
	python AFanalysis.py --data_file result_model_1_multimer_v2_pred_1.pkl
	python AFanalysis.py --data_file result_model_1_multimer_v2_pred_1.json
	python AFanalysis.py --pdb_file rank0.pdb
	python AFanalysis.py --data_file result_model_1_multimer_v2_pred_1.pkl --input_pdb rank0.pdb
	python AFanalysis.py --data_file result_model_1_multimer_v2_pred_1.json --input_pdb rank0.pdb
```

Recommended to use this script by creating a conda environment.
```
	conda create -n AFanalysis python=3.7
	conda activate AFanalysis
	conda install -c schrodinger pymol -y
	conda install matplotlib
```

Output figure examples are listed below: 

- PAE graphs (with and without residue numbers of monomers):

<img src="https://user-images.githubusercontent.com/62547137/230650827-6aecf698-285b-4fd7-b2e1-c33d4e4fc0fd.jpeg" width="370" height="400"><img src="https://user-images.githubusercontent.com/62547137/230651275-01160bd3-3372-4102-898b-68767532f450.jpeg" width="370" height="400">


- pLDDT graphs (with and without residue numbers of monomers):

<img src="https://user-images.githubusercontent.com/62547137/230650807-575f5178-f1af-4108-8545-43005fa545b9.jpeg" width="400" height="200"><img src="https://user-images.githubusercontent.com/62547137/230651186-4a51cd95-bc12-40f7-a24d-3a85be39e5f1.jpeg" width="400" height="200">

- PyMOL-made figures:

<img src="https://user-images.githubusercontent.com/62547137/230650419-808b7340-1d56-4c2b-bb09-004e66e49687.png" width="400" height="370"><img src="https://user-images.githubusercontent.com/62547137/230650456-964bcbc0-77e1-48c9-b26c-7ff726e01a55.png" width="400" height="370">



