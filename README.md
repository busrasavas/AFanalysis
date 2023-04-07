# AFanalysis

For a quick analysis of AlphaFold models using .pkl and/or .pdb files. 

	If only pkl file is provided, outputs:
	   - a json file containing plddt, pae, max_pae, ptm and iptm scores (%s.json)
	   - draws a PAE graph showing iptm and ptm values (%s_PAE.jpeg)
	   - draws a pLDDT graph showing the plddt score for each residue (%s_Plddt.jpeg)
	   
	If only pdb file is provided, outputs:
	- a PyMOL-made figure showing the structure colored by chains (%s_cbc.png)
	- a PyMOL-made figure showing the structure colored by AlphaFold coloring (%s_cbaf.png)
	
	If both pkl and pdb files are provided, outputs:
	- All mentioned above and updated versions of PAE and pLDDT graphs showing the residue numbers of each monomers as lines. 

Usage: 
```
	python AFanalysis.py --input_pkl <pkl file>
	python AFanalysis.py --input_pdb <pdb file>
	or 
	python AFanalysis.py --input_pkl <pkl file> --input_pdb <pdb file>
```
	
Example: 
```
	python AFanalysis.py --input_pkl result_model_5_multimer_v2_pred_4.pkl
	python AFanalysis.py --input_pdb rank0.pdb
	python AFanalysis.py --input_pkl result_model_5_multimer_v2_pred_4.pkl --input_pdb rank0.pdb
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



