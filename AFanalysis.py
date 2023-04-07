#!/usr/bin/env python
# Many thanks to ChatGPT for assisting me!

"""
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
	python AFanalysis.py --input_pkl <pkl file>
	python AFanalysis.py --input_pdb <pdb file>
	or 
	python AFanalysis.py --input_pkl <pkl file> --input_pdb <pdb file>

Example: 
	python AFanalysis.py --input_pkl result_model_5_multimer_v2_pred_4.pkl
	python AFanalysis.py --input_pdb rank0.pdb
	python AFanalysis.py --input_pkl result_model_5_multimer_v2_pred_4.pkl --input_pdb rank0.pdb


Recommended to use this script by creating a conda environment.
	conda create -n AFanalysis python=3.7
	conda activate AFanalysis
	conda install -c schrodinger pymol -y
	conda install matplotlib
	conda install biopython

"""

import argparse
import json
import pickle
import numpy as np
import os
import matplotlib
from matplotlib.pylab import *
import matplotlib.pyplot as plt
import pymol
from pymol import cmd, util
import threading
from Bio import PDB

# Parse command-line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--input_pkl", type=str, help="the input pickle file to convert json file")
parser.add_argument('--input_pdb', type=str, help='Input PDB file name', nargs='?')
args = parser.parse_args()

# Check if either argument is provided
if not args.input_pkl and not args.input_pdb:
    print("At least one argument (input_pkl or input_pdb) is required.")
    parser.print_help()
    exit()

# Check if the provided input_pkl is a pickle file
if args.input_pkl and not args.input_pkl.endswith(".pkl"):
    print("The input_pkl must be a pickle file (.pkl).")
    parser.print_help()
    exit()

# Check if the provided input_pdb is a PDB file
if args.input_pdb and not args.input_pdb.endswith(".pdb"):
    print("The input_pdb must be a PDB file (.pdb).")
    parser.print_help()
    exit()

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

if args.input_pkl:
    # Load the data from the input pickle file
    with open(args.input_pkl, 'rb') as f:
        data = pickle.load(f)

    # Rename the keys as "pae" and "max_pae"
    data["pae"] = data.pop("predicted_aligned_error")
    data["max_pae"] = data.pop("max_predicted_aligned_error")

    # Extract the desired keys from the data
    extracted_data = {}
    for key in ["plddt", "max_pae", "pae","ptm", "iptm"]:
        extracted_data[key] = data[key]

    # Generate the output file name based on the input file name
    output_file = os.path.splitext(args.input_pkl)[0]

    # Convert the extracted data to a JSON string
    json_data = json.dumps(extracted_data, cls=NumpyEncoder)

    # Save the JSON string to a file
    with open('%s.json' %output_file, 'w') as f:
        f.write(json_data)

    # Extract pae scores from data
    pae_data = []
    for key in ["pae"]:
        pae_data = data[key]

    # Extract iptm&ptm scores from data
    iptm = []
    for key in ["iptm"]:
        iptm= data[key]

    ptm = []
    for key in ["ptm"]:
        ptm= data[key]
        
    plddt = []
    for key in ["plddt"]:
        plddt= data[key]

    data= np.array(pae_data)
    fig, ax = plt.subplots(figsize=(7,9))
    
    if args.input_pdb:
        parser = PDB.PDBParser()
        structure = parser.get_structure("AFmodel", args.input_pdb)
        # Initialize a dictionary to store the residue counts for each chain
        residue_counts = {}       
        # Loop over the chains in the structure
        for chain in structure.get_chains():
            # Count the number of residues in the chain
            #residue_count = sum(1 for residue in chain.get_residues() if PDB.is_aa(residue))
            #residue_count = sum(1 for residue in chain.get_residues() if PDB.is_aa(residue))
            # Count the number of residues in the chain
            residue_count = np.sum(np.fromiter((1 for residue in chain.get_residues() if PDB.is_aa(residue)), dtype=int))
            # Store the residue count in the dictionary
            residue_counts[chain.get_id()] = residue_count
        residue_numbers = list(residue_counts.values())

        # Calculate cumulative sum for graphs
        cumulative_sum = [sum(residue_numbers[:i+1]) for i in range(len(residue_numbers)-1)]
        
        #Draw PAE graph
        im = ax.imshow(data, cmap=plt.cm.Greens_r, vmin=0, vmax=30)

        #Make adjustments
        cbar = ax.figure.colorbar(im, ax=ax, orientation = 'horizontal', pad=0.12)
        cbar.set_label(r"Expected position error (Ångströms)"
                       "\n"
                       "\n"
                       f"ptm=%.3f  iptm=%.3f" % (ptm, iptm)
                       ,fontsize=18)

        cbar.ax.tick_params(labelsize=18)
        plt.xlabel("Scored residue", fontsize=18)
        plt.ylabel("Aligned residue", fontsize=18)

        for tick in cbar.ax.get_yticklabels():
            ax.tick_params(labelsize=24)

        fig.tight_layout()
        for resnum in cumulative_sum:
            ax.axvline(x=resnum, color='darkred', linestyle='--', linewidth=2, alpha=0.5)
            ax.axhline(y=resnum, color='darkred', linestyle='--', linewidth=2, alpha=0.5)
        #Save graph
        fig.savefig('%s_PAE.jpeg' %output_file, dpi=500, bbox_inches='tight')
    
    
        #Draw plddt graph
        f = plt.figure(figsize=(20,10))
        ax=plt.gca()

        residues= list(range(1, len(plddt)+1))

        ax.set_xlabel("Residues", fontsize=46, labelpad=10)
        ax.set_ylabel("pLDDT", fontsize=46)
        ax.plot(residues,plddt)

        ax.tick_params(labelsize=28)
        for resnum in cumulative_sum:
            ax.axvline(x=resnum, color='darkred', linestyle='--', linewidth=2, alpha=0.5)
        ax.get_lines()[0].set_linewidth(5)
        plt.savefig('%s_Plddt.jpeg' %output_file, dpi=500, bbox_inches='tight')
            
    else: 
        #Draw PAE graph
        im = ax.imshow(data, cmap=plt.cm.Greens_r, vmin=0, vmax=30)

        #Make adjustments
        cbar = ax.figure.colorbar(im, ax=ax, orientation = 'horizontal', pad=0.12)
        cbar.set_label(r"Expected position error (Ångströms)"
                       "\n"
                       "\n"
                       f"ptm=%.3f  iptm=%.3f" % (ptm, iptm)
                       ,fontsize=18)

        cbar.ax.tick_params(labelsize=18)
        plt.xlabel("Scored residue", fontsize=18)
        plt.ylabel("Aligned residue", fontsize=18)

        for tick in cbar.ax.get_yticklabels():
            ax.tick_params(labelsize=24)

        fig.tight_layout()

        #Save graph
        fig.savefig('%s_PAE.jpeg' %output_file, dpi=500, bbox_inches='tight')
    
    
        #Draw plddt graph
        f = plt.figure(figsize=(20,10))
        ax=plt.gca()

        residues= list(range(1, len(plddt)+1))

        ax.set_xlabel("Residues", fontsize=46, labelpad=10)
        ax.set_ylabel("pLDDT", fontsize=46)
        ax.plot(residues,plddt)

        ax.tick_params(labelsize=28)

        ax.get_lines()[0].set_linewidth(5)
        plt.savefig('%s_Plddt.jpeg' %output_file, dpi=500, bbox_inches='tight')


if args.input_pdb:
    #Initiate a PyMOL session and load the provided structure. 
    pymol.finish_launching(['pymol', '-cq', '--gui'])
    pymol.cmd.load(args.input_pdb)

    # Generate the output file name based on the input file name
    fig_name = os.path.splitext(args.input_pdb)[0]

    #Make adjustements to show structure colored by chains
    cmd. bg_color('white')
    util.cbc()
    pymol.cmd.orient('all')
   #cmd.show("surface")
   #cmd.set('transparency', 0.85)

   #Save the image as a PNG file
    pymol.cmd.png('%s_cbc.png' %fig_name, width=1600, height=1200, dpi=500, ray=1)

   #Make adjustements to show structure colored by chains
    cmd.hide("surface")
    cmd.color("blue", f"(all) and b > 90")
    cmd.color("cyan", f"(all) and b < 90 and b > 70")
    cmd.color("yellow", f"(all) and b < 70 and b > 50")
    cmd.color("orange", f"(all) and b < 50")

  #Save the image as a PNG file
    pymol.cmd.png('%s_cbaf.png' %fig_name, width=1600, height=1200, dpi=500, ray=1)
  #Quit pymol
    pymol.cmd.quit()
