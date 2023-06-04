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
	python AFanalysis.py --data_file <pkl file>
	python AFanalysis.py --pdb_file <pdb file>
	or 
	python AFanalysis.py --data_file <pkl file> --pdb_file <pdb file>

Example: 
	python AFanalysis.py --data_file result_model_5_multimer_v2_pred_4.pkl
	python AFanalysis.py --pdb_file rank0.pdb
	python AFanalysis.py --data_file result_model_5_multimer_v2_pred_4.pkl --pdb_file rank0.pdb


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
parser.add_argument("--data_file", type=str, help="the input file (pickle or JSON)")
parser.add_argument("--pdb_file", type=str, help="Input PDB file name", nargs='?')
args = parser.parse_args()

# Check if at least one argument is provided
if not args.data_file and not args.pdb_file:
    print("At least one argument (data_file or pdb_file) is required.")
    parser.print_help()
    exit()

# Check if the provided data_file is a pickle or JSON file
if args.data_file:
    if not args.data_file.endswith(".pkl") and not args.data_file.endswith(".json"):
        print("The data_file must be a pickle file (.pkl) or a JSON file (.json).")
        parser.print_help()
        exit()

# Check if the provided pdb_file is a PDB file
if args.pdb_file and not args.pdb_file.endswith(".pdb"):
    print("The pdb_file must be a PDB file (.pdb).")
    parser.print_help()
    exit()

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

if args.data_file:
    # Load the data from the input file (pickle or JSON)
    if args.data_file.endswith(".pkl"):
        with open(args.data_file, 'rb') as f:
            data = pickle.load(f)
    elif args.data_file.endswith(".json"):
        with open(args.data_file, 'r') as f:
            data = json.load(f)

    # Rename the keys as "pae" and "max_pae" if loading from pickle
    if args.data_file.endswith(".pkl"):
        data["pae"] = data.pop("predicted_aligned_error")
        data["max_pae"] = data.pop("max_predicted_aligned_error")

    # Extract the desired keys from the data
    extracted_data = {}
    for key in ["plddt", "max_pae", "pae", "ptm", "iptm", "ranking_confidence"]:
        extracted_data[key] = data.get(key, [])

    # Generate the output file name based on the input file name
    output_file = os.path.splitext(args.data_file)[0]

    # Convert the extracted data to a JSON string
    json_data = json.dumps(extracted_data, cls=NumpyEncoder)

    # Save the JSON string to a file    
    if args.data_file.endswith(".pkl"):
        # Save the JSON string to a file
        with open('%s.json' % output_file, 'w') as f:
            f.write(json_data)
    

    # Extracting relevant data
    pae_data = data.get("pae", [])
    iptm = data.get("iptm", [])
    ptm = data.get("ptm", [])
    plddt = data.get("plddt", [])

    # Plotting pae_data
    data= np.array(pae_data)
    fig, ax = plt.subplots(figsize=(7,9))
    
    if args.pdb_file:
        parser = PDB.PDBParser()
        structure = parser.get_structure("AFmodel", args.pdb_file)
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


if args.pdb_file:
    #Initiate a PyMOL session and load the provided structure. 
    pymol.finish_launching(['pymol', '-cq', '--gui'])
    pymol.cmd.load(args.pdb_file)

    # Generate the output file name based on the input file name
    fig_name = os.path.splitext(args.pdb_file)[0]

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
