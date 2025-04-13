# -*- coding: utf-8 -*-
"""
rep_sample_generator.py input_folder

This script is a representative sample generator that selects certain structures from an input folder to maximise representation for sequence length and secondary structure proportion.
The global dataframe is also a useful summary for all the structures that will be tested.

A. Global dataframe generation
    1. Parse the input folder file by file
    2. For each file A. compute sequence length B. Compute DSSP string C. DSSP percentages and ¿D? Find PDB ID to track where they come from if possible
    3. Append each file to a dataframe with the following columns:
        [file_ID] ¿[PDB_ID]? [sequence] [sequence_length] [DSSP] [H%] [E%]
    4. Save the output to a csv file
B. Distributions of parameters
    1. Visualise the distributions of length and H% and E%
    2. Determine adequate tresholds for the following categories:
        a. Long length
        b. Short Length
        c. Medium length
        d. High alpha helix
        e. High beta strand
        f. Mixed SS
C. Selected representatives dataframe
    1. Select a representative for each category
    2. Generate a smaller dataframe with them with a new column [category]
    3. Save it to a csv file


NOTE: PDB estimate will be implemented in the future

07-04-2025
"""
#%% Import modules

#A.
from pathlib import Path
import argparse
import pandas as pd
from Bio.PDB import PDBParser, PPBuilder, DSSP

#B.
import numpy as np

#%%Functions


def sequence_calculator(input_file):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_file)

    # Extract sequences
    ppb = PPBuilder()
    sequences = []
    for pp in ppb.build_peptides(structure):
        seq = pp.get_sequence()
        sequences.append(str(seq))

    # Combine sequences for all peptides (if multiple chains)
    full_sequence = "".join(sequences)
    return full_sequence


def DSSP_string_generator(input_file):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_file)
    model = structure[0]  # DSSP works on one model, usually model 0

    # DSSP calculation
    # Ensure that you have DSSP installed and correctly configured
    dssp = DSSP(model, input_file)

    # Get DSSP codes (secondary structure assignments)
    dssp_codes = ''.join([dssp[key][2] for key in dssp.keys()])

    return dssp_codes

def DSSP_proportion_calculator(DSSP_string):
    H_count = DSSP_string.count("H")
    E_count = DSSP_string.count("E")
    length = len(DSSP_string)
    H_prop = (H_count /length)*100
    E_prop = (E_count /length)*100
    return([H_prop,E_prop])

def sample_lengths(df):
    
    #Select percentiles to compute and initialise df
    percentile_values = [0, 25, 50, 75, 100]
    new_df = pd.DataFrame(columns=list(df.columns) + ['category'])
    
    #Parse the percentiles, select the corresponding row and append it to the dataframe
    for i in percentile_values:
        category = f"Length percentile {i}"
        percentile = np.percentile(df['length'], i)
        idx = (df['length'] - percentile).abs().idxmin() # The minimum value of the substraction will be the closest to the percentile
        row = df.loc[[idx]].copy()
        row['category'] = category  # add new column
        new_df = pd.concat([new_df, row], ignore_index=True)

    return new_df

def sample_by_percentiles(df, col, percentiles, categories):
    rows = []
    for p, cat in zip(percentiles, categories):
        val = np.percentile(df[col], p)
        idx = (df[col] - val).abs().idxmin()
        row = df.loc[[idx]].copy()
        row['category'] = cat
        rows.append(row)
    return pd.concat(rows, ignore_index=True)

#%%Input check

parser = argparse.ArgumentParser(description="Generate summary dataframes for a project structures")
parser.add_argument('--folder', help="Folder that contains all the input structures", type=str, required= True)
args = parser.parse_args()
folder = Path(args.folder)
#folder = "inputs"

#%% Global dataframe generation

# Initialize pandas dataframe
global_df = pd.DataFrame(columns=['file_ID', 'sequence', 'length', 'DSSP', 'H%', 'E%'])

for file in folder.iterdir():
    if file.is_file():
        # Calculate length, DSSP string, proportions, etc.
        sequence = sequence_calculator(file)
        DSSP_string = DSSP_string_generator(file)
        DSSP_proportions = DSSP_proportion_calculator(DSSP_string)

        # Create a new row
        row = pd.Series([file.name, sequence, len(sequence), DSSP_string , DSSP_proportions[0], DSSP_proportions[1]], index=global_df.columns)

        # Append the row using pd.concat
        global_df = pd.concat([global_df, row.to_frame().T], ignore_index=True)

# Save the resulting DataFrame to CSV
global_df.to_csv('all_structures.csv', index=False)


#%%Selected representatives dataframe
# Step 1: Sample by length
length_df = sample_lengths(global_df)

# Step 2: Sample by alpha helix
global_df_minus_length = global_df[~global_df['file_ID'].isin(length_df['file_ID'])]
alpha_df = sample_by_percentiles(
    global_df_minus_length,
    'H%',
    [10, 50, 90],
    ['Low alpha helix', 'Medium alpha helix', 'High alpha helix']
)

# Step 3: Sample by beta sheet (excluding 0s)
global_df_minus_alpha = global_df_minus_length[~global_df_minus_length['file_ID'].isin(alpha_df['file_ID'])]
global_df_minus_alpha = global_df_minus_alpha[global_df_minus_alpha['E%'] > 0]
beta_df = sample_by_percentiles(
    global_df_minus_alpha,
    'E%',
    [10, 50, 100],
    ['Low beta sheet', 'Medium beta sheet', 'High beta sheet']
)

# Combine all and export
final_df = pd.concat([length_df, alpha_df, beta_df], ignore_index=True)
final_df.to_csv('representative_sample.csv', index=False)
