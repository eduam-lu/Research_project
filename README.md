Welcome to the github repo for my LU Msc in Bioinformatics' Research project titled "Enhancing shape directed design of protein containers’ subunits using the generative AI method RFdiffusion". Here I have made available all the code generated, all the input structures and data needed to replicate my analyses, the data produced as an output and used for visualization of results and ,lastly , all the environments, software versions and programs used in the project.

The project consists on the complementation of the Proteinshapedesign pipeline developed in "A dance with Protein Assemblies" Mads Jeppesen, 2025. This software is built to redesign an input pdb structure (identified with a potential to serve as a capsomer for self assembling capsids), optimising it for capsid formation. In this work, I have tested the viability of introducing an RF diffusion pretreatment to the input structures. RF diffusion is a denoising diffusion neural network model that is able to generate protein backbone variability from a given input structure. By generating and filtering structural variability over the inputs for ProteinShapedesign, we expected to obtained better capsomers from the original software. Results are discussed in the Report.pdf file.

Below, I break down the contents of each folder of the repo, explain how each script works briefly and how to run it; and finally, showcase external programs used in the project and how to download them
# 1_tuning_scripts
---
The scripts contained in this folder were used before the development of the improved pipeline in order to adjust some model parameters and select methods based on their performance
### rep_sample_generator.py (Representative sample generator)
The programs that require adjustments are computationally intensive and time consuming. Rather than assessing performance using all the inputs, this script selects a representative sample of a given folder of .pdb inputs and works as follows:
1. For each input file in the folder, sequence length, DSSP sequence (secondary structure prediction) and secondary structure proportions are calculated and stored in a pandas df
2. The following structures are selected: shortest and longest sequence, the sequences at the 25th, 50th, and 75th percentiles of sequence length, the structures with the lowest, highest, and 50th percentile α-helix content, and those with the lowest, highest, and 50th percentile β-sheet content.
3. The selected structures ID and their corresponding metrics are stored as a csv file

This way it is ensured that different sizes and folds will be considered when assessing performance. The script can be run as follows:
```bash
python /1_tuning_scripts/rep_sample_generator.py 
```
### partial_T_tuner.py
partial T is a parameter of the RF diffusion "partial diffusion" mode. It determines how much noise will be introduced in a given structure before RF diffusion denoises it, generating a variability of different backbones from the original one. Too little noise can lead to negligible changes, while too much can result in completely different folds. To adjust it, 3 values (10,20,30) were tested over the representative sample using the following process:
1. For each input structure, 3 backbones are generated with RF diffusion "partial diffusion", one per partial T value.
2. Sequence for each backbone is predicted using ProteinMPNN. 20 sequences are generated per backbone, lowest global_score is selected
3. Structure from sequence is predicted with Alphafold3
4. RF generated structures are aligned with their corresponding original structure and RMSD is calculated with Biopython's superimposer function
5. Metrics generated for the backbones are saved as a csv.file
Knowing the RMSD, alignments were visualised in Pymol to confirm fold similarity. After this analysis it was determined that a partial T of 20 would be used. The script can be run as follows
```bash
python /1_tuning_scripts/rep_sample_generator.py 
```
### prediction_performance.py 
Several protein structure prediction methods were considered for the project: Alphafold 3 (without MSA generation step), ESM Fold and Chai 1. To compare their performance, this script was developed and run to predict the structures of the backbones generated from the representative sample. For each input sequence structure is predicted and run time, pLDDT and RMSD with the original structure are calculated for each method.

At the light of the results, ESM Fold was selected for high throughput prediction for its outstanding speed and Chai 1 is selected for final predictions given its higher accuracy, being designed to work on single sequence as well, and for validation by predicting structure with two different architectures.

The script can be run as follows:
```bash
python /1_tuning_scripts/rep_sample_generator.py 
```
# 2_main_scripts 
---

### improved_pipeline.py and functions.py

### benchmark.py and functions_benchmark.py

### control_shapedesign.py 

# 3_auxiliary_scripts
---
- **folder_lower.py**
- **pdb_converter.py**
- **structure_relaxation.py**
- **delta_G_predictor**

# 4_visualisation scripts
---

# 5_inputs
---

# 6_outputs
---

# 7_envs
---

# External installations
---

# References
---
Mads thesis
RF diffusion
DSSP
Relaxation
