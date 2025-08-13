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
- **folder_lower.py** : this auxiliary script is designed to given a folder, it turns all the file names within it to all lowercase later. It is used within the main scripts to avoid case sensitivity issues when accesing folders or comparing structures. It can also be used as a standalone program
- **pdb_converter.py** : this script uses pymol to transform all the .cif files found within a folder into .pdb format. It is used several times within the main scripts as some prediction methods output cifs, others output directly pdbs. Given that ProteinShapedesign does not offer .cif support, .pdb format was selected as the standard.
- **structure_relaxation.py** : structures passed to Protein ShapeDesign must be in a minimal energy state, as the software evaluates the optimisation process with several rounds of minimisation. This scripts implements the Fast Relax protocol from the Rosetta documentation found here https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastRelaxMover and it needs to be run from the pyrosetta2 environmment (see section 7)
- **delta_G_predictor**: given a folder with protein structures, this script returns a dictionary with the global delta G (stability) in kcal/mol for each structure. This script is an adapted version of the .ipynb found in the repo for Cagiada et al.'s 2025 paper "Predicting absolute protein folding stability using generative models" (https://github.com/KULL-Centre/_2024_cagiada_stability/tree/main). Delta G is estimated from a regression model that correlates stability in kcal/mol with the sum of likelihood scores obtained from subjecting the sequence to the inverse folding model ESM-IF. This script needs to be run from the esm_stability environment (see section 7)

# 4_visualisation scripts
---
Two interactive python notebooks were developed to visualise benchmarking results and generate the plots for the report
**plot_monomer.ipynb** makes all the possible monomer condition comparison in terms of stability, solubility and hydropathy (and RMSD and pLDDT when possible) with the data contained in original_structures.csv (original structures), RF_only_metrics.csv (monomers treated with RF diffusion) and monomer_metrics.csv (capsomers generated with Proteinshapedesign with and without RF treatment). The notebook contains the paired plots, boxplots and statistical analyses for the comparisons. Also includes the code for the selected comparisons shown in the report.
**plot_oligomer.ipynb** compares oligomers and capsids with and without treatment. Dimers (dimer_metrics.csv), trimers (trimer_metrics.csv) and pentamers(pentamer_metrics.csv) with and without are compared in terms of pLDDT and TM score, while capsids with treatment (capsid_metrics_w.csv) and without (capsid_metrics_wo.csv) are compared in terms of interface delta delta G, SASA metrics and shape complementarity. Again, all possible boxplot, paired plots and statistical analysis of the comparisons are shown. 
The notebook also contains the selection process for "mixed oligomers", were the best predicted oligomer per sample is selected and plotted.

# 5_inputs
---
This folder contains all the data needed to replicate all the analyses.
Input pdbs contains 96 protein structures selected from the Protein Building Block database based on their previous performance in the ProteinShapeDesign pipeline, their potential for empirical testing and for their ability of predicting at least one oligomer properly with multimer prediction methods.

Symdefs contains the corresponding symmetry definition for each input pdb. A symmetry definition is the set of mathematical transformations needed to generate a whole capsid by repeating a single capsomer. Each input pdb has assigned to itself the symmetry definition of the natural capsid where the ideal shape with which they match was extracted.

Additionally, although it is technically an output from the improved shapedesign and the control, the capsids generated in the project have also been made available here: ; as they are needed to reproduce the benchmarking process
# 6_outputs
---
Here I break down each of the metrics generated at different stages of the project, in chronological order:
- all_structures.csv
- representative_sample.csv
- partial_T_tuning.csv
- prediction_methods_comparison.csv
- original_structures_def.csv
- RF_only_metrics.csv
- monomer_metrics.csv
- dimer_metrics.csv
- trimer_metrics.csv
- pentamer_metrics.csv
- capsid_metrics_w.csv
- capsid_metrics_wo.csv
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
