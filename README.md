Welcome to the github repo for my LU Msc in Bioinformatics' Research project titled "Enhancing shape directed design of protein containers’ subunits using the generative AI method RFdiffusion". Here I have made available all the code generated, all the input structures and data needed to replicate my analyses, the data produced as an output and used for visualization of results and ,lastly , all the environments, software versions and programs used in the project.

The project consists on the complementation of the Proteinshapedesign pipeline developed in "A dance with Protein Assemblies" Mads Jeppesen, 2025. This software is built to redesign an input pdb structure (identified with a potential to serve as a capsomer for self assembling capsids), optimising it for capsid formation. In this work, I have tested the viability of introducing an RF diffusion pretreatment to the input structures. RF diffusion is a denoising diffusion neural network model that is able to generate protein backbone variability from a given input structure. By generating and filtering structural variability over the inputs for ProteinShapedesign, we expected to obtained better capsomers from the original software. Results are discussed in the Report.pdf file.

Below, I break down the contents of each folder of the repo, explain how each script works briefly and how to run it; and finally, showcase external programs used in the project and how to download them
# 1_tuning_scripts
---

### rep_sample_generator.py (Representative sample generator)

### partial_T_tuner.py

### prediction_performance.py 

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
Relaxation
