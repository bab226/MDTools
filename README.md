# MDTools

Scripts using MDAnalysis and MDTraj to analyze MD simulations.

1. dssp_rg_and_ree_mdtraj.py has the code to calculate secondary structure, radius of gyration, and end-to-end-distance, as well as other biophysical metrics (meant to be run in interactive python env).
2. contact_prob_PDF.py has the code to calculate a probability density function of ligand-protein binding on an amino-acid level using the COM of the ligand and the COM of each residue on the protein. Contacts are signifigant when values drop below the null values (dotted lines in examples).
3. plotter_dssp_and_contacts.py plots and formats the outputs of the above scripts (1 and 2). Data is loaded from a compressed "pickle" format. Plots can be customized for different proteins. 
4. Mdtraj_helper.py has support functions for platting and calculating DSSP

For more information about use-case and applacations, please refer to: B. A. Bogin and Z. A. Levine, “Drugging Disordered Proteins by Conformational Selection to 
Inform Therapeutic Intervention” bioRxiv July 6, 2024, p 2024.07.03.601611 (In submission, Journal of Chemical Theory and Computation). https://doi.org/10.1101/2024.07.03.601611.



- 👋 Hi, I’m @bab226
- 👀 I am a computational biophysics researcher. I am interested in simulating molecular motions for single-molecule experiments, predicting various biophysical properties, and generating unique visualizations for complex data analysis.
- 🌱 I’m currently working on developing our in-house analysis scripts for MD simulations and FCS analysis.
- 📫 How to reach me: boginb1@gmail.com, https://www.linkedin.com/in/bryan-bogin150/ 

<!---
bab226/bab226 is a ✨ special ✨ repository because its `README.md` (this file) appears on your GitHub profile.
You can click the Preview link to take a look at your changes.
--->
