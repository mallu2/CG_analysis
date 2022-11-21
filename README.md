These are jupyter-notebooks used to analyse the trajectories of CG simulations,
this is the rackham branch.

pipeline_trajectory_analysis contains scripts for the trajectory analysis. The data obtained from these scripts can then be analysed and visualized with the jupyter-notebooks in this folder.

Processed trajetcory data (position Z along the DNA, the angle phi around the DNA and the distance d from the DNA) can the be plotted with the notebook plotting_CG_sim.ipynb.

The diffusion can be analysed and plotted with msd_diffusion_coefficient.ipynb.

The trajectory data can also be split into 1D and 3D diffusion and into groove tracking/sliding motions on the DNA with analysis_sliding_and_hopping.ipynb.

Interaction profiles of the protein on DNA (essentially d for the protein residues) can be plotted using  interaction_profiles.ipynb.

Finally different energies obtained from the simulation and bonds formed between protein and DNA can be analysed using CG_energies_analysis.ipynb based on ythe energies written out in the simulation and scripts that analyse specific interactions.
