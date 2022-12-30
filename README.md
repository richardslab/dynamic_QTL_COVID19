Author: Julian Daniel Sunday Willett, M.D.

Code contributing to dynamic eQTL investigation with respect to HGI COVID-19 outcomes.

Code involves identifying putatively causal variants with MR and then using colocalization (and other senstivity tests) to investigate assumptions of MR with the order of lines in the "main" file capturing the workflow. The process is as follows:

1. Gather all data used for MR. This was done prior to MR rather than during MR due to the long computation time and the large amount of data involved. This was a time-intensive step, particularly for the Soskic et al. data considering this was run on a number of cells and cell times across a large number of loci. I edited the code as possible to streamline it.
2. Run the MR for every saved series of data along with sensitivity testing. Tests with significant results that passed initial sensitivity testing were marked for colocalization.
3. Run the colocalization on those suspected causal MR results.
4. For strongly colocalizing (H4.PP >= 0.8) transcript:outcome analyses, run sensitivity testing to verify the single causal variant assumption.
5. Plot MR results for variants that passed MR and colocalization sensitivity testing
6. Collect data for all states for putatively causal variants (for making a heatmap plot)
7. Make the heatmap plot demonstrating the role of cell, cell stimulation, time after stimulation, and disease state on colocalization.

Package versions:
TwoSampleMR: v0.5.6
Coloc: v5.1.1
