Author: Julian Daniel Sunday Willett, M.D.

Code contributing to dynamic eQTL investigation with respect to HGI COVID-19 outcomes.

Code involves identifying putatively causal variants with MR and then using colocalization (and other senstivity tests) to investigate assumptions of MR with the order of lines in the "main" file capturing the workflow. The process is as follows:

1. Use the "Faster Pipeline" files to gather the necessary preliminary files that will be piped into MR. This was done (using a 'lapply' based method, run in parallel by locus using a batch scheduler) to greatly reduce the computational time (that took >3-4 days on all Soskic data when not done this way).
2. Run the MR for every saved series of data along with sensitivity testing. Tests with significant results that passed initial sensitivity testing were marked for colocalization.
3. Run the colocalization on those suspected causal MR results.
4. For strongly colocalizing (H4.PP >= 0.8) transcript:outcome analyses, run sensitivity testing to verify the single causal variant assumption.
5. Plot MR results for variants that passed MR and colocalization sensitivity testing
6. Collect data for all states for putatively causal variants (for making a heatmap plot)
7. Make the heatmap plot demonstrating the role of cell, cell stimulation, time after stimulation, and disease state on colocalization.

The code was organized to run these steps roughly in order. If you experience any challenges, or have questions, please file an "issue" on this github page and I will get back to you.

Package versions:
TwoSampleMR: v0.5.6
Coloc: v5.1.1
