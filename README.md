RMST_grp_seq_comp

This project contains the code necessary to produce the results in the manuscript titled, "Evaluating Restricted Mean Survival Time Methods in Group Sequential RCTs".

Below is a description of the programs that produce the figures, results and table used in the manuscript.

 - the program paper.figs.R is used to produce figures 1-4. First, run line 3 which sources the 'source.R' program that contains the necessary R functions for the manuscript. Then:

     - lines 7-10 produce figure 1. This figure is saved in the documentation folder and named scen1_fig.png.

     - Lines 14-48 produce figure 2. This figure is saved in the documentation folder and named scen2_fig.png.

     - Lines 53-88 produce figure 3. This figure is saved in the documentation folder and named scen3_fig.png.

      - Lines 93-128 produce figure 4. This figure is saved in the documentation folder and named scen4_fig.png.

Table 1-4 are produced by the analysis.R program. First, run the program. Then:

 - Table 1 is produced by lines 218-260.

 - Table 2 is produced by lines 300-347.

 - Table 3 is produced by lines 383-430.

 - Table 4 is produced by lines 467-514.

 To produce the results dataset that are stored in the "Results" folder, run the following programs:

  - sim1_grpseq.R
  - sim1_grpseq_type1.R
  - sim2_grpseq.R
  - sim2_grpseq_type1.R
  - sim3_grpseq.R
  - sim3_grpseq_type1.R
  - sim4_grpseq.R
  - sim4_grpseq_type1.R

The program source.R contains the necessary R function for the simulations. debug.R is an ancillary program used for debugging purposes and is not necessary for the results in the manuscript.
