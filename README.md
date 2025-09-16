# PenalizedLR_ConditionalLR

R code available for penalized conditional logistic regression for matched case-control data

- clogitf.R — Firth’s penalized conditional likelihood regression
- logFmatched_v3.R — log-F(m,m)-penalized conditional likelihood regression via data augmentation
- DES_application.R — Application on the DED dataset (see reference in DESdesc.txt)
- parent_child_trio_application.R — Application on a case-parent trio study
- Simulation studies
  - simdata.R
  - simstudyBinExp.R — Simulate binary exposure
  - simstudyCtsExp.R — Simulate continuous exposure
  - simSummary.R — Summarize simulation results
