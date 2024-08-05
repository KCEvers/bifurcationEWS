The HPC folder contains the scripts necessary to run the analysis on a high performance cluster (HPC). The master file is run_full_analysis_detGLV.R, which calls:

1. generate_full_GLV.R
2. generate_transitions_GLV.R
3. compute_EWS_GLV.R
4. eval_performance_EWS_GLV.R
5. summarise_EWS_GLV.R

with the support of visualise_helpers.R. The last file prep_GLV.R prepares the data sets for the package.

Please note that the analysis can be executed much more efficiently, but would need to be heavily edited.

