# Surrogacy_Codes

The scripts to be used:

- data_generate.R: generate Source and Target data with distribution shift
- sample_split.R: create the K-fold split
- nuisance_estimate.R: fit all nuisance estimators
- calculation_pq.R: compute P and Q (for solving beta)
- calculation_V.R: compute V
- calculation_gamma.R: compute xi and ha (for solving the NR-updated gamma)
- simulation.R: wrap all above together
 
Notes:

1. Simulation_only for reference.Rmd: notebook containing earlier setups (bootstrap, etc.), might not be helpful, kept for reference.

2. simulation(): recently drafted; still under construction. For real data runs, could use the body of the function directly.

3. Vab() and new_gamma(): newly constructed, may contain typos. Please let me know of any issues.




