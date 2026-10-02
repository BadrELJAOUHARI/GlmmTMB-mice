# glmmTMB function in mice simulation study

This project tests a custom `mice.impute.2l.glmmTMB` function for multilevel multiple imputation and compares it with the existing `mice` `2l.lmer` method.

The simulation focuses on two missing-data settings:
- age and study hours missing (predictors);
- score missing (outcome).

The main simulation is set to:

```r
nsim <- 100
runtime_reps <- 20
```

This is the version to use for the final results. It can take some time to run. If you only want to check that everything works first, you can temporarily use something smaller, for example:

```r
nsim <- 10
runtime_reps <- 2
```

Then change the values back before running the final simulation.

## Files

- `01_mice_impute_2l_glmmTMB.R` : the custom `2l.glmmTMB` imputation function.
- `02_simulation_study.R` : data generation, missingness, simulation, model fitting, tables, and runtime benchmark.
- `03_plots.R` :  creates and saves the final figures.
- `run_all.R` : sources the files in the correct order and runs the full project.

The results are saved in `/simulation_results`, including the raw simulation results, tables, runtime results, and graphs.
