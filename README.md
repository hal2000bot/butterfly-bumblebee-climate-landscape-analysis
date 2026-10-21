# Bumblebee and butterfly climate and community analyses

This repository contains a plain R script combining shared preparation with three final analyses:

1. A `glmmTMB` Gaussian model of **within-locality-standardized log(1 + total abundance)**, including `habitat_type`, taxon-by-year and taxon-by-climate interactions, quadratic summer/prior-summer temperatures, cubic spring temperature, and locality-level temporal AR(1) correlation.
2. A bumblebee fourth-corner negative-binomial GLLVM.
3. A butterfly fourth-corner negative-binomial GLLVM.

## Running

Provide R packages: `dplyr`, `tidyr`, `purrr`, `tibble`, `readxl`, `stringr`, `lubridate`, `janitor`, `bombLepiSurv`, `DBI`, `glmmTMB`, and `gllvm`. Configure authentication and a database connection `con` via `humlesommerfConnect()`.



