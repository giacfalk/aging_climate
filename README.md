# Replication code and data for: Aging in a warming world: global  projections of cumulative and acute heat exposure of older adults
By Giacomo Falchetta, Enrica De Cian, Ian Sue Wing and Deborah Carr

[![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

Software requirements:
- R v4.3+: https://cran.r-project.org/bin/windows/base/
- RStudio: v2023.06.0+: https://posit.co/download/rstudio-desktop/
- Package dependencies: raster, sf, tidyverse, rasterVis, rgdal, maptools, pbapply, terra, knitr, kableExtra, modelsummary, openxlsx, xtable, ggforce, maptools, weights, spatstat, rworldmap, scales, patchwork, stars, viridis

To replicate the analysis:
- Download input data from https://doi.org/10.5281/zenodo.8409700
- Download all the 1km age and gender-stratified global population counts rasters from the following WorldPop page https://hub.worldpop.org/geodata/summary?id=24798
- Clone this code repository
- Run the "project_pop.R" script to generate gridded age-stratified population data for each SSP scenario
- Run "projections_exposure_m.R" to quantify heat exposure and generate the figures and tables reported in the paper

To process the data and run succesfully, the script requires a computer with at least 32GB RAM. The running time varies based on CPU characteristics, but a runtime of at least 2 hours should be expected to generate all the output data, figures, and tables. All output files are saved in the working directory. 

Manuscript under peer review. Upon publication, a link to the paper will be made available at this repository.

___

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
