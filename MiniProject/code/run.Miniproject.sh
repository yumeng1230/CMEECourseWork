#!/bin/bash

# 1. Run data preparation script
Rscript 01_data_preprocessing.R

Rscript 02_define_and_fit_models.R

Rscript 03_compare_models_and_akaike_weights.R

Rscript 04_visualize_parameters.R


Rscript 05_plot_each_id.R

Rscript 06_compile_latex.R



# 3. Run Latex compilation script
chmod +x latex.sh
./latex.sh main.tex
