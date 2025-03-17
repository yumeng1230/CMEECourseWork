# Project README

## Overview
This project focuses on processing growth data, fitting various growth models, comparing them using information criteria, visualizing estimated parameters, and compiling a LaTeX report. The workflow is organized into six R scripts.

## Directory Structure
```
MiniProject/
├── code/
│   ├── 01_data_preprocessing.R              # Preprocess raw growth data
│   ├── 02_define_and_fit_models.R           # Define and fit nonlinear growth models
│   ├── 03_compare_models_and_akaike_weights.R # Compare models using AICc and BIC
│   ├── 04_visualize_parameters.R            # Visualize estimated model parameters
│   ├── 05_plot_each_id.R                    # Generate model plots for each dataset
│   ├── 06_compile_latex.R                   # Compile LaTeX report
│   ├── run_pipeline.sh                      # Script to run the entire workflow
│   ├── main.tex                             # LaTeX document for the report
│   ├── references.bib                        # Bibliography for LaTeX document
│   ├── wordcount.txt                         # Word count log for the report
├── data/
│   ├── LogisticGrowthData.csv               # Raw growth data
│   ├── modified_growth_data.csv             # Preprocessed growth data
├── results/
│   ├── *.rda                                # Model output files
│   ├── *.png                                # Visualization outputs
│   ├── .gitkeep                             # Placeholder file for results directory
├── readme.md                                # Project documentation
```

## Prerequisites
- **Operating System:** Windows, Linux, or macOS
- **Required Software:**
  - R (version 4.2.2 recommended)
  - LaTeX (for compiling the report)
  - Bash (for running shell scripts)
- **R Packages:**
  - dplyr, tidyr, ggplot2, minpack.lm, segmented, MuMIn, readr, broom

## How to Run
### Running the Full Pipeline
```bash
bash run_pipeline.sh
```

## Viewing Results
- **R Analysis:** Load `.rda` results in R using `load("file.rda")`.
- **Visualizations:** Generated `.png` files are available in the `results/` directory.
- **LaTeX Report:** The compiled PDF report summarizes the entire analysis.

## License
This project is released under [your preferred license]. Please see the LICENSE file for details.

## Contact
For any questions or issues, please contact the project maintainer.



Project Documentation: Growth Curve Model Analysis Pipeline

This project focuses on processing growth data, fitting various growth models, comparing them using information criteria, visualizing estimated parameters, and compiling a LaTeX report. The workflow is organized into six R scripts. Below is an overview of each script including its purpose, dependencies, runtime environment, and usage instructions.

Environment and Setup

    R Version: 4.2.2
    Operating System: Compatible with Windows, Linux, or macOS.
    Directory Structure:
        Raw and cleaned data are stored in ../data/.
        Results and plots are saved in ../results/.
        LaTeX document (main.tex) is located in the working directory for compilation.


***1. 01_data_preprocessing.R
Overview

    Purpose:
        Read raw data from LogisticGrowthData.csv.

        Clean the data by filtering records (e.g., Time >= 0 and PopBio > 0) and computing log(PopBio).

        Generate a unique ID for each dataset by concatenating key fields.

        Save the cleaned data to ../data/modified_growth_data.csv.

        (Optional) Fit quadratic and cubic polynomial models for each dataset and save model comparison metrics (AICc and BIC) 
        to ../results/model_comparison_linear.csv.

Dependencies

    dplyr, tidyr: Data manipulation and cleaning.
    ggplot2: Data visualization (for plotting and messages).
    minpack.lm: Nonlinear least squares fitting (nlsLM).
    segmented: (Optional) For segmented regression.
    MuMIn: AICc calculation.
    readr: Reading and writing CSV files.
    broom: Tidying model outputs.

Usage

    Place the raw data file (LogisticGrowthData.csv) in the ../data/ directory.
    Run the script. It generates the cleaned data file and, optionally, the polynomial model comparison file.

*****2. 02_define_and_fit_models.R
Overview

    Purpose:
        Read the cleaned data from ../data/modified_growth_data.csv.

        Define four nonlinear growth model formulas:
            Logistic Model
            Gompertz Model
            Baranyi Model
            Three-Phase Linear Model

        Implement multi-start fitting functions using nlsLM to search for optimal parameters with random initial values.

        These functions are used later for model comparison and visualization.

Dependencies

    dplyr, readr: Data handling.
    minpack.lm: Nonlinear model fitting.
    MuMIn: AICc calculation.
    broom: Tidying model outputs.

Usage

    Ensure that ../data/modified_growth_data.csv exists.
    Source or run this script to load the model definitions and fitting functions. It is typically sourced by subsequent scripts.

*******3. 03_compare_models_and_akaike_weights.R
Overview

    Purpose:
        Read the cleaned data and polynomial model comparison results.

        Strictly fit four nonlinear models (Logistic, Gompertz, Baranyi, and Three-Phase) for each dataset using the functions defined in the previous script.

        Merge the linear (polynomial) and nonlinear model results.

        Select the best model for each dataset based on AICc and BIC.

        Calculate Akaike weights for the fitted models.

        Generate comparison plots:
            A bar plot showing best model counts (based on AICc and BIC).
            A boxplot of the Akaike weight distributions.

Dependencies

    dplyr, readr: Data manipulation.
    MuMIn: AICc and BIC calculations.
    ggplot2: Plotting.
    broom: Tidying model outputs.

Usage

    Run this script after ensuring that the outputs from previous scripts are available. The script writes CSV files and saves plots (e.g., model_selection_comparison.png, akaike_weights_boxplot.png) in the ../results/ directory.

*****4. 04_visualize_parameters.R
Overview

    Purpose:
        Extract parameter estimates from the four nonlinear models fitted for each dataset.

        Tidy and combine these parameters into one data frame.

        Visualize the distribution of each parameter (e.g., r_max, K, N0, t_lag) across different models using boxplots.

Dependencies

    dplyr, readr: Data handling.
    minpack.lm, MuMIn, broom: Model fitting and tidying.
    ggplot2: Plotting.

Usage

    Run this script to generate a CSV file (nonlinear_model_parameters.csv) containing parameter estimates.
    Parameter-specific boxplots are saved in the ../results/ directory (e.g., nonlinear_parameter_comparison_r_max.png).

*****5. 05_plot_each_id.R
Overview

    Purpose:
        For each dataset (ID), fit six models (Quadratic, Cubic, Logistic, Gompertz, Baranyi, and Three-Phase) while allowing partial model fitting (i.e., if some models fail, plot those that succeed).

        If no model fits are successful, plot only the raw data.

        Otherwise, plot the raw data along with fitted lines from all successful models.

        Save the output plots for each ID to the ../results/ directory.

Dependencies

    dplyr, ggplot2, readr: Data manipulation and visualization.
    02_define_and_fit_models.R: Provides model definitions and fitting functions.

Usage

    Run the script after ensuring that ../data/modified_growth_data.csv is available.
    The script iterates over each unique ID and saves corresponding plots (e.g., 6models_ID_1.png) in the results folder.

*****6. 06_compile_latex.R
Overview

    Purpose:
        Check for an existing main.log file and remove it if present.

        Ensure the current directory is writable.

        Compile the LaTeX document main.tex using the pdflatex command.

        Provide feedback on whether the compilation succeeded or failed.

Dependencies

    Base R functions: For file handling and system commands.

Usage

    Make sure that the LaTeX document main.tex is in the current working directory.
    Run this script. It will remove any old log file and call pdflatex to compile the document.
    Review the console messages to confirm success or troubleshoot any errors.

*****7. LaTeX Document: main.tex
Overview

    Purpose:
        The main.tex file contains the LaTeX code for compiling the final report of the analysis.

        It integrates figures, tables, and summaries generated from the R scripts.

*****8. Pipeline Script: run_pipeline.sh
Overview

    Purpose:
        The run_pipeline.sh script is intended to automate the execution of the entire analysis pipeline.

