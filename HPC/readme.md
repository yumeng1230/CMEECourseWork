# Project README

## Overview
This project consists of scripts and data files for clustering analysis and visualization. The primary focus is on demographic and neutral clustering, with results stored as R data files and visualizations in PNG format.

## Directory Structure
```
code/
├── demographic_cluster_script.sh        # Script for demographic clustering
├── neutral_cluster.sh                   # Script for neutral clustering
├── *.rda                                # Output data files
├── *.png                                # Visualizations
├── *.sh.e* / *.sh.o*                    # Execution logs
```

## Prerequisites
- **Operating System:** Linux/Unix-based OS with Bash support
- **Required Software:**
  - Bash
  - R (for handling .rda files)

## How to Run
### Running Demographic Clustering
```bash
bash demographic_cluster_script.sh
```
### Running Neutral Clustering
```bash
bash neutral_cluster.sh
```
Results will be saved in `.rda` files.

## Viewing Results
- **R Analysis:** Load results in R using `load("file.rda")`.
- **Visualizations:** Check PNG images in the `code/` directory.

## License
This project is released under [your preferred license]. Please see the LICENSE file for details.

## Contact
For any questions or issues, please contact the project maintainer.

