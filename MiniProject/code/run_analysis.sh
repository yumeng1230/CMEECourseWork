#!/bin/bash

# Set the working directory
WORKDIR="$(dirname "$0")"
cd "$WORKDIR" || exit 1

# Create results directory if it doesn't exist
mkdir -p results

# Run DataClean.R for data cleaning
echo "Running DataClean.R..."
Rscript DataClean.R
if [ $? -ne 0 ]; then
    echo "Error running DataClean.R"
    exit 1
fi

# Run Main.R (make sure this file exists)
if [ -f Main.R ]; then
    echo "Running Main.R..."
    Rscript Main.R
    if [ $? -ne 0 ]; then
        echo "Error running Main.R"
        exit 1
    fi
else
    echo "Warning: Main.R not found, skipping..."
fi

# Run Parameter.R for extracting parameters
echo "Running Parameter.R..."
Rscript Parameter.R
if [ $? -ne 0 ]; then
    echo "Error running Parameter.R"
    exit 1
fi

# Run Parameter_visualization.R for parameter visualization
echo "Running Parameter_visualization.R..."
Rscript Parameter_visualization.R
if [ $? -ne 0 ]; then
    echo "Error running Parameter_visualization.R"
    exit 1
fi

# Generate LaTeX report (assuming main.tex is the main file)
if [ -f Main.tex ]; then
    echo "Compiling LaTeX report..."
    pdflatex -output-directory=results Main.tex
    bibtex results/Main  # Run BibTeX to process references
    pdflatex -output-directory=results Main.tex
    pdflatex -output-directory=results Main.tex  # Run twice to update references
    if [ $? -ne 0 ]; then
        echo "Error compiling LaTeX report"
        exit 1
    fi
else
    echo "Warning: main.tex not found, skipping LaTeX compilation..."
fi
