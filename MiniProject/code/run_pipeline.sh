#!/bin/bash
# run_pipeline.sh
#!/bin/bash
# run_pipeline.sh
# This script runs a sequence of R scripts, compiles a LaTeX document,
# prints the reference list (if available), and performs a word count.
#
# The R scripts executed are:
#   01_data_preprocessing.R
#   02_define_and_fit_models.R
#   03_compare_models_and_akaike_weights.R
#   04_visualize_parameters.R
#   05_plot_each_id.R
#   06_compile_latex.R

# Set the working directory to the directory containing this script
WORKDIR="$(dirname "$0")"
cd "$WORKDIR" || { echo "Error: Unable to change directory to $WORKDIR"; exit 1; }

# Define the list of R script files to run
FILES=(
 "01_data_preprocessing.R"
 "02_define_and_fit_models.R"
 "03_compare_models_and_akaike_weights.R"
 "04_visualize_parameters.R"
 "05_plot_each_id.R"
 "06_compile_latex.R"
)

# Run each R script in order
for file in "${FILES[@]}"; do
  if [ -f "$file" ]; then
    echo "Running $file ..."
    Rscript "$file"
    if [ $? -ne 0 ]; then
      echo "Error running $file"
      exit 1
    fi
  else
    echo "File $file not found, skipping..."
  fi
done

# If main.tex exists, compile the LaTeX document
if [ -f main.tex ]; then
  echo "Compiling LaTeX document..."
  pdflatex -output-directory=results main.tex
  bibtex main
  pdflatex -output-directory=results main.tex
  pdflatex -output-directory=results main.tex
  if [ $? -ne 0 ]; then
    echo "Error compiling LaTeX document"
    exit 1
  fi
else
  echo "main.tex not found, skipping LaTeX compilation."
fi

# Display the reference list if the bbl file exists
if [ -f results/main.bbl ]; then
  echo "Reference list:"
  cat results/main.bbl
else
  echo "Reference list not found (results/main.bbl not present)."
fi

# Perform word count using texcount if available
if command -v texcount >/dev/null 2>&1; then
  echo "Word count for main.tex:"
  texcount -inc -total main.tex
else
  echo "texcount command not found. Please install texcount for word count statistics."
fi

echo "Pipeline executed successfully."
