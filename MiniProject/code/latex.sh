#!/bin/bash

# Define output directory



# Compile LaTeX document to PDF with shell escape for texcount
pdflatex -shell-escape main.tex

# Run BibTeX for citations
bibtex main 

# Run pdflatex twice more to resolve references
pdflatex -shell-escape main.tex
pdflatex -shell-escape main.tex

# Move the output PDF to the output directory

# Clean auxiliary files
rm *.aux
rm *.log
rm *.bbl
rm *.blg
rm *.out
rm *.txt

echo "Compilation completed. Check $OUTPUT_DIR/$PDF_NAME for the final PDF."


