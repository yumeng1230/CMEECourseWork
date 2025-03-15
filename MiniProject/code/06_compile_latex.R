# compile_latex.R
# This script checks for the existence of a "main.log" file,
# removes it if present, ensures the current directory is writable,
# and then calls pdflatex to compile the LaTeX document "main.tex".

# Check if the "main.log" file exists. If it does, remove it.
if (file.exists("main.log")) {
  file.remove("main.log")
  message("Old main.log file removed.")
} else {
  message("No main.log file found.")
}

# Ensure the current directory has write permission.
# (This command may require appropriate permissions on your system.)
system("chmod +w .")

# Call pdflatex to compile the LaTeX document "main.tex".
# Make sure that "main.tex" exists in the current working directory.
compile_status <- system("pdflatex main.tex")

# Check the exit status. 0 means success.
if (compile_status == 0) {
  message("LaTeX compilation succeeded.")
} else {
  message("LaTeX compilation failed. Please check the log for details.")
}
