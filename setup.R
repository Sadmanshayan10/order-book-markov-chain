# setup.R: install the packages this project needs.
pkgs <- c("data.table", "Matrix", "rmarkdown", "tinytex", "knitr", "R.matlab",
          "formatR")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing)

# pdf_document output needs a LaTeX distribution. Install TinyTeX if none found.
if (!tinytex::is_tinytex() && Sys.which("pdflatex") == "") tinytex::install_tinytex()
