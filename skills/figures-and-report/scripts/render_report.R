# =============================================================================
# render_report.R  —  Build the PDF report from report_template.Rmd
# Author: Yue Hao
# -----------------------------------------------------------------------------
# Requires pandoc + LaTeX. In RStudio, pandoc is bundled (this script auto-detects
# it on macOS). Install LaTeX once with:  Rscript -e 'tinytex::install_tinytex()'
#
# Run:  Rscript render_report.R
# =============================================================================

source("config.R")

# Locate pandoc (RStudio bundles it; detect common macOS locations).
if (!rmarkdown::pandoc_available()) {
  cand <- c("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64",
            "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64",
            "/Applications/RStudio.app/Contents/Resources/app/bin/pandoc",
            Sys.getenv("RSTUDIO_PANDOC"))
  hit <- cand[file.exists(file.path(cand, "pandoc"))][1]
  if (!is.na(hit)) Sys.setenv(RSTUDIO_PANDOC = hit)
  rmarkdown::find_pandoc()
}
if (!rmarkdown::pandoc_available())
  stop("pandoc not found. Open report_template.Rmd in RStudio and Knit, or install pandoc.")

# Make a user TinyTeX visible if present.
tt <- file.path(Sys.getenv("HOME"), "Library/TinyTeX/bin")
if (dir.exists(tt)) {
  bin <- list.dirs(tt, recursive = FALSE)[1]
  if (!is.na(bin)) Sys.setenv(PATH = paste(bin, Sys.getenv("PATH"), sep = ":"))
}

rmarkdown::render("report_template.Rmd",
                  output_file = paste0(config$project_name, "_report.pdf"),
                  output_dir  = config$output_dir, quiet = FALSE)
message("Wrote report: ", file.path(config$output_dir, paste0(config$project_name, "_report.pdf")))
