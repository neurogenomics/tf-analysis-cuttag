req_packages <- c(
  "Rsamtools",
  "MACSr",
  "tidyverse",
  "MotifPeeker",
  "EpiCompare"
)


install.packages("BiocManager")
BiocManager::install(req_packages)
