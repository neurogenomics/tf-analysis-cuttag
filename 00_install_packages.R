req_packages <- c(
  "remotes",
  "devtools",
  "Rsamtools",
  "MACSr",
  "tidyverse"
)


install.packages("BiocManager")
BiocManager::install(req_packages, update = TRUE)

remotes::install_github("neurogenomics/EpiCompare", dependencies = TRUE)
remotes::install_github("neurogenomics/MotifPeeker", dependencies = TRUE)
