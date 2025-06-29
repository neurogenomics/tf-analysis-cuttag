
# nanoCUT&Tag - 17TFs - K562 cells - Analysis

**README Last Updated:** May 2025

## Authors
- **Experiment:** Isidora Gocmanac
- **Scripts:** Hiranyamaya (Hiru) Dash

This project is part of Isidora's PhD.

## Notes
- **Alignment data source:** Imperial HPC `/rds/general/user/$USER/projects/neurogenomics-lab/live/Data/cutandtag/processed_data/cuttag_18_04_2025_Isidora`
- **5'-shift peak calling method used for CUT&Tag:** [biorxiv - Rational design of peak calling parameters for TIP-seq based on pA-Tn5 insertion patterns improves predictive power](https://www.biorxiv.org/content/10.1101/2024.10.08.617149v1) 
- Default settings with control files were used for Chip-seq peak calling.  
- **Motif analysis:** [MotifPeeker R Package](https://github.com/neurogenomics/MotifPeeker)
- For our analysis, we use ChIP-seq peak files from ENCODE,
`01_peakcall_chipseq.R` was NOT used.