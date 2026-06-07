# Overview
HOMED (Hierarchically Optimized Methylation Deconvolution) is a reference-based framework for estimating cell-type composition from bulk DNA methylation (DNAm) data. Unlike conventional methylation deconvolution approaches that treat cell types as flat entities, HOMED explicitly incorporates hierarchical relationships among related cell populations to improve resolution of closely related subtypes.

The framework combines purified FACS/MACS-derived methylation references with pseudo-ground-truth cell proportions generated through scRNA-seq-guided deconvolution of paired bulk RNA-seq samples. HOMED then performs hierarchical differential methylation analysis and iterative IDOL-based optimization to identify CpG libraries that maximize deconvolution accuracy. The resulting hierarchical reference matrices can be applied to new bulk methylation datasets to estimate both major cell lineages and finer cellular subtypes while maintaining biological consistency across hierarchical levels.

HOMED is designed for tissue-specific methylation deconvolution and can be applied to a wide range of tissues where purified methylation references and single-cell transcriptomic atlases are available.

A complete PBMC tutorial is available here: 

- [PBMC Vignette](https://github.com/yhdu36/HOMED/blob/main/vignettes/HOMED_PBMC.html)

# Installation

```r
library(devtools)

install_github("yhdu36/HOMED")
```

# References
Andrews SV, Ladd-Acosta C, Feinberg AP, Hansen KD, Fallin MD (2016). 'Gap hunting' to characterize clustered probe signals in Illumina methylation array data. Epigenetics & Chromatin, 9, 56. doi:10.1186/s13072-016-0107-z.

Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, Irizarry RA (2014). Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA Methylation microarrays. Bioinformatics, 30(10), 13631369. doi:10.1093/bioinformatics/btu049.

Du Y, Benny PA, Lahiri S, AlAkwaa FM, Huang Q, Liu Y, Lassiter CB, Astern J, Riel J, Garmire LX (2026). Placental Molecular Subtypes of Severe Preeclampsia Reveal Divergent Aging Trajectories and Fetal Growth Outcomes. medRxiv. https://doi.org/10.64898/2026.06.02.26354756.

EpiDISH. n.d. Bioconductor. Accessed June 7, 2026. http://bioconductor.org/packages/EpiDISH/.

FlowSorted.BloodExtended.EPIC: A New Extended Cell Deconvolution for Peripheral Blood. n.d. GitHub. Accessed June 7, 2026. https://github.com/immunomethylomics/FlowSorted.BloodExtended.EPIC.

Fortin J, Hansen KD (2015). Reconstructing A/B compartments as revealed by Hi-C using long-range correlations in epigenetic data. Genome Biology, 16, 180. doi:10.1186/s13059-015-0741-y.

Fortin J, Labbe A, Lemire M, Zanke BW, Hudson TJ, Fertig EJ, Greenwood CM, Hansen KD (2014). Functional normalization of 450k methylation array data improves replication in large cancer studies. Genome Biology, 15(12), 503. doi:10.1186/s13059-014-0503-2.

Fortin J, Triche TJ, Hansen KD (2017). Preprocessing, normalization and integration of the Illumina HumanMethylationEPIC array with minfi. Bioinformatics, 33(4). doi:10.1093/bioinformatics/btw691.

IDOL: IDentifying Optimal DNA Methylation Libraries (IDOL). n.d. GitHub. Accessed June 7, 2026. https://github.com/immunomethylomics/IDOL.

Immunomethylomics/IDOL Source: R/IDOLoptimize.R. n.d. Accessed June 7, 2026. https://rdrr.io/github/immunomethylomics/IDOL/src/R/IDOLoptimize.R.

Maksimovic J, Gordon L, Oshlack A (2012). SWAN: Subset quantile Within-Array Normalization for Illumina Infinium HumanMethylation450 BeadChips. Genome Biology, 13(6), R44. doi:10.1186/gb-2012-13-6-r44.

R/estimateCellCounts2.R at Devel · Immunomethylomics/FlowSorted.Blood.EPIC. n.d. GitHub. Accessed June 7, 2026. https://github.com/immunomethylomics/FlowSorted.Blood.EPIC/blob/devel/R/estimateCellCounts2.R.

Salas LA, Zhang Z, Koestler DC, Butler RA, Hansen HM, Molinaro AM, Wiencke JK, Kelsey KT, Christensen BC (2022). Enhanced Cell Deconvolution of Peripheral Blood Using DNA Methylation for High-Resolution Immune Profiling. Nature Communications, 13(1), 761.

Triche TJ, Weisenberger DJ, Van Den Berg D, Laird PW, Siegmund KD (2013). Low-level processing of Illumina Infinium DNA Methylation BeadArrays. Nucleic Acids Research, 41(7), e90. doi:10.1093/nar/gkt090.



