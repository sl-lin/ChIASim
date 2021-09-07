# ChIASim
ChIA-Sim (Lou and Lin, 2021) is an in silico experimental protocol built on a realistic simulation scheme to mimic the biological experimental protocols of ChIA-PET, Hi-ChIP, or PLAC-seq, so that data simulated represent realistic specific-protein-mediated long-range interaction data.

# Installation 
ChIASim can be installed from GitHub using following commands:
```
install.packages("devtools")
devtools::install_git("https://github.com/sl-lin/ChIASim")
devtools::install_git("https://github.com/sy-lou/ChIASim")
```
ChIASim depends on other packages, i.e., parallel, BSgenome.Hsapiens.UCSC.hg19, and GenomicFeatures. They need to be installed, one may use the commands below, in advance to run ChIASim.  
```
  if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) BiocManager::install("GenomicFeatures")
```  

# Example
```
library(ChIASim)
data("POL2")
data("TSS")
output.base<- chiaSim(n.cells = 500, N.E=100, N.P=100, N.nE=100, N.nP=100,  mc.cores.user = 1)
head(output.base)
````

You may refer to the vignette for more details.

