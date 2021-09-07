\name{ChIASim-package}
\alias{ChIASim-package}
\alias{ChIASim}
\docType{package}
\title{
\packageTitle{ChIASim}
}
\description{
\packageDescription{ChIASim}
}
\details{
Chromosomal 3D interactions have been long studied to investigate the genome structure and functions. The understanding of the 3D structure of the genome is much deepened by the advances of biological technologies. ChIP-PET, PLAC-seq and Hi-ChIP are three methods that provide chromatin interactions that are associated with a specific protein. Over the past decade, many statistical tools have been developed to analyze ChIP-PET, PLAC-seq and Hi-ChIP data. While most of these methods claim better performance than the others, fair simulation studies is lacking, primarily due to the lacking of bench marking simulated data sets for fair comparisons. This package is designed to fill this gap.

In this package, we provide the following two functions:  

 - chiaSim: This function simulates a set of chromatin interactions mediated by a specific protein, which may be views are the in silico counterparts of the output data from ChIA-PET, PLAC-seq or Hi-ChIP. 
 
 - convertfmt: This function allows user to siwtch the output data format between "base" and "BEDPE" without re-running the whole simulation program.
 
 
Part of the code is this package is adopted from the open source package MDM (Niu and Lin, 2015).  
  
}
\author{
\packageAuthor{ChIASim}

Maintainer: \packageMaintainer{ChIASim}
}
\references{
Niu, L and Lin, S. A Bayesian mixture model for chromatin interaction data. Statistical Applications in Genetics and Molecular Biology, vol. 14, no. 1, (2015), pp. 53-64. https://doi.org/10.1515/sagmb-2014-0029
}
