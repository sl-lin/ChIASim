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
Chromosomal 3D interactions have been long studied to investigate the genome structure
and function. The understanding of 3D structure of genome is much deepened by the
advances of biological technologies that are designed to depict that. ChIP-PET and Hi-ChIP are two
methods that provide the chromatin interactions that are associated with a specific protein.
As researches on analyzing ChIA-PET and Hi-ChIP data accumulated over the past
decade, many corresponding statistical tools have been developed. While most of these
methods claiming better performance than the rest, a fair simulation study is lacking. In
particular, a bench marking simulated data for comparison is lacking. This package is designed to fill such blank that generates 
such simulated data.

In this package, we provide the following function:  

 - chiaSim: This function simulates a set of data of chromosomal 3D interactions mediated by a specific protein which meant to mimicing the output data of ChIA-PET or Hi-Chip.  
 - convertfmt: This function that allows user to siwtch the output data format between "base" and "BEDPE" without re-running the whole simulation program.
  
}
\author{
\packageAuthor{ChIASim}

Maintainer: \packageMaintainer{ChIASim}
}
\references{
Mumbach, M., Rubin, A., Flynn, R. et al. HiChIP: efficient and sensitive analysis of protein-directed genome architecture. Nat Methods 13, 919-922 (2016). 

Fudenberg, Geoffrey, and Leonid A Mirny. "Higher-order chromatin structure: bridging physics and biology." Current opinion in genetics & development vol. 22,2 (2012).

}
\note{
We thank Niu Liang for his work in R package MDM. Our work is inspaired by his work and also adopted part of his code. 
}
