\name{chiaSim}
\alias{chiaSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{chiaSim
}
\description{
This function simulates a set of chromatin interactions mediated by a specific protein, which may be views are the in silico counterparts of the output data from ChIA-PET, PLAC-seq or Hi-ChIP.   }
\usage{
chiaSim(TFBSfile="ERa", TSSfile=NULL, mean.frag.length = 250, 
	n.cells = 1E5, lp = 0.8, seqlength = 50, N.E=400, N.P=400, N.nE=400, N.nP=400,	REp = "AAGCTT", mc.cores.user = 1, beta = 1, gmclass = 1, outputformat= "base", seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TFBSfile}{
TFBSfile is a dataframe containing transcription factor binding sites (TFBS) information for a specific protein of interest. It should contain at least three columns with the first three columns providing information on the chromosome, start site, and end site of each TFBS for this protein. The ChIA-Sim package has two built-in TFBS data sets available, "ERa" and "POL2". If interested in a different protein, a user-supplied file may be used as long as it contains the required information.    }
  \item{TSSfile}{
TSSfile is a dataframe containing transcription start sites (TSS) information. Users do not need to provide a TSS file as the package has the TSS information of hg19 from Refseq (Pruitt et al, 2014). However, if users do provide a TSS file, it should follow the required format. }
  \item{mean.frag.length}{
mean.frag.length is the mean chromosome fragment length for sonication simulation. To simulate the sonication step in a real ChIA-PET experiment that cuts choromosomes into pieces, we need to specify the average length of fragments. The default is set to be 250bp.     }
  \item{n.cells}{
n.cells is the simulated number of cells which may be interpreted as the sequencing depth. By adjusting n.cells, different numbers of pairs can result. The default is set to be 10^5.}
  \item{lp}{
lp is the ligation probability. This is the probability that two selected segments to ligate. It is set to 0.8 by default. }
  \item{seqlength}{
seqlength is the length of the fragments to be seqenced. ChIP-PET seqences the two end segments of an interacting pair.  The typical length of PET seqencing is 50-100bp (Fullwood et al, 2009), hence it is set to be 50bp as the default seqencing length.  }
 
\item{N.E}{ is the number of TFBS sites randomly selected from the provided TFBS file for forming true pairs for simulation. The selected TFBS sites are allocated to each choromosome proportional to their lengths. It is set to 400 by default.}
  
\item{N.P}{ is the number of TSS sites randomly selected from the provided TSS file for forming true pairs for simulation. The selected TSS sites are allocated to each choromosome proportional to their lengths. It is set to 400 by default.} 

\item{N.nE}{ is the number of non-TFBS sites randomly selected across the whole genome for forming the false pairs for simulation. The selected non-TFBS sites are allocated to each chromosome proportional to their lengths. It is set to 400 by default.}
  
\item{N.nP}{ is the number of non-TSS sites randomly selected across the whole genome for forming the false pairs for simulation. The selected non-TSS sites are allocated to each choromosome proportional to their lengths. It is set to 400 by default.} 


 \item{REp}{ is the pattern of restriction enzyme. It is set to AAGCTT by default, which is the pattern of HindIII. }
  \item{mc.cores.user}{
is the number of cores used to do parallel computing when multi cores are available at the user's end.   }
  \item{beta}{
beta is the power parameter. Power parameter is the power based on the widely believed power law between interaction frequency and 1/log(genomic distance). It is usually set to 1 for human (Fudenberg and Mirny, 2012), and is the default value.}
  \item{gmclass}{
gmclass is a flagging variable on whether to identify self-ligation pairs. When it is set to 2, all resulting interaction pairs will be clustered into 2 groups according to their genomic distances; the cluster with smaller genomic distance will be treated as self-ligation and removed. When self-ligation is not a major issue, it is set to 1, the default value. 
}
  \item{outputformat}{
is the format to output the simulated data. Users can choose either "base" or "BEDPE". We provide these two formats becasue they are the two most popular data formats for real data. For example, the K562 CTCF ChIA-PET data set with accession number GSM970216 uses the BEDPE format, while the K562 POL2 ChIA-PET data set with accession number GSM832465 uses a 7-column base format. Defualt is set to the base format. We provide two examples in the following, one with the base format and the other with the BEDPE format.   
}
  \item{seed}{
is the random seed number. It may be set for reproducibility purpose. Otherwise, it should be set as null.  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A simulated ChIA-PET data set in either the "base" or "BEDPE" format.
}
\references{
%% ~put references to the literature/web site here ~
Mumbach, M., Rubin, A., Flynn, R. et al. HiChIP: efficient and sensitive analysis of protein-directed genome architecture. Nat Methods 13, 919-922 (2016). 

Fudenberg, Geoffrey, and Leonid A Mirny. "Higher-order chromatin structure: bridging physics and biology." Current opinion in genetics & development vol. 22,2 (2012).

Pruitt KD, Brown GR, Hiatt SM, Thibaud-Nissen F, Astashyn A, Ermolaeva O, Farrell CM, Hart J, Landrum MJ, McGarvey KM, Murphy MR, O'Leary NA, Pujar S, Rajput B, Rangwala SH, Riddick LD, Shkeda A, Sun H, Tamez P, Tully RE, Wallin C, Webb D, Weber J, Wu W, DiCuccio M, Kitts P, Maglott DR, Murphy TD, Ostell JM. RefSeq: an update on mammalian reference sequences. Nucleic Acids Res. 2014 Jan;42(Database issue):D756-63. doi: 10.1093/nar/gkt1114. Epub 2013 Nov 19. PMID: 24259432; PMCID: PMC3965018.

Pruitt, K. D., G. R. Brown, S. M. Hiatt, F. Thibaud-Nissen, A. Astashyn, O. Ermolaeva, C. M. Farrell, J. Hart, M. J. Landrum, K. M. McGarvey, M. R. Murphy, N. A. O'Leary, S. Pujar, B. Rajput, S. H. Rangwala, L. D. Riddick, A. Shkeda, H. Sun, P. Tamez, R. E. Tully, C. Wallin, D. Webb, J. Weber, W. Wu, M. DiCuccio, P. Kitts, D. R. Maglott, T. D. Murphy and J. M. Ostell (2014): "Refseq: an update on mammalian reference sequences," Nucleic Acids Res., 42, D756-D763. 

Fullwood MJ, Wei CL, Liu ET, Ruan Y. Next-generation DNA sequencing of paired-end tags (PET) for transcriptome and genome analyses. Genome Res. 2009;19(4):521-532. doi:10.1101/gr.074906.107
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("POL2")
data("TSS")
output.base<- chiaSim(n.cells = 1E3, N.E=100, N.P=100, N.nE=100, N.nP=100, mc.cores.user = 48)
head(output.base)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
