#' Simulate ChIA-PET Data
#'
#' This function simulates a set of data of chromosomal 3D interactions mediated by a specific protein which meant
#' to mimicing the output data of ChIA-PET or Hi-Chip.
#'
#' @usage data(ERa) A TFBS data of ERa.
#' @usage data(POL2) A TFBS data of POL2.
#' @usage data(TSS) A TSS data of hg19 obtained from Refseq.
#' @param TFBSfile is a dataframe contain TFBS information for a specific protein of interest. It should contain at least three columns with first three columns tell the chromosome, start site and end site of each TFBS for this protein (see example below). Package have two TFBS data sets available, "ERa" or "POL2". Otherwise, users need to provide TFBS file with required information.
#' @param TSSfile is a dataframe contain TSS information. Users do not have to provide TSS file as the package have the TSS information of hg19 that obtained from Refseq. However if users do provide a TSS file, it should follow the format as show in the example.
#' @param mean.frag.length is the mean chromosome fragment length for sonication simulation. To simulate the sonication step in the real ChIA-PET experiment that cut choromosomes into pieces, the mean chromosome fragment length need to be set to determine the average length of fragments. It is set to 250bp by defaut.
#' @param n.cells is the simulated number of cells which can be seen as the size of sampling pool for interactions or sequencing depth. ChIA-PET usually sequcence 10^5 cells and Hi-ChIP usually sequence 10^6 cells.
#' @param lp is the ligation probability which is the probability of two paired fragment to ligate and form an interaction.
#' @param seqlength is the length of seqencing read. ChIP-PET can only seqence two segments on the ends of an interacting pair. Usually the length of such seqencing is about 50bp, hence it is set to be the default seqencing length.
#' @param totalnumSites is the total number of TFBS and TSS sites selected, respectively, for interacton. Seleted TFBS and TSS are allocated to each choromosome according to their legnths. Sites for forming the true pairs are randomly selected from the provided TFBS and TSS file for true pairs, while sites for forming the true pairs are randomly selected across the whole genome.
#' @param mc.cores.user Number of cores used to do parallel computing
#' @param beta is the power parameter. Power parameter is the power based on the widely believed power law between interaction frequency and 1/log(genomic distance). Usually set to 1/2 for human, i.e., interaction frequency = (1/log(genomic distance))^(1/2).
#' @param gmclass is the variable of whether to identify self-ligation pairs. When it is set to 2, the pairs will be clusted into 2 clusters according to genomic distance and remove the cluster with smaller genomic distance (self-ligation).   and remove the cluster with smaller genomic distance (self-ligation)
#' @param outputformat Format of output simulated data. Users can choose either "base" or "BEDPE".
#' @import parallel
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import GenomicFeatures
#' @return A set of simulated ChIA-PET data.
#' @export


#Changed the sampling prob to according to distance
#Use hindIII to find the cutting segment
###### load parallel library
#library("parallel")


###### load BSgenome.Hsapiens package
#library("BSgenome.Hsapiens.UCSC.hg19")


##### set number of samples and other settings

#Randomly select sites from a set of TFBS and a set of TSS
chiaSim <- function(TFBSfile="ERa", TSSfile=NULL, mean.frag.length = 250, n.cells = 1E3, lp = 0.8,
<<<<<<< HEAD
                    seqlength = 50, N.E = 100, N.P = 100, N.nE = 100, N.nP = 100, REp = "AAGCTT", mc.cores.user = 1,
=======
                    seqlength = 50, N.E = 100, N.P = 100, REp = "AAGCTT", mc.cores.user = 1,
>>>>>>> cae7d715e79b7e5ec5fa6cb14f9865ddc43f3b44
                    beta = 1, gmclass = 1, outputformat= "base", seed=2021){

  if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) BiocManager::install("GenomicFeatures")
  require("parallel", quietly = TRUE)

  ###### load BSgenome.Hsapiens package
  require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  require("GenomicFeatures", quietly = TRUE)

  if(is.null(TSSfile)) {
    data("TSS",envir= environment())
    TSS<- TSS
  }

  if(TFBSfile=="ERa"){
    data("ERa",envir= environment())
    TFBS<-ERa
  }

  if(TFBSfile=="POL2"){
    data("POL2",envir= environment())
    TFBS<-POL2
  }

  colnames(TFBS)[1:3] <- c("chrom","chromStart","chromEnd")
  

  print("Reading TFBS and TSS")
  ###### get chromosome lengths
  #chr.length<-sapply(1:23,function(i){length(Hsapiens[[i]])})
  chr.length <- getChromInfoFromUCSC("hg19")$size[1:23]
  ###### get TFBS information
  TFBS<-TFBS[TFBS$chrom!="chrY",]
  TFBS$midpoint<-(TFBS$chromStart+TFBS$chromEnd)/2

  #Determine the portion of total number of pairs allocate to each chr
  allocateNumPairs.TFBS <- vector()
  for(i in 1:length(unique(TFBS$chrom))){
    allocateNumPairs.TFBS[i] <- length(which(TFBS$chrom == (unique(TFBS$chrom)[i])))/length(TFBS$chrom)
  }

  ######get TSS information####
  # read in annotation, remove isoforms with multiple gene symbols or no gene symbols
  temp<-strsplit(TSS[,6],",")
  TSS<-TSS[(sapply(temp,length)==1) & (sapply(temp,"[",1)!="n/a"),] #Remove isoforms with multiple gene symbols or no gene symbols
  TSS[,6]<-gsub(",","",TSS[,6])

  # remove "complicated" genes
  exclude.gene.1<-unique(TSS[grepl("_|M|Un|Y",TSS[,2]),6])
  temp<-by(TSS[,3],TSS[,6],function(x){length(unique(x))}) #Check if a gene appear in two different chrs
  exclude.gene.2<-names(temp)[temp>1]

  TSS<-TSS[!(TSS[,6] %in% union(exclude.gene.1,exclude.gene.2)),]

  # get TSS information
  TSS[,7]<-ifelse(TSS[,3]=="+",TSS[,4]+1,TSS[,5]) #get a starting position (TSS)
  TSS<-unique(TSS[,c(2,6,7)])
  colnames(TSS)<-c("chr","gene","TSS")

  #Determine the portion of total number of pairs allocate to each chr
  allocateNumPairs.TSS <- vector()
  for(i in 1:length(unique(TSS$chr))){
    allocateNumPairs.TSS[i] <- length(which(TSS$chr == (unique(TSS$chr)[i])))/length(TSS$chr)
  }

  print("Obtain cut loci")
  ###Get cut set of certain Restriction Enzyme, HindIII in this case####
  hind3 <- DNAString(REp)
  hind3 <- DNAStringSet(hind3)
  names(hind3)<-1:length(hind3)
  seqnames<-seqnames(Hsapiens)[1:23]
  matchHindIII<-function(dict0){

    dict<-PDict(dict0)  #The PDict class is a container for storing a preprocessed dictionary
    #of DNA patterns that can later be passed to the matchPDict
    #function for fast matching against a reference sequence (the subject).

    list.table<-lapply(seqnames,function(seqname){
      subject<-Hsapiens[[seqname]]
      m<-extractAllMatches(subject, matchPDict(dict, subject))
      return(data.frame(as.integer(names(m)),rep(seqname,length(m)),start(m),end(m),stringsAsFactors=FALSE))
    })
    table<-data.frame(do.call(rbind,list.table),stringsAsFactors=FALSE)
    colnames(table)<-c("PatternID","chromosome","start","end")
    return(table)
  }

  hind3cutset <- matchHindIII(hind3) # Take a couple of mins

  ###### select interested interaction sites

  print("Allocate interaction sites")
  numSites.TFBS <- round(N.E*allocateNumPairs.TFBS)# 3 - 17
  numSites.TSS <- round(N.P*allocateNumPairs.TSS) # 3 - 20
<<<<<<< HEAD
  numSites.nTFBS <- round(N.nE*allocateNumPairs.TFBS)# 3 - 17
  numSites.nTSS <- round(N.nP*allocateNumPairs.TSS) # 3 - 20



  #TFBS.selected<-mat.or.vec(23,2*max(numSites.TFBS)); TFBS.selected.id<-mat.or.vec(23, 2*max(numSites.TFBS)) #23*34
  #TSS.selected<-mat.or.vec(23,2*max(numSites.TSS)); TSS.selected.id<-mat.or.vec(23,2*max(numSites.TSS)) #23*40

   TFBS.selected<-mat.or.vec(23,max(numSites.TFBS)+max(numSites.nTFBS)); TFBS.selected.id<-mat.or.vec(23, max(numSites.TFBS)+max(numSites.nTFBS)) #23*34
   TSS.selected<-mat.or.vec(23,max(numSites.TSS)+max(numSites.nTSS)); TSS.selected.id<-mat.or.vec(23,max(numSites.TSS)+max(numSites.nTSS)) #23*40
=======
>>>>>>> cae7d715e79b7e5ec5fa6cb14f9865ddc43f3b44



  set.seed(seed=seed)
  for (i in 1:23){
    chr.idx<-ifelse(i==23,"X",i)
    temp.hind <- hind3cutset[which(hind3cutset$chromosome == paste("chr",chr.idx,sep="")),]

    temp.start <- c(1, temp.hind$end)
    temp.end <- c(temp.hind$start, chr.length[i])

    temp.seg <- cbind(temp.start, temp.end)
    temp.seg <- as.data.frame(temp.seg)
    colnames(temp.seg) <- c("start","end")

    get.start <- function(tempi){

      diff = tempi - temp.seg$start
      tempi.s <- temp.seg$start[which(diff == (min(diff[which(diff > 0)])))] #Start position of hindIII cutted TFBS segment

      return(tempi.s)
    }

    get.end <- function(tempi){

      diff =   temp.seg$end - tempi
      tempi.e <- temp.seg$end[which(diff == (min(diff[which(diff > 0)])))] #end position of hindIII cutted TFBS segment

      return(tempi.e)
    }


    temp<-sample(which(TFBS$chrom==paste("chr",chr.idx,sep="")),numSites.TFBS[i]) #randomly select n_i TFBS from a chrom
    extend.temp <- as.matrix(TFBS$midpoint[temp])

    extend.temp.TFBS <- cbind(temp, extend.temp, sapply(extend.temp, get.start), sapply(extend.temp, get.end))
    temp <- unique(extend.temp.TFBS)[,1]
    numSites.TFBS[i] <- length(temp)
    TFBS.selected[i,1:numSites.TFBS[i]]<-TFBS$midpoint[temp] #23xni to form the selected TFBS matrix
    TFBS.selected.id[i,1:numSites.TFBS[i]]<-temp #Also store their IDs

    temp.F <- sample.int(nrow(temp.seg),numSites.TFBS[i])

    temp.TFBS.F <- sapply(temp.F, function(rF){sample.int(temp.seg[rF,1]:temp.seg[rF,2],1)})
    TFBS.selected[i,(numSites.TFBS[i]+1):(numSites.TFBS[i]+numSites.nTFBS[i])]<- temp.TFBS.F
    TFBS.selected.id[i,(numSites.TFBS[i]+1):(numSites.TFBS[i]+numSites.nTFBS[i])]<-rep("No_gene",numSites.nTFBS[i])

    temp<-sample(which(TSS$chr==paste("chr",chr.idx,sep="")),numSites.TSS[i]) #randomly select m_i TFBS from a chrom
    extend.temp <- as.matrix(TSS$TSS[temp])
    extend.temp.TSS <- cbind(temp, extend.temp, sapply(extend.temp, get.start), sapply(extend.temp, get.end))
    temp <- unique(extend.temp.TSS)[,1]
    numSites.TSS[i] <- length(temp)
    TSS.selected[i,1:numSites.TSS[i]]<-TSS$TSS[temp]
    TSS.selected.id[i,1:numSites.TSS[i]]<-TSS$gene[temp]

    temp.F <- sample.int(nrow(temp.seg),numSites.TSS[i])

    temp.TSS.F <- sapply(temp.F, function(rF){sample.int(temp.seg[rF,1]:temp.seg[rF,2],1)})

    TSS.selected[i,(numSites.TSS[i]+1):(numSites.TSS[i]+numSites.nTSS[i])]<-temp.TSS.F
    TSS.selected.id[i,(numSites.TSS[i]+1):(numSites.TSS[i]+numSites.nTSS[i])]<-rep("No_gene",numSites.nTSS[i])
  }

  print("Make pairings")
  ######## make pairings
  # Set inter.prop to be the propotion of inter-chr pairs
  inter.prop <- 0.1
  z=1
  temp.1 <- list()
  for(i in 1:23){
    for(j in 1:23){
      temp.1[[z]] <- expand.grid(i,j,1:numSites.TFBS[i],1:numSites.TSS[j])#chr1, chr2, posID1,posID2
      z=z+1
    }
  }
  temp.1 <- do.call("rbind",temp.1) #sum(numSites.TFBS)*sum(numSites.TSS) = 40000 combinations (true pairs)

  temp.1[,5]<-sapply(1:nrow(temp.1),function(i){return(TFBS.selected[temp.1[i,1],temp.1[i,3]])}) #1:23, chr, id
  temp.1[,6]<-sapply(1:nrow(temp.1),function(i){return(TSS.selected[temp.1[i,2],temp.1[i,4]])})


  temp.size <- round(inter.prop*length(which(temp.1[,1]==temp.1[,2]&abs(temp.1[,5]-temp.1[,6])<=10^7))/(1-inter.prop))
  inter.1 <- sample(which(temp.1[,1]!=temp.1[,2]), temp.size)
  keep.1 <- c(which(temp.1[,1]==temp.1[,2]&abs(temp.1[,5]-temp.1[,6])<=10^7), inter.1)
  keep.1 <- sample(keep.1, length(keep.1)) #Mix intra and inter, so that TH TL could have similar proportion of intra/inter
  temp.1 <- temp.1[keep.1,] #Now the ratio of true inter to true total is inter.prop, 2120 combinations are kept
  # check length(which(temp.1[,1]!=temp.1[,2]))/length(which(temp.1[,1]==temp.1[,2])) = 0.1111 #inter/#intra = 1/9

  #summary(abs(temp.1[which(temp.1[,1]==temp.1[,2]),5]-temp.1[which(temp.1[,1]==temp.1[,2]),6]))

  z=1
  temp.2 <- list()
  for(i in 1:23){
    for(j in 1:23){
      temp.2[[z]] <- expand.grid(i,j,(numSites.TFBS[i]+1):(numSites.TFBS[i]+numSites.nTFBS[i]),(numSites.TSS[j]+1):(numSites.TSS[j]+numSites.nTSS[j])) #Non-specific TFBS and TSS
      z=z+1
    }
  }
  temp.2 <- do.call("rbind",temp.2) # sum(numSites.TFBS)*sum(numSites.TFBS) = 40000
  temp.2[,5]<-sapply(1:nrow(temp.2),function(i){return(TFBS.selected[temp.2[i,1],temp.2[i,3]])}) #1:23, chr, id
  temp.2[,6]<-sapply(1:nrow(temp.2),function(i){return(TSS.selected[temp.2[i,2],temp.2[i,4]])})

  temp.size.2 <- round(inter.prop*length(which(temp.2[,1]==temp.2[,2]&abs(temp.2[,5]-temp.2[,6])<=10^7))/(1-inter.prop))
  inter.2 <- sample(which(temp.2[,1]!=temp.2[,2]), temp.size.2)
  keep.2 <- c(which(temp.2[,1]==temp.2[,2]&abs(temp.2[,5]-temp.2[,6])<=10^7), inter.2)
  keep.2 <- sample(keep.2, length(keep.2)) #Mix intra and inter, so that F could have similar proportion of intra/inter
  temp.2 <- temp.2[keep.2,]

  #Now the ratio of true inter to true total is inter.prop, 2120 combinations are kept
  #check length(which(temp.2[,1]!=temp.2[,2]))/length(which(temp.2[,1]==temp.2[,2])) = 0.1111 #inter/#intra = 1/9
  #length(abs(temp.2[which(temp.2[,1]==temp.2[,2]),5]-temp.2[which(temp.2[,1]==temp.2[,2]),6]))


  pairing<-rbind(temp.1,temp.2) #chr#1, chr#2, id1, id2, [5]TFBS.pos, [6]TSS.pos  #7311

  bisect.temp <- round(nrow(temp.1)/2)
  PT<-rep("TH", bisect.temp)
  PT<-c(PT,rep("TL",(nrow(temp.1) - bisect.temp)))
  #PT<-sapply(1:nrow(temp.1),function(i){ifelse(temp.1[i,3]==temp.1[i,4],"TH","TL")}) #Group one consists the pair with same
  PT<-c(PT,rep("F",nrow(temp.2)))
  ######Inmitate Loop Ligation#######

  ###Extend the TFBS/TSS site to the segment of hindIII cut.
  for(i in 1:nrow(pairing)){


    #TFBS
    chr.idx<-ifelse(pairing[i,1]==23,"X", pairing[i,1])
    temp <- pairing[i,5] #5th column is the mid-point of TFBS
    temp.hind <- hind3cutset[which(hind3cutset$chromosome == paste("chr",chr.idx,sep="")),]

    diff = temp - temp.hind$end
    if(length(which(diff > 0)) == 0) {pairing[i,7] <- 1}
    else{
      pairing[i,7] <- temp.hind$end[which(diff == (min(diff[which(diff > 0)])))] #Start position of hindIII cutted TFBS segment
    }

    diff = temp.hind$start - temp
    if(length(which(diff > 0)) == 0) {pairing[i,8] <- chr.length[pairing[i,1]]}
    else{
      pairing[i,8] <- temp.hind$start[which(diff == (min(diff[which(diff > 0)])))] #end position of hindIII cutted TFBS segment
    }

    #TSS
    chr.idx<-ifelse(pairing[i,2]==23,"X", pairing[i,2])
    temp <- pairing[i,6] #6th column is the position of TSS
    temp.hind <- hind3cutset[which(hind3cutset$chromosome == paste("chr",chr.idx,sep="")),]

    diff = temp - temp.hind$end
    if(length(which(diff > 0)) == 0) {pairing[i,9] <- 1}
    else{
      pairing[i,9] <- temp.hind$end[which(diff == (min(diff[which(diff > 0)])))]#Start position of hindIII cutted TSS segment
    }

    diff = temp.hind$start - temp
    if(length(which(diff > 0)) == 0) {pairing[i,10] <- chr.length[pairing[i,2]]}
    else{
      pairing[i,10] <- temp.hind$start[which(diff == (min(diff[which(diff > 0)])))] #end position of hindIII cutted TSS segment
    }

  }

  ####### set sampling probabilities so that TH:TL:F=6:3:1
  same.cut <- which(pairing[,1]==pairing[,2]&pairing[,7]==pairing[,9]&pairing[,8]==pairing[,10]) #4 site
  if(length(same.cut)!=0){
    pairing <- pairing[-same.cut,] #chr#1, chr#2, id1, id2, [5]TFBS.pos, [6]TSS.pos
    PT <- PT[-same.cut]
    #pairing.gene <- pairing.gene[-same.cut]
  }

  ###Check how many segment of "F" are selected repetitively
  temp.F <- pairing[PT=="F",]

  temp.unique.F <- temp.F[!duplicated(temp.F[,c(1:2,7:10)]), ] #Select unique "F" pairings
  #Keep unique combinations for "F" group
  pairing <- rbind(pairing[PT!="F",], temp.unique.F)
  pairing.gene<-sapply(1:nrow(pairing),function(i){return(TSS.selected.id[pairing[i,2],pairing[i,4]])}) #Column 1&2 are true
  #TSS.selelceted matrix, last 4 are non-specific

  PT<-c(PT[PT!="F"],rep("F",nrow(temp.unique.F)))     #TH = Group 1, TL= group 2                  #locus index



  ####### set sampling probabilities so that TH:TL:F=6:3:1
  #beta <- 1

  temp <- vector() # a vector of 1/distance
  for(i in 1:nrow(pairing)){
    if(pairing[i,1] == pairing[i,2]){
      temp[i] <- (1/abs(pairing[i,5]-pairing[i,6]))^beta #Intra pair depend on distance
    }
    else {
      temp[i] <- (1/10^7)^beta
    }
  }


  TH.prob <- temp[which(PT=="TH")]/sum(temp[which(PT=="TH")])
  TL.prob <- temp[which(PT=="TL")]/sum(temp[which(PT=="TL")])
  F.prob <- temp[which(PT=="F")]/sum(temp[which(PT=="F")])

  sampling.prob<-rep(NA, nrow(pairing))
  sampling.prob[which(PT=="TH")] <- 0.6*TH.prob
  sampling.prob[which(PT=="TL")] <- 0.3*TL.prob
  sampling.prob[which(PT=="F")] <- 0.1*F.prob


  ###### sample the pairs
  sample.index<-sample(1:nrow(pairing),n.cells,sampling.prob,replace=T)
  #sample.ligation<-ligation[sample(1:7,n.cells,p.ligation,replace=T)]

  ###### perform sonications to get the PET reads

  ligation.function<-function(chr.1,chr.2,frag.1,frag.2, seqlength){
    if (frag.1 != "NA" && frag.2!="NA"){
      a<-toString(reverseComplement(subseq(Hsapiens[[chr.1]],frag.1[1],frag.1[1]+seqlength-1)))
      d<-toString(subseq(Hsapiens[[chr.2]],frag.1[2]-seqlength+1,frag.1[2]))

      b<-toString(subseq(Hsapiens[[chr.1]],frag.2[1],frag.2[1]+seqlength-1))
      c<-toString(reverseComplement(subseq(Hsapiens[[chr.2]],frag.2[2]-seqlength+1,frag.2[2])))


      temp.1 <- strsplit("ad",split="")[[1]]
      temp.2 <- strsplit("bc", split = "")[[1]]

      return(c(get(temp.1[1]),get(temp.1[2]), get(temp.2[1]), get(temp.2[2])))
    }
    else{
      return(NA)
    }
  }

  print("Sonication")
  temp.sindex <- sample.index
  ###### perform sonications to get the PET reads
  sonication.result<-mclapply(1:n.cells,function(idx.cell){
    i<-sample.index[idx.cell]

    chr.1<-pairing[i,1]
    chr.2<-pairing[i,2]

    a <-  pairing[i,7]
    b <- pairing[i,8]
    c <- pairing[i,9]
    d <- pairing[i,10]

    loop.seq <- c(a:b,c:d)
    ran.cut <- sample(1:length(loop.seq), 1)

    num.frag.1<-rpois(1,length(loop.seq)/mean.frag.length)

    sonicate.pt.1<-c(0,sample.int(length(loop.seq), num.frag.1), length(loop.seq)) + ran.cut

    sonicate.pt.1 <- sonicate.pt.1 %% length(loop.seq)

    #get the fragment that contains this position
    if(min(sonicate.pt.1)>=50 & (length(loop.seq)-max(sonicate.pt.1)) >=50){
      frag.1 <- c(a + min(sonicate.pt.1), d - (length(loop.seq)-max(sonicate.pt.1)))
    }
    else{frag.1 = "NA"}

    temp.diff <- sonicate.pt.1 - length(a:b)

    if(length(temp.diff[temp.diff > 0]) > 0 & length(temp.diff[temp.diff < 0])>0){
      temp.cut.1 <- sonicate.pt.1[which(temp.diff == min(temp.diff[temp.diff > 0]))]
      temp.cut.2 <- sonicate.pt.1[which(temp.diff == max(temp.diff[temp.diff < 0]))]
      temp.cut <- c(temp.cut.1, temp.cut.2)

      if((length(a:b) - min(temp.cut))>=50 & (max(temp.cut) - (length(a:b)+1)) >= 50){

        frag.2 <- c( b - (length(a:b) - min(temp.cut)), c + (max(temp.cut) - (length(a:b)+1)))

      }
      else{frag.2 ="NA"}
    }
    else{frag.2 ="NA"}
    #frag.1 <- c(a,b)
    #frag.2 <- c(c,d)
    return(ligation.function(chr.1,chr.2,frag.1,frag.2,seqlength))
  },mc.cores = mc.cores.user)


  ###### remove redundant sonication results
  temp.sindex <- temp.sindex[!sapply(sonication.result,function(i){any(is.na(i))})]
  sonication.result<-sonication.result[!sapply(sonication.result,function(i){any(is.na(i))})]

  PET.seqa<-sapply(sonication.result,"[",1)
  PET.seqc<-sapply(sonication.result,"[",4)
  PET.seqb<-sapply(sonication.result,"[",3)
  PET.seqd<-sapply(sonication.result,"[",2)

  idx<-!(grepl("N",PET.seqa) | grepl("N",PET.seqb)|grepl("N",PET.seqc)|grepl("N",PET.seqd))
  PET.seqa<-PET.seqa[idx]
  PET.seqc<-PET.seqc[idx]
  PET.seqb<-PET.seqb[idx]
  PET.seqd<-PET.seqd[idx]
  temp.sindex <- temp.sindex[idx]
  ###### make two DNAStringSet objects to contain the PET's
  PET.seqa<-DNAStringSet(PET.seqa)
  PET.seqc<-DNAStringSet(PET.seqc)
  PET.seqb<-DNAStringSet(PET.seqb)
  PET.seqd<-DNAStringSet(PET.seqd)

  names(PET.seqa)<-1:length(PET.seqa)
  names(PET.seqc)<-1:length(PET.seqc)
  names(PET.seqb)<-1:length(PET.seqb)
  names(PET.seqd)<-1:length(PET.seqd)

  print("Map the PET's to genome")
  ###### map the PET's to genome
  seqnames<-seqnames(Hsapiens)[1:23]
  new.matchpdict<-function(dict0,strand,mc.cores.user){
    if (strand == "-"){
      dict0<-reverseComplement(dict0)
    }
    dict<-PDict(dict0)

    list.table<-mclapply(seqnames,function(seqname){
      subject<-Hsapiens[[seqname]]
      m<-extractAllMatches(subject, matchPDict(dict, subject))
      return(data.frame(as.integer(names(m)),rep(seqname,length(m)),start(m),end(m),stringsAsFactors=FALSE))
    },mc.cores=mc.cores.user)
    table<-data.frame(do.call(rbind,list.table),strand,stringsAsFactors=FALSE)
    colnames(table)<-c("PatternID","chromosome","start","end","strand")
    return(table)
  }


  match.seqa<-new.matchpdict(PET.seqa,strand="-",mc.cores.user=mc.cores.user)
  match.seqc<-new.matchpdict(PET.seqc,strand="-",mc.cores.user=mc.cores.user)
  match.seqb<-new.matchpdict(PET.seqb,strand="+",mc.cores.user=mc.cores.user)
  match.seqd<-new.matchpdict(PET.seqd,strand="+",mc.cores.user=mc.cores.user)
  ###### remove non-unique alignments
  temp.a<-table(match.seqa$PatternID)
  temp.b<-table(match.seqb$PatternID)
  temp.c<-table(match.seqc$PatternID)
  temp.d<-table(match.seqd$PatternID)

  include.pattern1<-as.integer(intersect(names(temp.a)[temp.a==1],names(temp.c)[temp.c==1]))
  include.pattern2<-as.integer(intersect(names(temp.b)[temp.b==1],names(temp.d)[temp.d==1]))

  include.pattern <- as.integer(intersect(include.pattern1,include.pattern2))
  temp.sindex <- temp.sindex[include.pattern]
  match.seqa<-match.seqa[match.seqa$PatternID %in% include.pattern,]
  match.seqa<-match.seqa[order(match.seqa$PatternID),]
  match.seqc<-match.seqc[match.seqc$PatternID %in% include.pattern,]
  match.seqc<-match.seqc[order(match.seqc$PatternID),]


  match.seqb<-match.seqb[match.seqb$PatternID %in% include.pattern,]
  match.seqb<-match.seqb[order(match.seqb$PatternID),]
  match.seqd<-match.seqd[match.seqd$PatternID %in% include.pattern,]
  match.seqd<-match.seqd[order(match.seqd$PatternID),]

  match.seq1<-rbind(match.seqa,match.seqb)
  match.seq2<-rbind(match.seqd,match.seqc)
  temp.sindex <- c(temp.sindex,temp.sindex)


  if(gmclass == 2){
    ##### distinguish self-loops from inter-loops
    dist<-unlist(mclapply(1:nrow(match.seq1),function(i){
      if (match.seq1$chromosome[i]==match.seq2$chromosome[i]){
        return(abs(match.seq1$start[i]-match.seq2$start[i])+seqlength)
      }
      else {
        return(0)
      }
    },mc.cores=mc.cores.user))
    temp<-kmeans(log(dist[dist>0]),2)
    self.idx<-ifelse(temp$centers[1]>temp$centers[2],as.integer(2),as.integer(1))

    group<-ifelse(dist==0,"inter.diff","same")
    group[group=="same"]<-ifelse(temp$cluster==self.idx,"self","inter.same")

    ###### create GRanges
    match.seq1<-match.seq1[group!="self",]
    match.seq2<-match.seq2[group!="self",]
    temp.sindex<-temp.sindex[group!="self"]

  }

  names(chr.length)<-paste("chr",c(1:22,"X"),sep="")

  gr.seq1<-GRanges(seqnames=match.seq1$chromosome,ranges=IRanges(start=match.seq1$start,end=match.seq1$end),strand=match.seq1$strand,seqlengths=chr.length)
  gr.seq2<-GRanges(seqnames=match.seq2$chromosome,ranges=IRanges(start=match.seq2$start,end=match.seq2$end),strand=match.seq2$strand,seqlengths=chr.length)

  ##### extend to mean.frag.length
  gr.seq1<-resize(gr.seq1,mean.frag.length,fix="end")
  gr.seq2<-resize(gr.seq2,mean.frag.length,fix="end")


  ##### find interaction anchors, counts and marginal counts (mc)
  all.gr<-reduce(c(gr.seq1,gr.seq2),ignore.strand=TRUE)

  hit.seq1<-findOverlaps(gr.seq1,all.gr)
  hit.seq2<-findOverlaps(gr.seq2,all.gr)

  temp<-cbind(subjectHits(hit.seq1),subjectHits(hit.seq2))
  temp[subjectHits(hit.seq2)<subjectHits(hit.seq1),] <- temp[subjectHits(hit.seq2)<subjectHits(hit.seq1),2:1]


  raw.table<-temp
  temp.rt <- paste0(raw.table[,1],":", raw.table[,2])
  ones<-rep(1,nrow(raw.table))
  temp.ct<-aggregate(ones,by=list(temp.rt),sum) #count interactions between pairs
  temp.ctid <- do.call(rbind, strsplit(temp.ct[,1], ":"))
  count.table <- as.data.frame(cbind(as.numeric(temp.ctid[,1]), as.numeric(temp.ctid[,2]), temp.ct[,2]))

  temp.LT.index <- match(temp.ct[,1], temp.rt)
  LT.index <- temp.sindex[temp.LT.index]#label each pairing to the original LT group,i.e., TH,TL,F

  temp<-paste(seqnames(all.gr),":",start(all.gr),"..",end(all.gr),sep="")
  count.table[,4]<-temp[count.table[,1]]
  count.table[,5]<-temp[count.table[,2]]

  temp.reference<-table(as.vector(raw.table))
  count.table[,6]<-unlist(mclapply(1:nrow(count.table),function(i){
    return(temp.reference[as.character(count.table[i,1])]+temp.reference[as.character(count.table[i,2])])
  },mc.cores=mc.cores.user))
  count.table[,7] <- PT[LT.index]
  colnames(count.table)<-c("frag.1.id","frag.2.id","count","frag.1","frag.2","MC","LT")

  print("Calculate Mdists")
  ##### get distances
  #all.gr, TFBS, TSS
  all.gr.mid<-(start(all.gr)+end(all.gr))/2
  all.gr.chr<-as.character(seqnames(all.gr))
  all.gr.start<-as.numeric(start(all.gr)); all.gr.end<-as.numeric(end(all.gr))

  temp.TFBS<-mclapply(1:length(all.gr),function(i){
    temp<-TFBS[TFBS$chrom==all.gr.chr[i],]
    dist<-abs(temp$midpoint-all.gr.mid[i])

    dist.i<-min(dist)
    idx.i<-which.min(dist)

    return(c(dist.i,temp$midpoint[idx.i]))
  },mc.cores=mc.cores.user)

  temp.TSS<-mclapply(1:length(all.gr),function(i){
    temp<-TSS[TSS$chr==all.gr.chr[i],]
    dist<-abs(temp$TSS-all.gr.mid[i])

    dist.i<-min(dist)
    idx.i<-which.min(dist)

    return(list(dist.i,temp$gene[idx.i]))
  },mc.cores=mc.cores.user)

  pairing<-data.frame(pairing)
  pairing[,1]<-paste("chr",ifelse(pairing[,1]==23,"X",pairing[,1]),sep="")
  pairing[,2]<-paste("chr",ifelse(pairing[,2]==23,"X",pairing[,2]),sep="")

  mini.distance<-mclapply(1:nrow(count.table),function(i){
    index.gr.1<-count.table$frag.1.id[i]; index.gr.2<-count.table$frag.2.id[i];
    a<-(temp.TFBS[[index.gr.1]][1]+temp.TSS[[index.gr.2]][[1]])
    b<-(temp.TSS[[index.gr.1]][[1]]+temp.TFBS[[index.gr.2]][1])

    if (a>b){return(b)}
    else {return(a)}
  },mc.cores=mc.cores.user)

  count.table$Mdist<-as.numeric(sapply(mini.distance,"[",1))

    simu <- count.table[,c(3,6,7,8)]
    temp1 <- strsplit(as.character(count.table$frag.1), ":")
    temp2 <- strsplit(as.character(count.table$frag.2), ":")
    positionA <-  sapply(temp1,"[",2)
    simu$startA <- as.numeric(sapply(strsplit(positionA, "[..]"),"[",1))
    simu$endA <- as.numeric(sapply(strsplit(positionA, "[..]"),"[",3))
    positionB <-  sapply(temp2,"[",2)
    simu$startB <- as.numeric(sapply(strsplit(positionB, "[..]"),"[",1))
    simu$endB <- as.numeric(sapply(strsplit(positionB, "[..]"),"[",3))
    simu$chrA <- sapply(temp1,"[",1)
    simu$chrB <- sapply(temp2,"[",1)



    simu$Gdist <- abs((simu$startB + simu$endB - simu$startA - simu$endA)/2)
    if(length(which(simu$chrA != simu$chrB))>0){
      simu$Gdist[which(simu$chrA != simu$chrB)] <- 2*max(simu$Gdist[which(simu$chrA == simu$chrB)])
    }
    

    simu <- simu[,c(9,5,6,10,7,8,1,3,2,4,11)]

  if(outputformat=="BEDPE"){
    tempA<-paste(simu$chrA,":",simu$startA,"..",simu$endA,sep="")
    tempB<-paste(simu$chrB,":",simu$startB,"..",simu$endB,sep="")
    temp<- paste(tempA,"-",tempB,",",simu$count, sep="")

    simu$score <- "."
    simu$strand <- "."
    simu <- cbind(simu[,c(1,2,6)],temp,simu[,c(12,13,8:11)])
    colnames(simu)[1:4] <- c("chrom",	"chromStart",	"chromEnd",	"name")

    temp  <- list()
    temp$BEDPE <- simu[,1:6]
    temp$AddInfo <- simu[,7:10]
    simu <- temp	
  }


  return(simu)

}
