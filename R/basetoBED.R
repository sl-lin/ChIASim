basetoBED<-function(base){
  tempA<-paste(base$chrA,":",base$startA,"..",base$endA,sep="")
  tempB<-paste(base$chrB,":",base$startB,"..",base$endB,sep="")
  temp<- paste(tempA,"-",tempB,",",base$count, sep="")
  
  base$score <- "."
  base$strand <- "."
  base <- cbind(base[,c(1,2,6)],temp,base[,c(11,12,8:10)])
  colnames(base)[1:4] <- c("chrom",	"chromStart",	"chromEnd",	"name")
  
  return(base)
}
