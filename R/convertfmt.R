convertfmt<-function(input, inputformat=c("base","BEDPE")){
  if(inputformat=="base"){
  tempA<-paste(input$chrA,":",input$startA,"..",input$endA,sep="")
  tempB<-paste(input$chrB,":",input$startB,"..",input$endB,sep="")
  temp<- paste(tempA,"-",tempB,",",input$count, sep="")
  
  input$score <- "."
  input$strand <- "."
  input <- cbind(input[,c(1,2,6)],temp,input[,c(12,13,8:11)])
  colnames(input)[1:4] <- c("chrom",	"chromStart",	"chromEnd",	"name")
   
    temp  <- list()
    temp$BEDPE <- input[,1:6]
    temp$AddInfo <- input[,7:10]
    input <- temp
  }
  
  if(inputformat=="BEDPE"){
    input <- cbind(input[[1]],input[[2]])
    temp <- as.data.frame(do.call(rbind,strsplit(as.character(input[,4]),split = "[:.,-]+")))
    
    colnames(temp)<-c("chrA","startA","endA","chrB","startB","endB","count")
    temp$chrA <-as.character(temp$chrA)
    temp$startA <- as.numeric(as.character(temp$startA))
    temp$endA <- as.numeric(as.character(temp$endA))
    temp$chrB <-as.character(temp$chrB)
    temp$startB <- as.numeric(as.character(temp$startB))
    temp$endB <- as.numeric(as.character(temp$endB))
    temp$count <- as.numeric(as.character(temp$count))
    
    temp <- cbind(temp, input[,7:10])
    input <- temp
  }
  
  return(input)
}
