

rm(list=ls())

#####################################################
# DEBUG DEBUG DEBUG DEBUG DEBUG
# this disables warnings !! comment out for debugging
#####################################################
options(warn=-1)
#####################################################

# supported output extensions
outExt = data.frame(
  mode = c("png","pdf","svg"),
  regex = c("^.+\\.[pP][nN][gG]$","^.+\\.[pP][dD][fF]$","^.+\\.[sS][vV][gG]$")
)

#### dependencies ###

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(ggforce)
  library(ggthemes)
})

#### input ###

args <- commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
args <- c("RNA_dp.ps"
         ,"RNA-MUTANT_dp.ps"
         ,"RNA-MUTANT_dp_diff_weakened.dp"
         ,"RNA-MUTANT_dp_diff_increased.dp"
         ,"tmp.png")
}
# check arguments
if (length(args)!=5) {
  stop(paste("\nCall Rscript --vanilla MutaRNA_arcplot.R"
             ," <viennaDot_WT.ps>"
             ," <viennaDot_MUT.ps>"
             ," <viennaDot_INCREASED.ps>"
             ," <viennaDot_WEAKENED.ps>"
             ," <output.[",paste(outExt$mode,collapse="|"),"]>\n"
             ,sep=""), call.=FALSE)
}

# distribute input arguments
psFileW = args[1]
psFileM = args[2]
psFileD = args[3]
psFileI = args[4]
outFile = args[5]

#### parsing ###

# parse output mode via file ending
outMode = NA
for (i in 1:nrow(outExt) ) {
  if (1==sum(grep(outExt$regex[i],outFile))) {
    outMode = outExt$mode[i]
    break;
  }
}
if (is.na(outMode)) {
  stop(paste("output mode (defined by outFile extension) not supported. Has to be one of [",outExt$mode,"]",sep=""))
}

parsePS <- function( psFile ) {

  # parse dotplot PS file
  if (!file.exists(psFile)) {
    stop(paste("Cannot open input file",psFile,"for reading!"))
  }
  FileInput = readLines(psFile)

  # parse sequence
  seq = gsub("\\\\","",paste(FileInput[grep("^[acgutnACGUTN\\\\\\s]+$",FileInput)],collapse=""))

  # parse base pair probabilities
  # create empty data frame with column names i, j,
  bpProb = data.frame(i=numeric(0),j=numeric(0),prob=numeric(0),r=numeric(0),c=numeric(0))
  # if base pairs are present (not necessarily the case for subset visualization)
  if(length(str_subset(FileInput,"^(\\d+.+ubox)$"))>0) {
    # extract base pair probabilities
    bpProb =
      dplyr::mutate(
        dplyr::select(
          as_tibble(
            str_split( FileInput[grep("^(\\d+.+ubox)$",FileInput)], " ", simplify=T)
            , .name_repair=function(x){c("i","j","prob","ubox")})
          ,-ubox)
        , dplyr::across(dplyr::everything(), as.numeric))
    # compute additional values
    bpProb[,"r"] = (bpProb$j - bpProb$i)/2
    bpProb[,"c"] = bpProb$i + bpProb$r
    # sort by probability
    bpProb <- bpProb[sort.list(bpProb$prob,decreasing = F),]
  } else {
    bpProb <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(bpProb)=c("i","j","prob","r","c")
    bpProb[1,] <- NA_real_
  }
  return(list(bpProb=bpProb, seq=seq))
}

dW = parsePS(psFileW)
seqLength = str_length(dW$seq)
dM = parsePS(psFileM)
if(seqLength != str_length(dM$seq)) {
  stop("sequences of WT and MUT differ in lengths")
}
dI = parsePS(psFileI)
if(seqLength != str_length(dI$seq)) {
  stop("sequences of WT and INCREASED differ in lengths")
}
dD = parsePS(psFileD)
if(seqLength != str_length(dD$seq)) {
  stop("sequences of WT and DECREASED differ in lengths")
}

idxM = (1:seqLength)[unlist(strsplit(dW$seq,"")) != unlist(strsplit(dM$seq,""))]

# set output theme settings
old <- theme_set(theme_clean())
theme_update(axis.title.y=element_blank(),
             axis.text.y=element_blank(),
             axis.ticks.y=element_blank(),
             axis.line.y=element_blank(),
             legend.background = element_rect(fill="white",color="white"))

plotWidth = 1+ seqLength/100*6
plotHeight = 1 + (max(dW$bpProb$r)+max(dM$bpProb$r))/100*4
plotRes = min( 300, floor(32500/plotWidth) )

########## full length bp prob  #################

if (outMode=="pdf") {
  pdf(outFile, width= plotWidth, height = plotHeight)
} else if (outMode=="png") {
  png(outFile, units="in", res=plotRes, width= plotWidth, height = plotHeight)
} else if (outMode=="svg") {
  svg(outFile, bg="white", width= plotWidth, height = plotHeight)
}

# plot base pair probs
(p <- ggplot(dW$bpProb) +
    geom_vline(xintercept=idxM,col="red", linetype="dashed" ) +
    geom_arc(aes(x0 = c, y0 = 2, r=r,
                 start = -pi/2, end = pi/2,
                 col = prob), size=0.6)  +
    geom_arc(aes(x0 = c, y0 = -2, r=r,
                 start = pi/2, end = 1.5*pi,
                 col = prob), size=0.6, data = dM$bpProb) +
    scale_colour_gradient(low = "yellow", high = "red") +
    ylab(NULL) +
    xlim(1,seqLength) +
    labs(col="Pr(bp)",x="sequence position (nt)",y="") +
    annotate("text", x=1,y= max(dW$bpProb$r), label="WT",size=4,col="Darkgray",hjust=0)+
    annotate("text", x=1,y= -max(dM$bpProb$r), label="Mutant",size=4,col="Darkgray",hjust=0)+
    annotate("text", x=1:seqLength,y= 1, label=unlist(strsplit(dW$seq,"")),size=1.5,vjust=0.5)+
    annotate("text", x=idxM,y= -1, label=unlist(strsplit(dM$seq,""))[idxM],size=1.5,col="red",vjust=0.5)
)

# close output file
dev0 = dev.off()

########## full length diff prob  #################

if (outMode=="pdf") {
  pdf(str_replace(outFile,".[pP][dD][fF]$",".diff.pdf"), width= plotWidth, height = plotHeight)
} else if (outMode=="png") {
  png(str_replace(outFile,".[pP][nN][gG]$",".diff.png"), units="in", res=plotRes, width= plotWidth, height = plotHeight)
} else if (outMode=="svg") {
  svg(str_replace(outFile,".[sS][vV][gG]$",".diff.svg"), bg="white", width= plotWidth, height = plotHeight)
}

# plot base pair probs
(p <- ggplot(dD$bpProb) +
    geom_vline(xintercept=idxM,col="red", linetype="dashed" ) +
    geom_arc(aes(x0 = c, y0 = 2, r=r,
                 start = -pi/2, end = pi/2,
                 col = prob), size=0.6)  +
    geom_arc(aes(x0 = c, y0 = -2, r=r,
                 start = pi/2, end = 1.5*pi,
                 col = prob), size=0.6, data = dI$bpProb) +
    scale_colour_gradient(low = "yellow", high = "red") +
    ylab(NULL) +
    xlim(1,seqLength) +
    labs(col="Pr(bp)",x="sequence position (nt)",y="") +
    annotate("text", x=1,y= max(dW$bpProb$r), label="weakened (WT-MUT)[>0]",size=4,col="Darkgray",hjust=0)+
    annotate("text", x=1,y= -max(dM$bpProb$r), label="strengthened (MUT-WT)[>0]",size=4,col="Darkgray",hjust=0)+
    annotate("text", x=1:seqLength,y= 1, label=unlist(strsplit(dW$seq,"")),size=1.5,vjust=0.5)+
    annotate("text", x=idxM,y= -1, label=unlist(strsplit(dM$seq,""))[idxM],size=1.5,col="red",vjust=0.5)
)


# close output file
dev0 = dev.off()

########## cut-out diff prob  #################

if (outMode=="pdf") {
  pdf(str_replace(outFile,".[pP][dD][fF]$",".diff-cut.pdf"), width= plotWidth, height = plotHeight)
} else if (outMode=="png") {
  png(str_replace(outFile,".[pP][nN][gG]$",".diff-cut.png"), units="in", res=plotRes, width= plotWidth, height = plotHeight)
} else if (outMode=="svg") {
  svg(str_replace(outFile,".[sS][vV][gG]$",".diff-cut.svg"), bg="white", width= plotWidth, height = plotHeight)
}


# plot base pair probs
(p <- ggplot(dD$bpProb) +
    geom_vline(xintercept=idxM,col="red", linetype="dashed" ) +
    geom_arc(aes(x0 = c, y0 = 2, r=r,
                 start = -pi/2, end = pi/2,
                 col = prob), size=0.6)  +
    geom_arc(aes(x0 = c, y0 = -2, r=r,
                 start = pi/2, end = 1.5*pi,
                 col = prob), size=0.6, data = dI$bpProb) +
    scale_colour_gradient(low = "yellow", high = "red") +
    ylab(NULL) +
    labs(col="Pr(bp)",x="sequence position (nt)",y="") +
    annotate("text", x=1,y= max(dW$bpProb$r), label="weakened (WT-MUT)[>0]",size=4,col="Darkgray",hjust=0)+
    annotate("text", x=1,y= -max(dM$bpProb$r), label="strengthened (MUT-WT)[>0]",size=4,col="Darkgray",hjust=0)+
    annotate("text", x=1:seqLength,y= 1, label=unlist(strsplit(dW$seq,"")),size=1.5,vjust=0.5)+
    annotate("text", x=idxM,y= -1, label=unlist(strsplit(dM$seq,""))[idxM],size=1.5,col="red",vjust=0.5)
)

# close output file
dev0 = dev.off()


