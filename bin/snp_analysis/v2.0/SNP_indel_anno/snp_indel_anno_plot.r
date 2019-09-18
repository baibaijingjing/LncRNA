#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
  print("-------------------------------------------------------------------------------")
  print("Usage: Rscript bin/cog_anno_plot.r in.stat out.png")
  print("1) in.stat: the stat file for SNP/InDel Annoation")
  print("2) outdir: the dirname for output")
  print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 2 ) {
  print(args)
  usage()
  stop("the length of args != 2")
}



# load library
library(ggplot2)
library(grid)


# get args
args <-commandArgs(TRUE)


# reading data
data <- read.delim(args[1], header=TRUE,check.names =F, sep="\t")
data<-data[,-1]
head(data)

Sample<-colnames(data[,-1])

# plot
for (i in 1:length(Sample))
{

  df <- data.frame(Number=data[,i+1], type=factor(data$Type,levels=factor(data$Type[nrow(data):1])))
  p <- ggplot(data=df, aes(x=type, y=Number)) + geom_bar(aes(fill=type), stat="identity")
  p <- p + theme_classic()
  p<- p + theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  p<-p+coord_flip()+guides(fill=FALSE) #+ theme_bw()
  p<-p+geom_text(label=data[,i+1], hjust=-0.5, size=3.5)+ylim(c(0,max(data[,i+1]*1.1)))
  p<-p+labs(x="Type", title="Number of effects by type") 
  p <- p+ theme(axis.text.y = element_text(colour ="black"))
  
  # output plot	
  png(filename=paste(args[2],paste(Sample[i],"anno.stat.png",sep="."),sep='/'), height = 3000, width = 5000, res = 500, units = "px")
  print(p)
  dev.off()
}









