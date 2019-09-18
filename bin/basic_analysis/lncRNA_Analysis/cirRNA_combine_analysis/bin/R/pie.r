#!/share/nas2/genome/biosoft/R/current/bin/Rscript

# usage function
usage <- function() {
  print("-------------------------------------------------------------------------------")
  print("Usage: Rscript rpkm_density_plot_func.r read_count.txt outdir")
  print("1) ciriRNA.outfile")
  print("2) outdir: the dir for output")
  print("3) file type, 1: circRNA result; 0: statistics file")
  print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)

# check args length
if( length(args) != 3 ) {
  print(args)
  usage()
  stop("the length of args != 3")
}


# load library
#library(plotrix)

#Hex1<-rgb(red=190,green=186,blue=218,max = 255)
#Hex2<-rgb(red=128,green=177,blue=211,max = 255)
#Hex3<-rgb(red=179,green=222,blue=105,max = 255)
#Hex4<-rgb(red=252,green=205,blue=229,max = 255)
#Hex5<-rgb(red=141,green=211,blue=199,max = 255)
#Hex6<-rgb(red=188,green=128,blue=189,max = 255)
#Hex7<-rgb(red=253,green=192,blue=134,max = 255)
#Hex8<-rgb(red=217,green=217,blue=217,max = 255)
#colours1 <- c(Hex2,Hex3,Hex4,Hex5,Hex6,Hex7,Hex8,Hex1)
#other1<-rgb(red=130,green=57,blue=53,max = 255)
#other2<-rgb(red=137,green=190,blue=178,max = 255)
#other3<-rgb(red=201,green=186,blue=131,max = 255) 
#other4<-rgb(red=222,green=211,blue=140,max = 255)
#other5<-rgb(red=222,green=156,blue=83,max = 255) 
#colours2 <- c(other1,other2,other3,other4,other5)
# read count data
print("read count data ...")
if(args[3] == 1)
{
	count_data <- read.delim(args[1],header=TRUE,sep="\t")
	slices <- count_data[,"circRNA_type"]
	mytable <- table(slices)
	lbls <- paste(names(mytable), " : ", mytable, sep="")
	maintxt<-"Pie Chart of circRNA type"
}else
{
	count_data <- read.delim(args[1],header=FALSE,sep="\t")
	mytable <- count_data[,2]
	lbls <- paste(count_data[,1], " : ", mytable, sep="")
	maintxt<-"Pie Chart of known and new circRNA Statistics"
}
head(count_data)
print("read count data is over")
#colours <- rainbow(10)
colours <- c("#8dd3c7", "#ffffb3", "#bebada", "#80b1d3", "#fb8072", "#fdb462", "#b3de69", "#d9d9d9", "#fccde5")
# output
file <- paste(args[2], ".pie.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
pie(mytable, labels = lbls, main=maintxt, col = colours)
dev.off()







