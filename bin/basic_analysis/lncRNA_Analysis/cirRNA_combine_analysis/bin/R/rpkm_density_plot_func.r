#!/share/nas2/genome/biosoft/R/current/bin/Rscript

# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript rpkm_density_plot_func.r read_count.txt outdir")
	print("1) read_count.txt: the read count file for RNA-SEQ")
	print("2) outdir: the dir for output")
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
require(ggplot2)


# read count data
print("read count data ...")
count_data <- read.delim(args[1], row.names = 1, header=TRUE)
head(count_data)
print("read count data is over")




sample <- colnames(count_data)
rpkm <- count_data[ ,sample]
log10rpkm <- rpkm
for ( i in 1:(dim(rpkm)[2]) ){
	log10rpkm[,i] <- log10(rpkm[,i]+1)
}


finaldata<-data.frame(log10rpkm) 
# output
file <- paste(args[2], ".tpm_box.png", sep="")
png(filename=file, height = 3000, width = 3400, res = 500, units = "px")
boxplot(finaldata,aes(factor(sample)),col=c("steelblue","mediumturquoise","hotpink"),ylab="log10tmp+1",las=1, font.lab=2)
dev.off()







