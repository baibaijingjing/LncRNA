#!/share/nas2/genome/biosoft/R/3.1.1/lib64/R/bin/Rscript


# usage function
usage <- function() {
	print("-------------------------------------------------------------------------------")
	print("Usage: Rscript fpkm_cor_plot.r read_count.txt outdir")
	print("1) str:names of treated samples")
	print("2) str:names of control samples")
	print("3) infile:All genes expression list of Group")
	print("4) outfile: the name of output")
	print("-------------------------------------------------------------------------------")
}


# get args
args <-commandArgs(TRUE)


# check args length
if( length(args) != 4 ) {
	print(args)
	usage()
	stop("the length of args != 4")
}


require(ggplot2)
ALL<-read.delim(args[3], row.names = 1, header=TRUE,check.names =F)
v1<-strsplit(args[1],"_")
v2<-strsplit(args[2],"_")
v1<-as.vector(v1[[1]])
v2<-as.vector(v2[[1]])
	if( length(v1)== 1 ){
		lev1 <- ALL[,v1]
	}else{
		lev1 <- rowMeans( ALL[,v1] )
	}
	if( length(v2) == 1 ){
		lev2 <- ALL[,v2]
	}else{
		lev2 <- rowMeans( ALL[,v2] )
	}
p<-qplot(log10( lev2 ), log10( lev1 ), xlab=args[2], ylab=args[1], main="Correlation Plot", size=I(0.5))
p <- p + geom_abline(slope=1, size=0.5, colour=6, alpha=0.5)+theme_bw()
png(filename=args[4], height = 3000, width = 3000, res = 500, units = "px")
print(p)
dev.off()



