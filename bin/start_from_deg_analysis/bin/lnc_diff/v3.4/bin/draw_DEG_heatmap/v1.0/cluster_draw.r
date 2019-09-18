#by luml@biomarker.com.cn
#time 2015/09/09
# load library

library('getopt');
library("pheatmap")

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'outfile' , 'o', 1, "character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

#-----------------------------------------------------------------
# define usage function
#-----------------------------------------------------------------
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
------------------------------------------------------------------------------------------
      Usage example: 
      Rscript cluster.r --infile in.xls --outfile out\\
      Options: 
      --help      NULL          get this help
      --infile    character     the input file [forced]
      --outfile   character     the filename for output file[forced]
      --skip      integer       the number of line for skipping [optional, default: 0]
-------------------------------------------------------------------------------------------
      \n")
  q(status=1);
}


# if help was asked for print a friendly message
if ( !is.null(opt$help) ) { print_usage(spec) }

# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
  
  
x=read.table(file=opt$infile,header =T,comment.char = " ")
row.names(x)=x[,1]
x2=x[,-1]#####去掉第一列
png(file = paste0(opt$outfile), width = 900, height = 1200);
par(cex=0.8)
par(mar = c(7, 5, 5, 6));
pheatmap(x2,cex=1.3,cluster_rows = T,na.rm=T,display_numbers=F, cluster_cols = T,scale ="row",
          color=colorRampPalette(c("green","black","red"))(100),
          show_rownames = F, ##2015-11-09  by linhj
          show_colnames = T,
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          clustering_method = "ward.D2",
          main=" ")
dev.off()

