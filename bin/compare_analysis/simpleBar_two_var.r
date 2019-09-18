#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript


# load library
library('getopt');
require(ggplot2)


#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'outfile' , 'o', 1, "character",
  'x.col' , 'x', 1, "integer",
  'y.col' , 'y', 1, "integer",
  'z.col' , 'z', 1,"integer",
  'height' , 'H', 1, "integer",
  'width' , 'W', 1, "integer",
  'x.lab' , 'X', 1, "character",
  'y.lab' , 'Y', 1, "character",
  'lab.size' , 'l', 1, "integer",
  'number_label' ,'n',0,"logical",
  'axis.size' , 's', 1, "integer",
  'no.grid' , 'r', 0, "logical",
  'skip' , 'k', 1, "integer"
), byrow=TRUE, ncol=4);
opt = getopt(spec);


# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      1) Rscript simpleBar.r --infile in_simpleBar.data --outfile out_simpleBar.png \\
      --x.col 1 --y.col 2 --z.col 3 --x.lab \"x lab\" --y.lab \"y lab\"
      2) Rscript simpleBar.r --infile in_simpleBar.data --outfile out_simpleBar.png \\
      --x.col 1 --y.col 2 --z.col 3 --x.lab \"x lab\" --y.lab \"y lab\" \\
      -height 3000 --width 4000
      3) Rscript simpleBar.r --infile in_simpleBar.data --outfile out_simpleBar.png \\
      --x.col 1 --y.col 2 --z.col 3 --x.lab \"x lab\" --y.lab \"y lab\" \\
      --height 3000 --width 4000 --number_label
      
      Options: 
      --help		NULL 		get this help
      --infile 	character 	the input file [forced]
      --outfile 	character 	the filename for output graph,don't need expanded-name [forced]
      --x.col 	integer 	the col for x value [forced]
      --y.col 	integer 	the col for y value [forced]
      --z.col 	integer 	the col for z value [forced]
      --height 	integer 	the height of graph [optional, default: 3000]
      --width 	integer 	the width of graph [optional, default: 4000]
      --x.lab 	character 	the lab for x [forced]
      --y.lab 	character 	the lab for y [forced]
      --lab.size 	integer 	the font size of lab [optional, default: 15]
      --axis.size 	integer 	the font size of text for axis [optional, default: 10]
      --no.grid	NULL 		Do not drawing grid
      --number_label	NULL 		drawing value label on the bar
      --skip 		integer 	the number of line for skipping [optional, default: 0]
      \n")
  q(status=1);
}



# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) ) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )	{ print_usage(spec) }
if ( is.null(opt$outfile) )	{ print_usage(spec) }
if ( is.null(opt$x.col) )	{ print_usage(spec) }
if ( is.null(opt$y.col) )	{ print_usage(spec) }
if ( is.null(opt$z.col) )	{ print_usage(spec) }
if ( is.null(opt$x.lab) )	{ print_usage(spec) }
if ( is.null(opt$y.lab) )	{ print_usage(spec) }


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$axis.size ) )		{ opt$axis.size = 10 }
if ( is.null(opt$lab.size ) )		{ opt$lab.size = 15 }


#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile, skip=opt$skip,header = T)
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
  cat("Final Error: the format of infile is error, dim(data) is NULL\n")
  print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$x.col, opt$y.col,opt$z.col) ){
  cat("Final Error: max(x.col, y.col,z.col) > the col of infile\n")
  print_usage(spec)
}
# create df
df <- data.frame(x=as.factor(data[,opt$x.col]), y=data[,opt$y.col],z=as.factor(data[,opt$z.col]))



#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
p <-ggplot(data=df,aes(x,y,fill=z))


#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
#set bar
p <- p + geom_bar(stat="identity", position=position_dodge(),width=1.6)
# set lab
p <- p + xlab(opt$x.lab) + ylab(opt$y.lab) 
# set lab and axis test size
p <- p + theme(title = element_text(face="bold", size=opt$lab.size), 
               axis.text = element_text(face="bold", size=opt$axis.size))

#set number label
if ( !is.null(opt$number_label) ) {
  p <- p + geom_text(aes(label=round(y,2)),hjust=0.5, vjust=-0.5,size=4 )  #保留2位小数
}
# remove legend
p <- p + theme(legend.position = "none")

# set order
p <- p + scale_x_discrete(limits=factor(data[,opt$x.col])) 
# grid and background
if ( !is.null(opt$no.grid) ) {
  p <- p + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
                  panel.grid = element_blank())
}



#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
print(p)
dev.off()
png(filename=paste(opt$outfile,".png"), height=opt$height, width=opt$width, res=500, units="px")
print(p)
dev.off()







