#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript


# load library
library('getopt');
library("ggplot2");

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
  'height' , 'H', 1, "integer",
  'width' , 'W', 1, "integer",
  'legend.title' , 'l', 0, "character",
  'legend_title.size' ,'B' ,1,"integer",
  'legend_text.size' , 'C', 1, "integer",
  'no.grid' , 'r', 0, "logical",
  'percent_label.size' , 'D', 1, "integer",
  'percent_label.color','E',1,"character",
  'skip' , 'k', 1, "integer"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      1) Rscript simplePie.r --infile in_simplePie.data --outfile out_simplePie.png \\
         --x.col 1 --y.col 2 \\
      2) Rscript simplePie.r --infile in_simplePie.data --outfile out_simplePie.png \\
       --x.col 1 --y.col 2  --legend.title \"legend.title\" \\
      --title.size 10 --legend_title.size 20 --legend_text.size 13 --height 3000 --width 4000 \\
      3) Rscript simplePie.r --infile in_simplePie.data --outfile out_simplePie.png \\
      --x.col 1 --y.col 2  --legend.title \"legend.title\" \\
       --legend_title.size 20 --legend_text.size 13 --height 3000 --width 4000 \\
      --percent_label.size 28 --percent_label.color \"black\" \\
      
      Options: 
      --help		NULL 		    get this help
      --infile 	character 	the input file [forced]
      --outfile character 	the filename for output graph [forced]
      --x.col 	integer 	the col for x value [forced]
      --y.col 	integer 	the col for y value [forced]
      --height 	integer   	the height of graph [optional, default: 3000]
      --width 	integer 	  the width of graph [optional, default: 4000]
      --legend.title         character  the title of legend
      --legend_title.size    integer    the font size of legend.title[optional, default: 20]
      --legend_text.size     integer    the font size of[optional, default: 13]
      --no.grid 	NULL 	    	 Do not drawing grid
      --percent_label_size	 integer 		 the font size of percent label on the pie [optional, default: 10] 
      --percent_label.color  character   the color of percent label [optional, default: \"black\"]
      --skip 		             integer 	   the number of line for skipping [optional, default: 0]
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


#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$skip ) )		{ opt$skip = 0 }
if ( is.null(opt$height ) )		{ opt$height = 3000 }
if ( is.null(opt$width ) )		{ opt$width = 4000 }
if ( is.null(opt$legend.title ) ){opt$legend.title ="Type"}
if ( is.null(opt$legend_title.size ) )		{ opt$legend_title.size = 18 }
if ( is.null(opt$legend_text.size ) )		{ opt$legend_text.size = 10 }
if ( is.null(opt$percent_label.color) )		{ opt$percent_label.color = "black" }
if ( is.null(opt$percent_label.size ) )		{ opt$percent_label.size = 8 }


#-----------------------------------------------------------------
# reading data
#-----------------------------------------------------------------
# reading data
data <- read.table(opt$infile, skip=opt$skip)
# check dim
data.dim <- dim(data)
if ( is.null(data.dim) ){
  cat("Final Error: the format of infile is error, dim(data) is NULL\n")
  print_usage(spec)
}
# check col size
if ( data.dim[2] < max(opt$x.col, opt$y.col) ){
  cat("Final Error: max(x.col, y.col) > the col of infile\n")
  print_usage(spec)
}
# create df
percent_str<-paste(round(data[,opt$y.col]/sum(data[,opt$y.col])*100,1),"%", sep=" ")
percentage<-round(data[,opt$y.col]/sum(data[,opt$y.col]) * 100,1)
df <- data.frame(x=as.factor(data[,opt$x.col]), y=percentage,percent=percent_str)



#-----------------------------------------------------------------
# plot
#-----------------------------------------------------------------
# mian plot
pie <- ggplot(df, aes(x='',y=percentage,fill=as.factor(data[,opt$x.col]))) 

#-----------------------------------------------------------------
# theme
#-----------------------------------------------------------------
#remove lab
pie <- pie + xlab('') + ylab('')    

#set bar
pie <- pie + geom_bar(width = 1,stat="identity",alpha = 0.7)

# set legend
pie<- pie + coord_polar(theta = "y",start = pi/2) 
#
pie<- pie + scale_fill_discrete(name=opt$legend.title,
                      breaks=data[,opt$x.col],
                      labels=paste(data[,opt$x.col],"[",data[,opt$y.col],"]"))  
pie<- pie +theme_bw()                                                            #背景设置为空  20150730_1.png
pie<- pie + theme_linedraw()                                                     #标签设置阴影线 20150730_2.png
pie<- pie + theme_light()

pie <- pie +theme(legend.title = element_text(colour="black", size=opt$legend_title.size, face=2),     #图例标题字体设置
                 legend.text = element_text(colour="black", size = opt$legend_text.size, face =2))     #图例标题字体设置

# grid and background

pie <- pie +theme(panel.grid.major= element_blank(),                                  #消除主网格线                           
                  panel.grid.minor = element_blank(),                                 #消除次要网格线
                  panel.grid=element_blank(),                                         #将垂直的2条背景网格线线消除
                  panel.border=element_blank(),                                       #将绘图区边界去掉
                 # plot.background = element_rect(colour = "white") ,                 #绘制背景色
                 axis.ticks = element_blank(),                                        #清空坐标轴，去掉小胡子
                 axis.text.x = element_blank()                                        #去掉坐标轴
)
if ( !is.null(opt$no.grid) ) {
  pie <- pie + theme( panel.background = element_rect(colour="black", size=1, fill="white"),
                      panel.grid = element_blank())
}
#set title
pie<-pie + ggtitle("The proportion of each lncRNA")                             
pie<-pie + theme(plot.title = element_text(lineheight=.5,size = 20, face="bold"))


#set percent
pie<- pie +  geom_text(aes(y = percentage/2 + c(0, cumsum(percentage)[-length(percentage)]), label =percent),size=opt$percent_label.size,colour=opt$percent_label.color)




#-----------------------------------------------------------------
# output plot
#-----------------------------------------------------------------
pdf(file=paste(opt$outfile,".pdf",sep=""), height=opt$height*2/1000, width=opt$width*2/1000)
print(pie)
dev.off()
png(filename=(opt$outfile), height=opt$height, width=opt$width, res=500, units="px")
print(pie)
dev.off()
