#!/share/nas2/genome/biosoft/R/2.15.1/lib64/R/bin/Rscript
#####################################################################
# Copyright 2015, BMK
#
# Author:tengh <tengh@biomarker.com.cn>
#
# Function: draw genomewide cytosine coverage distribution map
#
# Modify date: 20150819
# Note: delete group label
#       reset opt$color="#263C8B,#4E74A6,#BDBF78,#BFA524"
#####################################################################
library("grid")
library("RColorBrewer")
library("scales")
library("gtable")
lo = function(rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, treeheight_col, treeheight_row, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, gaps_row, gaps_col, ...){
    # Get height of colnames and length of rownames
    if(!is.null(coln[1])){
        t = c(coln, colnames(annotation_row))
        longest_coln = which.max(strwidth(t, units = 'in'))
        gp = list(fontsize = fontsize_col, ...)
        coln_height = unit(1, "grobheight", textGrob(t[longest_coln], rot = 90, gp = do.call(gpar, gp))) + unit(10, "bigpts")
    }else{
        coln_height = unit(5, "bigpts")
    }
    
    if(!is.null(rown[1])){
        #t = c(rown, colnames(annotation_col))
	t = c(rown, "") #20150819
        longest_rown = which.max(strwidth(t, units = 'in'))
        gp = list(fontsize = fontsize_row, ...)
        rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], gp = do.call(gpar, gp))) + unit(10, "bigpts")
    }else{
        rown_width = unit(5, "bigpts")
    }
    
    gp = list(fontsize = fontsize, ...)
    # Legend position
    if(!is.na(legend[1])){
        longest_break = which.max(nchar(names(legend)))
        longest_break = unit(1.1, "grobwidth", textGrob(as.character(names(legend))[longest_break], gp = do.call(gpar, gp)))
        title_length = unit(1.1, "grobwidth", textGrob("Scale", gp = gpar(fontface = "bold", ...)))
        legend_width = unit(12, "bigpts") + longest_break * 1.2
        legend_width = max(title_length, legend_width)
    }else{
        legend_width = unit(0, "bigpts")
    }
    
    # Set main title height
    if(is.na(main)){
        main_height = unit(0, "npc")
    }else{
        main_height = unit(1.5, "grobheight", textGrob(main, gp = gpar(fontsize = 1.3 * fontsize, ...)))
    }
    
    # Column annotations
    textheight = unit(fontsize, "bigpts")
    
    if(!is.na(annotation_col[[1]][1])){
        # Column annotation height 
        annot_col_height = ncol(annotation_col) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")
        
        # Width of the correponding legend
        #t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col)) 
	t = c(as.vector(as.matrix(annotation_col)),"") #20150819
        annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
        if(!annotation_legend){
            annot_col_legend_width = unit(0, "npc")
        }
    }else{
        annot_col_height = unit(0, "bigpts")
        annot_col_legend_width = unit(0, "bigpts")
    }
    
    # Row annotations
    if(!is.na(annotation_row[[1]][1])){
        # Row annotation width 
        annot_row_width = ncol(annotation_row) * (textheight + unit(2, "bigpts")) + unit(2, "bigpts")
        
        # Width of the correponding legend
        t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row)) 
        annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], gp = gpar(...))) + unit(12, "bigpts")
        if(!annotation_legend){
            annot_row_legend_width = unit(0, "npc")
        }
    }else{
        annot_row_width = unit(0, "bigpts")
        annot_row_legend_width = unit(0, "bigpts")
    }
    
    annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
    
    # Tree height
    treeheight_col = unit(treeheight_col, "bigpts") + unit(5, "bigpts")
    treeheight_row = unit(treeheight_row, "bigpts") + unit(5, "bigpts") 
    
    # Set cell sizes
    if(is.na(cellwidth)){
        mat_width = unit(1, "npc") - rown_width - legend_width - treeheight_row - annot_row_width - annot_legend_width 
    }else{
        mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * unit(4, "bigpts")
    }
    
    if(is.na(cellheight)){
        mat_height = unit(1, "npc") - main_height - coln_height - treeheight_col - annot_col_height
    }else{
        mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * unit(4, "bigpts")
    }    
    
    # Produce gtable
    gt = gtable(widths = unit.c(treeheight_row, annot_row_width, mat_width, rown_width, legend_width, annot_legend_width), heights = unit.c(main_height, treeheight_col, annot_col_height, mat_height, coln_height), vp = viewport(gp = do.call(gpar, gp)))
    
    cw = convertWidth(mat_width - (length(gaps_col) * unit(4, "bigpts")), "bigpts", valueOnly = T) / ncol
    ch = convertHeight(mat_height - (length(gaps_row) * unit(4, "bigpts")), "bigpts", valueOnly = T) / nrow
    
    # Return minimal cell dimension in bigpts to decide if borders are drawn
    mindim = min(cw, ch) 
    
    res = list(gt = gt, mindim = mindim)
    
    return(res)
}

find_coordinates = function(n, gaps, m = 1:n){
    if(length(gaps) == 0){
        return(list(coord = unit(m / n, "npc"), size = unit(1 / n, "npc") ))
    }
    
    if(max(gaps) > n){
        stop("Gaps do not match with matrix size")
    }
    
    size = (1 / n) * (unit(1, "npc") - length(gaps) * unit("4", "bigpts"))
    
    gaps2 = apply(sapply(gaps, function(gap, x){x > gap}, m), 1, sum) 
    coord = m * size + (gaps2 * unit("4", "bigpts"))
    
    return(list(coord = coord, size = size))
}

draw_dendrogram = function(hc, gaps, horizontal = T){
    h = hc$height / max(hc$height) / 1.05
    m = hc$merge
    o = hc$order
    n = length(o)

    m[m > 0] = n + m[m > 0] 
    m[m < 0] = abs(m[m < 0])

    dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL, c("x", "y"))) 
    dist[1:n, 1] = 1 / n / 2 + (1 / n) * (match(1:n, o) - 1)

    for(i in 1:nrow(m)){
        dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1]) / 2
        dist[n + i, 2] = h[i]
    }
    
    draw_connection = function(x1, x2, y1, y2, y){
        res = list(
            x = c(x1, x1, x2, x2),
            y = c(y1, y, y, y2)
        )
        
        return(res)
    }
    
    x = rep(NA, nrow(m) * 4)
    y = rep(NA, nrow(m) * 4)
    id = rep(1:nrow(m), rep(4, nrow(m)))
    
    for(i in 1:nrow(m)){
        c = draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1], dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
        k = (i - 1) * 4 + 1
        x[k : (k + 3)] = c$x
        y[k : (k + 3)] = c$y
    }
    
    x = find_coordinates(n, gaps, x * n)$coord
    y = unit(y, "npc")
    
    if(!horizontal){
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a
    }
    res = polylineGrob(x = x, y = y, id = id)
    
    return(res)
}

draw_matrix = function(matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color){
    n = nrow(matrix)
    m = ncol(matrix)
    
    coord_x = find_coordinates(m, gaps_cols)
    coord_y = find_coordinates(n, gaps_rows)
    
    x = coord_x$coord - 0.5 * coord_x$size
    y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)
    
    coord = expand.grid(y = y, x = x)
    
    res = gList()
    
    res[["rect"]] = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = coord_y$size, gp = gpar(fill = matrix, col = border_color))
    
    if(attr(fmat, "draw")){
        res[["text"]] = textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
    }
    
    res = gTree(children = res)
    
    return(res)
}

draw_colnames = function(coln, gaps, ...){
    coord = find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 0, rot = 270, gp = gpar(...))
    
    return(res)
}

draw_rownames = function(rown, gaps, ...){
    coord = find_coordinates(length(rown), gaps)
    y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
    
    res = textGrob(rown, x = unit(3, "bigpts"), y = y, vjust = 0.5, hjust = 0, gp = gpar(...))
    
    return(res)
}

draw_legend = function(color, breaks, legend, ...){
    height = min(unit(1, "npc"), unit(150, "bigpts"))
    #message(paste(c("legend=",legend),collapse = "\t"))
    #message(paste(c("min(breaks)=",min(breaks)),collapse = "\t"))
    legend_pos = (legend - min(breaks)) / (max(breaks) - min(breaks))
    legend_pos = height * legend_pos + (unit(1, "npc") - height)
    breaks = (breaks - min(breaks)) / (max(breaks) - min(breaks))
    breaks = height * breaks + (unit(1, "npc") - height)
    
    h = breaks[-1] - breaks[-length(breaks)]
    
    rect = rectGrob(x = 0, y = breaks[-length(breaks)], width = unit(10, "bigpts"), height = h, hjust = 0, vjust = 0, gp = gpar(fill = color, col = "#FFFFFF00"))
    text = textGrob(names(legend), x = unit(14, "bigpts"), y = legend_pos, hjust = 0, gp = gpar(...))
    
    res = grobTree(rect, text)
    
    return(res)
}

convert_annotations = function(annotation, annotation_colors){
    new = annotation
    for(i in 1:ncol(annotation)){
        a = annotation[, i]
        b = annotation_colors[[colnames(annotation)[i]]]
        if(is.character(a) | is.factor(a)){
            a = as.character(a)
            if(length(setdiff(a, names(b))) > 0){
                stop(sprintf("Factor levels on variable %s do not match with annotation_colors", colnames(annotation)[i]))
            }
            new[, i] = b[a]
        }else{
            a = cut(a, breaks = 100)
            new[, i] = colorRampPalette(b)(100)[a]
        }
    }
    return(as.matrix(new))
}

draw_annotations = function(converted_annotations, border_color, gaps, fontsize, horizontal){
    n = ncol(converted_annotations)
    m = nrow(converted_annotations)
    
    coord_x = find_coordinates(m, gaps)
    
    x = coord_x$coord - 0.5 * coord_x$size
    
    # y = cumsum(rep(fontsize, n)) - 4 + cumsum(rep(2, n))
    y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1 
    y = unit(y, "bigpts")
    
    if(horizontal){
        coord = expand.grid(x = x, y = y)
        res = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = unit(fontsize, "bigpts"), gp = gpar(fill = converted_annotations, col = border_color))
    }else{
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a
        
        coord = expand.grid(y = y, x = x)
        res = rectGrob(x = coord$x, y = coord$y, width = unit(fontsize, "bigpts"), height = coord_x$size, gp = gpar(fill = converted_annotations, col = border_color))
    }
    
    return(res)
}

draw_annotation_names = function(annotations, fontsize, horizontal){
    n = ncol(annotations)
    
    x = unit(3, "bigpts")
    
    y = cumsum(rep(fontsize, n)) + cumsum(rep(2, n)) - fontsize / 2 + 1 
    y = unit(y, "bigpts")
    
    if(horizontal){
        res = textGrob(colnames(annotations), x = x, y = y, hjust = 0, gp = gpar(fontsize = fontsize, fontface = 2))
    }else{
        a = x
        x = unit(1, "npc") - y
        y = unit(1, "npc") - a
        
        res = textGrob(colnames(annotations), x = x, y = y, vjust = 0.5, hjust = 0, rot = 270, gp = gpar(fontsize = fontsize, fontface = 2))
    }
    
    return(res)
}

draw_annotation_legend = function(annotation, annotation_colors, border_color, ...){
    y = unit(1, "npc")
    text_height = unit(1, "grobheight", textGrob("FGH", gp = gpar(...)))
    
    res = gList()
    for(i in names(annotation)){
        res[[i]] = textGrob(i, x = 0, y = y, vjust = 1, hjust = 0, gp = gpar(fontface = "bold", ...))
        
        y = y - 1.5 * text_height
        if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
            n = length(annotation_colors[[i]])
            yy = y - (1:n - 1) * 2 * text_height
            
            res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = 2 * text_height, width = 2 * text_height, gp = gpar(col = border_color, fill = annotation_colors[[i]]))
            res[[paste(i, "t")]] = textGrob(names(annotation_colors[[i]]), x = text_height * 2.4, y = yy - text_height, hjust = 0, vjust = 0.5, gp = gpar(...))
            
            y = y - n * 2 * text_height
            
        }else{
            yy = y - 8 * text_height + seq(0, 1, 0.25)[-1] * 8 * text_height
            h = 8 * text_height * 0.25
            
            res[[paste(i, "r")]] = rectGrob(x = unit(0, "npc"), y = yy, hjust = 0, vjust = 1, height = h, width = 2 * text_height, gp = gpar(col = NA, fill = colorRampPalette(annotation_colors[[i]])(4)))
            res[[paste(i, "r2")]] = rectGrob(x = unit(0, "npc"), y = y, hjust = 0, vjust = 1, height = 8 * text_height, width = 2 * text_height, gp = gpar(col = border_color))
            
            txt = rev(range(grid.pretty(range(annotation[[i]], na.rm = TRUE))))
            yy = y - c(1, 7) * text_height
            res[[paste(i, "t")]]  = textGrob(txt, x = text_height * 2.4, y = yy, hjust = 0, vjust = 0.5, gp = gpar(...))
            y = y - 8 * text_height
        }
        y = y - 1.5 * text_height
    }
    
    res = gTree(children = res)
    
    return(res)
}

draw_main = function(text, ...){
    res = textGrob(text, gp = gpar(fontface = "bold", ...))
    
    return(res)
}

vplayout = function(x, y){
    return(viewport(layout.pos.row = x, layout.pos.col = y))
}

heatmap_motor = function(matrix, border_color, cellwidth, cellheight, tree_col, tree_row, treeheight_col, treeheight_row, filename, width, height, breaks, color, legend, annotation_row, annotation_col, annotation_colors, annotation_legend, main, fontsize, fontsize_row, fontsize_col, fmat, fontsize_number, number_color, gaps_col, gaps_row, labels_row, labels_col, ...){
    # Set layout
    lo = lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, gaps_row = gaps_row, gaps_col = gaps_col,  ...)
    
    res = lo$gt
    mindim = lo$mindim
    
    if(!is.na(filename)){
        if(is.na(height)){
            height = convertHeight(gtable_height(res), "inches", valueOnly = T)
        }
        if(is.na(width)){
            width = convertWidth(gtable_width(res), "inches", valueOnly = T)
        }
        
        # Get file type
        r = regexpr("\\.[a-zA-Z]*$", filename)
        if(r == -1) stop("Improper filename")
        ending = substr(filename, r + 1, r + attr(r, "match.length"))

        f = switch(ending,
            pdf = function(x, ...) pdf(x, ...),
            png = function(x, ...) png(x, units = "in", res = 500, ...),
            jpeg = function(x, ...) jpeg(x, units = "in", res = 500, ...),
            jpg = function(x, ...) jpeg(x, units = "in", res = 500, ...),
            tiff = function(x, ...) tiff(x, units = "in", res = 500, compression = "lzw", ...),
            bmp = function(x, ...) bmp(x, units = "in", res = 500, ...),
            stop("File type should be: pdf, png, bmp, jpg, tiff")
        )
        
        # print(sprintf("height:%f width:%f", height, width))
        
        # gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)

        f(filename, height = height, width = width)
        gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, border_color = border_color, tree_col = tree_col, tree_row = tree_row, treeheight_col = treeheight_col, treeheight_row = treeheight_row, breaks = breaks, color = color, legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, annotation_colors = annotation_colors, annotation_legend = annotation_legend, filename = NA, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number =  fontsize_number, number_color = number_color, labels_row = labels_row, labels_col = labels_col, gaps_col = gaps_col, gaps_row = gaps_row, ...)
        grid.draw(gt)
        dev.off()
        
        return(NULL)
    }
    
    # Omit border color if cell size is too small 
    if(mindim < 3) border_color = NA
    
    # Draw title
    if(!is.na(main)){
        elem = draw_main(main, fontsize = 1.3 * fontsize, ...)
        res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main")
    }
    
    # Draw tree for the columns
    if(!is.na(tree_col[[1]][1]) & treeheight_col != 0){
        elem = draw_dendrogram(tree_col, gaps_col, horizontal = T)
        res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
    }
    
    # Draw tree for the rows
    if(!is.na(tree_row[[1]][1]) & treeheight_row != 0){
        elem = draw_dendrogram(tree_row, gaps_row, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
    }
    
    # Draw matrix
    elem = draw_matrix(matrix, border_color, gaps_row, gaps_col, fmat, fontsize_number, number_color)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", name = "matrix")
    
    # Draw colnames
    if(length(labels_col) != 0){
        pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, ...)
        elem = do.call(draw_colnames, pars)
        res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", name = "col_names")
    }
    
    # Draw rownames
    if(length(labels_row) != 0){
        pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, ...)
        elem = do.call(draw_rownames, pars)
        res = gtable_add_grob(res, elem, t = 4, l = 4, clip = "off", name = "row_names")
    }
    
    # Draw annotation tracks on cols
    if(!is.na(annotation_col[[1]][1])){
        # Draw tracks
        converted_annotation = convert_annotations(annotation_col, annotation_colors)

        elem = draw_annotations(converted_annotation, border_color, gaps_col, fontsize, horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", name = "col_annotation")
        
        # Draw names
	annotation_col.tmp<-annotation_col
        colnames(annotation_col.tmp)<-""
        elem = draw_annotation_names(annotation_col.tmp, fontsize, horizontal = T)
        res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", name = "row_annotation_names")
        
    }
    
    # Draw annotation tracks on rows
    if(!is.na(annotation_row[[1]][1])){
        # Draw tracks
        converted_annotation = convert_annotations(annotation_row, annotation_colors)
        elem = draw_annotations(converted_annotation, border_color, gaps_row, fontsize, horizontal = F)
        res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", name = "row_annotation")
        
        # Draw names
        elem = draw_annotation_names(annotation_row, fontsize, horizontal = F)
        res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", name = "row_annotation_names")
    }
    
    # Draw annotation legend
    annotation = c(annotation_col[length(annotation_col):1], annotation_row[length(annotation_row):1])
    annotation = annotation[unlist(lapply(annotation, function(x) !is.na(x[1])))]
    
    if(length(annotation) > 0 & annotation_legend){
        elem = draw_annotation_legend(annotation, annotation_colors, border_color, fontsize = fontsize, ...)
        
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, clip = "off", name = "annotation_legend")
    }
    
    # Draw legend
    if(!is.na(legend[1])){
        elem = draw_legend(color, breaks, legend, fontsize = fontsize, ...)
        
        t = ifelse(is.null(labels_row), 4, 3)
        res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, clip = "off", name = "legend")
    }
    
    return(res)
}

generate_breaks = function(x, n, center = F){
    if(center){
        m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
        res = seq(-m, m, length.out = n + 1)
    }else{
        res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
    }
    
    return(res)
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
    return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

scale_colours = function(mat, col = rainbow(10), breaks = NA){
    mat = as.matrix(mat)
    return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

cluster_mat = function(mat, distance, method){
    if(!(method %in% c("ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid"))){
        stop("clustering method has to one form the list: 'ward', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
    }
    if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
        stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
    }
    if(distance[1] == "correlation"){
        d = as.dist(1 - cor(t(mat)))
    }else{
        if(class(distance) == "dist"){
            d = distance
        }else{
            d = dist(mat, method = distance)
        }
    }
    
    return(hclust(d, method = method))
}

scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

scale_mat = function(mat, scale){
    if(!(scale %in% c("none", "row", "column"))){
        stop("scale argument shoud take values: 'none', 'row' or 'column'")
    }
    mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
    return(mat)
}

generate_annotation_colours = function(annotation, annotation_colors, drop){
    if(is.na(annotation_colors)[[1]][1]){
        annotation_colors = list()
    }
    count = 0
    for(i in 1:length(annotation)){
        if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
            if (is.factor(annotation[[i]]) & !drop){
                count = count + length(levels(annotation[[i]]))
            }else{
                count = count + length(unique(annotation[[i]]))
            }
        }
    }
    
    factor_colors = dscale(factor(1:count), hue_pal(l = 75))
    
    set.seed(3453)
    
    cont_counter = 2
    for(i in 1:length(annotation)){
        if(!(names(annotation)[i] %in% names(annotation_colors))){
            if(is.character(annotation[[i]]) | is.factor(annotation[[i]])){
                n = length(unique(annotation[[i]]))
                if (is.factor(annotation[[i]]) & !drop){
                    n = length(levels(annotation[[i]]))
                }
                ind = sample(1:length(factor_colors), n)
                annotation_colors[[names(annotation)[i]]] = factor_colors[ind]
                l = levels(as.factor(annotation[[i]]))
                l = l[l %in% unique(annotation[[i]])]
                if (is.factor(annotation[[i]]) & !drop){
                    l = levels(annotation[[i]])
                }
                
                names(annotation_colors[[names(annotation)[i]]]) = l
                factor_colors = factor_colors[-ind]
            }else{
                annotation_colors[[names(annotation)[i]]] = brewer_pal("seq", cont_counter)(5)[1:4]
                cont_counter = cont_counter + 1
            }
        }
    }
    return(annotation_colors)
}

kmeans_pheatmap = function(mat, k = min(nrow(mat), 150), sd_limit = NA, ...){
    # Filter data
    if(!is.na(sd_limit)){
        s = apply(mat, 1, sd)
        mat = mat[s > sd_limit, ]    
    }
    
    # Cluster data
    set.seed(1245678)
    km = kmeans(mat, k, iter.max = 100)
    mat2 = km$centers
    
    # Compose rownames
    t = table(km$cluster)
    rownames(mat2) = sprintf("cl%s_size_%d", names(t), t)
    
    # Draw heatmap
    pheatmap2(mat2, ...)
}

find_gaps = function(tree, cutree_n){
    v = cutree(tree, cutree_n)[tree$order]
    gaps = which((v[-1] - v[-length(v)]) != 0)
    
}
 
#' A function to draw clustered heatmaps.
#' 
#' A function to draw clustered heatmaps where one has better control over some graphical 
#' parameters such as cell size, etc. 
#' 
#' The function also allows to aggregate the rows using kmeans clustering. This is 
#' advisable if number of rows is so big that R cannot handle their hierarchical 
#' clustering anymore, roughly more than 1000. Instead of showing all the rows 
#' separately one can cluster the rows in advance and show only the cluster centers. 
#' The number of clusters can be tuned with parameter kmeans_k.
#'
#' @param mat numeric matrix of the values to be plotted.
#' @param color vector of colors used in heatmap.
#' @param kmeans_k the number of kmeans clusters to make, if we want to agggregate the 
#' rows before drawing heatmap. If NA then the rows are not aggregated.
#' @param breaks a sequence of numbers that covers the range of values in mat and is one 
#' element longer than color vector. Used for mapping values to colors. Useful, if needed 
#' to map certain values to certain colors, to certain values. If value is NA then the 
#' breaks are calculated automatically.
#' @param border_color color of cell borders on heatmap, use NA if no border should be 
#' drawn.
#' @param cellwidth individual cell width in points. If left as NA, then the values 
#' depend on the size of plotting window.
#' @param cellheight individual cell height in points. If left as NA, 
#' then the values depend on the size of plotting window.
#' @param scale character indicating if the values should be centered and scaled in 
#' either the row direction or the column direction, or none. Corresponding values are 
#' \code{"row"}, \code{"column"} and \code{"none"}
#' @param cluster_rows boolean values determining if rows should be clustered,
#' @param cluster_cols boolean values determining if columns should be clustered.
#' @param clustering_distance_rows distance measure used in clustering rows. Possible 
#' values are \code{"correlation"} for Pearson correlation and all the distances 
#' supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. If the value is none 
#' of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols distance measure used in clustering columns. Possible 
#' values the same as for clustering_distance_rows.
#' @param clustering_method clustering method used. Accepts the same values as 
#' \code{\link{hclust}}.
#' @param cutree_rows number of clusters the rows are divided into, based on the
#'  hierarchical clustering (using cutree), if rows are not clustered, the 
#' argument is ignored
#' @param cutree_cols similar to \code{cutree_rows}, but for columns
#' @param treeheight_row the height of a tree for rows, if these are clustered. 
#' Default value 50 points.
#' @param treeheight_col the height of a tree for columns, if these are clustered. 
#' Default value 50 points.
#' @param legend logical to determine if legend should be drawn or not.
#' @param legend_breaks vector of breakpoints for the legend.
#' @param legend_labels vector of labels for the \code{legend_breaks}.
#' @param annotation_row data frame that specifies the annotations shown on left
#'  side of the heatmap. Each row defines the features for a specific row. The 
#' rows in the data and in the annotation are matched using corresponding row
#'  names. Note that color schemes takes into account if variable is continuous
#'  or discrete.
#' @param annotation_col similar to annotation_row, but for columns. 
#' @param annotation deprecated parameter that currently sets the annotation_col if it is missing
#' @param annotation_colors list for specifying annotation_row and 
#' annotation_col track colors manually. It is  possible to define the colors 
#' for only some of the features. Check examples for  details.
#' @param annotation_legend boolean value showing if the legend for annotation 
#' tracks should be drawn. 
#' @param drop_levels logical to determine if unused levels are also shown in 
#' the legend
#' @param show_rownames boolean specifying if column names are be shown.
#' @param show_colnames boolean specifying if column names are be shown.
#' @param main the title of the plot
#' @param fontsize base fontsize for the plot 
#' @param fontsize_row fontsize for rownames (Default: fontsize) 
#' @param fontsize_col fontsize for colnames (Default: fontsize) 
#' @param display_numbers logical determining if the numeric values are also printed to 
#' the cells. If this is a matrix (with same dimensions as original matrix), the contents
#' of the matrix are shown instead of original values.
#' @param number_format format strings (C printf style) of the numbers shown in cells. 
#' For example "\code{\%.2f}" shows 2 decimal places and "\code{\%.1e}" shows exponential 
#' notation (see more in \code{\link{sprintf}}).
#' @param number_color color of the text    
#' @param fontsize_number fontsize of the numbers displayed in cells
#' @param gaps_row vector of row indices that show shere to put gaps into
#'  heatmap. Used only if the rows are not clustered. See \code{cutree_row}
#'  to see how to introduce gaps to clustered rows. 
#' @param gaps_col similar to gaps_row, but for columns.
#' @param labels_row custom labels for rows that are used instead of rownames.
#' @param labels_col similar to labels_row, but for columns.
#' @param filename file path where to save the picture. Filetype is decided by 
#' the extension in the path. Currently following formats are supported: png, pdf, tiff,
#'  bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is 
#' calculated so that the plot would fit there, unless specified otherwise.
#' @param width manual option for determining the output file width in inches.
#' @param height manual option for determining the output file height in inches.
#' @param silent do not draw the plot (useful when using the gtable output)
#' @param \dots graphical parameters for the text used in plot. Parameters passed to 
#' \code{\link{grid.text}}, see \code{\link{gpar}}. 
#' 
#' @return 
#' Invisibly a list of components 
#' \itemize{
#'     \item \code{tree_row} the clustering of rows as \code{\link{hclust}} object 
#'     \item \code{tree_col} the clustering of columns as \code{\link{hclust}} object
#'     \item \code{kmeans} the kmeans clustering of rows if parameter \code{kmeans_k} was 
#' specified 
#' }
#' 
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#' # Create test matrix
#' test = matrix(rnorm(200), 20, 10)
#' test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
#' test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
#' test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
#' colnames(test) = paste("Test", 1:10, sep = "")
#' rownames(test) = paste("Gene", 1:20, sep = "")
#' 
#' # Draw heatmaps
#' pheatmap2(test)
#' pheatmap2(test, kmeans_k = 2)
#' pheatmap2(test, scale = "row", clustering_distance_rows = "correlation")
#' pheatmap2(test, color = colorRampPalette(c("navy", "white", "firebrick3"))(50))
#' pheatmap2(test, cluster_row = FALSE)
#' pheatmap2(test, legend = FALSE)
#' 
#' # Show text within cells
#' pheatmap2(test, display_numbers = TRUE)
#' pheatmap2(test, display_numbers = TRUE, number_format = "\%.1e")
#' pheatmap2(test, display_numbers = matrix(ifelse(test > 5, "*", ""), nrow(test)))
#' pheatmap2(test, cluster_row = FALSE, legend_breaks = -1:4, legend_labels = c("0",
#' "1e-4", "1e-3", "1e-2", "1e-1", "1"))
#' 
#' # Fix cell sizes and save to file with correct size
#' pheatmap2(test, cellwidth = 15, cellheight = 12, main = "Example heatmap")
#' pheatmap2(test, cellwidth = 15, cellheight = 12, fontsize = 8, filename = "test.pdf")
#' 
#' # Generate annotations for rows and columns
#' annotation_col = data.frame(
#'                     CellType = factor(rep(c("CT1", "CT2"), 5)), 
#'                     Time = 1:5
#'                 )
#' rownames(annotation_col) = paste("Test", 1:10, sep = "")
#' 
#' annotation_row = data.frame(
#'                     GeneClass = factor(rep(c("Path1", "Path2", "Path3"), c(10, 4, 6)))
#'                 )
#' rownames(annotation_row) = paste("Gene", 1:20, sep = "")
#' 
#' # Display row and color annotations
#' pheatmap2(test, annotation_col = annotation_col)
#' pheatmap2(test, annotation_col = annotation_col, annotation_legend = FALSE)
#' pheatmap2(test, annotation_col = annotation_col, annotation_row = annotation_row)
#' 
#' 
#' # Specify colors
#' ann_colors = list(
#'     Time = c("white", "firebrick"),
#'     CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
#'     GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E")
#' )
#' 
#' pheatmap2(test, annotation_col = annotation_col, annotation_colors = ann_colors, main = "Title")
#' pheatmap2(test, annotation_col = annotation_col, annotation_row = annotation_row, 
#'          annotation_colors = ann_colors)
#' pheatmap2(test, annotation_col = annotation_col, annotation_colors = ann_colors[2]) 
#' 
#' # Gaps in heatmaps
#' pheatmap2(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14))
#' pheatmap2(test, annotation_col = annotation_col, cluster_rows = FALSE, gaps_row = c(10, 14), 
#'          cutree_col = 2)
#' 
#' # Show custom strings as row/col names
#' labels_row = c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", 
#' "", "", "Il10", "Il15", "Il1b")
#' 
#' pheatmap2(test, annotation_col = annotation_col, labels_row = labels_row)
#' 
#' # Specifying clustering from distance matrix
#' drows = dist(test, method = "minkowski")
#' dcols = dist(t(test), method = "minkowski")
#' pheatmap2(test, clustering_distance_rows = drows, clustering_distance_cols = dcols)
#' 
#' @export
pheatmap2 = function(mat,color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, border_color = "grey60", cellwidth = NA, cellheight = NA, scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete", cutree_rows = NA, cutree_cols = NA,  treeheight_row = ifelse(cluster_rows, 50, 0), treeheight_col = ifelse(cluster_cols, 50, 0), legend = TRUE, legend_breaks = NA, legend_labels = NA, annotation_row = NA, annotation_col = NA, annotation = NA, annotation_colors = NA, annotation_legend =FALSE, drop_levels = TRUE, show_rownames = T, show_colnames = T, main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, display_numbers = F, number_format = "%.2f", number_color = "grey30", fontsize_number = 0.8 * fontsize, gaps_row = NULL, gaps_col = NULL, labels_row = NULL, labels_col = NULL, filename = NA, width = NA, height = NA, silent = FALSE, ...){
    
    # Set labels
    if(is.null(labels_row)){
        labels_row = rownames(mat)
    }
    if(is.null(labels_col)){
        labels_col = colnames(mat)
    }
    
    # Preprocess matrix
    mat = as.matrix(mat)
    if(scale != "none"){
        mat = scale_mat(mat, scale)
        if(is.na(breaks)){
            breaks = generate_breaks(mat, length(color), center = T)
        }
    }
    
    
    # Kmeans
    if(!is.na(kmeans_k)){
        # Cluster data
        km = kmeans(mat, kmeans_k, iter.max = 100)
        mat = km$centers
        
        # Compose rownames
        t = table(km$cluster)
        labels_row = sprintf("Cluster: %s Size: %d", names(t), t)
    }else{
        km = NA
    }
    
    # Format numbers to be displayed in cells
    if(is.matrix(display_numbers) | is.data.frame(display_numbers)){
        if(nrow(display_numbers) != nrow(mat) | ncol(display_numbers) != ncol(mat)){
            stop("If display_numbers provided as matrix, its dimensions have to match with mat")
        }
        
        display_numbers = as.matrix(display_numbers)
        fmat = matrix(as.character(display_numbers), nrow = nrow(display_numbers), ncol = ncol(display_numbers))
        fmat_draw = TRUE
        
    }else{
        if(display_numbers){
            fmat = matrix(sprintf(number_format, mat), nrow = nrow(mat), ncol = ncol(mat))
            fmat_draw = TRUE
        }else{
            fmat = matrix(NA, nrow = nrow(mat), ncol = ncol(mat))
            fmat_draw = FALSE
        }
    }
    
    # Do clustering
    if(cluster_rows){
        tree_row = cluster_mat(mat, distance = clustering_distance_rows, method = clustering_method)
        mat = mat[tree_row$order, , drop = FALSE]
        fmat = fmat[tree_row$order, , drop = FALSE]
        labels_row = labels_row[tree_row$order]
        if(!is.na(cutree_rows)){
            gaps_row = find_gaps(tree_row, cutree_rows)
        }else{
            gaps_row = NULL
        }
    }else{
        tree_row = NA
        treeheight_row = 0
    }
    
    if(cluster_cols){
        tree_col = cluster_mat(t(mat), distance = clustering_distance_cols, method = clustering_method)
        mat = mat[, tree_col$order, drop = FALSE]
        fmat = fmat[, tree_col$order, drop = FALSE]
        labels_col = labels_col[tree_col$order]
        if(!is.na(cutree_cols)){
            gaps_col = find_gaps(tree_col, cutree_cols)
        } else{
            gaps_col = NULL
        }
    }else{
        tree_col = NA
        treeheight_col = 0
    }
    
    attr(fmat, "draw") = fmat_draw
    
    # Colors and scales
    if(!is.na(legend_breaks[1]) & !is.na(legend_labels[1])){
        if(length(legend_breaks) != length(legend_labels)){
            stop("Lengths of legend_breaks and legend_labels must be the same")
        }
    }
    
    
    if(is.na(breaks[1])){
        breaks = generate_breaks(as.vector(mat), length(color))
    }
    if (legend & is.na(legend_breaks[1])) {
        legend = grid.pretty(range(as.vector(breaks)))
        names(legend) = legend
    }else if(legend & !is.na(legend_breaks[1])){
        legend = legend_breaks[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
        
        if(!is.na(legend_labels[1])){
            legend_labels = legend_labels[legend_breaks >= min(breaks) & legend_breaks <= max(breaks)]
            names(legend) = legend_labels
        }else{
            names(legend) = legend
        }
    }else {
        legend = NA
    }
    mat = scale_colours(mat, col = color, breaks = breaks)
    
    # Preparing annotations
    if(is.na(annotation_col[[1]][1]) & !is.na(annotation[[1]][1])){
        annotation_col = annotation
    }
    # Select only the ones present in the matrix
    if(!is.na(annotation_col[[1]][1])){
        annotation_col = annotation_col[colnames(mat), , drop = F]
    }
    
    if(!is.na(annotation_row[[1]][1])){
        annotation_row = annotation_row[rownames(mat), , drop = F]
    }
    
    annotation = c(annotation_row, annotation_col)
    annotation = annotation[unlist(lapply(annotation, function(x) !is.na(x[1])))]
    if(length(annotation) != 0){
        annotation_colors = generate_annotation_colours(annotation, annotation_colors, drop = drop_levels)
    } else{
        annotation_colors = NA
    }
    
    if(!show_rownames){
        labels_row = NULL
    }
    
    if(!show_colnames){
        labels_col = NULL
    }
    
    # Draw heatmap
    gt = heatmap_motor(mat, border_color = border_color, cellwidth = cellwidth, cellheight = cellheight, treeheight_col = treeheight_col, treeheight_row = treeheight_row, tree_col = tree_col, tree_row = tree_row, filename = filename, width = width, height = height, breaks = breaks, color = color, legend = legend, annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors, annotation_legend = annotation_legend, main = main, fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col, fmat = fmat, fontsize_number = fontsize_number, number_color = number_color, gaps_row = gaps_row, gaps_col = gaps_col, labels_row = labels_row, labels_col = labels_col)
    
    if(is.na(filename) & !silent){
        grid.newpage()
        grid.draw(gt)
    }
    
    invisible(list(tree_row = tree_row, tree_col = tree_col, kmeans = km, gtable = gt))
}


# load library
library('getopt');
#opt<-data.frame(infile="E:/R_workplace/20150626heatmap/T1_T2_vs_T3_T4.DEG.final.cluster",groupfile="E:/R_workplace/20150626heatmap/groupfile.heatmap")
.sourcePath<-"/share/nas1/tengh/research/Rsource/"

#-----------------------------------------------------------------
# getting parameters
#-----------------------------------------------------------------
#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
  'help' , 'h', 0, "logical",
  'infile' , 'i', 1, "character",
  'groupfile' , 'G', 2, "character",
  'outfile','o',2,"character",
  'cell.width' , 'w', 2, "double",
  'cell.height','e',2,"double",
  'title','t',2,"character",
  'width','W',2,"integer",
  'height','H',2,"integer",
  'size','s',2,"double",
  'rowname','R',2,"logical",
  'colname','C',2,"logical",
  'color','c',2,"character" ,
  'zero','z',2,"double",
  'log','l',2,"character",
  'scale','S',2,"character"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

#遗传图与基因组的共线性分析
# define usage function
print_usage <- function(spec=NULL){
  cat(getopt(spec, usage=TRUE));
  cat("Usage example: \n")
  cat("
      Usage example: 
      1) Rscript heatmap.R --infile in.heatmap --outfile heatmap --color BrBG
      2) Rscript heatmap.R --infile in.heatmap --outfile heatmap --groupfile group.heatmap --title heatmap --size 10 --rownames F
      3) Rscript heatmap.R --infile in.heatmap --outfile heatmap --title heatmap --size 10 --cell.width 7 --cell.height 7



      
      Options: 
      --help          -h  NULL        get this help
      --infile        -i  character   the tab delimited input file saving numeric matrix of the values to be plotted.[forced]
      --outfile       -o  character   file path where to save the picture. Filetype is decided by the extension in the path. [optional,heatmap in current working directory]
      --groupfile     -G  character   the tab delimited input file saving data frame that specifies the annotations shown on top side of the heatmap [optional, default:NA]
      --cell.width    -w  double      individual cell width in points[optional, default: 7]
      --cell.height   -e  double      individual cell height in points[optional, default: 7]
      --size          -s  double      base fontsize for the plot[optional, default: 10]
      --width         -W  double      manual option for determining the output file width in pixel.[optional, default: NA]
      --heigth        -H  double      manual option for determining the output file height in pixel.[optional, default:NA]
      --title         -t  character   a title for the plot[optional, default: ]
      --rowname       -R  logical     boolean specifying if row names are be shown.[optional, default: TRUE]
      --colname       -C  logical     boolean specifying if column names are be shown.[optional, default:NA]
      --color         -c  character   choose the colour set(redgreen BrBG PiYG PRGn PuOr RdBu RdGy RdYlBu RdYlGn Spectral)or set colour splited by , .[optional, default: BrBG]
      --zero          -z  double      Set the minima vlaue: set mat values less than minima to minima.[optional, default:1]
      --log           -l  character   a logarithmic log scale is in use.[optional, default:log2]
      --scale         -S  character   character indicating if the values should be centered and scaled in either the row direction or the column direction, or none..[optional, default:none]
      \n")
  q(status=1);
}

#if(file.exists(paste(.sourcePath,"heatmap/pheatmap2.r",sep="")))source(paste(.sourcePath,"heatmap/pheatmap2.r",sep=""))else stop(paste(.sourcePath,"heatmap/pheatmap2.r does not exist!",sep=""))

# if help was asked for print a friendly message
# and exit with a non-zero error code
if (!is.null(opt$help)) { print_usage(spec) }


# check non-null args
if ( is.null(opt$infile) )  { print_usage(spec) }else {opt$infile<-gsub("\\\\",replacement = "/",x = opt$infile)}

#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$groupfile) )  { opt$groupfile=NA }else {opt$groupfile<-gsub("\\\\",replacement = "/",x = opt$groupfile)}

if( is.null(opt$outfile))opt$outfile="heatmap"
if(is.null(opt$title))opt$title=""

if(is.null(opt$width)){
  opt$width=NA  
}else if(!(is.numeric(opt$width)&&opt$width>0)){
  stop("Parameter Error：outfile width must be positive integer")  
}else{
  opt$width=opt$width/500
}

if(is.null(opt$height)){
  opt$height=NA
}else if(!(is.numeric(opt$height)&&opt$height>0)){
  stop("Parameter Error：outfile height must be positive integer")  
}else{
  opt$height=opt$height/500
}


if(is.null(opt$cell.width)){
  opt$cell.width=ifelse(is.na(opt$width),7,NA)
}else if(!(is.numeric(opt$cell.width)&&opt$cell.width>0)){
  stop("Parameter Error：cell width must be positive integer")
}

if(is.null(opt$cell.height)){
  opt$cell.height=ifelse(is.na(opt$height),7,NA )  
}else if(!(is.numeric(opt$cell.height)&&opt$cell.height>0)){
  stop("Parameter Error：cell height must be positive integer")
}

if(is.null(opt$rowname))opt$rowname=T
if(is.null(opt$colname))opt$colname=T
#if(is.null(opt$color))opt$color="RdYlGn"
if(is.null(opt$color))opt$color="#263C8B,#4E74A6,#BDBF78,#BFA524"
if(is.null(opt$zero))opt$zero=0
if(is.null(opt$size))opt$size=10
if(is.null(opt$log))opt$log=NA
if(is.null(opt$scale))opt$scale="none"

##import data
rawdat<-read.table(opt$infile,head=T,sep="\t",comment.char = "",check.names =F,quote="")
message(nrow(rawdat))
n<-nrow(rawdat)
n
if(nrow(rawdat)>30){
#rawdat<-read.table(as.vector(opt$infile),header=T,sep="\t",comment.char = "")
rawdat <- rawdat[1:30,]
}else if(nrow(rawdat) >0 && nrow(rawdat) < 30){
    rawdat <- rawdat
}
rownames(rawdat)<-as.matrix(rawdat)[,1]
#rownames(rawdat)
#rawdat=as.matrix(rawdat[,grepl("[0-9]+$",colnames(rawdat))])
#rawdat<-as.matrix(rawdat[,2:(ncol(rawdat)-3)])
rawdat<-as.matrix(rawdat[,-1])
rawdat<-rawdat+opt$zero
if(is.na(opt$log)){
  rawdat=rawdat
}else if(opt$log=="log2"){  
    rawdat<-(rawdat)
}else if(opt$log=="log10"){
  rawdat<-log10(rawdat)
}else{
  stop("Paramter error: a logarithmic scale parameter log can only be NA log10 or log2!")
}


#
if(is.na(opt$groupfile)){
  anColor = NA
  colGroup =NA
  heat.dat=rawdat
  
}else{
  groupdat<-read.table(as.vector(opt$groupfile),header=F,sep="\t",comment.char = "")
  group<-as.vector(groupdat[,2])
  names(group)<-as.vector(groupdat[,1])
  if(sum(!is.element(names(group),colnames(rawdat)))>0){
    stop(paste(c("the following samples in group file not exist:",setdiff(names(group),colnames(rawdat)),"please check your groupfile!"),sep="\n"))
  }
  if(sum(!is.element(colnames(rawdat),names(group)))>0){
    warning(paste(c("the following samples in infile will not be ploted:",setdiff(names(group),colnames(rawdat))),sep="\n"))
  }
  #多类样品热图添加分类条
  heat.dat<-rawdat[,names(group)]
  colGroup<-data.frame(Group=group)
  colGroup$Group= factor(colGroup$Group, levels = c(unique(group), "other"))
  row.names(colGroup)<-names(group)#设置样品颜色类
  gColor<-c( "#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666", "#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC")
  gColor=gColor[1:length(unique(group))]
  names(gColor)<-unique(group)
  anColor<-list(Group=gColor)
}


if(length(opt$color)==1&&is.element(opt$color,c("BrBG","PiYG","PRGn","PuOr","RdBu","RdGy","RdYlBu","RdYlGn","Spectral"))){  
  require(RColorBrewer)
  hmColors=colorRampPalette(rev(brewer.pal(n = 7, name = opt$color)))(100)
}else if(length(opt$color)==1&&(opt$color=="redgreen")){
  library(gplots) 
  message(paste("color=",opt$color,sep=""))
  hmColors=redgreen(255)
}else{
  hmColors<-strsplit(opt$color,split = ",")[[1]]
  hmColors=colorRampPalette(hmColors)(256)    
}
hl<-hclust(dist(heat.dat))
capture.output(str(as.dendrogram(hclust(dist(heat.dat)))),file =paste(c(opt$outfile,".txt"),collapse =""))
#message(c("width",opt$width,"height",opt$height))
pheatmap2(filename =paste(c(opt$outfile,".png"),collapse =""),width = opt$width,height = opt$height,mat=heat.dat,cellwidth=opt$cell.width,color = hmColors,cellheight=opt$cell.height,main=opt$title,cluster_rows=F,cluster_cols=F,annotation_col = colGroup,annotation = colGroup,annotation_colors = anColor,fontsize=opt$size,col=hmColors,show_rownames=opt$rowname,show_colnames=opt$colname,fontsize_col=ifelse(is.na(opt$cell.width),opt$size,min(opt$size,opt$cell.width)),fontsize_row=ifelse(is.na(opt$cell.height),opt$size,min(opt$size,opt$cell.height)),scale=opt$scale)
dev.off()
pheatmap2(filename =paste(c(opt$outfile,".pdf"),collapse =""),width = opt$width,height = opt$height,mat=heat.dat,cellwidth=opt$cell.width,color = hmColors,cellheight=opt$cell.height,main=opt$title,cluster_rows=F,cluster_cols=F,annotation_col = colGroup,annotation = colGroup,annotation_colors = anColor,fontsize=opt$size,col=hmColors,show_rownames=opt$rowname,show_colnames=opt$colname,fontsize_col=ifelse(is.na(opt$cell.width),opt$size,min(opt$size,opt$cell.width)),fontsize_row=ifelse(is.na(opt$cell.height),opt$size,min(opt$size,opt$cell.height)),scale=opt$scale)
dev.off()



