#!/usr/bin/env Rscript
# Script adapted from: https://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/topology-weighting/ and https://github.com/simonhmartin/twisst/blob/master/example_plot.R


# color to use
cols = c(
  "#0075DC", #Blue
  "#2BCE48", #Green
  "#FFA405", #Orpiment
  "#5EF1F2", #Sky
  "#FF5005", #Zinnia
  "#005C31", #Forest
  "#00998F", #Turquoise
  "#FF0010", #Red
  "#9DCC00", #Lime
  "#003380", #Navy
  "#F0A3FF", #Amethyst
  "#740AFF", #Violet
  "#426600", #Quagmire
  "#C20088", #Mallow
  "#94FFB5") #Jade

#Import ape and arguments (i.e. start and end of plot, plus level of smoothing)
library(ape)
args <- commandArgs(trailingOnly = TRUE)

# List of functions used to plot twisst results. Copied and pasted from https://github.com/simonhmartin/twisst/blob/master/plot_twisst.R for conveniency.
simple.loess.predict <- function(x, y, span, new_x=NULL, weights = NULL, max = NULL, min = NULL, family=NULL){
    y.loess <- loess(y ~ x, span = span, weights = weights, family=family)
    if (is.null(new_x)) {y.predict <- predict(y.loess,x)}
    else {y.predict <- predict(y.loess,new_x)}
    if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
    if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
    y.predict
    }

smooth.df <- function(x, df, span, new_x = NULL, col.names=NULL, weights=NULL, min=NULL, max=NULL, family=NULL){
    if (is.null(new_x)) {smoothed <- df}
    else smoothed = df[1:length(new_x),]
    if (is.null(col.names)){col.names=colnames(df)}
    for (col.name in col.names){
        print(paste("smoothing",col.name))
        smoothed[,col.name] <- simple.loess.predict(x,df[,col.name],span = span, new_x = new_x, max = max, min = min, weights = weights, family=family)
        }
    smoothed
    }

smooth.weights <- function(window_positions, weights_dataframe, span, new_positions=NULL, window_sites=NULL){
    weights_smooth <- smooth.df(x=window_positions,df=weights_dataframe,
                                span=span, new_x=new_positions, min=0, max=1, weights=window_sites)

    #return rescaled to sum to 1
    weights_smooth / apply(weights_smooth, 1, sum)
    }


stack <- function(mat){
    upper <- t(apply(mat, 1, cumsum))
    lower <- upper - mat
    list(upper=upper,lower=lower)
    }

interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }


sum_df_columns <- function(df, columns_list){
    new_df <- df[,0]
    for (x in 1:length(columns_list)){
        if (length(columns_list[[x]]) > 1) new_df[,x] <- apply(df[,columns_list[[x]]], 1, sum, na.rm=T)
        else new_df[,x] <- df[,columns_list[[x]]]
        if (is.null(names(columns_list)[x]) == FALSE) names(new_df)[x] <- names(columns_list)[x]
        }
    new_df
    }


plot.weights <- function(weights_dataframe,positions=NULL,line_cols=NULL,fill_cols=NULL,density=NULL,lwd=1,xlim=NULL,ylim=c(0,1),stacked=FALSE,
                                        ylab="Weighting", xlab = "Position", main="",xaxt=NULL,yaxt=NULL,bty="n", add=FALSE){
    #get x axis
    x = positions
    #if a two-column matrix is given - plot step-like weights with start and end of each window    
    if (dim(as.matrix(x))[2]==2) {
        x = interleave(positions[,1],positions[,2])
        yreps=2
        }
    else {
        if (is.null(x)==FALSE) x = positions
        else x = 1:nrow(weights_dataframe)
        yreps=1
        }
    
    #set x limits
    if(is.null(xlim)) xlim = c(min(x), max(x))
    
    #if not adding to an old plot, make a new plot
    if (add==FALSE) plot(0, pch = "", xlim = xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main,xaxt=xaxt,yaxt=yaxt,bty=bty)
    
    if (stacked == TRUE){
        y_stacked <- stack(weights_dataframe)
        for (n in 1:ncol(weights_dataframe)){
            y_upper = rep(y_stacked[["upper"]][,n],each=yreps)
            y_lower = rep(y_stacked[["lower"]][,n],each = yreps)
            polygon(c(x,rev(x)),c(y_upper, rev(y_lower)), col = fill_cols[n], density=density[n], border=NA)
            }
        }
    else{
        for (n in 1:ncol(weights_dataframe)){
            y = rep(weights_dataframe[,n],each=yreps)
            polygon(c(x,rev(x)),c(y, rep(0,length(y))), col=fill_cols[n], border=NA,density=density[n])
            lines(x,y, type = "l", col = line_cols[n],lwd=lwd)
            }
        }
    }




#read topologies file
topos = read.tree(file="Phylogeny.topos")

#weights file with a column for each topology
weights_file = 'Phylogeny.weights.tsv.gz'

weights = read.table(weights_file, header = T)

#normalise rows so weights sum to 1
weights = weights / apply(weights, 1, sum)
weights[is.na(weights)] <- 0


tot_weight = sort(apply(weights, 2, mean, na.rm=TRUE),decreasing=T)
Contribution = tot_weight/sum(tot_weight)*100

# Assign color to each topology
corresp <- as.data.frame(cbind(names(Contribution),cols[1:length(topos)]))

# with group names in 'V1' and color hex codes in 'V2'
group_names <- names(weights)  # Extract the group names from 'weights'

# Create a named vector of colors using the 'corresp' data frame
color_mapping <- setNames(corresp$V2, corresp$V1)

# Map the group names in 'weights' to their corresponding colors
group_colors <- color_mapping[group_names]

# Plot each topology contribution in a barplot
pdf("Contribution.pdf",10,10)
names(Contribution) = gsub("topo", "topo_", names(Contribution), ignore.case = FALSE, perl = TRUE,fixed = FALSE, useBytes = FALSE)
barplot(Contribution[1:length(topos)], col = corresp$V2,las=2)
dev.off()

# Extract ranking of the tree with a weight 
good_trees = as.numeric(as.character(gsub("topo", "", names(tot_weight), ignore.case = FALSE, perl = TRUE,fixed = FALSE, useBytes = FALSE)))

# Plot each topology ranked by frequency
pdf("Principal_topology.pdf",12,8)
par(mfrow = c(3,5), mar = c(1,1,2,1), xpd=NA)
a = 0
for (n in good_trees[1:length(topos)]){
  a = a + 1
  plot.phylo(topos[[n]], type = "clad", edge.color=corresp$V2[a], edge.width=5, label.offset=.4, cex=1)
  mtext(side=3,text=paste0("topo_",n," (",round(Contribution[a],digits=1)," %)"))
}
dev.off()

# Extract coordinate for each window and prepare it in the format expected by twisst
window_data_file = 'positions.txt'
window_data = read.table(window_data_file, header = T)
colnames(window_data) = c("start","end")
window_data$mid = (window_data$start + window_data$end) / 2
window_data$sites <- window_data$end - window_data$start + 1

# Plot stacked barplot without smoothing
pdf("twisst_barplot_no_smoothing.pdf",12,8)
par(mfrow =c(1,1), mar = c(3,3,1,1), xpd = FALSE)
plot.weights(weights_dataframe=weights, positions=cbind(window_data$start,window_data$end),
             line_cols=NA, lwd= 0, fill_cols= group_colors,stacked=TRUE,xlim =c(as.numeric(args[1]),as.numeric(args[2])))
dev.off()


# Plot stacked barplot with smoothing
weights_smooth = smooth.weights(window_positions=window_data$mid, weights_dataframe = weights,
                                span = as.numeric(args[3]), window_sites=window_data$sites)
pdf("twisst_barplot_with_smoothing.pdf",12,8)
plot.weights(weights_dataframe=weights_smooth, positions=window_data$mid,
             line_cols=cols, fill_cols=group_colors,stacked=TRUE,xlim = c(as.numeric(args[1]), as.numeric(args[2])))
dev.off()
