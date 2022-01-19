
#' Plot a summary with boxes, each box is an alignment between a read and a construct using the output of dotplot.pl
#'
#' @param fileout a path to a pdf or a png
#' @param seq.lengths a named vector with the length of each read, the names should be the names of the reads
#' @param myConstructsLengths as seq.lengths but for reference.
#' @param inchPerBase the number of each per base to plot
#' @param annotation.for.reads a named vector with the subtitle for each read, the names should be the names of the reads
#' @param myConstructsOffsets a named vector with the offset for the construct, the names should be the names of the constructs (default is NULL)
#' @param myAnnotations a data.frame with at least 4 columns (start, end, name, construct) will be used to draw color rectangles (default is NULL)
#' @param myAnnotationsCols a named vector with colors for all the annotations, the names should match the name of `myAnnotations` (default is NULL)
#' @param colorBorder logical whether the border of the annotations should be the color of the annotation (without alpha) (default is FALSE)
#' @param myGuides a data.frame like myAnnotations but it will be used to draw lines (default is NULL)
#' @param filterHitsAbove integer maximum number of hit allowed for a sequence in the read on the construct (default is 1e15)
#' @param keepClosestMatch logical whether only one hit for a sequence in the read on the construct should be kept and only the one closest to 'unique hits' (this slows the process).
#' @param decreasingX a named vector of logical whether the x-axis should be plotted decreasing. The names should be the names of the reads (default is NULL = all FALSE)
#' 
#' @details assume there is in the same folder a file called shortName_construct.txt with shortName beeing the first part of the read name before the `-` which corresponds to the output of dot_plot.pl
#' 
plot_summary <- function(fileout, seq.lengths, myConstructsLengths, inchPerBase,
                         annotation.for.reads, 
                         myConstructsOffset = NULL, colorBorder = FALSE,
                         myAnnotations = NULL, myGuides = NULL, myAnnotationsCols = NULL,
                         filterHitsAbove = 1e15, keepClosestMatch = FALSE, decreasingX = NULL,
                         keepClosestToLongestMatch = FALSE, roundSize = 1000, onlyHighlight = F){
  if (keepClosestMatch & keepClosestToLongestMatch){
    stop("Incompatible keepClosest option.")
  }
  if (is.null(myConstructsOffset)){
    myConstructsOffset = rep(0, length(myConstructsLengths))
    names(myConstructsOffset) = names(myConstructsLengths)
  }
  if (is.null(decreasingX)){
    decreasingX <- rep(FALSE, length(seq.lengths))
    names(decreasingX) <- names(seq.lengths)
  } else if (! all(names(seq.lengths) %in% names(decreasingX))){
    warning("Some values of decreasingX are missing will be set to FALSE")
    decreasingX[setdiff(names(seq.lengths), names(decreasingX))] <- FALSE
  }
  if (grepl("pdf$", fileout)){
    pdf(fileout, title = basename(fileout), 
        width = 2 * 2 * length(seq.lengths) # margins
        + sum(seq.lengths) * inchPerBase, # plots
        height = 2 * 2 * (length(myConstructsLengths)) # margins
        + sum(myConstructsLengths) * inchPerBase # plots
        + 4) # I do not know why, probably for the read name
  } else {
    png(fileout,
        width = 2 * 2 * length(seq.lengths) # margins
        + sum(seq.lengths) * inchPerBase # plots
        + 1, # With png you need to be slightly above
        height = 2 * 2 * (length(myConstructsLengths)) # margins
        + sum(myConstructsLengths) * inchPerBase # plots
        + 4 # I do not know why, probably for the read name
        + 1, # With png you need to be slightly above
        units = "in",
        res = 96)
  }
  layout(matrix(1:((length(myConstructsLengths)) * length(seq.lengths)),
                ncol = length(seq.lengths)), 
         widths = lcm((4 + seq.lengths * inchPerBase) * 2.54),
         heights = lcm((4 + c(myConstructsLengths * inchPerBase) + c(rep(0, length(myConstructsLengths) - 1), 4)) * 2.54))
  par(mai = rep(2, 4),
      cex = 4,
      cex.axis = 1.2)
  for(read in names(seq.lengths)){
    shortName <- strsplit(read, "-")[[1]][1]
    for(construct in names(myConstructsLengths)){
      myTable <- tryCatch(read.delim(paste0(shortName, "_", construct, ".txt"),
                                     header = F),
                          error = function(e){data.frame(V1 = -1000000, V2 = -1000000)})
      if(construct == tail(names(myConstructsLengths), 1)){
        par(mai = rep(2, 4) + c(2, 0, 0, 0))
      }
      myTable$V2 <- myTable$V2 + myConstructsOffset[construct]
      nb.map <- table(myTable$V1)
      if (onlyHighlight){
        myOriginalTable <- myTable
      }
      myTable <- unique(subset(myTable, ! V1 %in% names(nb.map[nb.map > filterHitsAbove])))
      if (keepClosestMatch){
        print(paste0("Computing ", shortName, "_", construct, ".txt"))
        myTable.unique <- subset(myTable, V1 %in% names(nb.map[nb.map == 1]))
        print(paste0("Computing distance"))
        all.d <- as.matrix(dist(myTable))
        print(paste0("Done"))
        i.unique <- which(myTable$V1 %in% myTable.unique$V1)
        print(paste0("Objective: ", length(unique(myTable$V1))))
        while(nrow(myTable.unique) < length(unique(myTable$V1))){
          print(nrow(myTable.unique))
          i.to.check <- which(! myTable$V1 %in% myTable.unique$V1)
          # cat(paste("Get min of", length(i.to.check),"x", length(i.unique), "..."))
          min.value <- min(all.d[i.to.check, i.unique])
          # cat("Get match")
          i.to.add <- i.to.check[which(all.d[i.to.check, i.unique] == min.value, arr.ind=T)[, 1]]
          # cat("Done.\n")
          # print(i.to.add)
          # print(myTable[i.to.add, ])
          myTable.unique <- rbind(myTable.unique, unique(myTable[i.to.add, ]))
          i.unique <- c(i.unique, i.to.add)
        }
        myTable <- myTable.unique
        # minDist <- rep(0, nrow(myTable))
        # 
        # minDist[non.unique] <- apply(myTable[non.unique, ], 1, min_dist, df = myTable.unique)
        # myTable <- myTable[order(minDist), ]
        # myTable.filtered <- myTable[! duplicated(myTable$V1), ]
        # # I do a regression:
        # 
        # minDist <- rep(0, nrow(myTable))
        # non.unique <- which(! myTable$V1 %in% myTable.unique$V1)
        # minDist[non.unique] <- apply(myTable[non.unique, ], 1, min_dist, df = myTable.filtered)
        # myTable <- myTable[order(minDist), ]
        # myTable <- myTable[! duplicated(myTable$V1), ]
      }
      if (keepClosestToLongestMatch){
        myTable.unique <- subset(myTable, V1 %in% names(nb.map[nb.map == 1]))
        current.round <- roundSize
        myTable$diff <- myTable$V2 - myTable$V1
        myTable$sum <- myTable$V2 + myTable$V1
        # while(current.round >= 1 & nrow(myTable.unique) < length(unique(myTable$V1))){
        myWorkingTable <- myTable.unique
        myTable$roundDiff <- round(myTable$diff / current.round)
        myTable$roundSum <- round(myTable$sum / current.round)
        print(paste0("Objective: ", length(unique(myTable$V1))))
        while(length(unique(myWorkingTable$V1)) < length(unique(myTable$V1))){
          print(length(unique(myWorkingTable$V1)))
          # print(current.round)
          i.to.check <- which(! myTable$V1 %in% myWorkingTable$V1)
          t1 <- sort(table(myTable$roundDiff[i.to.check]))
          roundedFrequentValues <- names(t1)
          names(roundedFrequentValues) <- names(t1)
          for (v in names(t1)){
            roundedFrequentValues[as.character(as.integer(v) + -1:1)] <- v
          }
          myTable$roundDiffFreq <- roundedFrequentValues[as.character(myTable$roundDiff)]
          t2 <- sort(table(unique(myTable[i.to.check, c("V1", "roundDiffFreq")])$roundDiffFreq), decreasing = T)
          best.coeff.diff <- names(t2[1])
          nb.with.best.diff <- t2[1]
          t1 <- sort(table(myTable$roundSum[i.to.check]))
          roundedFrequentValues <- names(t1)
          names(roundedFrequentValues) <- names(t1)
          for (v in names(t1)){
            roundedFrequentValues[as.character(as.integer(v) + -1:1)] <- v
          }
          myTable$roundSumFreq <- roundedFrequentValues[as.character(myTable$roundSum)]
          t2 <- sort(table(unique(myTable[i.to.check, c("V1", "roundSumFreq")])$roundSumFreq), decreasing = T)
          best.coeff.sum <- names(t2[1])
          nb.with.best.sum <- t2[1]
          if (nb.with.best.diff > nb.with.best.sum){
            myWorkingTable <- rbind(myWorkingTable, myTable[i.to.check, ][myTable$roundDiffFreq[i.to.check] == best.coeff.diff, colnames(myWorkingTable)])
          } else {
            myWorkingTable <- rbind(myWorkingTable, myTable[i.to.check, ][myTable$roundSumFreq[i.to.check] == best.coeff.sum, colnames(myWorkingTable)])
          }
        }
        my.nb.map <- table(myWorkingTable$V1)
        myTable.unique <- subset(myWorkingTable, V1 %in% names(my.nb.map[my.nb.map == 1]))
        # current.round <- round(current.round / 2)
        # }
        i.to.check <- which(! myWorkingTable$V1 %in% myTable.unique$V1)
        print(paste0("Computing distance"))
        temp.df <- data.frame(dist = sapply(i.to.check, function(i){min(sqrt((myTable.unique$V1 - myWorkingTable$V1[i])^2 + (myTable.unique$V2 - myWorkingTable$V2[i])^2))}),
                              i = i.to.check,
                              V1 = myWorkingTable$V1[i.to.check])
        temp.df <- temp.df[order(temp.df$dist), ]
        temp.df.u <- temp.df[! duplicated(temp.df$V1), ]
        myTable <- rbind(myTable.unique, myWorkingTable[temp.df.u$i, colnames(myTable.unique)])
      }
      if (decreasingX[read]){
        xlim <- c(seq.lengths[read], 0)
      } else {
        xlim <- c(0, seq.lengths[read])
      }
      if (onlyHighlight && nrow(myOriginalTable) > nrow(myTable)){
        plot(myOriginalTable, pch = 19, xlab = shortName, ylab = construct,
             xlim = xlim, ylim = myConstructsOffset[construct] + c(0, myConstructsLengths[construct]),
             sub = annotation.for.reads[read], col = "grey")
        points(myTable, pch = 19)
      } else {
        plot(myTable, pch = 19, xlab = shortName, ylab = construct,
             xlim = xlim, ylim = myConstructsOffset[construct] + c(0, myConstructsLengths[construct]),
             sub = annotation.for.reads[read])
      }
      if(! is.null(myAnnotations)){
        annotations.df <- myAnnotations[myAnnotations$construct == construct, ]
        if(nrow(annotations.df) > 0){
          if(is.null(myAnnotationsCols)){
            myAnnotationsCols <- head(rainbow(length(unique(myAnnotations$name)) + 1), - 1)
            names(myAnnotationsCols) <- unique(myAnnotations$name)
          }
          if (colorBorder){
            border_color = myAnnotationsCols[annotations.df$name]
          } else {
            border_color = NULL
          }
          rect(xleft = 0, xright = seq.lengths[read],
               ybottom = annotations.df$start, ytop = annotations.df$end,
               col = alpha(myAnnotationsCols[annotations.df$name], 0.2),
               border = border_color)
        }
      }
      if(! is.null(myGuides)){
        guides.df <- myGuides[myGuides$construct == construct, ]
        if(nrow(guides.df) > 0){
          abline(h = (guides.df$start + guides.df$end) / 2)
        }
      }
      if(construct == tail(names(myConstructsLengths), 1)){
        par(mai = rep(2, 4))
      }
    }
  }
  dev.off()
}


#' Give mean, mean - sd, mean + sd to use it in ggplot stat_summary
#'
#' @param x a vector with numerical data
#' 
#' @details compute mean, mean - sd, mean + sd
#' 
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


#' Give minimum euclidean distance strictly positive between one point and a set of points
#'
#' @param x coordinates of a point
#' @param df a data.frame with coordinates of a set of points (each row is a point)
#' 
#' @details return the minimum distance strictly positive
#' 
min_dist <- function(x, df) {
  d <- apply(df, 1, function(v){sqrt(sum((x - v) ^ 2))})
  return(min(d[d > 0]))
}
