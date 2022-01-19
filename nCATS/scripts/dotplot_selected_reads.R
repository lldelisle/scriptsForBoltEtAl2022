options(stringsAsFactors = F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("GenomicRanges", cranRep="https://cloud.r-project.org/")
safelyLoadAPackageInCRANorBioconductor("BiocGenerics", cranRep="https://cloud.r-project.org/")
safelyLoadAPackageInCRANorBioconductor("rtracklayer", cranRep="https://cloud.r-project.org/")
devtools::install_github("lldelisle/usefulLDfunctionsGR")
library(usefulLDfunctionsGR)
# To work with fasta:
safelyLoadAPackageInCRANorBioconductor("seqinr", cranRep="https://cloud.r-project.org/")
# To add alpha
safelyLoadAPackageInCRANorBioconductor("scales", cranRep="https://cloud.r-project.org/")

initial.options <- commandArgs(trailingOnly=FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
# Source the functions which are in the other file:
other.name <- file.path(script.basename, "nanopore_functions.R")
source(other.name)

sample <- commandArgs(trailingOnly=TRUE)[1]

# For the pdf plot
inchPerBase <- 0.0005
if (grepl("542", sample)){
  myConstructs <- c("542_around_Tg")
  good.construct <- "542_around_Tg"
  my.order <- c("075732f4")
} else {
  myConstructs <- c("542_around_Tg", "320_guess")
  good.construct <- "320_guess"
  my.order <- c("255ebdba", "62e2c507", "919cd74e", "b954ccaf", "ee36c6a1")
}

myConstructs_files <- sapply(myConstructs, function(cn){paste0(cn, ".fa")})
names(myConstructs_files) <- myConstructs
myConstructs_seq <- lapply(myConstructs_files, function(cf){read.fasta(cf, seqtype = "DNA")[[1]]})
myConstructsLengths <- sapply(myConstructs_seq, length)
myConstructsName <- sapply(myConstructs_seq, attr, which="name")
myConstructsOffset <- sapply(sapply(myConstructs_seq, attr, which="Annot"), function(x){
  m <- regexec(" range=([^:]+):([0-9, ]+)-([0-9, ]+) ", x)
  split <- regmatches(x, m)
  if (m[[1]][1] == -1){
    m <- regexec(">([^:]+):([0-9, ]+)-([0-9, ]+)$", x)
    split <- regmatches(x, m)
    if (m[[1]][1] == -1){
      return(0)
    } else {
      return(as.numeric(split[[1]][3]))
    }
  } else {
    return(as.numeric(split[[1]][3]))
  }
})

# Get the length of the interesting reads from the fasta:
interesting.reads.files <- list.files(pattern = ".fa$", path = "targetted_reads.faCONTIGS/")
names(interesting.reads.files) <- sapply(interesting.reads.files, function(f){strsplit(f, "-")[[1]][1]})
interesting.reads.files <- interesting.reads.files[my.order]
interesting.reads <- lapply(interesting.reads.files, function(f){read.fasta(file.path("targetted_reads.faCONTIGS", f), seqtype = "DNA")[[1]]})
seq.lengths <- sapply(interesting.reads, length)

# Run the dotplot
# dot_plot.pl was downloaded from jura.wi.mit.edu/page/papers/Hughes_et_al_2005/tables/dot_plot.pl
# Copyright (c) 2008 Helen Skaletsky and Whitehead Institute
for (shortName in names(interesting.reads)){
  f <- file.path("targetted_reads.faCONTIGS", interesting.reads.files[shortName])
  for (shortName2 in names(myConstructs_files)){
    f2 <- myConstructs_files[shortName2]
    basenameOutput <- paste0(shortName, "_", shortName2)
    if (! file.exists(paste0(basenameOutput, ".txt"))){
      system(paste0("perl ", script.basename, "/dot_plot.pl -w 20 -s 1 -1 ",
                    f, " -2 ", f2, " -o \"", basenameOutput, ".png\" -d \"", basenameOutput,
                    ".txt\" -t \"", shortName, " against ", shortName2, "\""))
    }
  }
}

# Load the annotations to put on the plots
myAnnotations <- readBed(file.path(script.basename, "annotations.bed"))
myAnnotations$construct <- myAnnotations$chr
myAnnotationsNames <- unique(myAnnotations$name)
# Give them colors
myAnnotationsCols <- c("#FF0000", "#00FF00", "#0000FF", "black")
names(myAnnotationsCols) <- c("LHA", "RHA", "LacZ", "pUC_back-bone")
# myAnnotationsCols <- sapply(strsplit(myAnnotations$itemRgb, ","), function(v){rgb(as.numeric(v[1]), as.numeric(v[2]), as.numeric(v[3]), maxColorValue=256)})
# names(myAnnotationsCols) <- myAnnotations$name

# Get the guides info
myGuides <- readBed(file.path(script.basename, "guides.bed"))
print(myGuides)
myGuides$construct <- myGuides$chr
myGuides$pretty_name <- myGuides$name


# Plot the legend of each contruct
for(construct in myConstructs){
  pdf(paste0(construct, "_legend.pdf"), title = paste0(construct, "_legend"),
      width = 2 * 2 # margins
      + myConstructsLengths[construct] * inchPerBase, # plot
      height = 2 * 2 # margins
      + 8)
  plot(1, type = "n", xlab = construct, xlim = myConstructsOffset[construct] + c(0, myConstructsLengths[construct]), ylim = c(0, 1.2), ylab = '', yaxt = 'n')
  annotations.df <- myAnnotations[myAnnotations$construct == construct, ]
  annotations.df <- annotations.df[order(annotations.df$start), ]
  if(nrow(annotations.df) > 0){
    rect(ybottom = 0, ytop = 1,
         xleft = annotations.df$start, xright = annotations.df$end,
         col = alpha(myAnnotationsCols[annotations.df$name], 0.2),
         border = myAnnotationsCols[annotations.df$name])
    annotations.used <- unique(annotations.df$name)
    legend("topleft", legend = annotations.used,
           fill = alpha(myAnnotationsCols[annotations.used], 0.2),
           border = myAnnotationsCols[annotations.used],
           bg = "white")
  }
  guides.df <- myGuides[myGuides$construct == construct, ]
  guides.df <- guides.df[order(guides.df$start), ]
  abline(v = (guides.df$start + guides.df$end) / 2)
  text(guides.df$start, y = seq(0, 1, length = nrow(guides.df)), labels = guides.df$pretty_name)
  dev.off()
  pdf(paste0(construct, "_legend_v.pdf"), title = paste0(construct, "_legend_v"),
      height = 2 * 2 # margins
      + myConstructsLengths[construct] * inchPerBase, # plot
      width = 2 * 2 # margins
      + 8)
  plot(1, type = "n", ylab = construct,
       ylim = myConstructsOffset[construct] + c(0, myConstructsLengths[construct]),
       xlim = c(0, 1), xlab = '', xaxt = 'n')
  annotations.df <- myAnnotations[myAnnotations$construct == construct, ]
  annotations.df <- annotations.df[order(annotations.df$start), ]
  if(nrow(annotations.df) > 0){
    rect(xleft = 0, xright = 1,
         ybottom = annotations.df$start, ytop = annotations.df$end,
         col = alpha(myAnnotationsCols[annotations.df$name], 0.2),
         border = myAnnotationsCols[annotations.df$name])
    annotations.used <- unique(annotations.df$name)
    legend("topleft", legend = annotations.used,
           fill = alpha(myAnnotationsCols[annotations.used], 0.2),
           border = myAnnotationsCols[annotations.used],
           bg = "white")
  }
  guides.df <- myGuides[myGuides$construct == construct, ]
  guides.df <- guides.df[order(guides.df$start), ]
  abline(h = (guides.df$start + guides.df$end) / 2)
  text(y = guides.df$start, x = seq(0, 1, length = nrow(guides.df)),
       labels = guides.df$pretty_name, adj = c(0.5, 1))
  dev.off()
}

if (length(myConstructs) > 1){
  plot_summary(fileout = paste0(sample, "_bothConstructs.pdf"),
               seq.lengths = seq.lengths, myConstructsLengths = myConstructsLengths,
               inchPerBase = inchPerBase, annotation.for.reads = NULL,
               myConstructsOffset = myConstructsOffset, colorBorder = TRUE,
               myAnnotations = myAnnotations, myGuides = myGuides,
               myAnnotationsCols = myAnnotationsCols, # decreasingX = decreasing.x,
               keepClosestMatch = FALSE)
}
plot_summary(fileout = paste0(sample, "_goodConstruct.pdf"),
             seq.lengths = seq.lengths, myConstructsLengths = myConstructsLengths[good.construct],
             inchPerBase = inchPerBase, annotation.for.reads = NULL,
             myConstructsOffset = myConstructsOffset[good.construct], colorBorder = TRUE,
             myAnnotations = myAnnotations, myGuides = myGuides,
             myAnnotationsCols = myAnnotationsCols, # decreasingX = decreasing.x,
             keepClosestToLongestMatch = TRUE, roundSize = 2500)

plot_summary(fileout = paste0(sample, "_goodConstruct_highlight.pdf"),
             seq.lengths = seq.lengths, myConstructsLengths = myConstructsLengths[good.construct],
             inchPerBase = inchPerBase, annotation.for.reads = NULL,
             myConstructsOffset = myConstructsOffset[good.construct], colorBorder = TRUE,
             myAnnotations = myAnnotations, myGuides = myGuides,
             myAnnotationsCols = myAnnotationsCols, # decreasingX = decreasing.x,
             keepClosestToLongestMatch = TRUE, roundSize = 2500, onlyHighlight = TRUE)

plot_summary(fileout = paste0(sample, "_goodConstruct_highlight.png"),
             seq.lengths = seq.lengths, myConstructsLengths = myConstructsLengths[good.construct],
             inchPerBase = inchPerBase, annotation.for.reads = NULL,
             myConstructsOffset = myConstructsOffset[good.construct], colorBorder = TRUE,
             myAnnotations = myAnnotations, myGuides = myGuides,
             myAnnotationsCols = myAnnotationsCols, # decreasingX = decreasing.x,
             keepClosestToLongestMatch = TRUE, roundSize = 2500, onlyHighlight = TRUE)
