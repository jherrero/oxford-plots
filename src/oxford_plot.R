#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

parser$add_argument("--genomes", default="genomes/",
        help="Specifies the directory containing the genome sequences (and the chr lengths)");
parser$add_argument("--alignments", default="alignments/callithrix_jacchus/",
        help="Specifies the directory containing the main alignments");
parser$add_argument("--colour", default="blue",
        help="Specifies the colour for drawing the main alignments");
parser$add_argument("--alignments2",
        help="Specifies the directory containing the main alignments");
parser$add_argument("--colour2", default="black",
        help="Specifies the colour for drawing the main alignments");
parser$add_argument("--pdf",
        help="Specifies the output filename");
parser$add_argument("--swap", action="store_true", default=FALSE,
        help="Swap the species in the plot");

args <- parser$parse_args()

genomes.dir <- args$genomes;
alignments.dir <- args$alignments;
user.col <- args$colour;
alignments.dir2 <- args$alignments2;
user.col2 <- args$colour2;
pdf.filename <- args$pdf;

cat("Preparing a new Oxford plot.\n");
cat("Files with chromosome lengths are in", genomes.dir,"\n");
cat("Alignment data are in", alignments.dir,"\n");


## FUNCTION read.chromosome.lengths

read.chromosome.lengths <- function(genomes.dir) {
    files <- list.files(genomes.dir, ".txt");
    chr.lengths <- list();
    for (this.file in files) {
        these.lengths <- read.table(paste(genomes.dir, this.file, sep="/"), row.names=1);
        chr.lengths[this.file] = list(these.lengths);
    }
    return(chr.lengths);
}


plot.chains <- function(filename, my.lengths.1, my.length.2, lw, col=red) {

    these.chain.lengths.1 <- c();
    these.chain.lengths.2 <- c();
    

    if (!is.na(file.info(filename)$size) & file.info(filename)$size > 0) {
        # cat(species.1, assembly.1, chr.1, species.2, assembly.2, chr.2, "\n");
        chains <- read.table(filename, header=T);

        if (dim(chains)[1] > 0) {
            chains[,1] <- chains[,1] + my.lengths.1[i];
            chains[,2] <- chains[,2] + my.lengths.2[j];
            if (args$swap) {
                chains[,c(1,2)] <- chains[,c(2,1)];
            }
            lines(chains, type="l", lw=lw, col=col, cex=0.3);
    
            for (c in 1:(dim(chains)[1]/3)) {
                # cat("c is ", c,"\n")
                chain.length.1 = abs(chains[(c-1)*3+2,1] - chains[(c-1)*3+1,1]);
                # cat(chains[(c-1)*3+1,1], chains[(c-1)*3+2,1], chain.length, "\n")
                these.chain.lengths.1 <- c(these.chain.lengths.1, chain.length.1);
            }
    
            for (c in 1:(dim(chains)[1]/3)) {
                # cat("c is ", c,"\n")
                chain.length.2 = abs(chains[(c-1)*3+2,2] - chains[(c-1)*3+1,2]);
                # cat(chains[(c-1)*3+1,2], chains[(c-1)*3+2,2], chain.length, "\n")
                these.chain.lengths.2 <- c(these.chain.lengths.2, chain.length.2);
            }
        }
    }
    return(list(these.chain.lengths.1, these.chain.lengths.2));
}





if (length(grep("^alignments.\\d+[MK]\\.\\d+[MK]/", alignments.dir, perl=T)) > 0) {
    min.length <- alignments.dir
    min.length <- sub("alignments.", "", alignments.dir)
    min.length <- sub("/.+", "", min.length)
    max.gap <- min.length
    min.length <- sub("\\..+", "", min.length)
    max.gap <- sub(".+\\.", "", max.gap)
}

chr.lengths <- read.chromosome.lengths(genomes.dir);

## Reads the list of files in this directory. Actually, this is only done to extract some information from the file names.
# cat("Reading files from ", alignments.dir);
files <- list.files(alignments.dir, ".rdotplot");
if (length(alignments.dir2) != 0) {
    files <- c(files, list.files(alignments.dir2, ".rdotplot"));
}

## Species and assembly names (the species names are used to label the axes)
species.1 <- sub(".dna_sm.chromosome.+", "", files[1]);
species.2 <- sub(".dna_sm.chromosome.+", "", sub(".+.fa.vs.", "", files[1]));
assembly.1 <- sub("\\w+\\.", "", species.1);
assembly.2 <- sub("\\w+\\.", "", species.2);
species.1 <- sub("\\..+", "", species.1);
species.2 <- sub("\\..+", "", species.2);

cat("Species 1: ", species.1, " ", assembly.1, "\n", sep="");
cat("Species 2: ", species.2, " ", assembly.2, "\n", sep="");
# print(chr.lengths);

## Extract a sorted list of unique chromosome names. This might be better extracted from the file with chr sizes in the future.
chromosomes.1 <- unique(sub("\\..+", "", sub(paste0(species.1, ".", assembly.1, ".dna_sm.chromosome."), "", files)));
chromosomes.1 <- grep("MT", chromosomes.1, invert=T, value=T)
##sorted.chromosomes.1 <- ("X");
sorted.chromosomes.1 <- c(sort(grep("^\\d\\D*$", chromosomes.1, value=T)), sort(grep("^\\d\\d", chromosomes.1, value=T)), sort(grep("^\\D*$", chromosomes.1, value=T)));

chromosomes.2 <- unique(sub(".fa(.chain)?(.\\d+[MK].\\d+[MK])?(.\\d+)?.rdotplot", "", sub(paste0(".+.fa.vs.", species.2, ".", assembly.2, ".dna_sm.chromosome."), "", files)));
chromosomes.2 <- grep("MT", chromosomes.2, invert=T, value=T)
##sorted.chromosomes.2 <- ("X");
sorted.chromosomes.2 <- c(sort(grep("^\\d\\D*$", chromosomes.2, value=T)), sort(grep("^\\d\\d", chromosomes.2, value=T)), sort(grep("^\\D*$", chromosomes.2, value=T)));

## Check we haven't missed any chromosome
if (length(chromosomes.1) != length(sorted.chromosomes.1)) {
    cat("Error while sorting chromsomes for", species.1, assembly.1);
    quit("no", 1);
}
if (length(chromosomes.2) != length(sorted.chromosomes.2)) {
    cat("Error while sorting chromosomes for", species.2, assembly.2, "\n");
    cat(chromosomes.2, "\n")
    cat(sorted.chromosomes.2, "\n")
    quit("no", 1);
}

cat("Chromosomes found for each species:\n");
cat("-", species.1, assembly.1, ":", sorted.chromosomes.1, "\n");
cat("-", species.2, assembly.2, ":", sorted.chromosomes.2, "\n");


my.lengths.1 <- c(0);
these.lengths <- data.frame(chr.lengths[tolower(paste0(species.1,".txt"))]);
if (length(these.lengths) == 0) {
    these.lengths <- data.frame(chr.lengths[tolower(paste0(species.1,"_",assembly.1, ".txt"))]);
}

for (chr.1 in sorted.chromosomes.1) {
    if (is.na(these.lengths[chr.1,])) {
        this.length <- these.lengths[paste0("chr",chr.1),];
    } else {
        this.length <- these.lengths[chr.1,];
    }
    my.lengths.1 <- c(my.lengths.1, tail(my.lengths.1, n=1) + this.length)
}
max.x <- tail(my.lengths.1, n=1);
#cat("max-x is", max.x, "\n");

my.lengths.2 <- c(0);
these.lengths <- data.frame(chr.lengths[tolower(paste0(species.2,".txt"))]);
if (length(these.lengths) == 0) {
    these.lengths <- data.frame(chr.lengths[tolower(paste0(species.2,"_",assembly.2, ".txt"))]);
}
for (chr.2 in sorted.chromosomes.2) {
    if (is.na(these.lengths[chr.2,])) {
        this.length <- these.lengths[paste0("chr",chr.2),];
    } else {
        this.length <- these.lengths[chr.2,];
    }
    my.lengths.2 <- c(my.lengths.2, tail(my.lengths.2, n=1) + this.length)
}
max.y <- tail(my.lengths.2, n=1);
#cat("max-y is", max.y, "\n");

if (length(pdf.filename) == 0) {
    if (exists("min.length") & exists("max.gap")) {
        pdf.filename <- paste("new", species.1, min.length, max.gap, "pdf", sep=".");
    } else {
        pdf.filename <- paste0(species.1, ".pdf");
    }
}
cat("Storing figure in", pdf.filename, "\n");
pdf(pdf.filename);

if (args$swap) {
    plot(NULL, xlim=c(1,max.y), ylim=c(1,max.x),
        xlab=paste0(sub("_", " ", species.2), " (", assembly.2, ")"),
        ylab=paste0(sub("_", " ", species.1), " (", assembly.1, ")"),
        xaxt="n", yaxt="n", mgp=c(0.5,0,0), bty="n");
    for (l in my.lengths.1) lines(c(0,max.y), c(l, l), col="grey");
    for (l in my.lengths.2) lines(c(l,l), c(0,max.x), col="grey");
    for (i in 1:length(sorted.chromosomes.2)) {
        text(sum(my.lengths.2[c(i, i+1)])/2, 1, sorted.chromosomes.2[i], pos=1, cex=0.6)
    }
    for (i in 1:length(sorted.chromosomes.1)) {
        text(1, sum(my.lengths.1[c(i, i+1)])/2, sorted.chromosomes.1[i], pos=2, cex=0.6)
    }
} else {
    plot(NULL, xlim=c(1,max.x), ylim=c(1,max.y),
        xlab=paste0(sub("_", " ", species.1), " (", assembly.1, ")"),
        ylab=paste0(sub("_", " ", species.2), " (", assembly.2, ")"),
        xaxt="n", yaxt="n", mgp=c(0.5,0,0), bty="n");
    for (l in my.lengths.2) lines(c(0,max.x), c(l, l), col="grey");
    for (l in my.lengths.1) lines(c(l,l), c(0,max.y), col="grey");
    for (i in 1:length(sorted.chromosomes.1)) {
        text(sum(my.lengths.1[c(i, i+1)])/2, 1, sorted.chromosomes.1[i], pos=1, cex=0.6)
    }
    for (i in 1:length(sorted.chromosomes.2)) {
        text(1, sum(my.lengths.2[c(i, i+1)])/2, sorted.chromosomes.2[i], pos=2, cex=0.6)
    }
}

all.chain.lengths.1 <- c();
all.chain.lengths.2 <- c();

# rect(0, data$Start, max.y, data$Stop, col="grey", border=NA)

for (i in 1:length(sorted.chromosomes.1)) {
    
    chr.1 <- sorted.chromosomes.1[i];

    for (j in 1:length(sorted.chromosomes.2)) {

        chr.2 <- sorted.chromosomes.2[j];

        if (exists("min.length") & exists("max.gap")) {
            filename <- paste(species.1, assembly.1, "dna_sm", "chromosome", chr.1, "fa", "vs", species.2, assembly.2, "dna_sm", "chromosome", chr.2, "fa", "chain", min.length, max.gap, "rdotplot", sep=".");
        } else {
            filename <- paste(species.1, assembly.1, "dna_sm", "chromosome", chr.1, "fa", "vs", species.2, assembly.2, "dna_sm", "chromosome", chr.2, "fa", "chain", "rdotplot", sep=".");
        }
        
        if (length(alignments.dir2)!=0) {
#            filename <- paste(species.1, assembly.1, "dna_sm", "chromosome", chr.1, "fa", "vs", species.2, assembly.2, "dna_sm", "chromosome", chr.2, "fa", "chain", "rdotplot", sep=".");
            plot.chains(paste0(alignments.dir2, "/", filename), my.lengths.1, my.lengths.2, 2, user.col2);
        }
#            filename <- paste(species.1, assembly.1, "dna_sm", "chromosome", chr.1, "fa", "vs", species.2, assembly.2, "dna_sm", "chromosome", chr.2, "fa", "chain", "2", "rdotplot", sep=".");
        these.chain.lengths <- plot.chains(paste0(alignments.dir, "/", filename), my.lengths.1, my.lengths.2, 5, user.col);
        all.chain.lengths.1 <- c(all.chain.lengths.1, these.chain.lengths[[1]]);
        all.chain.lengths.2 <- c(all.chain.lengths.2, these.chain.lengths[[2]]);

    }
}

sum.all.chain.lengths.1 = sum(all.chain.lengths.1, na.rm=T)
sum.all.chain.lengths.2 = sum(all.chain.lengths.2, na.rm=T)

cat("Total coverage:", sum.all.chain.lengths.1, sum.all.chain.lengths.2, "\n");

n50.1 <- 0;
for (x in rev(sort(all.chain.lengths.1))) {
    n50.1 <- n50.1 + x;
    if (n50.1 > sum.all.chain.lengths.1 / 2) {
        n50.1 <- x;
    }
}

n50.2 <- 0;
for (x in rev(sort(all.chain.lengths.2))) {
    n50.2 <- n50.2 + x;
    if (n50.2 > sum.all.chain.lengths.2 / 2) {
        n50.2 <- x;
    }
}
cat("N50:", n50.1, n50.2, "\n");


n50.1 <- 0;
for (x in sort(all.chain.lengths.1)) {
    n50.1 <- n50.1 + x;
    if (n50.1 > sum.all.chain.lengths.1 / 2) {
        n50.1 <- x;
    }
}

n50.2 <- 0;
for (x in sort(all.chain.lengths.2)) {
    n50.2 <- n50.2 + x;
    if (n50.2 > sum.all.chain.lengths.2 / 2) {
        n50.2 <- x;
    }
}
cat("N50:", n50.1, n50.2, "\n");

if (exists("min.length") & exists("max.gap")) {
    title(main = paste0("Collinear blocks (>", min.length, "bp long; gaps <", max.gap, "bp)"), sub = paste(length(all.chain.lengths.2), "blocks; length:", format(sum.all.chain.lengths.2/1000000000, digits=3), "Gbp; N50:", format(n50.2/1000000, digits=3), "Mbp"));
} else {
#    title(main = "Collinear blocks (>10Mbp long; gaps <10Mbp)", sub = paste(length(all.chain.lengths.2), "blocks; length:", format(sum.all.chain.lengths.2/1000000000, digits=3), "Gbp; N50:", format(n50.2/1000000, digits=3), "Mbp"));
}
cat("Done.\n");

