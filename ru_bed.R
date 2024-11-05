#!/usr/bin/env Rscript --vanilla

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))

option_list <- list(make_option(c("-c", "--controls"), action="store", help="Controls genes to output, separated by '-'. default COL1A, NDN and FMR.", default="COL1A1-FMR1-NDN"),
                    make_option(c("-b", "--buffer"), action="store", type="integer", default=100000, help="buffer to add to each side of each target. default is 100kb."),
                    make_option(c("-g", "--nonamesave"), action="store_true", help="Omit saving a copy of bed file with gene names, default false.", default=FALSE),
                    make_option(c("-e", "--ensembl"), action="store", help="Ensembl library file (csv). Including one increases performance speed. default included.", default="resources/ensemble.library.csv"),
                    make_option(c("-B", "--controlBuffer"), action="store", help="Specify control buffer length, if different from gene buffer. default is 50kb.", default=50000),
                    make_option(c("-n", "--naked"), action="store_true", help="Do not add any buffer or round targets, do not include controls", default=FALSE)
                    )

parser <- OptionParser(usage="%prog [options] genelist prefix", option_list=option_list, description="\nSupply a list of genes separated by '-' and a prefix for the name of your target files.\n
By default outputs a bedfile for use with adaptive sampling and a tab separated bed file with gene names for downstream use with samtools.\n
Output is also printed to stdout.\n
Example: ./makeRUTarget.R F8-F9-NDN testcase     -->     outputs testcase.targets.bed and testcase.named.targets.bed")

arguments <- parse_args(parser, positional_arguments=2)
opt <- arguments$options
genestring <- arguments$args[1]
outname <- arguments$args[2]

chrom.ok <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
"chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM")

geneList <- stringr::str_split(genestring, "-")[[1]]
if(length(geneList) < 1){
  stop(sprintf("ru_bed did not understand gene list input: %s", genestring))
}

controlGenes <- stringr::str_split(opt$controls, pattern="-")[[1]]
if(length(controlGenes) < 1){
  stop(sprintf("ru_bed did not understand control list input: %s", opt$controls))
}

buffer=opt$buffer
controlBuffer=opt$controlBuffer
ensemblpath=opt$ensembl
naked=opt$naked

namesaveskip=opt$nonamesave

fixOverlapRecur <- function(df, index=1){
  df <- df %>% dplyr::arrange(start_position)
  i = index
  while(i < dim(df)[1]){
    #test for overlap
    if(df[i+1, "start_position"] <= df[i, "end_position"]){
      startpos <- df[i, "start_position"]
      gene_name <- paste(df[i,"external_gene_name"], df[i+1, "external_gene_name"], sep="-")
      chrname <- df[i, "chromosome_name"]
      if(df[i, "end_position"] > df[i+1, "end_position"]){
        endpos <- df[i, "end_position"]
      } else {
        endpos <- df[i+1, "end_position"]
      }
      subdf <- df[-c(i, i+1),]
      smalldf <- data.frame(gene_name, chrname, startpos, endpos)
      names(smalldf) <- colnames(subdf)
      df <- rbind(subdf, smalldf)
      return(fixOverlapRecur(df, i))
    } else {
      i = i+1
      return(fixOverlapRecur(df, i))
    }
  }
  return(df)
}

makeTargetedBed <- function(geneList, ensembleLibraryFile=NA, buffer=100000, controlBuffer=50000, controlList=c("COL1A1", "FMR1"), returnGeneName=FALSE, round=TRUE){
  if(is.na(ensembleLibraryFile)){
    biolist <- as.data.frame(biomaRt::listMarts())
    ensemble <- biomaRt::useMart("ensembl")
    ensemble.list <- as.data.frame(biomaRt::listDatasets(ensemble))
    ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart=ensemble)
    attributes <- biomaRt::listAttributes(ensembl)
    ens.subset <- biomaRt::getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'ensembl_gene_id_version', 'chromosome_name', 'start_position', 'end_position'), mart=ensembl)
  } else {
    ens.subset <- read.csv(ensembleLibraryFile, header=TRUE, row.names=1)
  }
  
  ags.df <- ens.subset[which(ens.subset$external_gene_name %in% geneList),]
  control.df <- ens.subset[which(ens.subset$external_gene_name %in% controlList),]
  column.names <- c("chromosome_name", "start_position", "end_position")
  if(returnGeneName){
    column.names <- c("external_gene_name", column.names)
  }
  bed.df <- ags.df[,column.names] %>% dplyr::arrange(chromosome_name)
  cont.bed.df <- control.df[,column.names] %>% dplyr::arrange(chromosome_name)
  bed.df$chromosome_name <- paste("chr", bed.df$chromosome_name, sep="")
  cont.bed.df$chromosome_name <- paste("chr", cont.bed.df$chromosome_name, sep="")
  if(round==TRUE){
    bed.df$start_position <- sapply(bed.df$start_position, function(x) round(plyr::round_any(x-buffer, accuracy=50000, f=floor)))
    bed.df$end_position <- sapply(bed.df$end_position, function(x) round(plyr::round_any(x+buffer, accuracy=50000, f=ceiling)))
    cont.bed.df$start_position <- sapply(cont.bed.df$start_position, function(x) round(plyr::round_any(x-controlBuffer, accuracy=50000, f=floor)))
    cont.bed.df$end_position <- sapply(cont.bed.df$end_position, function(x) round(plyr::round_any(x+controlBuffer, accuracy=50000, f=ceiling)))
  }
  
  bed.df <- rbind(bed.df, cont.bed.df)
  
  bed.df$start_position <- format(bed.df$start_position, scientific=FALSE)
  bed.df$end_position <- format(bed.df$end_position, scientific=FALSE)
  bed.df <- bed.df %>% distinct(start_position, end_position, .keep_all=TRUE)
  
  #if(saveToFile){
  #  write.table(bed.df[,c("chromosome_name", "start_position", "end_position")], file=outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  #}
  
  chrs <- unique(bed.df$chromosome_name)
  chdf1 <- bed.df[which(bed.df$chromosome_name==chrs[1]),]
  if(naked==TRUE){
    mergedf <- chdf1 %>% dplyr::distinct(external_gene_name, .keep_all=TRUE)
  }else{
    mergedf <- fixOverlapRecur(chdf1)
  }
  for (chr in chrs[2:length(chrs)]){
    df <- bed.df[which(bed.df$chromosome_name==chr),]
    if(naked==TRUE){
      newdf <- df %>% dplyr::distinct(external_gene_name, .keep_all=TRUE)
    }else{
     newdf <- fixOverlapRecur(df)
    }
    mergedf <- rbind(mergedf, newdf)
  }
  
  mergedf <- mergedf %>% dplyr::distinct()
  return(mergedf)
}

# runtime
if(naked==TRUE){
  gene.bed.df <- makeTargetedBed(geneList[2:length(geneList)], ensembleLibraryFile=ensemblpath, buffer=0, controlBuffer=0, controlList=geneList[1], returnGeneName=TRUE, round=FALSE)
  gene.bed.df <- gene.bed.df %>% filter(chromosome_name %in% chrom.ok)
} else {
gene.bed.df <- makeTargetedBed(geneList, ensembleLibraryFile=ensemblpath, buffer=buffer, controlBuffer=controlBuffer, controlList=controlGenes, returnGeneName=TRUE)
gene.bed.df <- gene.bed.df %>% filter(chromosome_name %in% chrom.ok)
}
print(gene.bed.df)


#write.csv2(gene.bed.df, file=paste(outname, ".named.targets.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
#write.csv2(subset(gene.bed.df, select=-c(external_gene_name)), file=paste(outname, ".targets.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(gene.bed.df, file=paste(outname, ".named.targets.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")
write.table(subset(gene.bed.df, select=-c(external_gene_name)), file=paste(outname, ".targets.bed", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep=",")

# print out the total % of genome covered by these targets

gene.bed.df$span <- as.numeric(gene.bed.df$end_position) - as.numeric(gene.bed.df$start_position)
totalLength <- sum(gene.bed.df$span)

sprintf("Total target size is %2.2f MB", totalLength/1000000)
sprintf("Total percent of diploid human genome is %2.2f percent", totalLength/31000000)