#!/usr/bin/env Rscript --vanilla

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))

option_list <- list(make_option(c("-n", "--network"), action="store_true", help="Output .bed and .tsv to network locations. default true.", default=FALSE),
                    make_option(c("-c", "--controls"), action="store", help="Controls genes to output, separated by '-'. default COL1A and FMR.", default="COL1A1-FMR1"),
                    make_option(c("-b", "--buffer"), action="store", type="integer", default=100000, help="buffer to add to each side of each target. default is 10kb."),
                    make_option(c("-g", "--namesave"), action="store_true", help="Save gene names in alignment target file, default true.", default=TRUE),
                    make_option(c("-e", "--ensembl"), action="store", help="Ensembl library file (csv). Including one increases performance speed. default included.", default="/Users/mirandagaley/Documents/Computational/R_scripts/ensemble.library.csv"),
                    make_option(c("--alignmentOut"), action="store", help="Directory to write target file to (gene names included)", default=""),
                    make_option(c("--gridIONOut"), action="store", help="Directory to write bed file to (gene names not included)", default=""),
                    make_option(c("--controlBuffer"), action="store", help="Specify control buffer length, if different from gene buffer. default is 5kb.", default=5000)
                    )

parser <- OptionParser(usage="%prog [options] genelist prefix", option_list=option_list, description="\nSupply a list of genes separated by '-' and a prefix for the name of your target files.\n
By default outputs a bedfile for use with adaptive sampling on the GridION and a tab separated region file with gene names for downstream use in samtools or minimap2.\n
Example: ./makeRUTarget.R F8-F9-NDN testcase     -->     outputs targets.testcase.tsv and testcase.targets.bed")

arguments <- parse_args(parser, positional_arguments=2)
opt <- arguments$options
genestring <- arguments$args[1]
outname <- arguments$args[2]


geneList <- stringr::str_split(genestring, "-")[[1]]
if(length(geneList) < 1){
  stop(sprintf("TargetMaker did not understand gene list input: %s", genestring))
}
if(opt$network){
  alignOut="/Volumes/eichler-vol27/projects/clinicalLongReadSeq/nobackups/targetsForGridION/"
  bedOut="/Volumes/eichler-vol27/projects/clinicalLongReadSeq/nobackups/alignments/targetFiles/"
}
if(opt$alignmentOut != ""){
  alignOut=opt$alignmentOut
}
if(opt$namesave & opt$gridIONOut != ""){
  bedOut=opt$gridIONOut
} else {
  stop(sprintf("no output specified"))
}
if(opt$gridIONOut=="" & opt$network != TRUE){
  bedOut=""
  saveToFile=FALSE
} else {
  saveToFile=TRUE
}

controlGenes <- stringr::str_split(opt$controls, pattern="-")[[1]]
if(length(controlGenes) < 1){
  stop(sprintf("Target Maker did not understand control list input: %s", opt$controls))
}


fixOverlapRecur <- function(df, index=1){
  df <- df %>% dplyr::arrange(start_position)
  i = index
  while(i < dim(df)[1]){
    #test for overlap
    if(df[i+1, "start_position"] <= df[i, "end_position"]){
      print(paste("adjusting", i))
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
      print(i)
      return(fixOverlapRecur(df, i))
    }
  }
  return(df)
}


makeTargetedBed <- function(geneList, outfile, ensembleLibraryFile=NA, buffer=100000, controlBuffer=20000, controlList=c("COL1A1", "FMR1"), saveToFile=TRUE, returnGeneName=FALSE, round=TRUE){
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
  
  if(saveToFile){
    write.table(bed.df[,c("chromosome_name", "start_position", "end_position")], file=outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
  }
  
  chrs <- unique(bed.df$chromosome_name)
  chdf1 <- bed.df[which(bed.df$chromosome_name==chrs[1]),]
  mergedf <- fixOverlapRecur(chdf1)
  
  for (chr in chrs[2:length(chrs)]){
    df <- bed.df[which(bed.df$chromosome_name==chr),]
    newdf <- fixOverlapRecur(df)
    mergedf <- rbind(mergedf, newdf)
  }
  
  mergedf <- mergedf %>% dplyr::distinct()
  return(mergedf)
}

# runtime
gene.bed.df <- makeTargetedBed(geneList, paste(alignOut, outname, ".targets.bed", sep=""), ensembleLibraryFile=opt$ensembl, buffer=opt$buffer, controlBuffer=opt$controlBuffer, controlList=controlGenes, saveToFile=saveToFile, returnGeneName=TRUE)

print(gene.bed.df)

write.table(gene.bed.df, file=paste(bedOut, "targets.", outname, ".tsv", sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

