#!/usr/bin/env Rscript

####Load libraries and functions:
library("optparse")
library("deepSNV")
library("cowplot")
library("reshape2")
library("dplyr")
library("rlist")
library("tidyverse")
library("abind")
library("readr")
library("stringr")
# library("ggplot2")
source("functions/modified_deepsnv_function_bbb.R")

option_list = list(
    make_option(c("-m", "--mc_cores"), type="integer", default=6, 
              help="number of cores to use", metavar="character"),
    make_option(c("-g", "--genome_ver"), type="character", default="GRCh37",
              help="BS-genome to use", metavar="character"),
    make_option(c("-p", "--project_name"), type="character", default="project", 
              help="project_name", metavar="character"),
    make_option(c("-o", "--opt_q"), type="integer", default=25, 
              help="base quality for bam", metavar="character"),
    make_option(c("-q", "--opt_mq"), type="integer", default= 0, 
              help="map quality for bam", metavar="character"),
    make_option(c("-c", "--contig_id"), type="character", default="chr1", 
              help="contig id", metavar="character"),
    make_option(c("-d", "--contig_details"), type="character", default='length=249250621,file="Homo_sapiens.GRCh37.dna.primary_assembly.fa",species="Homo sapiens"', 
              help="contig details", metavar="character"),
    make_option(c("-v", "--minCoverage"), type="integer", default=200, 
              help="minimum coverage", metavar="character"),
    make_option(c("-b", "--postProbCutoff"), type="double", default=0.05, 
              help="minimum coverage", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

####Options:
cores          <- opt$mc_cores
if(opt$genome_ver == "GRCh37"){
  genome_ver <- "BSgenome.Hsapiens.UCSC.hg19"
} else if (opt$genome_ver == "GRCh38"){
  genome_ver <- "BSgenome.Hsapiens.UCSC.hg38"
}
project        <- opt$project_name
opt_q          <- opt$opt_q     # base quality filter for loading of bam files
opt_mq         <- opt$opt_mq      # map quality filter for loading of bam files 
contig_id      <- opt$contig_id # contig id
contig_details <- opt$contig_details
minCoverage    <- opt$minCoverage    # Minimum depth to use locus for mutation calls 
postProbCutoff <- opt$postProbCutoff   # Cutoff for the posterior prob to call mutation

# load appropriate library
library(genome_ver,character.only = TRUE)

targetFile2Contexts <- function(bedFile, genome=genome_ver,chrom, chr_prefix=NULL){
  require(GenomicRanges)
  if (is.character(genome)) { genome = get(genome) }
  stopifnot(file.exists(bedFile))
  stopifnot("BSgenome" %in% class(genome))
  # Load the target list:
  targets =  read_tsv(bedFile) %>% as("GRanges") 
  targets = targets[seqnames(targets) == chrom,]
  if(!grepl("chr", seqlevels(targets)[1])) {
    chr_prefix="chr"
  }
  result_data = NULL
  for (i in seq_along(targets)) {        # for each target interval 
    pos = tile(targets[i], width=1)[[1]] # data for each position
    seqlevels(pos) = paste0(chr_prefix, seqlevels(pos))
    cdata = 
      data.frame(
        chr=as.character(seqnames(targets[i])), 
        pos=start(pos), 
        target=as.character(targets[i]), 
        context=getSeq(genome, pos+1),   # gets context one either side of ref
        ref=getSeq(genome, pos),
        row.names=as.character(pos)
      )  
    result_data = rbind(result_data, cdata)
  }
  return(result_data)
}
bams_n <- list.files(path = 'normal/', pattern = "bam", full.names = TRUE)
bams_t <- list.files(path = 'tumor/', pattern = "bam", full.names = TRUE)

# chr <- str_split(bedfile,"_|\\.",simplify = T)[2]
contig <- opt$contig_id
targets <- list.files(path = '.', pattern = "bed") # bed file of targets
print(targets)
# Load target data:
target_data <- read_tsv(targets, col_names = TRUE, 
                        col_types = cols(.default = col_character())) %>% 
  dplyr::select(chr,start,end,annotation) %>% as("GRanges")
# subset by chromosome input
target_data <- target_data[seqnames(target_data) == contig]

# targetFile2Contexts produces this
#id      chr   pos   target    context                           ref
#chr-pos chr   pos   bed       trinucleotide around position     #ref-nucleotide

#eg
#              chr      pos              target context ref
#chr1:7202129    1  7202129   1:7202129-7202249     TGC   G
#chr1:7202130    1  7202130   1:7202129-7202249     GCT   C

####Determine context of all mutations:

# Expand the target object and determine mutation context for all variants:
contexts <- targetFile2Contexts(targets,chrom=contig)

####Extract count data for normal (i.e. BUFFY coat) samples:
# Determine all available reference bam files 

# This function uses the parallel package and the bam2R interface to read the
#nucleotide counts on each position of a .bam  alignment. The counts of both 
# strands are reported separately and nucleotides below a quality cutoff are
# masked. 

# It is called by deepSNV to parse the alignments of the test and
# control experiments, respectively.
countsN <- loadAllData(bams_n,
              target_data,
              q = opt_q,
              mq = opt_mq,
              mc.cores = cores,
              verbose = FALSE
              )

# A named matrix with rows corresponding to genomic positions and columns for the
# nucleotide counts (A, T, C, G, -), masked nucleotides (N) (INS)ertions, 
# (DEL)etions, (HEAD)s and (TAIL)s that count how often a read begins and ends at 
# the given position,  respectively,  and the sum of alignment  (QUAL)ities,  which 
# can  be  indicative  of  alignment  problems. 
# Counts  from  matches on the reference strand (s=0) are uppercase, counts on the
# complement (s=1) are lowercase.  The returned matrix has 11 * 2 (strands) = 22
# columns and (stop - start + 1) rows.

dimnames(countsN) <- list(basename(bams_n),
                          contexts$context,
                          c("A", "T", "C", "G", "-", "a", "t", "c", "g", "_"))



# output:
# array of matrices for each bam, each matrix rows = position (named ref 
# nucleotide, in context)
# columns  are "A", "T", "C", "G", "-", "a", "t", "c", "g", "_"



####Extract count data for other (tumour and matched normal) samples:

countsT <-
  loadAllData(bams_t,
              target_data,
              q = opt_q,
              mq = opt_mq,
              mc.cores = cores,
              verbose = FALSE)

dimnames(countsT) <- list(basename(bams_t),
                          contexts$context,
                          c("A", "T", "C", "G", "-", "a", "t", "c", "g", "_"))

## DeepSNV mutation calling
# Create composite reference for shearwater algorithm:
bgrNormal <- countsN
for (i in seq_len(dim(bgrNormal)[1])) {
  # Add the counts from the forward and reverse reads:
  comp <- bgrNormal[i, , 1:5] + bgrNormal[i, , 6:10]
  # Determine reference base
  refs <- substr(dimnames(bgrNormal)[[2]], 2, 2)
  # Get reference and total counts:
  refCount <-
    sapply(seq_len(nrow(comp)), function(j) {
      comp[j, refs[j]]
    })
  
  totCount <- apply(comp, 1, sum)

  # Blacklist positions with low coverage or high VAF:
  blackL <- (refCount / totCount) < 0.9 | totCount < minCoverage
  bgrNormal[i, blackL, ] <- 0
}

####Calling of mutations in tumour set with shearwater algorithm:

prior <- 0.01 
odds  <- prior / (1 - prior)

ppAllT <- lapply(seq_len(dim(countsT)[1]), function(i) {
  # binds the background counts
  comp <- abind(countsT[i, , , drop = FALSE], bgrNormal, along = 1)
  bf <- bbb_mod(
    comp,
    truncate = 0.1,
    model = "AND",
    mu.min = 1e-5,
    rho.min = 1e-5,
    min.depth = minCoverage
  )[1, ,]
  pp <- bf / (bf + odds)
  
  return(pp)
})

# returns a list of matrices with rownames = context triplicate and colnames the
# pp <- bf / (bf + odds) of either (A, T, C, G or -)
ppAllT <- abind(ppAllT, along = 3)
# creates an array of the matrices and names the matrix with the tumour sample 
# name
dimnames(ppAllT)[[3]] <- dimnames(countsT)[[1]]

# Mask sites where only one reference sample was covered
goodSites <-
  apply(apply(bgrNormal, c(1, 2), sum) > minCoverage, 2, sum) >= 1
ppAllT[!goodSites, ,] <- NA

# Determine reference base
bgrTumour <- countsT
refs <- substr(dimnames(bgrTumour)[[2]], 2, 2)
count_list <- list()
# Get SNV counts
for (i in seq_len(dim(bgrTumour)[1])) {
  # get the total counts
  total_count <- apply(bgrTumour[i, , ],1, sum)
  
  # get the count matrix for each sample
  count_matrix <- bgrTumour[i,,]
  # get the non ref count for each sample
  alt_count <- sapply(seq_len(nrow(count_matrix)), function(j) {
        sum(count_matrix[j,
            -which((colnames(count_matrix) == refs[j])
                   | (colnames(count_matrix) == tolower(refs[j])) )])
        })
  
  count_matrix <- cbind(count_matrix, total_count,alt_count)
  count_list <- list.append(count_list,count_matrix)
  #bgrTumour[i, , ] <- count_matrix
}

snv_array <- abind(count_list,along = 3)
dimnames(snv_array) <- list(contexts$context,
                            c("A", "T", "C", "G", "-", "a", "t", "c", "g", "_","total_count","alt_count"),
                            dimnames(bgrTumour)[[1]])
wh <- which(ppAllT < 0.01, arr.ind = TRUE)
snv_df <- data.frame(contexts$chr[wh[,1]],contexts$pos[wh[,1]],contexts$ref[wh[,1]],colnames(ppAllT)[wh[,2]],dimnames(ppAllT)[[3]][wh[,3]],
                     stringsAsFactors = FALSE )

colnames(snv_df) <- c("chr","pos","ref","alt","sample_id")
total_counts <- sapply( seq_len(nrow(wh)), function(i) {
    snv_array[wh[i,1],c(11),wh[i,3]]
})
alt_counts <- sapply( seq_len(nrow(wh)), function(i) {
    snv_array[wh[i,1],c(12),wh[i,3]]
})
snv_df <- add_column(snv_df, total_counts, alt_counts, .before = "sample_id")
snv_df <- snv_df %>% unite("alt/tot", alt_counts:total_counts, sep = "/") %>% pivot_wider(id_cols = c(chr,pos,ref,alt), values_from = "alt/tot", names_from = "sample_id")
snv_df <- snv_df[order(snv_df$pos),]
snv_df <- mutate(snv_df, ref = ifelse(ref == "-", ".", ref), alt = ifelse(alt  == "-", ".", alt ))

# Start writing to an output file
sink(paste0(project,"_",contig,".vcf"))
cat("##fileformat=VCFv4.2\n")
cat(paste0("##fileDate=",str_split(Sys.time()," ")[[1]][1] %>% str_remove_all("-"),"\n"))
cat(paste0("##source=",package.version("DeepSNV"),"\n"))
cat(paste0('##reference=',providerVersion(BSgenome.Hsapiens.UCSC.hg19),"\n"))
cat(paste0('##contig=<ID=',contig_id,",",contig_details,">\n"))
cat('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
cat('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
cat('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
cat(paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", paste0(colnames(snv_df)[-(1:4)],collapse = "\t"),"\n"))

get_row_info <- function(read_count,ref,alt){
  gt <- ''
  ad <- '.,.'
  dp <- '.'
  if(is.na(read_count)){gt <- './.'} else {
  read_count_split <- str_split(read_count, "/",simplify = TRUE)
  if(read_count_split[1] == read_count_split[2]){
    gt <- '1/1'
    } else {
    gt <- '0/1'
    }
  }
  if (is.na(read_count)) {ad <- '.,.'} else {
    read_count_split <- as.numeric(str_split(read_count, "/",simplify = TRUE))
    ref_count <- read_count_split[2] - read_count_split[1]
    alt_count <- read_count_split[1]
    ifelse(ref == "." | alt == ".", ad <- '.,.', ad <- paste0(ref_count,",", alt_count) )
    dp <- read_count_split[2]
  }
  return(paste0(gt,":",ad,":",dp))
 }

vcf_lines <- function(df){
paste0(df[1],"\t",as.integer(df[2]),"\t.\t",df[3],"\t",df[4],"\t",".","\t",".","\t",
                   paste0('NS=',sum(!is.na(df[-(1:4)]))),"\t",'GT:AD:DP',"\t",
                  str_c(sapply(df[-(1:4)], get_row_info, df[3],df[4]),collapse = "\t"),"\n")
}

cat(apply(snv_df,1,vcf_lines),sep = '')
sink()

