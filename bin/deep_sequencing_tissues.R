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

option_list = list(
    make_option(c("-m", "--mc_cores"), type="integer", default=6, 
              help="number of cores to use", metavar="character"),
    make_option(c("-g", "--genome_ver"), type="character", default="GRCh37",
              help="BS-genome to use", metavar="character"),
    make_option(c("-o", "--opt_q"), type="integer", default=25, 
              help="base quality for bam", metavar="character"),
    make_option(c("-q", "--opt_mq"), type="integer", default= 0, 
              help="map quality for bam", metavar="character"),
    make_option(c("-c", "--contig_id"), type="character", default="chr1", 
              help="contig id", metavar="character"),
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
opt_q          <- opt$opt_q     # base quality filter for loading of bam files
opt_mq         <- opt$opt_mq      # map quality filter for loading of bam files 
contig_id      <- opt$contig_id # contig id
minCoverage    <- opt$minCoverage    # Minimum depth to use locus for mutation calls 
postProbCutoff <- opt$postProbCutoff   # Cutoff for the posterior prob to call mutation

# load appropriate library
library(genome_ver,character.only = TRUE)

estimateRho <- deepSNV:::estimateRho

bbb_mod <-
  function (counts,
            rho = NULL,
            alternative = "greater",
            truncate = 0.1,
            rho.min = 1e-04,
            rho.max = 0.1,
            pseudo = .Machine$double.eps, # the smallest positive floating-point number x
            # such that 1 + x != 1. 
            # It equals double.base ^ ulp.digits if either 
            # double.base is 2 or double.rounding is 0;
            # otherwise, 
            # it is (double.base ^ double.ulp.digits)/ 2. Normally 2.220446e-16.
            return.value = c("BF", "P0", "err"),
            model = c("OR", "AND",
                      "adaptive"),
            min.cov = NULL,
            max.odds = 10,
            mu.min = 1e-06,
            mu.max = 1 - mu.min,
            min.depth = 1000)
  {
    pseudo.rho = .Machine$double.eps
    model = match.arg(model)
    return.value = match.arg(return.value)
    ncol = dim(counts)[3] / 2
    x.fw = counts[, , 1:ncol, drop = FALSE]
    x.bw = counts[, , 1:ncol + ncol, drop = FALSE]
    n.fw = rep(rowSums(x.fw, dims = 2), dim(x.fw)[3])
    n.bw = rep(rowSums(x.bw, dims = 2), dim(x.bw)[3])
    x <- x.fw + x.bw
    n = array(n.fw + n.bw, dim = dim(x)[1:2])
    mu = (x + pseudo.rho) / (rep(n + ncol * pseudo.rho, dim(x)[3]))
    ix = (mu < truncate)
    ix[1, ,] <- FALSE
    for (i in seq_len(ncol)) {
      ix[, , i] <- (n > min.depth) & ix[, , i]
    }
    
    if (is.null(rho)) {
      rho = estimateRho(x, mu, ix)
      rho = pmin(pmax(rho, rho.min), rho.max)
      rho[is.na(rho)] = rho.min
    }
    X = colSums(x, dims = 1)
    bound = function(x, xmin, xmax) {
      x = pmax(x, xmin)
      x = pmin(x, xmax)
      return(x)
    }
    disp = (1 - rho) / rho
    rdisp <- rep(disp, each = nrow(counts))
    mu = (x + pseudo) / (rep(n + ncol * pseudo, dim(x)[3]))
    mu = bound(mu, mu.min, mu.max) * rdisp
    tr.fw = x.fw * ix
    X.fw = rep(colSums(tr.fw, dims = 1), each = nrow(counts)) - tr.fw
    N.fw = rep(colSums(n.fw * ix), each = nrow(counts)) - n.fw *  ix
    nu0.fw <- (X.fw + x.fw + pseudo) / (N.fw + n.fw + ncol * pseudo)
    nu0.fw <- bound(nu0.fw, mu.min, mu.max) * rdisp
    mu0.bw <- (x.bw + pseudo) / (n.bw + ncol * pseudo)
    mu0.bw <- bound(mu0.bw, mu.min, mu.max) * rdisp
    nu.fw <- (X.fw + pseudo) / (N.fw + ncol * pseudo)
    nu.fw <- bound(nu.fw, mu.min, mu.max) * rdisp
    tr.bw = x.bw * ix
    X.bw = rep(colSums(tr.bw, dims = 1), each = nrow(counts)) -  tr.bw
    N.bw = rep(colSums(n.bw * ix), each = nrow(counts)) - n.bw *  ix
    nu0.bw <- (X.bw + x.bw + pseudo) / (N.bw + n.bw + ncol * pseudo)
    nu0.bw <- bound(nu0.bw, mu.min, mu.max) * rdisp
    mu0.fw <- (x.fw + pseudo) / (n.fw + ncol * pseudo)
    mu0.fw <- bound(mu0.fw, mu.min, mu.max) * rdisp
    nu.bw <- (X.bw + pseudo) / (N.bw + ncol * pseudo)
    nu.bw <- bound(nu.bw, mu.min, mu.max) * rdisp
    if (return.value == "err") {
      nu0 <- (X.bw + tr.fw + X.fw + tr.bw + pseudo) / (N.bw + n.bw + N.fw + n.fw + ncol * pseudo)
      nu0 <- bound(nu0, mu.min, mu.max)
      return(list(
        nu = nu0[1, , ],
        nu.fw = (nu0.fw / rdisp)[1, , ],
        nu.bw = (nu0.bw / rdisp)[1, , ],
        rho = rho
      ))
    }
    rm(tr.fw)
    rm(tr.bw)
    mu = pmax(mu, nu0.fw)
    mu = pmax(mu, nu0.bw)
    mu0.fw = pmax(mu0.fw, nu0.fw)
    mu0.bw = pmax(mu0.bw, nu0.bw)
    if (model %in% c("OR", "adaptive")) {
      Bf.fw <- deepSNV:::logbb(x.fw, n.fw, nu0.fw, rdisp) + deepSNV:::logbb(x.bw, n.bw, mu0.bw, rdisp) + deepSNV:::logbb(X.fw, N.fw, nu0.fw, rdisp) - deepSNV:::logbb(x.fw, n.fw, mu, rdisp) - deepSNV:::logbb(x.bw, n.bw, mu, rdisp) - deepSNV:::logbb(X.fw, N.fw, nu.fw, rdisp)
      Bf.fw = exp(Bf.fw)
      Bf.both = deepSNV:::logbb(x.fw, n.fw, nu0.fw, rdisp) + deepSNV:::logbb(X.fw,N.fw, nu0.fw, rdisp) - deepSNV:::logbb(x.fw, n.fw, mu, rdisp) - deepSNV:::logbb(X.fw, N.fw, nu.fw, rdisp)
      rm(X.fw, N.fw, mu0.bw, nu.fw)
      Bf.bw <- deepSNV:::logbb(x.fw, n.fw, mu0.fw, rdisp) + deepSNV:::logbb(x.bw, n.bw, nu0.bw, rdisp) + deepSNV:::logbb(X.bw, N.bw, nu0.bw,rdisp) - deepSNV:::logbb(x.fw, n.fw, mu, rdisp) - deepSNV:::logbb(x.bw, n.bw, mu, rdisp) - deepSNV:::logbb(X.bw, N.bw, nu.bw, rdisp)
      Bf.bw = exp(Bf.bw)
      Bf.both = Bf.both + deepSNV:::logbb(x.bw, n.bw, nu0.bw, rdisp) + deepSNV:::logbb(X.bw, N.bw, nu0.bw, rdisp) - deepSNV:::logbb(x.bw, n.bw, mu, rdisp) - deepSNV:::logbb(X.bw, N.bw, nu.bw, rdisp)
      Bf.both = exp(Bf.both)
      rm(X.bw, N.bw, mu0.fw, nu.bw)
      rm(mu, nu0.fw, nu0.bw)
      Bf = Bf.fw + Bf.bw - Bf.both + .Machine$double.xmin
    }
    else {
      Bf.both = deepSNV:::logbb(x.fw, n.fw, nu0.fw, rdisp) + deepSNV:::logbb(X.fw,  N.fw, nu0.fw, rdisp) - deepSNV:::logbb(x.fw, n.fw, mu, rdisp) - deepSNV:::logbb(X.fw, N.fw, nu.fw, rdisp) + deepSNV:::logbb(x.bw, n.bw, nu0.bw, rdisp) + deepSNV:::logbb(X.bw, N.bw, nu0.bw, rdisp) - deepSNV:::logbb(x.bw, n.bw, mu, rdisp) - deepSNV:::logbb(X.bw, N.bw, nu.bw, rdisp)
      Bf = exp(Bf.both)
    }
    if (model == "adaptive") {
      if (!is.null(min.cov))
        ix <- n.fw < min.cov | n.bw < min.cov
      else
        ix <- na.omit(abs(log10(n.fw / n.bw)) > log10(max.odds))
      Bf[ix] <- Bf.both[ix]
    }
    cons = apply(X, 1, which.max)
    for (i in 1:ncol(Bf))
      Bf[, i, cons[i]] = NA
    if (return.value == "P0") {
      return(Bf / (1 + Bf))
    }
    else {
      return(Bf)
    }
  }

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

write_vcf_header <- function(contig){
  sink(paste0(contig,".vcf"))
  cat("##fileformat=VCFv4.2\n")
  cat(paste0("##fileDate=",str_split(Sys.time()," ")[[1]][1] %>% str_remove_all("-"),"\n"))
  cat(paste0("##source=",package.version("DeepSNV"),"\n"))
  cat(paste0('##reference=',providerVersion(BSgenome.Hsapiens.UCSC.hg19),"\n"))
  cat(paste0('##contig=<ID=',contig,">\n"))
  cat('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n')
  cat('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
  cat('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n')
  cat('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
  cat(paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", paste0(basename(bams_t),collapse = "\t"),"\n"))
  sink()
}

# chr <- str_split(bedfile,"_|\\.",simplify = T)[2]
contig <- opt$contig_id
targets <- list.files(path = '.', pattern = "bed") # bed file of targets

# Load target data:
target_data <- read_tsv(targets, col_names = TRUE, 
                        col_types = cols(.default = col_character())) %>% 
  dplyr::select(chr,start,end,annotation) %>% as("GRanges")
# subset by chromosome input
target_data <- target_data[seqnames(target_data) == contig]

# targetFile2Contexts produces this
#              chr      pos              target context ref
#chr1:7202129    1  7202129   1:7202129-7202249     TGC   G
#chr1:7202130    1  7202130   1:7202129-7202249     GCT   C

####Determine context of all mutations:

# Expand the target object and determine mutation context for all variants:
contexts <- targetFile2Contexts(targets,chrom=contig)
if(is.null(contexts)){
  write_vcf_header(contig = contig)
  quit(status = 0)
}
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
write_vcf_header(contig)
sink(paste0(contig,".vcf"),append = TRUE)
cat(apply(snv_df,1,vcf_lines),sep = '')
sink()
