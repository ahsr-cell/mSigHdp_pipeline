#!/usr/bin/env Rscript

### Script for mSigHdp 

application <- "mSigHdp mutational signature extraction pipeline."

message(paste("Welcome to", application, "For any questions, please consort the original paper (Liu et al. 2023, https://doi.org/10.1093/nargab/lqad005), GitHub/manual: (https://github.com/steverozen/mSigHdp/blob/v2.1.2-branch/mSigHdp_2.1.2.pdf), or the Stratton group. \n"))

### Load in required packages 
suppressPackageStartupMessages(require(mSigHdp))
suppressPackageStartupMessages(require(hdpx))
suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = F )

### Setting up CLI argument parsing 
# Create parser
parser = ArgumentParser(prog = 'mSigHdp', description='mSigHdp pipeline')
#Command line arguments
parser$add_argument("matrix_path", nargs = 1, help = "Specify path to input mutational matrix.") 

parser$add_argument("-a", "--analysis_type", type = "character", default = "Testing", help = "Specify type of analysis run. Options are [testing] or [analysis].", required=TRUE)

parser$add_argument("-c", "--mutational_context", type = 'character', default = "SBS96", help = "Specify context of mutational matrix; options are SBS96 (default), SBS288, SBS1536, DBS78, or ID83.", required = TRUE)

#mSigHdp analysis run options
parser$add_argument("-b", "--burnin_iterations", type = 'double', default = "5000", help = "Specify number of burn-in iterations. Default set to 5000.", required=FALSE) 
parser$add_argument("-x", "--burnin_multiplier", type = 'double', default = "10", help = "Specify burin-in iteration multiplier. Default set to 10.", required=FALSE) 
parser$add_argument("-o", "--posterior", type = 'double', default = "250", help = "Specify number of posterior samples to collect off each posterior sampling chain. Default set to 250.", required=FALSE) 
parser$add_argument("-i", "--posterior_iterations", type = 'double', default = "100", help = "Specify number of iterations collected between each posterior sampling chain. Default set to 100.", required=FALSE) 
parser$add_argument("-ch", "--chains", type='double', default="20", help = "Specify number of chains to run. Note that this is will also be fed to the number of CPUs as the number of chains must always equal the number of CPUs. Default set to 20.", required=FALSE)
parser$add_argument("-k", "--clusters", type='double', default="16", help = "Specify number of clusters. Default set to 16.", required=FALSE)
parser$add_argument("-ga", "--alpha", type='double', default="1", help = "Specify number of clusters. Default set to 1.", required=FALSE)
parser$add_argument("-gb", "--beta", type='double', default="20", help = "Specify number of clusters. Default set to 20.", required=FALSE)
parser$add_argument("-h", "--confidence", type='double', default="0.9", help = "Specify number of clusters. Default set to 0.9.", required=FALSE)

#Parse arguments
args <- parser$parse_args()

mutation_matrix_path <- args$matrix_path
if(!exists("mutation_matrix_path")) {
  stop(sprintf("Mutation matrix not provided. Please specify by providing path at end of command; Use -h for further information."))
}

if (!is.null(args$mutational_context)) {
  mut_context <- args$mutational_context  
}

if(!exists("mut_context")) {
  stop(sprintf("Mutational signature context not specified. Please specify using -c or --mutational_context; Use -h for further information."))
}

u.analysis.type <- args$analysis_type

if (u.analysis.type == 'analysis') {
  if (!is.null("args$burnin_iterations")) {
    u.burnin <- args$burnin_iterations
    }
  if (!is.null("args$burnin_multiplier")) {
    u.burnin.multip <- args$burnin_multiplier
  }
  if (!is.null("args$posterior")) {
    u.post <- args$posterior
  }
  if (!is.null("args$posterior_iterations")) {
    u.post.space <- args$posterior_iterations
  }
  if (!is.null("args$clusters")) {
    u.clusters <- args$clusters
  }
  if (!is.null("args$alpha")) {
    u.alpha <- args$alpha
  }
  if (!is.null("args$beta")) {
    u.beta <- args$beta
  }
  if (!is.null("args$confidence")) {
    u.confidence <- args$confidence
  }
  if (!is.null("args$chains")) {
    u.chains <- args$chains
  }
}

if(u.analysis.type == 'Analysis' | u.analysis.type == 'analysis') {
  message(paste("You have selected Analysis - please note that this is intended to be run on Lustre as it requires 20 CPUs, unless specified otherwise via --chains."))
}

if (mut_context == 'SBS96' | mut_context == 'SBS288' | mut_context == 'SBS1536') {
  u.mc <- 'SBS'
}
if (mut_context == 'DBS78') {
  u.mc = 'DBS'
}
  if (mut_context == 'ID83') {
    u.mc = 'ID'
}

##### Setting up mSigHdp

message("Importing user datasets and conducting necessary data wrangling.")

### Data wrangling
##Mutation matrix 
#change context to rownames
mutation_types <- mutation_matrix_path
mutation_types <- read.csv(file = mutation_types, 
                           header = TRUE, sep="\t")

if (ncol(mutation_types) == 1 ) {
  message("Incorrect input mutation matrix format, attempting comma-delimited import. \n")
  mutation_types <- read.table(mutation_matrix, header=T, sep = ",")
}

if ("MutationType" %in% colnames(mutation_types)) {
  mutation_types <- tibble::column_to_rownames(mutation_types, "MutationType")
} else {
  stop(sprintf("Error: Input mutation matrix does not provide mutations under a column labelled as 'MutationType'. Please conduct the necessary data wrangling to ensure your mutation matrix is compatible with the pipeline. Stopping mSigHdp pipeline."))
}

#if (ncol(mutation_types) != nrow(mutation_types)) {
#  if (mut_context == 'SBS96') {
#    if (ncol(mutation_types) == 96 | nrow(mutation_types) != 96 ) {
#    message("Input mutation matrix detected with columns as mutation type. Conducting data wrangling to make compatible with mSigHdp pipeline.")
#    mutation_types <- as.data.frame(t(mutation_types))
#  }
#  }
#}

#if (tibble::has_rownames(mutation_types)==TRUE) {
#  trinuc_code <- rownames(mutation_types)
#}

### User specification of options
message(paste0("Successfully imported datasets and completed data wrangling. Proceeding with mSigHdp ", u.analysis.type, " run.")) 

message(paste0("Creating output subdirectory for run"))  
  main_dir <- getwd()
  sub_dir <- paste0("deNovo_signatures")
  if (!file.exists(sub_dir)){
    dir.create(file.path(main_dir, sub_dir))
    u.work.dir <- file.path(main_dir,sub_dir)
    u.work.dir
  } else {
    u.work.dir <- file.path(main_dir,sub_dir)
    message(paste0("Work directory is ",u.work.dir))
  }

if (u.analysis.type == 'testing' | u.analysis.type == 'Testing' | u.analysis.type == 'test' | u.analysis.type == 'Test') {
    message(paste0("Executing mSigHdp with test settings: 100 burn-in iterations with 1x burn-in multiplier, collecting 5 posterior samples off each posterior sampling chain with 5 iterations between each."))
    results <- mSigHdp::RunHdpxParallel(
      input.catalog = mutation_types,
      seedNumber = 123,
      K.guess = 5,
      out.dir = u.work.dir,
      burnin = 100,
      burnin.multiplier = 1,
      post.n = 5,
      post.space = 5,
      num.child.process = 1,
      CPU.cores = 1,
      multi.types = FALSE, 
      overwrite = TRUE,
      high.confidence.prop = 0.9,  
      gamma.alpha = 1,   
      gamma.beta = 20,   
      checkpoint = TRUE,   
      verbose = TRUE
    )
    
    message(paste0("mSigHdp analysis complete. Preparing final output files."))
  } else if (u.analysis.type == 'analysis' | u.analysis.type == 'Analysis') {
    message(paste0("Executing mSigHdp with ", u.burnin, " burn-in iterations, using a ", u.burnin.multip, "x multiplier. Collecting ", u.post, " posterior samples off each posterior sampling chain. Collecting ", u.post.space, " iterations between each chain."))
    results <- mSigHdp::RunHdpxParallel(
      input.catalog        = mutation_types,
      out.dir              = u.work.dir, 
      num.child.process    = as.integeter(u.chains), 
      CPU.cores            = as.integeter(u.chains), 
      seedNumber           = 123,
      K.guess              = as.integer(u.clusters), 
      burnin               = as.integer(u.burnin),
      burnin.multiplier    = as.integer(u.burnin.multip),
      post.n               = as.integer(u.post), 
      post.space           = as.integer(u.post.space), 
      multi.types          = FALSE,  
      overwrite            = TRUE, 
      gamma.alpha          = as.integer(u.alpha), 
      gamma.beta           = as.integer(u.beta), 
      high.confidence.prop = as.integer(u.confidence), 
      checkpoint           = TRUE,
      verbose              = FALSE 
      )
      message(paste0("mSigHdp analysis complete. Preparing final output files."))
  }

mSigHdp_Extracted_Signatures <- read.csv(file = paste0(u.work.dir,"/extracted.signatures.csv"), 
                                           header = TRUE)

if (tibble::has_rownames(mutation_types)==TRUE) {
  mSigHdp_Extracted_Signatures$MutationType <- rownames(mutation_types)
}

mSigHdp_Extracted_Signatures <- mSigHdp_Extracted_Signatures %>% select(MutationType, everything())
  
#colnames(mSigHdp_Extracted_Signatures) = c('MutationType', paste0(u.mc,"_", LETTERS[1:ncol(mSigHdp_Extracted_Signatures) - 1]))
  
write.table(mSigHdp_Extracted_Signatures, file = paste0(u.work.dir,"/mSigHdp_deNovoSignatures.txt"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = TRUE)

mSigHdp_SigPA_ExtractedSigs <- mSigHdp_Extracted_Signatures

colnames(mSigHdp_SigPA_ExtractedSigs) = c('MutationType', paste0(mut_context, LETTERS[1:ncol(mSigHdp_SigPA_ExtractedSigs) - 1]))

write.table(mSigHdp_SigPA_ExtractedSigs, file = paste0(u.work.dir,"/mSigHdp_deNovoSigs_sigPADecomp.txt"), 
                quote = FALSE, row.names = FALSE, sep = '\t')

lowconfsigs_path <- paste0(u.work.dir,"/low.confidence.signatures.csv")

if (file.exists(lowconfsigs_path)) {
  mSigHdp_LowConf_Signatures <- read.csv(file = paste0(u.work.dir,"/low.confidence.signatures.csv"), 
                                           header = TRUE)
                                        
  if (tibble::has_rownames(mutation_types)==TRUE) {
    mSigHdp_LowConf_Signatures$MutationType <- rownames(mutation_types)
  }
  
  mSigHdp_LowConf_Signatures <- mSigHdp_LowConf_Signatures %>% select(MutationType, everything())
  
  #colnames(mSigHdp_Extracted_Signatures) = c('MutationType', paste0(u.mc,"_", LETTERS[1:ncol(mSigHdp_Extracted_Signatures) - 1]))
  
  write.table(mSigHdp_LowConf_Signatures, file = paste0(u.work.dir,"/mSigHdp_lowConfSignatures.txt"), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = TRUE)

}

message(paste0("mSigHdp run successfully executed. Output files (including mSigHdp_deNovoSignatures.txt) can be found in specified directory: ", u.work.dir))