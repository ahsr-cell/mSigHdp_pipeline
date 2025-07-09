#!/usr/bin/env Rscript

### Script for mSigHdp 

application <- "mSigHdp mutational signature extraction pipeline."

message(paste("Welcome to", application, "For any questions, please consort the original paper (Liu et al. 2023, https://doi.org/10.1093/nargab/lqad005), GitHub/manual: (https://github.com/steverozen/mSigHdp/blob/v2.1.2-branch/mSigHdp_2.1.2.pdf), or the Stratton group. \n"))

### Load in required packages 
suppressPackageStartupMessages(require(mSigHdp))
suppressPackageStartupMessages(require(hdpx))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(argparse))

options(stringsAsFactors = F )

### Setting up CLI argument parsing 
# Create parser
parser = ArgumentParser(prog = 'mSigHdp', description='mSigHdp pipeline')
#Command line arguments
parser$add_argument("matrix_path", nargs = 1, help = "Specify path to input mutational matrix.") 

parser$add_argument("-s","--sample_matrix", type = 'character', help = "If available, specify path to sample key.") 

parser$add_argument("-a", "--analysis_type", type = "character", default = "Testing", help = "Specify type of analysis run. Options are [testing] or [analysis].", required=TRUE)

parser$add_argument("-b", "--burnin_iterations", type = 'double', default = "5000", help = "Specify number of burn-in iterations. Default set to 5000.", required=FALSE) 
parser$add_argument("-x", "--burnin_multiplier", type = 'double', default = "10", help = "Specify burin-in iteration multiplier. Default set to 10.", required=FALSE) 
parser$add_argument("-o", "--posterior", type = 'double', default = "250", help = "Specify number of posterior samples to collect. Default set to 250.", required=FALSE) 
parser$add_argument("-i", "--posterior_iterations", type = 'double', default = "100", help = "Specify number of iterations collected between samples. Default set to 100.", required=FALSE) 

parser$add_argument("-c", "--mutational_context", type = 'character', default = "SBS96", help = "Specify context of mutational matrix; options are SBS96 (default), SBS288, SBS1536, DBS78, or ID83.", required = TRUE)

#Parse arguments
args <- parser$parse_args()

mutation_matrix_path <- args$matrix_path
if(!exists("mutation_matrix_path")) {
  stop(sprintf("Mutation matrix not provided. Please specify by providing path at end of command; Use -h for further information."))
}

if (!is.null(args$sample_path)) {
  sample_path <- args$sample_path
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
}

if(u.analysis.type == 'Analysis' | u.analysis.type == 'analysis') {
  message(paste("You have selected Analysis - please note that this is intended to be run on Lustre as it requires 20 CPUs."))
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
rownames(mutation_types) <- NULL
mutation_types <- tibble::column_to_rownames(mutation_types, "MutationType")

##Sample key  
if (exists("sample_path")) {
  message(paste("Sample key provided. Incorporating into mutational matrix."))
  sample_key <- read.csv(file = sample_path)
  samples <- colnames(mutation_types)
  sample_key_vector <- as.vector(sample_key$sample_type)
  sample_key_vector_colnames <- paste0(sample_key_vector,"::",samples)
  
  colnames(mutation_types) = sample_key_vector_colnames
} else {
  message(paste("No sample key provided, please note that mSigHdp will run assuming a single sample type."))
}

### User specification of options
message(paste0("Successfully imported datasets and completed data wrangling. Proceeding with mSigHdp ", u.analysis.type, " run.")) 

message(paste0("Creating output file subdirectory called mSigHdp_run"))  
  main_dir <- getwd()
  sub_dir <- paste0("mSigHdp")
  if (!file.exists(sub_dir)){
    dir.create(file.path(main_dir, sub_dir))
    u.work.dir <- file.path(main_dir,sub_dir)
  } else {
    u.work.dir <- file.path(main_dir,sub_dir)
  }

if (exists("sample_path")) {
  if (u.analysis.type == 'testing' | u.analysis.type == 'Testing' | u.analysis.type == 'test' | u.analysis.type == 'Test') {
    message(paste0("Executing mSigHdp with test settings: 1000 burn-in iterations with 1x burn-in multiplier, collecting 5 posterior samples with 5 iterations between samples."))
    results <- mSigHdp::RunHdpxParallel(
      input.catalog = mutation_types,
      seedNumber = 123,
      K.guess = 5,
      out.dir = u.work.dir,
      burnin = 1000,
      burnin.multiplier = 1,
      post.n = 5,
      post.space = 5,
      num.child.process = 1,
      CPU.cores = 1,
      multi.types = TRUE, 
      overwrite = TRUE,
      high.confidence.prop = 0.9,  
      gamma.alpha = 1,   
      gamma.beta = 20,   
      checkpoint = TRUE,   
      verbose = TRUE
    )
    mSigHdp_Extracted_Signatures <- read.csv(file = paste0(u.work.dir,"/extracted.signatures.csv"), 
                                             header = TRUE)
    mSigHdp_Extracted_Signatures$MutationTypes <- rownames(mutation_types)
    
    mSigHdp_Extracted_Signatures <- mSigHdp_Extracted_Signatures %>% select(MutationTypes, everything())
    
    colnames(mSigHdp_Extracted_Signatures) = c('MutationType', paste0(u.mc,"_", LETTERS[1:ncol(mSigHdp_Extracted_Signatures) - 1]))
    
    write.table(mSigHdp_Extracted_Signatures, file = "mSigHdp_deNovoSignatures.txt", sep = "\t",
                  row.names = TRUE, col.names = TRUE)
    
    message(paste0("Analysis run for mSigHdp successfully executed. Output files (including mSigHdp_deNovoSignatures.txt) can be found in specified directory: ", u.work.dir," . Thank you for using mSigHdp."))
    
  }
 if (u.analysis.type == 'analysis' | u.analysis.type == 'Analysis') {
  message(paste0("Executing mSigHdp with ", u.burnin, " burn-in iterations, using a ", u.burnin.multip, "x multiplier. Collecting ", u.post, " posterior samples. Collecting ", u.post.space, " iterations between samples."))
  results <- mSigHdp::RunHdpxParallel(
    input.catalog        = mutation_types,
    out.dir              = u.work.dir, 
    num.child.process    = 20, 
    CPU.cores            = 20, 
    seedNumber           = 123,
    K.guess              = 16,
    burnin               = as.integer(u.burnin),
    burnin.multiplier    = as.integer(u.burnin.multip),
    post.n               = as.integer(u.post), 
    post.space           = as.integer(u.post.space), 
    multi.types          = TRUE, 
    overwrite            = TRUE,
    gamma.alpha          = 1,
    gamma.beta           = 20, 
    high.confidence.prop = 0.9,
    checkpoint           = TRUE,
    verbose              = FALSE)
  
  mSigHdp_Extracted_Signatures <- read.csv(file = paste0(u.work.dir,"/extracted.signatures.csv"), 
                                           header = TRUE)
  mSigHdp_Extracted_Signatures$MutationTypes <- rownames(mutation_types)
  
  mSigHdp_Extracted_Signatures <- mSigHdp_Extracted_Signatures %>% select(MutationTypes, everything())
  
  colnames(mSigHdp_Extracted_Signatures) = c('MutationType', paste0(u.mc,"_", LETTERS[1:ncol(mSigHdp_Extracted_Signatures) - 1]))
  
  write.table(mSigHdp_Extracted_Signatures, file = "mSigHdp_deNovoSignatures.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE)
  
  message(paste0("Analysis run for mSigHdp successfully executed. Output files (including mSigHdp_deNovoSignatures.txt) can be found in specified directory: ", u.work.dir," . Thank you for using mSigHdp."))
}

if (!exists("sample_path")) { #Single sample type, therefore, multi.types option turned to FALSE
  if (u.analysis.type == 'testing' | u.analysis.type == 'Testing' | u.analysis.type == 'test' | u.analysis.type == 'Test'){
    message(paste0("Executing mSigHdp with test settings: 1000 burn-in iterations with 1x burn-in multiplier, collecting 5 posterior samples with 5 iterations between samples."))
    results <- mSigHdp::RunHdpxParallel(
      input.catalog = mutation_types,
      seedNumber = 123,
      K.guess = 5,
      out.dir = u.work.dir,
      burnin = 1000,
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
    
    mSigHdp_Extracted_Signatures <- read.csv(file = paste0(u.work.dir,"/extracted.signatures.csv"), 
                                             header = TRUE)
    mSigHdp_Extracted_Signatures$MutationTypes <- rownames(mutation_types)
    
    mSigHdp_Extracted_Signatures <- mSigHdp_Extracted_Signatures %>% select(MutationTypes, everything())
    
    colnames(mSigHdp_Extracted_Signatures) = c('MutationType', paste0(u.mc,"_", LETTERS[1:ncol(mSigHdp_Extracted_Signatures) - 1]))
    
    write.table(mSigHdp_Extracted_Signatures, file = "mSigHdp_deNovoSignatures.txt", sep = "\t",
                  row.names = TRUE, col.names = TRUE)
    
    message(paste0("Analysis run for mSigHdp successfully executed. Output files (including mSigHdp_deNovoSignatures.txt) can be found in specified directory: ", u.work.dir," . Thank you for using mSigHdp."))
    
  } else if (u.analysis.type == 'analysis' | u.analysis.type == 'Analysis') {
    message(paste0("Executing mSigHdp with ", u.burnin, " burn-in iterations, using a ", u.burnin.multip, "x multiplier. Collecting ", u.post, " posterior samples. Collecting ", u.post.space, " iterations between samples."))
    results <- mSigHdp::RunHdpxParallel(
      input.catalog        = mutation_types,
      out.dir              = u.work.dir, 
      num.child.process    = 20, 
      CPU.cores            = 20, 
      seedNumber           = 123,
      K.guess              = 16,
      burnin               = as.integer(u.burnin),
      burnin.multiplier    = as.integer(u.burnin.multip),
      post.n               = as.integer(u.post), 
      post.space           = as.integer(u.post.space), 
      multi.types          = FALSE, 
      overwrite            = TRUE,
      gamma.alpha          = 1,
      gamma.beta           = 20, 
      high.confidence.prop = 0.9,
      checkpoint           = TRUE,
      verbose              = FALSE)
    
    mSigHdp_Extracted_Signatures <- read.csv(file = paste0(u.work.dir,"/extracted.signatures.csv"), 
                                             header = TRUE)
    mSigHdp_Extracted_Signatures$MutationTypes <- rownames(mutation_types)
    
    mSigHdp_Extracted_Signatures <- mSigHdp_Extracted_Signatures %>% select(MutationTypes, everything())
    
    colnames(mSigHdp_Extracted_Signatures) = c('MutationType', paste0(u.mc,"_", LETTERS[1:ncol(mSigHdp_Extracted_Signatures) - 1]))

    write.table(mSigHdp_Extracted_Signatures, file = "mSigHdp_deNovoSignatures.txt", sep = "\t",
                  row.names = TRUE, col.names = TRUE)
    
    message(paste0("Analysis run for mSigHdp successfully executed. Output files (including mSigHdp_deNovoSignatures.txt) can be found in specified directory: ", u.work.dir," . Thank you for using mSigHdp."))
  }
}
