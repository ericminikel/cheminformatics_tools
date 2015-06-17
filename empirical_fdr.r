#!/usr/bin/env Rscript

# This is just a wrapper for the empirical_fdr function (located in associate_fragments.r) to make it easy
# to submit jobs to the LSF cluster to do a large number of iterations to determine empirical FDR cutoffs
# Example command line usage:
# bsub -q week -o emp1000.o -e emp1000.e "Rscript empirical_fdr.r -c data/dos_with_smiles.txt -x 124 -i 1000"

options(stringsAsFactors=F)
source('associate_fragments.r')

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-c", "--cpds"), action="store", default=NA, type='character',
              help="Path to table of compounds"),
  make_option(c("-x", "--n_hits"), action="store", default=NA, type='integer',
              help="Number of hits to simulate"),
  make_option(c("-i", "--n_iter"), action="store", default=NA, type='integer',
              help="Number of iterations for empirical FDR")
)
opt = parse_args(OptionParser(option_list=option_list))

cpds = load_compounds(opt$cpds)
mat = create_fragment_matrix(cpds$smiles)
empirical_result = empirical_fdr(mat, n_hits=opt$n_hits, n_iter=opt$n_iter)

# print out cutoffs
empirical_result[['fdr_cutoffs_all']]
empirical_result[['fdr_cutoffs_min']]


