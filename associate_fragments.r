options(stringsAsFactors=FALSE)
require(rcdk)

# minimum number of compounds in which a fragment must appear in order to be considered
minimum_appearances = 2

compounds = read.table('dos_with_smiles.txt',sep='\t',header=FALSE,comment.char='',quote='',na.string='')
colnames(compounds) = c('broad_id', 'pubchem_id', 'smiles')

# remove ones with no smiles
compounds = compounds[!is.na(compounds$smiles),]

obj = parse.smiles(compounds$smiles)
mur = get.murcko.fragments(obj)
tbl_frags = table(unlist(mur))
use_frags = names(tbl_frags)[tbl_frags >= minimum_appearances]
n_frags = length(use_frags)

mat = matrix(nrow=dim(compounds)[1],ncol=n_frags)
rownames(mat) = compounds$smiles
colnames(mat) = use_frags

for (i in 1:dim(compounds)[1]) {
	mat[i,] = use_frags %in% unlist(mur[[i]])
}

# make up a list of screening hits for an example
hits = sample(1:dim(compounds)[1],size=300,replace=FALSE)

# convert to logical array
is_hit = 1:dim(compounds)[1] %in% hits

association = data.frame(fragment=all_frags, n_compounds=integer(n_frags), n_hits=integer(n_frags), or=numeric(n_frags), p_val=numeric(n_frags))

for (j in 1:n_frags) {
	n_compounds = sum(mat[,j])
	n_hits = sum(mat[,j] & is_hit)
	ctable = table(data.frame(hit=is_hit, frag=as.logical(mat[,j])))
	fisher_result = fisher.test(ctable,alternative='two.sided')
	association$n_compounds[j] = n_compounds
	association$n_hits[j] = n_hits
	association$or[j] = fisher_result$estimate
	association$p_val[j] = fisher_result$p.value
}


associate_fragments = function(mat, is_hit) {
  # initialize data frame of results to return
  association = data.frame(fragment=colnames(mat), 
                           n_compounds=integer(dim(mat)[2]), 
                           n_hits=integer(dim(mat)[2]), 
                           or=numeric(dim(mat)[2]), 
                           p_val=numeric(dim(mat)[2]))
  # test each fragment for association with hit status
  for (j in 1:dim(mat)[2]) {
    n_compounds = sum(mat[,j])
    n_hits = sum(mat[,j] & is_hit)
    ctable = ftable(hit=is_hit, frag=as.logical(mat[,j])) # * see comments below
    fisher_result = fisher.test(ctable,alternative='two.sided')
    association$n_compounds[j] = n_compounds
    association$n_hits[j] = n_hits
    association$or[j] = fisher_result$estimate
    association$p_val[j] = fisher_result$p.value
  }
  return (association)
}

# * based on benchmarking, this seems to be the most computationally intensive step, at about .01 seconds per iteration.
# we have 1977 columns and the whole function takes ~35s (elapsed) so this line is >50% of the compute time here.
# it appears to be just slightly faster to use ftable than table, and to feed 'hit' and 'frag' in directly as vectors
# rather than first wrapping them into a data.frame; but beyond that I am at a loss for how to make this faster.
# to benchmark the whole function: require(rbenchmark); benchmark(associate_fragments(mat,is_hit),replications=1)
# but unlike in Python, benchmark doesn't break it down by calls; you have to wrap each line in benchmark separately

bootstrap = function(mat, n_hits, n_iter, fdr=c(.10,.05,.01)) {
  p_val_distributions = {}
  for (iter in 1:n_iter) {
    hits = sample(1:dim(mat)[1],size=n_hits,replace=FALSE)
    is_hit = 1:dim(mat)[1] %in% hits
    associations = associate_fragments(mat, is_hit)
    p_val_distributions[[iter]] = associations$p_val
  }
  fdr_cutoffs = quantile(unlist(p_val_distributions), fdr)
  return_value = {}
  return_value[['fdr_cutoffs']] = fdr_cutoffs
  return_value[['p_val_distributions']] = p_val_distributions
  return (return_value)
}


