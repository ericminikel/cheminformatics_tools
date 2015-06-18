options(stringsAsFactors=FALSE)
require(rcdk)

# Loads a data.frame with compound IDs and smiles, throws away any compounds without smiles
# This is hard-coded for a particular input data format at present, will need to be
# made more flexible in the future
load_compounds = function(path) {
  compounds = read.table(path,sep='\t',header=FALSE,comment.char='',quote='',na.string='')
  colnames(compounds) = c('broad_id', 'pubchem_id', 'smiles')
  # remove ones with no smiles
  compounds = compounds[!is.na(compounds$smiles),]
  return (compounds)
}

# Create a matrix of compounds (rows) vs. fragments in them (columns)
# smiles: vector of SMILES strings, cannot contain NA values
# minimum_appearances: minimum number of compounds in which a fragment must appear in order to be considered
# type: type of fragments to investigate. options are 'murcko' and 'exhaustive'
# note that sometimes Java heap will run out of memory on the parse.smiles step
# I don't yet have a solution for this other than re-starting R and trying again
create_fragment_matrix = function(smiles, minimum_appearances=2, type='murcko') {
  # check that inputs are valid
  stopifnot(sum(is.na(smiles))==0)
  stopifnot(type %in% c('murcko','exhaustive'))
  # parse smiles to get a list of rcdk objects representing molecules
  obj = parse.smiles(smiles)
  # get a list of fragments in these objects
  if (type == 'murcko') {
    frag = get.murcko.fragments(obj)
  } else if (type == 'exhaustive') {
    frag = get.exhaustive.fragments(obj)
  }
  tbl_frags = table(unlist(frag))
  use_frags = names(tbl_frags)[tbl_frags >= minimum_appearances]
  n_frags = length(use_frags)
  n_compounds = length(smiles)
  mat = matrix(nrow=n_compounds, ncol=n_frags)
  rownames(mat) = smiles
  colnames(mat) = use_frags
  for (i in 1:n_compounds) {
    mat[i,] = use_frags %in% unlist(frag[[i]])
  }
  return (mat)
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

empirical_fdr = function(mat, n_hits, n_iter, fdr=c(.10,.05,.01)) {
  p_val_distributions = {}
  for (iter in 1:n_iter) {
    hits = sample(1:dim(mat)[1],size=n_hits,replace=FALSE)
    is_hit = 1:dim(mat)[1] %in% hits
    associations = associate_fragments(mat, is_hit)
    p_val_distributions[[iter]] = associations$p_val
  }
  fdr_cutoffs_all = quantile(unlist(p_val_distributions), fdr)
  fdr_cutoffs_min = quantile(mapply(min,p_val_distributions), fdr)
  return_value = {}
  return_value[['fdr_cutoffs_all']] = fdr_cutoffs_all
  return_value[['fdr_cutoffs_min']] = fdr_cutoffs_min
  return_value[['p_val_distributions']] = p_val_distributions
  return (return_value)
}


# Create an image of a fragment and the top 5 compounds that are hits that contain this fragment
plot_lead = function(mat, is_hit, frag, max_hits=5, dest=NULL, width=1000, height=1000) {
  hits_in_class = head(rownames(mat)[is_hit & mat[,frag]], n = max_hits)
  if (!is.null(dest)) { # otherwise just plot to current graphics device
    pdf(dest, width=9, height=6)
  }
  par(mfrow=c(2,3), mar=c(2,2,2,2), oma=c(1,1,3,1))
  n_compounds = sum(mat[,frag])
  n_hits = sum(is_hit & mat[,frag])
  summary_text = paste('Fragment found in ',n_compounds,' compounds, ',n_hits,' of which are hits',sep='')
  temp = view.image.2d(parse.smiles(frag)[[1]],width,height) # get Java representation into an image matrix. set number of pixels you want horiz and vertical
  plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='') # create an empty plot
  rasterImage(temp,1,1,10,10) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries
  mtext(side=3, text='lead fragment', font=2)
  mtext(side=1, text=frag, cex=.5)
  for (hit in hits_in_class) {
    temp = view.image.2d(parse.smiles(hit)[[1]],width,height) # get Java representation into an image matrix. set number of pixels you want horiz and vertical
    plot(NA,NA,xlim=c(1,10),ylim=c(1,10),xaxt='n',yaxt='n',xlab='',ylab='') # create an empty plot
    rasterImage(temp,1,1,10,10) # boundaries of raster: xmin, ymin, xmax, ymax. here i set them equal to plot boundaries
    mtext(side=1, text=hit, cex=.5)
  }
  mtext(side=3, outer=T, text=summary_text)
  if (!is.null(dest)) { # otherwise just plot to current graphics device
    dev.off()
  }  
}

