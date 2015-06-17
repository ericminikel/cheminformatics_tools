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

# will need to do some bootstrapping or something to decide how to handle the multiple testing burden