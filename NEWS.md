#phylogram 2.1.0

* Addresses performance issues when converting between "dendrogram" and "phylo" objects
* More examples and detail added to vignette
* `as.cladogram` replaces `ultrametricize`
* `read.dendrogram` now wraps `ape::read.tree` and converts to dendrogram via `as.dendrogram.phylo`

#phylogram 2.0.1

Patch release addressing issue where as.phylo was masked from ape.

#phylogram 2.0.0

Major release retaining the newick parsing and tree manipulation 
functions but migrating the k-mer counting and clustering functions 
to the **kmer** package. 

Package no longer needs compilation as all Rcpp functions 
are migrated. 

All dependencies are eliminated with the exception of **ape**.


# phylogram 1.0.0

Released on CRAN 2017-06-12.
