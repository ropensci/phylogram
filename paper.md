---
title: "phylogram: an R package for phylogenetic analysis with dendrograms"
author: "Shaun P. Wilkinson *^1,2^* and Simon K. Davy *^1^*"
output:
  html_document:
    keep_md: true
bibliography: vignettes/phylogram.bib
vignette: >
  %\VignetteIndexEntry{Introduction to the phylogram package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



--------------------------------------------------------------------------------

*^1^* School of Biological Sciences, Victoria University of Wellington, Kelburn Parade, Wellington 6012, New Zealand.

*^2^* Correspondence author. E-mail: shaunpwilkinson\@gmail.com

## Abstract
The **phylogram** R package is a tool for for developing 
phylogenetic trees as deeply-nested lists known as "dendrogram" objects. 
It provides functions for importing and exporting trees in 
parenthetic text format, as well as several tools for command-line 
tree manipulation and conversion between "dendrogram" and "phylo" class objects.
This facilitates access to a comprehensive range of object-specific 
analytical and tree-visualization functions found across a wide array of 
bioinformatic R packages.
The **phylogram** package is released under the GPL-3 license, 
and is available for download from 
CRAN <https://CRAN.R-project.org/package=phylogram> and GitHub 
<https://github.com/shaunpwilkinson/phylogram>.  


## Introduction
The R environment continues to gain popularity as a platform for
bioinformatic analysis, due to the reproducible 
code-based workflow and the many powerful analytical tools 
available in a suite of open-source packages such as **ape** 
[@Paradis2004], **phangorn** [@Schliep2011] and **Phytools** [@Revell2012]. 
These packages typically employ a tree structure known as the 
"phylo" object, whose primary element is an integer matrix 
with one row for each edge in the graph, 
and two columns giving the indices of the connecting nodes. 
This is a highly versatile and memory-efficient format suitable 
for most applications encountered by evolutionary biologists.

An alternative tree structure is the "dendrogram" object, generated 
using the `as.dendrogram` function in the **stats** package [@RCoreTeam2015].
Rather than a matrix of edges, a dendrogram is a hierarchical list. 
These 'lists of lists' can be deeply nested, with the limit depending on 
the C stack size (settable *via* `options("expressions")`).
A useful feature of this representation is its modularity, whereby the 
sub tree of a tree is itself a tree - a dendrogram within a dendrogram. 
This means that dendrogram objects are subsettable in the same 
way that standard lists are, which in addition to the inbuilt
editing functions such as `cut` and `merge`, 
facilitates intuitive command-line tree manipulation. 
An especially powerful feature of this object type is that tree 
editing operations can be carried out recursively 
using fast inbuilt functions in the "apply" family such as `dendrapply` 
and `lapply`. 

Each node of a dendrogram object has the following mandatory attributes: 

* "height" the position of the node along the vertical axis 
   (assuming the graph is orientated vertically)
* "midpoint" the horizontal distance of the node from the left-most member 
   of the sub tree (where the horizontal distance between adjacent leaves is 1 unit)
* "members" the number of terminal leaf nodes belonging to the node 
* "class" all nodes have the class attribute "dendrogram" 

Rather than lists, terminal leaf nodes are length-1 integer vectors whose
values correspond to the indices of the members in the set. 
The "members" attributes of leaf nodes is always 1, 
the "midpoint" attribute is 0, and they have two additional attributes:    

* "leaf" TRUE for terminal nodes (NULL otherwise)
* "label" an optional character string giving the name of the taxon or group

Aside from those listed above, users may attach other objects 
as attributes to the dendrogram nodes. For example, "label" attributes
can be attached to inner nodes, and users can specify plotting parameters
for each node by setting the attributes "nodePar" and "edgePar".

While generally not as memory efficient as "phylo" objects, 
the flexibility, modularity and intuitive structure of dendrogram 
objects are appealing to many users, particularly where highly dynamic 
tree structures are required for applications such as 
machine learning clustering and classification.


## The 'phylogram' package
Here, we introduce **phylogram**, an R package for working with 
evolutionary trees as deeply-nested lists. 
The package contains functions for importing and exporting dendrogram 
objects to and from parenthetic text, 
and assembling, visualizing, manipulating and rendering trees for publication.
These functions are detailed below with examples of their utility.

#### Example 1: Building a dendrogram object manually
Consider the simple example of a tree with three members named 
"A", "B" and "C", where "B" and "C" are more closely related
to each other than they are to "A". 
An unweighted Newick string for this tree would be *(A,(B,C));*
We can manually create a dendrogram object for this basic phylogeny
and plot the tree as follows:


```r
x <- list(1, list(2, 3))
## attach "leaf" and "label" attributes to leaf nodes
attr(x[[1]], "leaf") <- TRUE
attr(x[[2]][[1]], "leaf") <- attr(x[[2]][[2]], "leaf") <- TRUE
attr(x[[1]], "label") <- "A"
attr(x[[2]][[1]], "label") <- "B"
attr(x[[2]][[2]], "label") <- "C"
## set "height" attributes for all nodes
attr(x, "height") <- 1
attr(x[[1]], "height") <- 0
attr(x[[2]], "height") <- 0.5
attr(x[[2]][[1]], "height") <- attr(x[[2]][[2]], "height") <- 0
## set "midpoints" attributes for all nodes
attr(x, "midpoint") <- 0.75
attr(x[[1]], "midpoint") <- 0
attr(x[[2]], "midpoint") <- 0.5
attr(x[[2]][[1]], "midpoint") <- attr(x[[2]][[2]], "midpoint") <- 0
## set "members" attributes for all nodes
attr(x, "members") <- 3
attr(x[[1]], "members") <- 1
attr(x[[2]], "members") <- 2
attr(x[[2]][[1]], "members") <- attr(x[[2]][[2]], "members") <- 1
## set class as "dendrogram" 
## Note that setting the class for the root node
## automatically sets the class of all nested subnodes
class(x) <- "dendrogram"

plot(x, yaxt = "n")
```

<img src="vignettes/figures/unnamed-chunk-3-1.png" width="500px" style="display: block; margin: auto auto auto 0;" />

\ \ \ \ **Figure 1:** A simple dendrogram with three terminal leaf nodes


As demonstrated in this example, manually setting attributes on dendrogram
objects can be rather tedious, motivating the development of functions 
to automate the generation and manipulation of these tree 
structures.

### Importing and exporting trees
The Newick (a.k.a. New Hampshire) parenthetic text 
format [@Felsenstein1986] is a universal phylogenetic tree 
representation that is compatible with most tree-editing software. 
The **phylogram** package features the text parser `read.dendrogram`
that reads a character string or text file in the Newick
format and creates a dendrogram object.
This function supports weighted edges, labels with special meta-characters 
(enclosed in single quotation marks), comments 
(enclosed in square brackets; ignored by the parser), 
multifuricating nodes, and both rooted and unrooted trees.
Inner-node labels are currently ignored; however, 
the inclusion of "label" attributes for non-leaf 
nodes will be available in a future version.
Objects of class "dendrogram" can be exported as 
Newick-style parenthetic text using the function 
`write.dendrogram`.


#### Example 2: Import and export a tree from a Newick string
The simple Newick string in Example 1 can be imported as a 
dendrogram object using the `read.dendrogram` function 
as follows:


```r
library(phylogram)
newick <- "(A,(B,C));"
x <- read.dendrogram(text = newick)
x
#> 'dendrogram' with 2 branches and 3 members total, at height 2
```

The following command writes the object back to the console in 
Newick format without edge weights:


```r
write.dendrogram(x, edges = FALSE)
#> [1] "(A,(B,C));"
```

The syntax is similar when reading and writing text files, 
except that the `text` argument is replaced by `file`, 
and a valid file path is passed to the function. 

If required, the dendrogram can be converted to an object of class "phylo"
using the `as.phylo.dendrogram` method, and converted back to a dendrogram
with `as.dendrogram.phylo`.


```r
y <- as.phylo(x)
z <- as.dendrogram(y)
identical(x, z)
#> [1] TRUE
```


### Tree editing/manipulation
The **phylogram** package features 
several additional functions to facilitate some of the more common 
manipulation operations.
Leaf nodes and internal branching nodes can be removed 
using the function `prune`, which identifies and 
recursively deletes nodes based on regular expression pattern 
matching of node "label" attributes.
To aid visualization, the function `ladder` rearranges
the tree, sorting nodes by the number of members (analogous to the
`ladderize ` function in the **ape** package). Another function 
aiding in tree visualization is `ultrametricize `, which
resets the "height" attributes of all terminal leaf nodes to zero 
(note that unlike `ape::chronos()` there is no mathematical basis 
for this operation; rather it is merely a visualization aid). 
The function `reposition` scales the height of all nodes in
a tree by a given constant (passed *via* the argument `shift`), 
and features the option to reset all node heights so that height of 
the farthest terminal leaf node from the root is zero (by specifying 
`shift = "reset"`). 
The function `remidpoint` recursively corrects all "midpoint", 
"members" and "leaf" attributes following manual editing of a tree.


### Tree vizualization 
Publication-quality trees can be generated from dendrogram objects 
using the **stats** plotting function `plot.dendrogram`, and the extensive
plotting functions available in dendrogram-enhancing packages such as 
**circlize** [@Gu2014] and **dendextend** [@Galili2015].
The latter also offers the facility to convert dendrograms to "ggdend" objects, 
for which many powerful 'grammar of graphics' plotting functions are available in 
the **ggplot2** [@Wickam2009] and **ggdendro** [@deVries2016] packages. 
Moreover, there are several advanced plotting options for "phylo" objects in 
the **ape** package [@Paradis2004], as well as the
Bioconductor package **ggtree** [@Guangchuang2017].
Given the extensive tree visualization options already available, 
we do not include any additional plotting functions in the **phylogram** package. 


### Summary
The **phylogram** package offers a high quality tree parser, 
several new dendrogram-editing functions,
and a bridge between the "dendrogram"" and "phylo" object types that 
facilitates access to a comprehensive number of object-specific 
functions across a suite of contributed packages.
Future versions of the package will aim to further expand the range of input 
formats and object types available, 
thereby helping to integrate the wide variety of 
bioinformatic applications implemented in the R programming language.
This software is still under active development, and will
continue to be upgraded and expanded as new applications arise.
Bug reports and other feedback are welcomed and can be directed to the GitHub issues
page at <http://github.com/shaunpwilkinson/phylogram/issues>, 
or the **phylogram** google group at <https://groups.google.com/group/phylogram>.

## Acknowledgements 
This software was developed with funding from a Rutherford Foundation Postdoctoral 
Research Fellowship from the Royal Society of New Zealand. The authors declare no 
competing interests.

## Author Contributions
Both authors conceived and designed the software. 
SPW wrote the package functions and documentation.
Both authors wrote the manuscript and gave final approval for publication.

## References 