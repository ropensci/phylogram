# phylogram

--------------------------------------------------------------------------------

`phylogram` is an R package for developing evolutionary trees 
as deeply-nested lists known as "dendrogram" objects. 
It provides functions for importing and exporting trees in the Newick 
parenthetic text format, as well as several functions for command-line 
tree manipulation.
With an emphasis on speed and computational efficiency, `phylogram` also 
includes a suite of tools for rapidly computing distance matrices and 
building large trees using fast alignment-free k-mer counting and 
divisive clustering techniques.
This package makes R's powerful nested-list architecture more 
accessible to evolutionary biologists, and facilitates the analysis 
of very large sequence datasets.


### Installation
To download `phylogram` from CRAN and load the package, run
```R
install.packages("phylogram")
library("phylogram")
```
To download the development version from 
GitHub you will first need to ensure you have a C/C++ compliler and the 
[devtools](https://github.com/hadley/devtools) R package installed. 
Linux users will generally have a compiler such as `gcc` installed by default; 
however Windows users will need to download 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) and Mac 
OSX users will need [Xcode](https://developer.apple.com/xcode) 
(note that Rtools and Xcode are not R packages). To download and install 
devtools, run 
```R
install.packages("devtools")
``` 
Then install and load the `phylogram` package by running 
```R
devtools::install_github("shaunpwilkinson/phylogram", build_vignettes = TRUE) 
library("phylogram")
```

### Use and Examples
#### Example 1: reading and writing trees
Consider the simple example of a tree with three members named 
"A", "B" and "C", where "B" and "C" are more closely related
to eachother than they are to "A". 
An unweighted Newick string for this tree would be "(A,(B,C));".
This text can be imported as a 
dendrogram object using the `read.dendrogram` function 
as follows:

```R
library("phylogram")
newick <- "(A,(B,C));"
x <- read.dendrogram(text = newick)
plot(x)
```

The following command writes the object back to the console in 
Newick format without edge weights:
```R
write.dendrogram(x, edges = FALSE)
```
The syntax is similar when reading and writing text files, 
except that the `text` argument is replaced by `file` and a 
valid file path is passed to the function.

#### Example 2: building and editing trees
The function `topdown` builds a tree by divisive clustering.
This is done by counting k-mers and recursively partitioning 
the sequence set using successive k-means clustering steps. 
No alignment is necessary and no distance matrix is computed,
making it possible to rapidly and efficiently build trees 
from very large sequence datasets.

This following code demonstrates how to build and plot a divisive 
tree using the `woodmouse` data from the ape package:

```R
library("phylogram")
library("ape")
data(woodmouse)
x <- topdown(woodmouse, k = 5, nstart = 20)
op <- par(no.readonly = TRUE)
par(mar = c(4, 4, 4, 5))
plot(x, horiz = TRUE)
par(op)
```
These and more examples are available in the package vignette.
If downloading the package from github, users will need to have 
LaTeX installed to build the vignette. RStudio recommends 
[MiKTeX Complete](http://miktex.org/2.9/setup) for Windows and
[TexLive 2013 Full](http://tug.org/) for Mac OS X and Linux. 
To view the vignette, run `vignette(package = "phylogram")`


### Help
An overview of the package with links to the function documentation can be found by running
```R
?phylogram
```
If you experience a problem using this package please
either raise it as an issue on [GitHub](http://github.com/shaunpwilkinson/phylogram/issues) 
or post it on the [phylogram google group](http://groups.google.com/group/phylogram).


### Acknowledgements
This software was developed at 
[Victoria University of Wellington](http://www.victoria.ac.nz/) 
with funding from a Rutherford Foundation Postdoctoral Research Fellowship 
award from the Royal Society of New Zealand.
