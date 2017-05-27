# phylogram

[![Build Status](https://travis-ci.org/shaunpwilkinson/phylogram.svg?branch=master)](https://travis-ci.org/shaunpwilkinson/phylogram)

--------------------------------------------------------------------------------

`phylogram` is an R package for viewing, editing and publishing phylogenetic trees. 
It differs from others in that trees are represented as deeply nested lists, which 
enables users to manipulate nodes recursively using fast inbuilt functions such as 
`dendrapply` and `lapply`. Trees can be ported between `phylogram`  and other 
software packages via the Newick (aka New Hampshire) text format, using the 
functions `read.dendrogram` and `write.dendrogram`.


### Installation
`phylogram` is currently available as a development version, with a stable
release available on CRAN shortly. To download the package from 
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
devtools::install_github("shaunpwilkinson/phylogram") 
library("phylogram")
```

### Help
An overview of the package and it's functions can be found by running
```R
?phylogram
```
and more detail on the individual functions can be found using the 
`?<function-name>` command.

To build the vignette users will need to have LaTeX installed. RStudio recommends 
[MiKTeX Complete](http://miktex.org/2.9/setup) for Windows and
[TexLive 2013 Full](http://tug.org/) for Mac OS X and Linux.

If you experience a problem using this package please
either raise it as an issue on [GitHub](http://github.com/shaunpwilkinson/phylogram/issues) 
or post it on the [phylogram google group](http://groups.google.com/group/phylogram).


### Acknowledgements
This software was developed at 
[Victoria University of Wellington](http://www.victoria.ac.nz/) 
with funding from a Rutherford Foundation Postdoctoral Research Fellowship 
award from the Royal Society of New Zealand.





