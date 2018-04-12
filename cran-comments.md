#phylogram version 2.0.0

This release follows the migration of several k-mer counting and clustering functions
to the new package **kmer**. The package no longer needs compilation since all C++
source code has been removed.

## Test environments

 * local ubuntu 16.04.2 x86_64-pc-linux-gnu; R version 3.4.0 
 * travis-ci ubuntu 14.04.5 x86_64-pc-linux-gnu; R 3.4.4
 * winbuilder devel R version 3.5.0 beta (2018-04-08 r74552)

## R CMD check results

There were no ERRORs WARNINGs or NOTEs.

## Downstream dependencies

kmer:  OK
