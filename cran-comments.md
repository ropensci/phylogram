This is a patch release addressing a solaris build error discovered by Prof. Ripley:
"error: call of overloaded ‘pow(int, int&)’ is ambiguous"
The patch replaces C++ calls to 'pow' with recursive multiplication.

## Test environments
 * local ubuntu 16.04 x86_64-pc-linux-gnu; R version 3.4.0 
 * travis-ci ubuntu 12.04.5 x86_64-pc-linux-gnu; R 3.4.0
 * winbuilder devel (2017-06-12 r72786)

## R CMD check results
There were no ERRORs or WARNINGs.
The winbuilder check returned one NOTE 
"Possibly mis-spelled words in DESCRIPTION:
  Dendrograms (3:8)"
This spelling is correct.

## Downstream dependencies
There were no issues.
