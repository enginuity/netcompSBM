##@S Loads all needed libraries and functions
##@S Sourcing this file (from this directory) will load all the functions into the current global environment (since this runs 'source' on all the codefiles)

## Load own library functions
allfiles = list.files("netcompSBM/R/", full.names = T)
for(s in allfiles) source(s)

## Source other libraries
require(abind)
