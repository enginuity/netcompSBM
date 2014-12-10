##@S Contains code to compile the netcomp_sbm package (to be run from R)


# Optional Steps ----------------------------------------------------------
## These only need to be run, if the source code was edited. 

## Load up needed libraries
require(codeProcessing) ## This is a package I've written to speed up my own coding
require(stringr)

## Update documentation -- this looks for new parameters & 
update_fx_documentation(FD = FilesDescription(dirlist = "netcompSBM/R/"), fill_emptyparam = FALSE)


# Compile Package ---------------------------------------------------------
## Both these steps need to be run in order to build the package properly

## Generate the documentation -- THIS MUST be run before building packages (since the documentation files are not version-controlled, as the version-controlled version is in the raw source code)
require(roxygen2)
roxygenise("netcompSBM/", clean = TRUE)

## Install the package
system("R CMD INSTALL netcompSBM")

