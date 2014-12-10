##@S Contains code to compile the netcomp_sbm package (to be run from R)

## Load up needed libraries
require(metacode)
require(stringr)
require(roxygen2)

## Update documentation, and re-run Roxygen on it
update_fx_documentation(FD = FilesDescription(dirlist = "netcomp_sbm/R/"), fill_emptyparam = FALSE)
roxygenise("netcomp_sbm/", clean = TRUE)

## Install the package
system("R CMD INSTALL netcomp_sbm")
