# First, Rtools must be installed if using a Windows system
# https://cran.r-project.org/bin/windows/Rtools
# Version 4.4. of Rtools is not currently supported, please use an older version

# On a Mac computer, install Xcode and check out the details about installing
# gfortran at:
# https://github.com/fxcoudert/gfortran-for-macOS/releases



#Finally, run the following lines in R:

#Install packages required for TrainSel if they are 
#not already installed:
list.of.packages <- c("doParallel", 
                      "RcppArmadillo",
                      "RcppEigen",
                      "foreach",
                      "RcppProgress",
                      "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

#Install TrainSel itself:
library(devtools)
install_github("TheRocinante-lab/TrainSel")



######################################
# Further instructions for Mac systems
######################################

#If you are having trouble with the installation in Mac, there may be 
#some problems with the gfortran installation. You may perform the 
#following checks: 

################################
# Checking gfortran Installation
################################

# To check where gfortran is installed, open a terminal and run the following
# command:
#   find /usr/local -name "libgfortran.*"
# 
# Your results should be similar to the following:
# /usr/local/gfortran/lib/libgfortran.a
# /usr/local/gfortran/lib/libgfortran.dylib
# /usr/local/gfortran/lib/libgfortran.la
# /usr/local/gfortran/lib/libgfortran.spec
# /usr/local/gfortran/lib/libgfortran.5.dylib

################################
# Updating the Makevars File
################################

# Now that we know the libgfortran libraries are located in
# /usr/local/gfortran/lib/, we need to update the Makevars file.

# Open the Makevars file by typing the following command in the terminal:
#   nano ~/.R/Makevars

# Add the following lines to the file:
#   FC = /usr/local/bin/gfortran
#   F77 = /usr/local/bin/gfortran
#   FLIBS = -L/usr/local/gfortran/lib -lgfortran
#   CC = /usr/bin/clang
#   CXX = /usr/bin/clang++

#   If your location of gfortran is different from /usr/local/bin, then use
# that path. For example, if the location is /usr/local/trend/gfortran/lib,
# then you will have to use:

#   FC = /usr/local/trend/gfortran
#   F77 = /usr/local/trend/gfortran
#   FLIBS = -L/usr/local/gfortran/lib -lgfortran
#   CC = /usr/bin/clang
#   CXX = /usr/bin/clang++


#   Save the changes by pressing Ctrl + O, then press Enter to confirm. Exit
# nano by pressing Ctrl + X.

