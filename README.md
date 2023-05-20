# TrainSel
Selection of training populations

To install the package in R use:

```
library(devtools)
install_github("TheRocinante-lab/TrainSel")
On MAC Install Xcode from apple, and gfortran6.3 or gfortran8.2 from 
https://github.com/rmacoslib/r-macos-rtools/
Please install also Rcpp and Rcpp-Armadillo packages.
```

Read LICENSE file for profit organizations.


If you are encountering issues related to gfortran6.3 while installing packages on a Silicon Mac, it could be due to compatibility issues with the ARM architecture. Currently, there might be limitations with certain packages that require gfortran on Silicon Macs.

One possible solution is to try installing the package using Rosetta 2, which allows you to run Intel-based software on Apple Silicon Macs. Follow these steps:

1. Open a Terminal window on your Mac.
2. Type the following command to open a new Terminal session running under Rosetta 2:

```
arch -x86_64 /bin/bash
```
3. Install R by running:

```
curl -O https://cran.rstudio.com/bin/macosx/base/R-4.1.0.pkg
sudo installer -pkg R-4.1.0.pkg -target /
```
4. Launch R under Rosetta 2 by running the command:

```
arch -x86_64 /usr/local/bin/R
```
5. Install the 'TrainSel' package using the devtools::install_github() 
6. Close Rstudio





