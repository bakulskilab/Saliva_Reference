.libPaths("C:/Users/HP/Documents/R/win-library/3.6")
library(devtools)
dev_mode(on = F)

if(any(endsWith(search(), "ewastools"))){detach("package:ewastools")}

install_github("hhhh5/ewastools", ref = "devel", dependencies = T)
library(ewastools)