# staCMRsetup 
# sources STACMR functions, loads required packages, and links to java runtime library
# execute this program before doing anything else
bookinfo = 'Dunn, J. C. & Kalish, M. L. (2018). State-Trace Analysis. Springer.'
cat ('STACMR program library Version 26.03.2020\n')
cat ('Utility programs for use with the book:\n')
cat (paste0(bookinfo,'\n\n'))
source ('staSTATS.R')
source ('staPLOT.R')
source ('staMR.R')
source ('staMRFIT.R')
source ('staCMR.R')
source ('staCMRFIT.R')
source ('staPLOTBN.R')
source ('staMRBN.R')
source ('staCMRBN.R')
source ('staMRFITBN.R')
source ('staCMRFITBN.R')

source ('Utility functions/staCMRx.R')
source ('Utility functions/gen2list.R')
source ('Utility functions/list2adj.R')
source ('Utility functions/adj2list.R')
source ('Utility functions/shrinkDiag.R')
source ('Utility functions/gen2listBN.R')
source ('Utility functions/staSTATSBN.R')
source ('Utility functions/BNframe2list.R')
source ('Utility functions/binSTATS.R')
source ('Utility functions/tiesort.R')
source ('Utility functions/LoftusMasson.R')

source ('java/jMR.R')
source ('java/jMRfits.R')
source ('java/jMRBNfits.R')
source ('java/jCMRxBNfits.R')
source ('java/jCMRx.R')
source ('java/jCMRBN.R')
source ('java/jCMRfitsx.R')
source ('java/jCMRBN.R')
if(require("magic")==F){install.packages ("magic"); library(magic)} 
if(require("ggplot2")==F){install.packages ("ggplot2"); library (ggplot2)}
if(require("RColorBrewer")==F){install.packages ("RColorBrewer"); library(RColorBrewer)}
if(require("rJava")==F){install.packages ("rJava"); library(rJava)}
j = list.files(path='java/',pattern = ".jar$");
if (length(j)==0) {print("Error: Java runtime library not found")
  } else {
  j=sort(j,decreasing=T); vm=paste0('java/',j[1])
  .jinit (classpath=vm) # initialize java VM
  cat(paste('STACMR linked to java library',j[1]),'\n')
}
