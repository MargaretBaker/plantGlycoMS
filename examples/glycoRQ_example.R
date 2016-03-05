#########################################################################
# 1. Prepare the environment
#########################################################################
setwd("<drive>:/pGlyco/")
library("protViz")
source("R/plantGlycoMSFunctions.R") 

###########################
# 1. make the config file for msaccess
###############################

RQ_table <- read.csv("data/RQ_input.csv", stringsAsFactors=F)

RQ_config <- paste("exec = sic mzCenter=", RQ_table$exact.precursor.mz, 
                   " radius=5 radiusUnits=ppm delimiter=comma", sep="")

write(RQ_config, "data/RQ_Data/config.txt", sep="/n")

### using ProteoWizard, run msaccess with working directory: /data/RQ_Data/
### The .mzML file must be filtered to contain only MS1 data!
# c:msaccess *.mzML -c config.txt

###########################
# 2. run glycoRQ
###############################

Table2 <- read.csv("data/Table2.csv", header=TRUE, stringsAsFactors=FALSE)

RQ <- Read.RQ(input="data/RQ_input.csv", dir="data/RQ_Data")

RQresults <- glycoRQ(RQ, Table2, dir="data/RQ_Data", rt.min.minus=1, rt.min.plus=3)

