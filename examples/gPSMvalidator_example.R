#gPSMvalidator to validate glycopeptide spectrum matches from a massRecon search

### With ProteoWizard
## use msconvert to convert MS files to .mzML and to filter for MS2 data
#       msconvert.exe data.RAW --filter "msLevel 2"

## use chainsaw for in silico digest of the N-glycoproteinSequence and of the
## same sequence with glycosylation sites having lowercase "n"
# chainsaw.exe WPP_sequence.fasta -c Chymotrypsin -n 2 -s semi
# chainsaw.exe WPP_sequence_glyco.fasta -c Chymotrypsin -n 2 -s semi

## use msaccess to obtain all of the MS2 data (mZ,intensity)
# msaccess.exe data.mzML -x "binary sn=1,19887" -o MS2Data

### with BumberShoot (see supplementary information for the cfg file)
# c:tagrecon.exe data.mzML -cfg massReconConfiguration.cfg
# open results in IDPicker and export the results (massRecon_results.tsv)

### in R
#########################################################################
# 1. Prepare the R environment
#########################################################################
setwd("<drive>:/pGlyco")
library("protViz")
source("R/plantGlycoMSFunctions.R") 

#########################################################################
# 2. Read in data PSMs & gPSMs, in silico digest, and MS2 binary data)
#########################################################################
# make a data.frame of glycopeptide masses

digest.N <- read.table(file="data/WPP.fasta_digestedPeptides_semi.tsv", 
                       stringsAsFactors= FALSE, header=TRUE, sep="\t")
digest.n <- read.table(file="data/WPP_glyco.fasta_digestedPeptides_semi.tsv",
                       stringsAsFactors= FALSE, header=TRUE, sep="\t")


WPP_gps <- glycoChainSaw(digest.N,digest.n,carbamidomethylation=T, glycoOnly=T)

# read in massRecon results

massRecon_results <- read.table(file="data/massRecon_results.tsv", 
                                sep="\t", header=T, stringsAsFactors=F)

IDPdb <- Read.IDPdb(IDPdb=massRecon_results, ChainSaw=WPP_gps,dir="data/MS2Data")

MS2Data <- Read.MS2Data(IDPdb=IDPdb, dir="data/MS2Data")

# Calculate ions (Y1, Y0, F1)

Ions <- calculateIons(IDPdb)

################################################################################
#   3. Make a list to use in the data argument of gPSMvalidator
###############################################################################

data <- compileData(analysis="massRecon", sampleName="results", 
                    IDPdb=IDPdb, MS2Data=MS2Data, Ions=Ions)

  #check the data
str(data[[400]], nchar.max = 30)	


#########################################################################
#   4. Prepare the data for the modifications argument of gPSMvalidator
#########################################################################

ptm.0<-cbind(AA="-",
             mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)

ptm.1<-cbind(AA='N',
             mono=0, avg=NA, desc="GlycanLoss",
             unimodAccID=1)

ptm.2<-cbind(AA='M',
             mono=15.994915, avg=NA, desc="Oxidation",
             unimodAccID=2)

ptm.3<-cbind(AA='K',
             mono=57.021464, avg=NA, desc="IAM artifact",
             unimodAccID=3)
ptm.4<-cbind(AA='H',
             mono=57.021464, avg=NA, desc="IAM artifact",
             unimodAccID=4)
ptm.5<-cbind(AA='D',
             mono=57.021464, avg=NA, desc="IAM artifact",
             unimodAccID=5)
ptm.6<-cbind(AA='E',
             mono=57.021464, avg=NA, desc="IAM artifact",
             unimodAccID=6)
ptm.7<-cbind(AA='S',
             mono=57.021464, avg=NA, desc="IAM artifact",
             unimodAccID=7)
ptm.8<-cbind(AA='T',
             mono=57.021464, avg=NA, desc="IAM artifact",
             unimodAccID=8)
ptm.9<-cbind(AA='Y',
             mono=57.021464, avg=NA, desc="IAM artifact",
             unimodAccID=9)



m<-as.data.frame(rbind(ptm.0, ptm.1, ptm.2,
                       ptm.3, ptm.4, ptm.5,
                       ptm.6, ptm.7, ptm.8,
                       ptm.9))


#########################################################################
#   5. Prepare the data for the marker ions argument
#########################################################################

otherOxonium <- c( 528.1923, 366.1400, 163.0601, 325.11292)
HexNAc_markerIons <- c(126.05495, 138.05495, 168.06552, 186.07608, 204.08665)

markerIons <- c(HexNAc_markerIons)

#########################################################################
#   6. gPSM validation
#########################################################################

pdf(file="output/annotatedSpectra.pdf", 11, 8.5)

s <- gPSMvalidator (data=data,
    modification=m$mono, modificationName=m$desc,
    mZmarkerIons=markerIons, minNumberIons=2, itol_ppm=5,
    minMarkerIntensityRatio=2, PEAKPLOT = TRUE, validate=FALSE)
 	
dev.off()

write.table(s, file="output/infoTable.csv", sep=",", row.names=FALSE,
            col.names=TRUE, quote=FALSE)        

##### The End ###########################################################