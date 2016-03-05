#gPSMvalidator to validate glycopeptide spectrum matches from a GlycoMod search

### With ProteoWizard
## use msconvert to convert MS files to .mzML and to filter for MS2 data
#       msconvert.exe data.RAW --filter "msLevel 2"

## use chainsaw for in silico digest of the N-glycoproteinSequence and of the
## same sequence with glycosylation sites having lowercase "n"
# chainsaw.exe WPP_sequence.fasta -c Chymotrypsin -n 2 -s semi
# chainsaw.exe WPP_sequence_glyco.fasta -c Chymotrypsin -n 2 -s semi

## use msaccess to obtain all of the MS2 data (mZ,intensity)
# msaccess.exe data.mzML -x "binary sn=1,19887" -o MS2Data

## use msaccess to make a spectrum table
#       msaccess.exe data.mzML -x spectrum_table

# calculate neutral precursor monoisotopic mass values in R
setwd("<drive>:/pGlyco/")
source("R/plantGlycoMSFunctions.R")
MS2.dataTable <- spectrumTable(input="data/data.mzML.spectrum_table.csv")
write(MS2.dataTable$MonoPrecursorMass, 
      file="output/experimentalMasses.txt", sep="\n")

### with GlycoMod (http://web.expasy.org/glycomod/) and spreadsheet program
## Run GlycoMod
# Enter a list of experimental masses: "experimentalMasses.txt"
# All mass values are: monoisotopic; Mass tolerance: +/- 15 ppm, neutral [M]
# N-linked oligosaccharides: Glycopeptides (motif N-X-S/T/C (X not P))
# Enter a set of unmodified peptide masses ([M]): "peptideMasses.tryp.txt"
# Hexose (e.g. Man, Gal): Yes 0-9; HexNAc (e.g. GlcNAc, GalNAc): Yes 1-4
# Deoxyhexose (e.g. fucose): Yes 0-3; NeuAc (e.g. sialic acid): No
# NeuGc: No; Pentose (e.g. xylose): Yes 0-1; Sulfate: No; Phosphate: No
# KDN: No; HexA (e.g. glucuronic acid): No
# (Yes) List compositions reported in UniCarbKB separately.
## Save GlycoMod results in a .csv file

#########################################################################
# 1. Prepare the environment
#########################################################################

setwd("<drive>:/pGlyco")
require(protViz)
source("R/plantGlycoMSFunctions.R") 

#########################################################################
# 2. Read data and run pGlycoFilter
#########################################################################
#read GlycoMod results
input <- read.csv("data/gmod_glycopeptideID.csv",
                header=TRUE, stringsAsFactors=F, check.names=F)

###RUN pGlycoFilter
structure = input[,3]
pGlycoFilter(structure, input) -> output

# make a list of peptide masses
digest.N <- read.table(file="data/WPP.fasta_digestedPeptides_semi.tsv", 
                       stringsAsFactors= FALSE, header=TRUE, sep="\t")
digest.n <- read.table(file="data/WPP_glyco.fasta_digestedPeptides_semi.tsv",
                       stringsAsFactors= FALSE, header=TRUE, sep="\t")

WPP_gps <- glycoChainSaw(digest.N, digest.n,carbamidomethylation=T, glycoOnly=T)

IDPdb.ALL <- Read.GlycoMod (input=output, 
                            spectrum.table="data/data.mzML.spectrum_table.csv",
                            ChainSaw=WPP_gps,
                            dir="data/MS2Data")

#####restricting the rt (how many minutes plus and minus of the rt?)

Table2 <- read.csv("data/Table2.csv")

IDPdb <- rt.restrict(IDPdb=IDPdb.ALL, Table2=Table2, rt.min.minus=5, rt.min.plus=5)

MS2Data <- Read.MS2Data(IDPdb=IDPdb, dir="data/MS2Data")

#Calculate ions (Y1, Y0)
Ions <- calculateIons(IDPdb)

###############################################################################################
##   3. Make a list to use in the data argument   
###############################################################################################

data <- compileData(analysis="glycoMod", sampleName="results", 
                    IDPdb=IDPdb, MS2Data=MS2Data, Ions=Ions)

#check the data
str(data[[1]], nchar.max = 30)        

#########################################################################
#   4. Prepare the data for the modifications argument
#########################################################################
ptm.0<-cbind(AA="-",
             mono=0.0, avg=0.0, desc="unmodified", unimodAccID=NA)

ptm.1<-cbind(AA='N',
             mono=0, avg=NA, desc="GlycanLoss",
             unimodAccID=1)

ptm.2<-cbind(AA='M',
             mono=15.994915, avg=NA, desc="Oxidation",
             unimodAccID=2)

m<-as.data.frame(rbind(ptm.0, ptm.1, ptm.2))

#########################################################################
#   5. Prepare the data for the marker ions argument
#########################################################################

otherOxonium <- c( 528.1923, 366.1400, 163.0601, 325.11292)
HexNAc_markerIons <- c(126.05495, 138.05495, 168.06552, 186.07608, 204.08665)

markerIons <- c(HexNAc_markerIons)


#########################################################################
#   6. gPSM validation
#########################################################################

pdf(file="output/gmod_annotatedSpectra.pdf", 11, 8.5)

s <- gPSMvalidator (data=data,
                    modification=m$mono, modificationName=m$desc,
                    mZmarkerIons=markerIons, minNumberIons=2, itol_ppm=5,
                    minMarkerIntensityRatio=2, PEAKPLOT = TRUE, validate=TRUE)

dev.off()

write.table(s, file="output/gmod_infoTable.csv", sep=",", row.names=FALSE,
            col.names=TRUE, quote=FALSE)        

##### The End ###########################################################