#pGlycoFilter to make a plant glycan database from MS2 data

### With ProteoWizard
## use msconvert to convert MS files to .mzML and to filter for MS2 data
#       msconvert.exe data.RAW --filter "msLevel 2"

## use msaccess to make a spectrum table
#       msaccess.exe data.mzML -x spectrum_table

## use chainsaw for in silico digest of the N-glycoproteinSequence and of the
## same sequence with glycosylation sites having lowercase "n"
# chainsaw.exe WPP_sequence.fasta -c Chymotrypsin -n 2 -s fully
# chainsaw.exe WPP_sequence_glyco.fasta -c Chymotrypsin -n 2 -s fully

### in R
setwd("<drive>:/pGlyco")
require(protViz)
source("R/plantGlycoMSFunctions.R") 

# calculate neutral precursor monoisotopic mass values
MS2.dataTable <- spectrumTable(input="data/data.mzML.spectrum_table.csv")
write(MS2.dataTable$MonoPrecursorMass, 
      file="output/step2_experimentalMasses.txt", sep="\n")

# make a list of peptide masses
digest.N <- read.table(file="data/WPP.fasta_digestedPeptides_fully.tsv", 
                       stringsAsFactors= FALSE, header=TRUE, sep="\t")
digest.n <- read.table(file="data/WPP_glyco.fasta_digestedPeptides_fully.tsv",
                       stringsAsFactors= FALSE, header=TRUE, sep="\t")


WPP_gps <- glycoChainSaw(digest.N,digest.n,carbamidomethylation=T, glycoOnly=T)
write(WPP_gps$mass, file="output/unmodifiedPeptideMasses.txt", sep="\n")

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

### in R
# read GlycoMod results
input <- read.csv("data/gmod_glycanID.csv", header=T, stringsAsFactors=F, check.names=F)
#RUN pGlycoFilter
structure = input[,3]
output <- pGlycoFilter(structure, input)

# make a list of 'dynamic delta masses' for MassRecon config file
PreferredDeltaMasses(input=output)
## The End ##