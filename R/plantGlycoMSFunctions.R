################################################################################
#   Function 1. Add in the spectrum table and calculate the MonoPrecursorMass
#
# input is the spectrum table obtained from an msaccess search (proteowizard)
# this function is needed to calculate the monoisotopic precursor mass from the
# precursor mz
################################################################################

spectrumTable <- function (input) {
        
spectrumTable <- read.csv(file=input, header=TRUE, skip=1)

        
spectrumTable$MonoPrecursorMass <- (spectrumTable$precursorMZ * 
                                            spectrumTable$charge - 
                                            (spectrumTable$charge * 1.007276 )) 
        
        return(spectrumTable)
        
        
}

###############################################################################################################
# Function 1: glycoChainSaw
# this function takes in results from chainsaw with glycosylation site labeled (digest.n)
# and unlabled (digest.N). It combines these and optionally adds the mass for carbamidomethylation
# the output is a dataframe that has a column identifying the glycosylation site in the sequence with "n"
# it can optionally return a dataframe containing only data for peptide with glycosylation sites
###############################################################################################################

glycoChainSaw <- function(digest.N, digest.n, carbamidomethylation=TRUE, glycoOnly =TRUE) {
        
        
        names(digest.n) <- c("n.sequence", "protein", "mass"  ,  "missedCleavages",  "specificity",     
                             "nTerminusIsSpecific" ,"cTerminusIsSpecific" )
        
        digest.n$sequence <- gsub("n", "N", digest.n$n.sequence)
        
        digest.Nn <- merge(digest.N, digest.n, by="sequence")
        
        digest.Nn <- digest.Nn[c("sequence" , "protein.x" ,  "mass.x" ,  "missedCleavages.x" ,   
                                 "specificity.x"   ,      "nTerminusIsSpecific.x", "cTerminusIsSpecific.x" ,"n.sequence"  )]         
        
        names(digest.Nn) <- c("sequence" ,  "protein" ,  "mass"    ,   "missedCleavages" ,   
                              "specificity"   ,      "nTerminusIsSpecific", "cTerminusIsSpecific" ,"n.sequence"  )
        digest.Nn <- digest.Nn[!duplicated(digest.Nn$sequence),]
        
        if (carbamidomethylation ==TRUE ) {
                modMass <- gsub( "C", "c", digest.Nn$sequence ,fixed=TRUE)
                modMass <- gsub( "[A-Z]", "0", modMass)
                modMass <- strsplit(modMass , split="")
                modMass <- lapply(modMass , function(x){gsub("c", "57.021464", x)})
                modMass <- lapply(modMass , function(x){as.numeric(x)})
                modMass <- sapply(modMass , function(x) {sum(x)})
                digest.Nn$mass <- mapply(sum, modMass , digest.Nn$mass)
        }
        
        if (glycoOnly==TRUE){
                glycoOnly <- grep("n", digest.Nn$n.sequence)
                digest.Nn <- digest.Nn[glycoOnly,]
                
        }
        
        return(digest.Nn)
}


################################################################################
#  3. GlycoMod Results (filtered) - PreferredDeltaMasses
# input for this function is the output of pGlycoFilter
################################################################################

PreferredDeltaMasses <- function (input)
{
        
        deltaMass<- paste("N!{P}[ST] * ", input$glycoform)
        deltaMass <- unique(deltaMass)
        cat(deltaMass)

}


##############################################################################
#     Function 5. Read in data 
# input: table of PSMs from IDPicker
# ChainSaw: table of peptides from digest software like chainsaw
# dir: a directory containing all of the MS2 data files
##############################################################################

Read.IDPdb <- function (IDPdb, ChainSaw,dir) {
        

        IDPdb <- IDPdb[, -c(2:6)]
       IDPdb <- IDPdb[-c(1,2),]
       
      # x <-  data.frame (IDPdb$Group.Source.Spectrum, IDPdb$Sequence)
     #  IDPdb <- IDPdb[!duplicated(x),]
        
        ppm.mass.error <- IDPdb$Mass.Error/IDPdb$Exact.Mass *10^6
        IDPdb$ppm.mass.error <- ppm.mass.error     
        
        
        # change: factor 0.1.4586  to: numeric 4586 
        sc <- IDPdb$Group.Source.Spectrum
        sc <- strsplit(sc, split=".", fixed=T)
        scans <- sapply(sc, function(x) x[3])
        scans <- as.numeric(scans)
        IDPdb$Group.Source.Spectrum <- scans
        
        # Make a data frame 'MS2.all' with columns 'title' and 'scans'; 
        # it contains this info for all MS2 binary files
        
        title <- list.files(dir)
        
        #chr "ljz_20131022_MR_Chym1_HILIC_MS2.mzML.binary.sn1830.txt" to: num  10011
        sc <- strsplit(title, split="sn")
        scans <- sapply(sc, function(x) x[2])
        scans <- gsub(pattern=".txt", replacement="", x=scans)
        Group.Source.Spectrum <- as.numeric(scans)
        
        # itle and scan number. data will be filled in after merging with IDPdb
        MS2.all <- data.frame(Group.Source.Spectrum, title) 
        IDPdb <- merge(  x = MS2.all,y = IDPdb,  by = "Group.Source.Spectrum", all.y = TRUE)
        IDPdb$title <- as.character(IDPdb$title)        
        
        
        #Modify the column contents of IDPdb
        
        # change: chr "QHGFTMM[16]NVYN[1170]STK" to: chr "QHGFTMMNVYNSTK"
        peptideSequence <- IDPdb$Sequence        
        peptideSequence <- gsub(pattern="[0123456789]", replacement="", peptideSequence)
        peptideSequence <- gsub(pattern="[]", replacement="", peptideSequence, 
                                fixed = TRUE)
     
        IDPdb$peptideSequence <- peptideSequence
        
        # change: chr "QHGFTMM[16]NVYN[1170]STK" to: chr "QHGFTMM16NVYNSTK"
        table2Sequence <- IDPdb$Sequence        
        table2Sequence <- gsub( "C[57]", "C", table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub( "M[16]", "$",  table2Sequence ,fixed=TRUE)
        
        table2Sequence <- gsub( "K[57]", "@",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub( "H[57]", "#",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub( "D[57]", "%",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub( "E[57]", "^",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub( "S[57]", "&",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub( "T[57]", "*",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub( "Y[57]", "(",  table2Sequence ,fixed=TRUE)
        
        
        table2Sequence <- gsub(pattern="[0123456789]", replacement="", table2Sequence)
        table2Sequence <- gsub(pattern="[]", replacement="", table2Sequence, 
                               fixed = TRUE)
        
        table2Sequence <- gsub(  "$", "M16",  table2Sequence ,fixed=TRUE)
        
        table2Sequence <- gsub(  "@", "K57",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub(  "#", "H57",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub(  "%", "D57",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub(  "^", "E57",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub(  "&", "S57",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub(  "*", "T57",  table2Sequence ,fixed=TRUE)
        table2Sequence <- gsub(  "(", "Y57",  table2Sequence ,fixed=TRUE)
        
        IDPdb$table2Sequence <- table2Sequence
       
        # change: chr "QHGFTMM[16]NVYN[1170]STK" to: chr "0000000200010000"
        modification <- IDPdb$Sequence
        modification<- gsub(pattern="[0123456789]", replacement="0", 
                            modification, fixed=FALSE)
        modification <- gsub( "C[00]", "C", modification ,fixed=TRUE)
        modification <- gsub( "M[00]", "2",  modification ,fixed=TRUE)
        
        modification <- gsub( "K[00]", "3",  modification ,fixed=TRUE)
        modification <- gsub( "H[00]", "4",  modification ,fixed=TRUE)
        modification <- gsub( "D[00]", "5",  modification ,fixed=TRUE)
        modification <- gsub( "E[00]", "6",  modification ,fixed=TRUE)
        modification <- gsub( "S[00]", "7",  modification ,fixed=TRUE)
        modification <- gsub( "T[00]", "8",  modification ,fixed=TRUE)
        modification <- gsub( "Y[00]", "9",  modification ,fixed=TRUE)
        
        modification<- gsub( "N[000]", "1",  modification,fixed=TRUE)
        modification<- gsub( "N[0000]", "1",  modification,fixed=TRUE)
        modification <- gsub( "[A-Z]", "0", modification)
        modification[1:length(modification)] <- paste("0", 
                                                      modification[1:length(modification)], "0", sep="")
        IDPdb$modification <- modification
        
        IDPdb$charge <- as.numeric(IDPdb$Charge)
        
        names(ChainSaw) <- c("peptideSequence", "protein", "mass", "missedCleavages","specificity",        
                             "nTerminusIsSpecific", "cTerminusIsSpecific", "n.sequence") 
        
        IDPdb <- merge( x = ChainSaw, y = IDPdb, by = "peptideSequence" )
        
        # make the GlycanMass column for IDPdb data frame 
        # change: chr "Q[-17]HGFTMM[16]NVYN[1170]STK"
        #     to: num 1170
        
        GlycanMass <- IDPdb$Sequence
        GlycanMass <- gsub( "C[57]", "C", GlycanMass ,fixed=TRUE)
        
        GlycanMass <- gsub( "M[16]", "M",  GlycanMass ,fixed=TRUE)
        
        GlycanMass <- gsub( "K[57]", "K",  GlycanMass ,fixed=TRUE)
        GlycanMass <- gsub( "H[57]", "H",  GlycanMass ,fixed=TRUE)
        GlycanMass <- gsub( "D[57]", "D",  GlycanMass ,fixed=TRUE)
        GlycanMass <- gsub( "E[57]", "E",  GlycanMass ,fixed=TRUE)
        GlycanMass <- gsub( "S[57]", "S",  GlycanMass ,fixed=TRUE)
        GlycanMass <- gsub( "T[57]", "T",  GlycanMass ,fixed=TRUE)
        GlycanMass <- gsub( "Y[57]", "Y",  GlycanMass ,fixed=TRUE)
        
        
        GlycanMass <- gsub( "[A-Z]", "", GlycanMass)
        
        GlycanMass <- strsplit(GlycanMass , split="][", fixed=TRUE)
        GlycanMass <- lapply(GlycanMass, function(x){
                gsub("[", "", x, fixed=TRUE)})
        GlycanMass <- lapply(GlycanMass, function(x){
                gsub("]", "", x, fixed=TRUE)})
        GlycanMass <- lapply(GlycanMass, function(x){
                as.numeric(x)})
        GlycanMass <- lapply(GlycanMass, function(x){
                sum(x)})
        GlycanMass <- unlist(GlycanMass)
        
        IDPdb$GlycanMass <- GlycanMass
        
        # Calculate Y1 values
        
        # make modMass column for IDPdb data.frame; use to calculate MonoisotopicY1mass
        # change: chr "QHGFTMM[16]NVYN[1170]STK"
        #     to: num [1:14] 0 0 0 0 0 0 16 0 0 0 203 0 0 0
        #     to: num 219
        
        modMass <- IDPdb$Sequence
        modMass <- gsub( "C[57]", "C", modMass ,fixed=TRUE)
        modMass <- gsub( "M[16]", "m",  modMass ,fixed=TRUE)
        modMass <- gsub( "K[57]", "k",  modMass ,fixed=TRUE)
        modMass <- gsub( "H[57]", "h",  modMass ,fixed=TRUE)
        modMass <- gsub( "D[57]", "d",  modMass ,fixed=TRUE)
        modMass <- gsub( "E[57]", "e",  modMass ,fixed=TRUE)
        modMass <- gsub( "S[57]", "s",  modMass ,fixed=TRUE)
        modMass <- gsub( "T[57]", "t",  modMass ,fixed=TRUE)
        modMass <- gsub( "Y[57]", "y",  modMass ,fixed=TRUE)
        modMass <- gsub(pattern="[0123456789]", replacement="0", modMass )
        modMass <- gsub( "N[000]", "n",  modMass ,fixed=TRUE)
        modMass <- gsub( "N[0000]", "n",  modMass ,fixed=TRUE)
        modMass <- gsub( "[A-Z]", "0", modMass)
        modMass <- strsplit(modMass , split="")
        modMass <- lapply(modMass, function(x){gsub("m", "15.994915", x)})
        modMass <- lapply(modMass, function(x){gsub("k", "57.021464", x)})
        modMass <- lapply(modMass, function(x){gsub("h", "57.021464", x)})
        modMass <- lapply(modMass, function(x){gsub("d", "57.021464", x)})
        modMass <- lapply(modMass, function(x){gsub("e", "57.021464", x)})
        modMass <- lapply(modMass, function(x){gsub("s", "57.021464", x)})
        modMass <- lapply(modMass, function(x){gsub("t", "57.021464", x)})
        modMass <- lapply(modMass, function(x){gsub("y", "57.021464", x)})
        modMass <- lapply(modMass, function(x){gsub("n", "0", x)})
        modMass <- lapply(modMass, function(x){as.numeric(x)})
        modMass <- sapply(modMass, function(x) {sum(x)})
        IDPdb$modMass <- modMass
        IDPdb$MonoisotopicPeptideMass <- mapply(sum, IDPdb$modMass , IDPdb$mass)
        
        glycoModMass <- IDPdb$Sequence
        glycoModMass <- gsub( "C[57]", "0", glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "M[16]", "0",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "K[57]", "0",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "H[57]", "0",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "D[57]", "0",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "E[57]", "0",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "S[57]", "0",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "T[57]", "0",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "Y[57]", "0",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub(pattern="[0123456789]", replacement="0", 
                             glycoModMass )
        glycoModMass <- gsub( "N[000]", "n",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "N[0000]", "n",  glycoModMass ,fixed=TRUE)
        glycoModMass <- gsub( "[A-Z]", "0", glycoModMass)
        glycoModMass <- strsplit(glycoModMass , split="")
        glycoModMass <- lapply(glycoModMass,function(x){gsub("n" , "203.079373", x, 
                                                             fixed=TRUE)})
        glycoModMass <- lapply(glycoModMass, function(x){as.numeric(x)})
        glycoModMass <- sapply(glycoModMass, function(x) {sum(x)})
        is.na(glycoModMass) <- glycoModMass ==0
        IDPdb$glycoModMass <- glycoModMass 
        IDPdb$MonoisotopicY1mass <- mapply(sum, IDPdb$glycoModMass , IDPdb$MonoisotopicPeptideMass)
        
        MonoisotopicY1mass <- IDPdb$MonoisotopicY1mass
        MonoisotopicY1mass[is.na(MonoisotopicY1mass)] <- 0
        IDPdb$MonoisotopicY1mass <- MonoisotopicY1mass
        
        
        
        
        return(IDPdb) 
        
}


###########################################################
#        Function 6. Read in MS2 data with IDPdb
###########################################################

Read.MS2Data <- function (IDPdb, dir="convertedData/MS2Data_Chym1_ELUTE") {
        
read.dat <- function(file="ljz_20131022_MR_Chym2_ELUTE.mzML.binary.sn1801.txt") 
        {				
                
        dat <- read.table(paste(dir, "/", file, sep=""), sep="\t", skip=18)
        names(dat) <- c("mZ", "intensity")
        return(dat)
        }
        
# Make a list 'MS2Data' containing all MS2 binary data
# these are all of the titles matching ids in IDPdb
        title.identified <- IDPdb$title						
        
        MS2Data <- vector( mode="list", length=length(title.identified))
        MS2Data <- lapply( title.identified, read.dat )
        names(MS2Data) <- title.identified
        
# remove all mZ, intensity pairs with intensity=0
        MS2Data <- lapply( MS2Data, function(x) 
        {subset(x, x$intensity >0)})
        
        return(MS2Data)
}

###############################################################################
#  3. Calculate the Y1, Y0, Y2, and Y3 masses and save each to a list
# These ions are specific to a particular glycopeptide and are key to 
# validating the gPSMs.
##############################################################################

calculateIons <- function (IDPdb) {
        
        Ions <- vector( mode="list", length=length(1:nrow(IDPdb)))
        
        ii <- 0        							
        
        for(i in 1:length(1:nrow(IDPdb))) {			
                
                ii = ii +1 

Ions[[ii]]$PepPlus <- 
        find_PepPlus(MonoisotopicPeptideMass=IDPdb$MonoisotopicPeptideMass[i],charge=1:IDPdb$charge[i])   

Ions[[ii]]$Y0NH3 <- 
        find_Y0NH3(MonoisotopicPeptideMass=IDPdb$MonoisotopicPeptideMass[i],charge=1:IDPdb$charge[i]) 
        
Ions[[ii]]$Y1 <- 
        find_Y1(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])

Ions[[ii]]$Y1Mox <- 
        find_Y1Mox(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])

Ions[[ii]]$Y1Ca <- 
        find_Y1Ca(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])

Ions[[ii]]$Y1F <- 
        find_Y1F(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])

Ions[[ii]]$Y2F <- 
        find_Y2F(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])

Ions[[ii]]$Y2 <- 
        find_Y2(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])

Ions[[ii]]$Y3 <- 
        find_Y3(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])

Ions[[ii]]$Y3X <- 
        find_Y3X(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])
Ions[[ii]]$Y3F <- 
        find_Y3F(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])
Ions[[ii]]$Y3FX <- 
        find_Y3FX(MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i],charge=1:IDPdb$charge[i])

}



return(Ions)

}
#########################################
#        Function 7. Calculate the Y1 masses
#########################################

# identify the Y1 (MonoisotopicY1mass) peak.

find_Y1 <- function(MonoisotopicY1mass, charge ) {
        
        Y1 <- ((MonoisotopicY1mass + ( charge * 1.007276 ) ) / charge)
        
        return(Y1)
        
}

#########################################
#        Function 7. Calculate the Y1Mox masses, ie the mass if lost methane sulfenic acid
#########################################


find_Y1Mox <- function(MonoisotopicY1mass, charge ) {
        
        Y1Mox <- (((MonoisotopicY1mass-64.10686) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y1Mox)
        
}
#########################################
#        Function 7. Calculate the Y1Ca masses, ie the mass if lost carbamidomethyl
#########################################


find_Y1Ca <- function(MonoisotopicY1mass, charge ) {
        
        Y1Ca <- (((MonoisotopicY1mass-57.021464) + ( charge * 1.007276 ) ) 
                  / charge)
        
        return(Y1Ca)
        
}


#########################################
#        Function 8. Calculate the Peptide masses 
#########################################

# identify the Peptide masses (of glycopeptides that lost the whole glycan)

find_PepPlus <- function(MonoisotopicPeptideMass, charge ) {
        
        PP <- ((MonoisotopicPeptideMass + ( charge * 1.007276 ) ) / charge)
        
        return(PP)
        
}

#########################################
#        Function 8. Calculate the Peptide masses with loss of ammonia
#########################################

# identify the Peptide masses (of glycopeptides that lost the whole glycan)

find_Y0NH3 <- function(MonoisotopicPeptideMass, charge ) {
        
        Y0NH3 <- (((MonoisotopicPeptideMass- 17.026549) + ( charge * 1.007276 ) ) 
                  / charge)
        
        return(Y0NH3)
        
}

#########################################
#        Function 7. Calculate the Y1F masses, ie the mass if core fucose present
#########################################


find_Y1F <- function(MonoisotopicY1mass, charge ) {
        
        Y1F <- (((MonoisotopicY1mass+146.057909) + ( charge * 1.007276 ) ) 
               / charge)
        
        return(Y1F)
        
}

#########################################
#        Function 7. Calculate the Y1F masses, ie the mass if core fucose present
#########################################


find_Y2F <- function(MonoisotopicY1mass, charge ) {
        
        Y2F <- (((MonoisotopicY1mass+146.057909+203.079373) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y2F)
        
}

#########################################
#        Function 7. Calculate the Y1F masses, ie the mass if core fucose present
#########################################


find_Y3F <- function(MonoisotopicY1mass, charge ) {
        
        Y3F <- (((MonoisotopicY1mass+146.057909+203.079373+162.052824) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y3F)
        
}

#########################################
#        Function 7. Calculate the Y1F masses, ie the mass if core fucose present
#########################################


find_Y3FX <- function(MonoisotopicY1mass, charge ) {
        
        Y3FX <- (((MonoisotopicY1mass+146.057909+203.079373+162.052824+132.042259) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y3FX)
        
}


#########################################
#        Function 7. Calculate the Y2 masses, ie the mass if core fucose present
#########################################


find_Y2 <- function(MonoisotopicY1mass, charge ) {
        
        Y2 <- (((MonoisotopicY1mass+203.079373 ) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y2)
        
}

#########################################
#        Function 7. Calculate the Y3 masses, ie the mass if core fucose present
#########################################


find_Y3 <- function(MonoisotopicY1mass, charge ) {
        
        Y3 <- (((MonoisotopicY1mass+ 203.079373 + 162.052824  ) + ( charge * 1.007276 ) ) 
               / charge)
        
        return(Y3)
        
}

#########################################
#        Function 7. Calculate the Y3X masses, ie the mass if core xylose present
#########################################


find_Y3X <- function(MonoisotopicY1mass, charge ) {
        
        Y3X <- (((MonoisotopicY1mass+ 203.079373 + 162.052824 +  132.042259 ) + ( charge * 1.007276 ) ) 
               / charge)
        
        return(Y3X)
        
}

###############################################################################
#  3. Compile all of the data into one list
##############################################################################

compileData <- function (analysis="analysis", sampleName="results", 
                         IDPdb=IDPdb, MS2Data=MS2Data, Ions=Ions) {
        
        data <- vector( mode="list", length=nrow(IDPdb))
        
        ii <- 0                						
        
        for(i in 1:nrow(IDPdb)) {			
                
                ii = ii +1 
                data[[ii]] <- (list( analysis=analysis, 
                                     sampleName=sampleName,
                                     proteinInformation=IDPdb$protein[i], 
                                     title=IDPdb$title[i],
                                     
                                     specificity = IDPdb$specificity[i], 
                                     q.value=IDPdb$Q.Value[i], 
                                     
                                     scans=IDPdb$Group.Source.Spectrum[i], 
                                     scan.time=IDPdb$Scan.Time[i], 
                                     
                                     sequence = IDPdb$Sequence[i], 
                                     peptideSequence=IDPdb$peptideSequence[i], 
                                     n.sequence=IDPdb$n.sequence[i],
                                     table2Sequence=IDPdb$table2Sequence[i],
                                     
                                     exact.mass=IDPdb$Exact.Mass[i],
                                     obs.mass=IDPdb$Observed.Mass[i], 
                                     precursorMZ=IDPdb$Precursor.m.z[i],
                                     charge=IDPdb$charge[i],
                                     
                                     modification=IDPdb$modification[i],
                                     GlycanMass=IDPdb$GlycanMass[i], 
                                     glycoform.mass=IDPdb$glycoform.mass[i],
                                     structure=IDPdb$structure[i],
                                     mZ=MS2Data[[i]]$mZ, 
                                     intensity=MS2Data[[i]]$intensity, 
                                     
                                     MonoisotopicY1mass=IDPdb$MonoisotopicY1mass[i], 
                                     MonoisotopicPeptideMass=IDPdb$MonoisotopicPeptideMass[i],
                                     PepPlus=Ions[[i]]$PepPlus,
                                     Y0NH3=Ions[[i]]$Y0NH3, 
                                     Y1=Ions[[i]]$Y1, 
                                     Y1Mox=Ions[[i]]$Y1Mox,
                                     Y1Ca=Ions[[i]]$Y1Ca,
                                     Y1F=Ions[[i]]$Y1F, 
                                     Y2F=Ions[[i]]$Y2F,
                                     Y3F=Ions[[i]]$Y3F,
                                     Y2=Ions[[i]]$Y2,
                                     Y3=Ions[[i]]$Y3,
                                     Y3X=Ions[[i]]$Y3X,
                                     Y3FX=Ions[[i]]$Y3FX
                ))
                
        }
        
        
        
        return(data)
        
}

#################################################################
#        Function 10. gPSMvalidation 
# validation and spectrum annotation
#################################################################

gPSMvalidator <-
        
        function (data, modification, modificationName, mZmarkerIons, 
                  minNumberIons = 2, itol_ppm = 10, minMarkerIntensityRatio = 2, 
                  PEAKPLOT=TRUE, validate=FALSE) 
        {
                
                query.idx <- 1:length(data)
                query.to.scan <- as.integer(as.character(lapply(data, 
                                                                function(x) {
                        if (length(x$scans) == 1) {
                                return(x$scans)
                        } else {
                                return(x$scans[1])
                        }
                })))
                scan.to.query <- rep(-1, max(query.to.scan, na.rm = TRUE))
                scan.to.query[query.to.scan[query.idx]] <- query.idx
                
                rr <- numeric()
                for (i in 1:length(data)) {
                        
        idx <- findNN(mZmarkerIons, data[[i]]$mZ)
        ppm.error <- 1e-06 * itol_ppm * data[[i]]$mZ[idx]
        b <- (abs(mZmarkerIons - data[[i]]$mZ[idx]) < ppm.error)
        sum.mZmarkerIons.intensity <- sum(data[[i]]$intensity[idx[b]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZmarkerIons.intensity
        percent.mZmarkerIons <- round(100 * (sum.mZmarkerIons.intensity/
                        (sum.mZmarkerIons.intensity + sum.intensity)), 1)
        idx.ppm.error <- (mZmarkerIons[b]-data[[i]]$mZ[idx][b])/mZmarkerIons[b]*1e06
        idx.mZ <- data[[i]]$mZ[idx][b]
                             
        idY1 <- findNN(data[[i]]$Y1 , data[[i]]$mZ)
        ppm.errorY1 <- 1e-06 * 20 * data[[i]]$mZ[idY1]
        c <- (abs(data[[i]]$Y1  - data[[i]]$mZ[idY1]) < ppm.errorY1)
        sum.mZY1Ions.intensity <- sum(data[[i]]$intensity[idY1[c]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY1Ions.intensity
        percent.mZY1Ions <- round(100 *(sum.mZY1Ions.intensity/
                        (sum.mZY1Ions.intensity +  sum.intensity)), 1)
        idY1.ppm.error <- (data[[i]]$Y1[c]-data[[i]]$mZ[idY1][c])/data[[i]]$Y1[c]*1e06
        idY1.mZ <- data[[i]]$mZ[idY1][c]
        
        idY1Mox <- findNN(data[[i]]$Y1Mox , data[[i]]$mZ)
        ppm.errorY1Mox <- 1e-06 * 20 * data[[i]]$mZ[idY1Mox]
        bb <- (abs(data[[i]]$Y1Mox  - data[[i]]$mZ[idY1Mox]) < ppm.errorY1Mox)
        sum.mZY1MoxIons.intensity <- sum(data[[i]]$intensity[idY1Mox[bb]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY1MoxIons.intensity
        percent.mZY1MoxIons <- round(100 *(sum.mZY1MoxIons.intensity/
                                                   (sum.mZY1MoxIons.intensity +  sum.intensity)), 1)
        idY1Mox.ppm.error <- (data[[i]]$Y1Mox[bb]-data[[i]]$mZ[idY1Mox][bb])/data[[i]]$Y1Mox[bb]*1e06
        idY1Mox.mZ <- data[[i]]$mZ[idY1Mox][bb]
        
        
        idY1Ca <- findNN(data[[i]]$Y1Ca , data[[i]]$mZ)
        ppm.errorY1Ca <- 1e-06 * 20 * data[[i]]$mZ[idY1Ca]
        ee <- (abs(data[[i]]$Y1Ca  - data[[i]]$mZ[idY1Ca]) < ppm.errorY1Ca)
        sum.mZY1CaIons.intensity <- sum(data[[i]]$intensity[idY1Ca[ee]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY1CaIons.intensity
        percent.mZY1CaIons <- round(100 *(sum.mZY1CaIons.intensity/
                                                  (sum.mZY1CaIons.intensity +  sum.intensity)), 1)
        idY1Ca.ppm.error <- (data[[i]]$Y1Ca[ee]-data[[i]]$mZ[idY1Ca][ee])/data[[i]]$Y1Ca[ee]*1e06
        idY1Ca.mZ <- data[[i]]$mZ[idY1Ca][ee]
        
        idY2 <- findNN(data[[i]]$Y2 , data[[i]]$mZ)
        ppm.errorY2 <- 1e-06 * 20 * data[[i]]$mZ[idY2]
        dd <- (abs(data[[i]]$Y2  - data[[i]]$mZ[idY2]) < ppm.errorY2)
        sum.mZY2Ions.intensity <- sum(data[[i]]$intensity[idY2[dd]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY2Ions.intensity
        percent.mZY2Ions <- round(100 *(sum.mZY2Ions.intensity/
                                (sum.mZY2Ions.intensity +  sum.intensity)), 1)
        idY2.ppm.error <- (data[[i]]$Y2[dd]-data[[i]]$mZ[idY2][dd])/data[[i]]$Y2[dd]*1e06
        idY2.mZ <- data[[i]]$mZ[idY2][dd]
        
        idY1F <- findNN(data[[i]]$Y1F , data[[i]]$mZ)
        ppm.errorY1F <- 1e-06 * 20 * data[[i]]$mZ[idY1F]
        ff <- (abs(data[[i]]$Y1F  - data[[i]]$mZ[idY1F]) < ppm.errorY1F)
        sum.mZY1FIons.intensity <- sum(data[[i]]$intensity[idY1F[ff]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY1FIons.intensity
        percent.mZY1FIons <- round(100 *(sum.mZY1FIons.intensity/
                                                (sum.mZY1FIons.intensity +  sum.intensity)), 1)
        idY1F.ppm.error <- (data[[i]]$Y1F[ff]-data[[i]]$mZ[idY1F][ff])/data[[i]]$Y1F[ff]*1e06
        idY1F.mZ <- data[[i]]$mZ[idY1F][ff]
        
        idY2F <- findNN(data[[i]]$Y2F , data[[i]]$mZ)
        ppm.errorY2F <- 1e-06 * 20 * data[[i]]$mZ[idY2F]
        gg <- (abs(data[[i]]$Y2F  - data[[i]]$mZ[idY2F]) < ppm.errorY2F)
        sum.mZY2FIons.intensity <- sum(data[[i]]$intensity[idY2F[gg]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY2FIons.intensity
        percent.mZY2FIons <- round(100 *(sum.mZY2FIons.intensity/
                                                 (sum.mZY2FIons.intensity +  sum.intensity)), 1)
        idY2F.ppm.error <- (data[[i]]$Y2F[gg]-data[[i]]$mZ[idY2F][gg])/data[[i]]$Y2F[gg]*1e06
        idY2F.mZ <- data[[i]]$mZ[idY2F][gg]
        
        idY3F <- findNN(data[[i]]$Y3F , data[[i]]$mZ)
        ppm.errorY3F <- 1e-06 * 20 * data[[i]]$mZ[idY3F]
        kk <- (abs(data[[i]]$Y3F  - data[[i]]$mZ[idY3F]) < ppm.errorY3F)
        sum.mZY3FIons.intensity <- sum(data[[i]]$intensity[idY3F[kk]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY3FIons.intensity
        percent.mZY3FIons <- round(100 *(sum.mZY3FIons.intensity/
                                                 (sum.mZY3FIons.intensity +  sum.intensity)), 1)
        idY3F.ppm.error <- (data[[i]]$Y3F[kk]-data[[i]]$mZ[idY3F][kk])/data[[i]]$Y3F[kk]*1e06
        idY3F.mZ <- data[[i]]$mZ[idY3F][kk]
        
        idY3FX <- findNN(data[[i]]$Y3FX , data[[i]]$mZ)
        ppm.errorY3FX <- 1e-06 * 20 * data[[i]]$mZ[idY3FX]
        ll <- (abs(data[[i]]$Y3FX  - data[[i]]$mZ[idY3FX]) < ppm.errorY3FX)
        sum.mZY3FXIons.intensity <- sum(data[[i]]$intensity[idY3FX[ll]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY3FXIons.intensity
        percent.mZY3FXIons <- round(100 *(sum.mZY3FXIons.intensity/
                                                  (sum.mZY3FXIons.intensity +  sum.intensity)), 1)
        idY3FX.ppm.error <- (data[[i]]$Y3FX[ll]-data[[i]]$mZ[idY3FX][ll])/data[[i]]$Y3FX[ll]*1e06
        idY3FX.mZ <- data[[i]]$mZ[idY3FX][ll]
        
        idY3 <- findNN(data[[i]]$Y3 , data[[i]]$mZ)
        ppm.errorY3 <- 1e-06 * 20 * data[[i]]$mZ[idY3]
        cc <- (abs(data[[i]]$Y3  - data[[i]]$mZ[idY3]) < ppm.errorY3)
        sum.mZY3Ions.intensity <- sum(data[[i]]$intensity[idY3[cc]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY3Ions.intensity
        percent.mZY3Ions <- round(100 *(sum.mZY3Ions.intensity/
                                (sum.mZY3Ions.intensity +  sum.intensity)), 1)
        idY3.ppm.error <- (data[[i]]$Y3[cc]-data[[i]]$mZ[idY3][cc])/data[[i]]$Y3[cc]*1e06
        idY3.mZ <- data[[i]]$mZ[idY3][cc]
        
        idY3X <- findNN(data[[i]]$Y3X , data[[i]]$mZ)
        ppm.errorY3X <- 1e-06 * 20 * data[[i]]$mZ[idY3X]
        hh <- (abs(data[[i]]$Y3X  - data[[i]]$mZ[idY3X]) < ppm.errorY3X)
        sum.mZY3XIons.intensity <- sum(data[[i]]$intensity[idY3X[hh]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY3XIons.intensity
        percent.mZY3XIons <- round(100 *(sum.mZY3XIons.intensity/
                                                 (sum.mZY3XIons.intensity +  sum.intensity)), 1)
        idY3X.ppm.error <- (data[[i]]$Y3X[hh]-data[[i]]$mZ[idY3X][hh])/data[[i]]$Y3X[hh]*1e06
        idY3X.mZ <- data[[i]]$mZ[idY3X][hh]
                      
        idPepPlus <- findNN(data[[i]]$PepPlus , data[[i]]$mZ)
        ppm.errorPepPlus <- 1e-06 * 20 * data[[i]]$mZ[idPepPlus]
        f <- (abs(data[[i]]$PepPlus  - data[[i]]$mZ[idPepPlus]) < 
                      ppm.errorPepPlus)
        sum.mZPepPlusIons.intensity <- sum(data[[i]]$intensity[idPepPlus[f]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZPepPlusIons.intensity
        percent.mZPepPlusIons <- round(100 * (sum.mZPepPlusIons.intensity/
                        (sum.mZPepPlusIons.intensity +  sum.intensity)), 1)
        idPepPlus.ppm.error <- (data[[i]]$PepPlus[f]-data[[i]]$mZ[idPepPlus][f])/data[[i]]$PepPlus[f]*1e06
        idPepPlus.mZ <- data[[i]]$mZ[idPepPlus][f]
         
        idY0NH3 <- findNN(data[[i]]$Y0NH3 , data[[i]]$mZ)
        ppm.errorY0NH3 <- 1e-06 * 20 * data[[i]]$mZ[idY0NH3]
        mm <- (abs(data[[i]]$Y0NH3  - data[[i]]$mZ[idY0NH3]) < ppm.errorY0NH3)
        sum.mZY0NH3Ions.intensity <- sum(data[[i]]$intensity[idY0NH3[mm]])
        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY0NH3Ions.intensity
        percent.mZY0NH3Ions <- round(100 *(sum.mZY0NH3Ions.intensity/
                                                   (sum.mZY0NH3Ions.intensity +  sum.intensity)), 1)
        idY0NH3.ppm.error <- (data[[i]]$Y0NH3[mm]-data[[i]]$mZ[idY0NH3][mm])/data[[i]]$Y0NH3[mm]*1e06
        idY0NH3.mZ <- data[[i]]$mZ[idY0NH3][mm]
                       
        ido <- findNN(otherOxonium, data[[i]]$mZ)
        ppm.errorOO <- 1e-06 * itol_ppm * data[[i]]$mZ[ido]
        e <- (abs(otherOxonium - data[[i]]$mZ[ido]) < ppm.errorOO)
        sum.otherOxonium.intensity <- sum(data[[i]]$intensity[ido[e]])
        ido.ppm.error <- (otherOxonium[e]-data[[i]]$mZ[ido][e])/otherOxonium[e]*1e06
         ido.mZ <-   data[[i]]$mZ[ido][e]     
        
   if(validate)  {   
                 
if ((length(data[[i]]$mZ[idx[b]]) >= minNumberIons) & 
        percent.mZmarkerIons > minMarkerIntensityRatio 
        & (length(data[[i]]$mZ[idY1[c]]) >= 1) &
           (length(data[[i]]$mZ[idPepPlus[f]]) >= 1))

        
        {
                                
        r <- cbind(scans = data[[i]]$scans, 
                   query = i, 
                   obs.mass = data[[i]]$obs.mass,
                   GlycanMass = data[[i]]$GlycanMass, 
                   glycoform.mass = data[[i]]$glycoform.mass,
                   MonoisotopicY1mass = data[[i]]$MonoisotopicY1mass,
                   MonoisotopicPeptideMass = data[[i]]$MonoisotopicPeptideMass,
                   pepmass_unmod = data[[i]]$pepmass_unmod,
                   pepmass_mod = data[[i]]$pepmass_mod,
                   scan.time = data[[i]]$scan.time,
                   structure = data[[i]]$structure,
                   precursorMZ = data[[i]]$precursorMZ,
                   exact.mass = data[[i]]$exact.mass,
                   precursorMassErrorPpm = data[[i]]$precursorMassErrorPpm,
                   charge=data[[i]]$charge, 
                   peptideSequence=data[[i]]$peptideSequence,
                   sequence=data[[i]]$sequence,
                   specificity=data[[i]]$specificity,
                   analysis=data[[i]]$analysis,
                   q.value=data[[i]]$q.value,
                   table2Sequence=data[[i]]$table2Sequence,
                   n.sequence=data[[i]]$n.sequence,
                   title=data[[i]]$title, modification=data[[i]]$modification)
                                
        rr <- rbind(rr, r)
        
        def.par <- par(no.readonly = TRUE) # save default, for resetting...
        
        layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
               widths=c(1,1), heights=c(1.25,1))
        
        
        if (!is.na(data[[i]]$q.value)) {
                
            
                
        fi <- fragmentIon(sequence = data[[i]]$peptideSequence, 
                           FUN = defaultIons, 
                           modified = substr(data[[i]]$modification, 2, 
                                             nchar(data[[i]]$modification) - 1), 
                           modification = modification)
        
        fi.by <- as.data.frame(cbind(b = fi[[1]]$b, y = fi[[1]]$y))
                           
       check<- peakplot(data[[i]]$peptideSequence, spec = data[[i]], 
                           fi = fi.by, ion.axes = F, 
                           main = list(paste(data[[i]]$sampleName, 
                                             data[[i]]$sequence, 
                                             data[[i]]$glycoform.mass, 
                                             data[[i]]$structure, sep=" : "), 
                                       cex = 1),
                           xlim = c(0, max(data[[i]]$mZ)))
      aa.mZ<- data[[i]]$mZ[check$idx][abs(check$mZ.Da.error)<=0.02]
      sum.aa.intensity <- sum(data[[i]]$intensity[check$idx]) 
      aa.Da.error <- check$mZ.Da.error[abs(check$mZ.Da.error)<=0.02]
      aa.ppm.error <- aa.Da.error/ aa.mZ * 1000000
                              }
         else {
         par(cex = 1)
         plot(data[[i]]$mZ, data[[i]]$intensity, type = "h", 
             xlab = "m/z", ylab = "Intensity", 
             main = list(paste(data[[i]]$sampleName, 
                               data[[i]]$sequence, sep=" : "), cex = 1), 
             xlim = c(0, max(data[[i]]$mZ)))
                                }
                   
        
          points(data[[i]]$mZ[idx[b]], data[[i]]$intensity[idx[b]],
                 pch = 22, col = "green", bg = "green",cex = 0.75)
                        
          points(data[[i]]$mZ[idY1[c]], data[[i]]$intensity[idY1[c]], 
                                       pch = 22, col = "red", bg = "red", 
                                       cex = 0.75)
     points(data[[i]]$mZ[idY1Mox[bb]], data[[i]]$intensity[idY1Mox[bb]], 
           pch = 22, col = "blue", bg = "blue", 
            cex = 0.75)
   #  points(data[[i]]$mZ[idY1Ca[ee]], data[[i]]$intensity[idY1Ca[ee]], 
    #         pch = 22, col = "lightblue", bg = "lightblue", 
     #        cex = 0.75)
        points(data[[i]]$mZ[idY2[dd]], data[[i]]$intensity[idY2[dd]], 
               pch = 22, col = "red", bg = "red", 
               cex = 0.75)
      points(data[[i]]$mZ[idY1F[ff]], data[[i]]$intensity[idY1F[ff]], 
             pch = 22, col = "purple", bg = "purple", 
             cex = 0.75)
      points(data[[i]]$mZ[idY2F[gg]], data[[i]]$intensity[idY2F[gg]], 
             pch = 22, col = "purple", bg = "purple", 
             cex = 0.75)
      points(data[[i]]$mZ[idY3F[kk]], data[[i]]$intensity[idY3F[kk]], 
             pch = 22, col = "purple", bg = "purple", 
             cex = 0.75)
      points(data[[i]]$mZ[idY3FX[ll]], data[[i]]$intensity[idY3FX[ll]], 
             pch = 22, col = "orange", bg = "orange", 
             cex = 0.75)
        points(data[[i]]$mZ[idY3[cc]], data[[i]]$intensity[idY3[cc]], 
              pch = 22, col = "red", bg = "red", 
              cex = 0.75)
      points(data[[i]]$mZ[idY3X[hh]], data[[i]]$intensity[idY3X[hh]], 
             pch = 22, col = "orange", bg = "orange", 
             cex = 0.75)
          points(data[[i]]$mZ[idPepPlus[f]], data[[i]]$intensity[idPepPlus[f]], 
                                       pch = 22, col = "pink", bg = "pink", 
                                       cex = 0.75)
    points(data[[i]]$mZ[idY0NH3[mm]], data[[i]]$intensity[idY0NH3[mm]], 
           pch = 22, col = "pink4", bg = "pink4", 
           cex = 0.75)
   
          points(data[[i]]$mZ[ido[e]], data[[i]]$intensity[ido[e]], 
                                       pch = 22, col = "black",bg = "black", 
                                       cex = 0.75)
    sum.ID.intensity <- sum(
           # sum.mZmarkerIons.intensity, 
            sum.mZY1Ions.intensity, 
                            sum.mZY1MoxIons.intensity, 
            #sum.mZY1CaIons.intensity, 
                            sum.mZY2Ions.intensity, sum.mZY1FIons.intensity,
                            sum.mZY2FIons.intensity, sum.mZY3FIons.intensity, sum.mZY3FXIons.intensity,
                            sum.mZY3Ions.intensity, sum.mZY3XIons.intensity,sum.mZPepPlusIons.intensity, 
                            sum.mZY0NH3Ions.intensity, sum.otherOxonium.intensity,
                            sum.aa.intensity)
    percent.ID <- round(100*sum.ID.intensity/sum(data[[i]]$intensity))
    #sum.intensity <- sum(data[[i]]$intensity) - sum.ID.intensity
   # percent.ID <- round(100 *(sum.ID.intensity/
       #                               (sum.ID.intensity +  sum.intensity)), 1)

          
          legend("topright", paste(c( "m/z", "charge", "scan",
                                      "query"
                                      #, "% ID intensity"
                                      ),
          c(round(data[[i]]$precursorMZ,3), data[[i]]$charge,  
          data[[i]]$scans, i
          #, percent.ID
          )),cex = 1)
        
  
        
        
        ppm.error <- c(idx.ppm.error, idY1.ppm.error, 
                       idY1Mox.ppm.error, 
                       #idY1Ca.ppm.error, 
                       idY2.ppm.error, idY1F.ppm.error,
                       idY2F.ppm.error, idY3F.ppm.error, idY3FX.ppm.error,
                       idY3.ppm.error, idY3X.ppm.error,idPepPlus.ppm.error, idY0NH3.ppm.error, ido.ppm.error,
                       aa.ppm.error)
        mZ <- c(idx.mZ,  idY1.mZ, 
                idY1Mox.mZ,
                #idY1Ca.mZ, 
                idY2.mZ,idY1F.mZ, idY2F.mZ,idY3F.mZ,idY3FX.mZ,
                idY3.mZ, idY3X.mZ, idPepPlus.mZ, idY0NH3.mZ, ido.mZ,
                aa.mZ)
        plot(mZ, ppm.error,
             main="Error Plot",
             pch='o',
             cex=0.5, 
             ylim = c(-10* itol_ppm, 10* itol_ppm)
             )
      
     
       
       plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
       
       
       legend("center", paste(c( "HexNAc+", "Oxonium", 
                                                  "Y0", "Y0-NH3","Y1,2,3", "YF1,2,3", "YX3,F"
                                 , "Y1-Mox"
                                 #, "Y1-Ca"
                                 )),
              fill=c( "green", 
                      "black", 
                      "pink", "pink4", "red", "purple", "orange"
                      , "blue"
                      #, "lightblue"
                      ), bty="n",cex=1, ncol=2)
       
       
       
       
       par(def.par)  #- reset to default
                                
                        
                                         
                        }}
else {
        r <- cbind(scans = data[[i]]$scans, 
                   query = i, 
                   obs.mass = data[[i]]$obs.mass,
                   GlycanMass = data[[i]]$GlycanMass, 
                   glycoform.mass = data[[i]]$glycoform.mass,
                   MonoisotopicY1mass = data[[i]]$MonoisotopicY1mass,
                   MonoisotopicPeptideMass = data[[i]]$MonoisotopicPeptideMass,
                   pepmass_unmod = data[[i]]$pepmass_unmod,
                   pepmass_mod = data[[i]]$pepmass_mod,
                   scan.time = data[[i]]$scan.time,
                   structure = data[[i]]$structure,
                   precursorMZ = data[[i]]$precursorMZ,
                   exact.mass = data[[i]]$exact.mass,
                   precursorMassErrorPpm = data[[i]]$precursorMassErrorPpm,
                   charge=data[[i]]$charge, 
                   peptideSequence=data[[i]]$peptideSequence,
                   sequence=data[[i]]$sequence,
                   specificity=data[[i]]$specificity,
                   analysis=data[[i]]$analysis,
                   q.value=data[[i]]$q.value,
                   table2Sequence=data[[i]]$table2Sequence,
                   n.sequence=data[[i]]$n.sequence,
                   title=data[[i]]$title, modification=data[[i]]$modification)
        
        rr <- rbind(rr, r)
        
        def.par <- par(no.readonly = TRUE) # save default, for resetting...
        
        layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
               widths=c(1,1), heights=c(1.25,1))
        
        
        if (!is.na(data[[i]]$q.value)) {
                
                
                
                fi <- fragmentIon(sequence = data[[i]]$peptideSequence, 
                                  FUN = defaultIons, 
                                  modified = substr(data[[i]]$modification, 2, 
                                                    nchar(data[[i]]$modification) - 1), 
                                  modification = modification)
                
                fi.by <- as.data.frame(cbind(b = fi[[1]]$b, y = fi[[1]]$y))
                
                check<- peakplot(data[[i]]$peptideSequence, spec = data[[i]], 
                                 fi = fi.by, ion.axes = F, 
                                 main = list(paste(data[[i]]$sampleName, 
                                                   data[[i]]$sequence, 
                                                   data[[i]]$glycoform.mass, 
                                                   data[[i]]$structure, sep=" : "), 
                                             cex = 1),
                                 xlim = c(0, max(data[[i]]$mZ)))
                aa.mZ<- data[[i]]$mZ[check$idx][abs(check$mZ.Da.error)<=0.02]
                sum.aa.intensity <- sum(data[[i]]$intensity[check$idx]) 
                aa.Da.error <- check$mZ.Da.error[abs(check$mZ.Da.error)<=0.02]
                aa.ppm.error <- aa.Da.error/ aa.mZ * 1000000
        }
        else {
                par(cex = 1)
                plot(data[[i]]$mZ, data[[i]]$intensity, type = "h", 
                     xlab = "m/z", ylab = "Intensity", 
                     main = list(paste(data[[i]]$sampleName, 
                                       data[[i]]$sequence, sep=" : "), cex = 1), 
                     xlim = c(0, max(data[[i]]$mZ)))
        }
        
        
        points(data[[i]]$mZ[idx[b]], data[[i]]$intensity[idx[b]],
               pch = 22, col = "green", bg = "green",cex = 0.75)
        
        points(data[[i]]$mZ[idY1[c]], data[[i]]$intensity[idY1[c]], 
               pch = 22, col = "red", bg = "red", 
               cex = 0.75)
        points(data[[i]]$mZ[idY1Mox[bb]], data[[i]]$intensity[idY1Mox[bb]], 
               pch = 22, col = "blue", bg = "blue", 
               cex = 0.75)
       # points(data[[i]]$mZ[idY1Ca[ee]], data[[i]]$intensity[idY1Ca[ee]], 
             #  pch = 22, col = "lightblue", bg = "lightblue", 
            #   cex = 0.75)
        points(data[[i]]$mZ[idY2[dd]], data[[i]]$intensity[idY2[dd]], 
               pch = 22, col = "red", bg = "red", 
               cex = 0.75)
        points(data[[i]]$mZ[idY1F[ff]], data[[i]]$intensity[idY1F[ff]], 
               pch = 22, col = "purple", bg = "purple", 
               cex = 0.75)
        points(data[[i]]$mZ[idY2F[gg]], data[[i]]$intensity[idY2F[gg]], 
               pch = 22, col = "purple", bg = "purple", 
               cex = 0.75)
        points(data[[i]]$mZ[idY3F[kk]], data[[i]]$intensity[idY3F[kk]], 
               pch = 22, col = "purple", bg = "purple", 
               cex = 0.75)
        points(data[[i]]$mZ[idY3FX[ll]], data[[i]]$intensity[idY3FX[ll]], 
               pch = 22, col = "orange", bg = "orange", 
               cex = 0.75)
        points(data[[i]]$mZ[idY3[cc]], data[[i]]$intensity[idY3[cc]], 
               pch = 22, col = "red", bg = "red", 
               cex = 0.75)
        points(data[[i]]$mZ[idY3X[hh]], data[[i]]$intensity[idY3X[hh]], 
               pch = 22, col = "orange", bg = "orange", 
               cex = 0.75)
        points(data[[i]]$mZ[idPepPlus[f]], data[[i]]$intensity[idPepPlus[f]], 
               pch = 22, col = "pink", bg = "pink", 
               cex = 0.75)
        points(data[[i]]$mZ[idY0NH3[mm]], data[[i]]$intensity[idY0NH3[mm]], 
               pch = 22, col = "pink4", bg = "pink4", 
               cex = 0.75)
        
        points(data[[i]]$mZ[ido[e]], data[[i]]$intensity[ido[e]], 
               pch = 22, col = "black",bg = "black", 
               cex = 0.75)
        sum.ID.intensity <- sum(
                # sum.mZmarkerIons.intensity, 
                sum.mZY1Ions.intensity, 
                sum.mZY1MoxIons.intensity, 
                #sum.mZY1CaIons.intensity, 
                sum.mZY2Ions.intensity, sum.mZY1FIons.intensity,
                sum.mZY2FIons.intensity, sum.mZY3FIons.intensity, sum.mZY3FXIons.intensity,
                sum.mZY3Ions.intensity, sum.mZY3XIons.intensity,sum.mZPepPlusIons.intensity, 
                sum.mZY0NH3Ions.intensity, sum.otherOxonium.intensity,
                sum.aa.intensity)
        percent.ID <- round(100*sum.ID.intensity/sum(data[[i]]$intensity))
        #sum.intensity <- sum(data[[i]]$intensity) - sum.ID.intensity
        # percent.ID <- round(100 *(sum.ID.intensity/
        #                               (sum.ID.intensity +  sum.intensity)), 1)
        
        
        legend("topright", paste(c( "m/z", "charge", "scan",
                                    "query"
                                    #, "% ID intensity"
        ),
        c(round(data[[i]]$precursorMZ,3), data[[i]]$charge,  
          data[[i]]$scans, i
          #, percent.ID
        )),cex = 1)
        
        
        
        
        ppm.error <- c(idx.ppm.error, idY1.ppm.error, 
                       idY1Mox.ppm.error, 
                       #idY1Ca.ppm.error, 
                       idY2.ppm.error, idY1F.ppm.error,
                       idY2F.ppm.error, idY3F.ppm.error, idY3FX.ppm.error,
                       idY3.ppm.error, idY3X.ppm.error,idPepPlus.ppm.error, idY0NH3.ppm.error, ido.ppm.error,
                       aa.ppm.error)
        mZ <- c(idx.mZ,  idY1.mZ, 
                idY1Mox.mZ, 
                #idY1Ca.mZ, 
                idY2.mZ,idY1F.mZ, idY2F.mZ,idY3F.mZ,idY3FX.mZ,
                idY3.mZ, idY3X.mZ, idPepPlus.mZ, idY0NH3.mZ, ido.mZ,
                aa.mZ)
        plot(mZ, ppm.error,
             main="Error Plot",
             pch='o',
             cex=0.5, 
             ylim = c(-10* itol_ppm, 10* itol_ppm)
        )
        
        
        
        plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
        
        
        legend("center", paste(c( "HexNAc+", "Oxonium", 
                                  "Y0", "Y0-NH3","Y1,2,3", "YF1,2,3", "YX3,F"
                                  , "Y1-Mox"
                                  #, "Y1-Ca"
        )),
        fill=c( "green", 
                "black", 
                "pink", "pink4", "red", "purple", "orange"
                , "blue"
                #, "lightblue"
        ), bty="n",cex=1, ncol=2)
        
        
        
        
        par(def.par)  #- reset to default
        
}
                }
                close.screen(all.screens = TRUE)
                return(as.data.frame(rr, stringsAsFactors=FALSE))
        }


##############################################################################
#  Function 14. #default ions
##############################################################################

defaultIons <-
        function (b, y)
        {
                Hydrogen <- 1.007825
                Oxygen <- 15.994915
                Nitrogen <- 14.003074
                c <- b + (Nitrogen + (3 * Hydrogen))
                z <- y - (Nitrogen + (3 * Hydrogen))
                return(cbind(b, y, c, z))
        }


##############################################################################
#   Function 19. Read.RQ. Read in quantitation information 
# and associate it with the SIC peaks file
###############################################################################

Read.RQ <- function (input="Output/quantitation_Chym1_ELUTE_211.csv",
                     dir="quantitation_Chym1_ELUTE") {
      
        
RQ <- read.csv(file=input, header=TRUE)
        
# Make a data frame 'SIC.all' with columns 'title' and 'exact.precursor.mz'; 
# it is a list of the SIC files and the exact mz
        
title <- list.files(dir)
peaks <- grep("peaks", title)
        
data <- title[peaks]
title <- title[peaks]
        
# change: chr "ljz_20131022_MR_Chym1_HILIC_MS2.mzML.binary.sn1830.txt" 
#to: num  10011
        
data <- strsplit(data, split="sic.")
data <- sapply(data, function(x) x[2])
data <- strsplit(data, split=".peaks")
data <- sapply(data, function(x) x[1])
exact.precursor.mz <- as.numeric(data)
        
# this is the title and scan number. 
#data will be filled in after merging with IDPdb
SIC.all <- data.frame(exact.precursor.mz, title) 
RQ <- merge(  x = SIC.all, y = RQ,  by = "exact.precursor.mz", all.y = TRUE)
RQ$title <- as.character(RQ$title)
return(RQ)
}

##############################################################################
#   Function. 20. Read.SICData   Read in SIC data 
##############################################################################


Read.SICData <- function (dir="quantitation_Chym1_ELUTE") {
        
read.dat <- function(file="ljz_20131022_MR_Chym1_ELUTE.mzML.sic.
                     519.5880.peaks.csv") 
        {                		
                
                dat <- read.csv(paste(dir, "/", file, sep=""), skip=1)
                dat <- dat[c("rt", "sumIntensity","peakMZ")]
                dat$rt.min <- dat$rt/60
                
                tt <- dat$peakMZ
                tt <- as.character(tt)
                tt <- strsplit(tt, split=".", fixed=T)
                tag <- sapply(tt, function(x) x[1])
                dat$tag <- as.numeric(tag)
                
                return(dat)
        }
        
        # Make a list 'MS2Data' containing all MS2 binary data
# these are all of the titles matching ids in RQ        
title.identified <- RQ$title							
        
        SICData <- vector( mode="list", length=length(title.identified))
        SICData <- lapply( title.identified, read.dat )
        names(SICData) <- title.identified
        
        
        return(SICData)
}





###############################################################################
#   #peakplot
##############################################################################

peakplot <- 
        function (peptideSequence, spec, FUN = defaultIons, 
                  fi = fragmentIon(peptideSequence, FUN = FUN)[[1]], 
                  main = NULL, 
                  sub = NULL, 
                  xlim = range(spec$mZ, na.rm = TRUE), 
                  ylim = range(spec$intensity, na.rm = TRUE), 
                  itol = 0.02, pattern.abc = "[abc].*", 
                  pattern.xyz = "[xyz].*", ion.axes = TRUE) 
        {
                n <- nchar(peptideSequence)
                m <- psm(peptideSequence, spec, FUN, fi = fi, plot = FALSE)
                max.intensity <- max(spec$intensity, na.rm = TRUE)
                yMax <- 1 * max.intensity
                
                
                op <- par(mar = c(5, 5, 5, 5), cex = 0.75)
                
                
                plot(spec$mZ, spec$intensity, xlab = "m/z", ylab = "Intensity", 
                     type = "h", main = main, xlim = c(min(spec$mZ), max(spec$mZ)), 
                     ylim = c(0, 1.2 * yMax), sub = sub, axes = "T")
                
                LABEL.abc <- (abs(m$mZ.Da.error) < itol) & (regexpr(pattern.abc, m$label) > 0)
                LABEL.xyz <- (abs(m$mZ.Da.error) < itol) & (regexpr(pattern.xyz, m$label) > 0)
                if (length(m$idx[LABEL.abc]) > 0) {
                for (i in 1:length(spec$mZ[m$idx[LABEL.abc]])) {
                        lines(spec$mZ[m$idx[LABEL.abc]][i], spec$intensity[m$idx[LABEL.abc]][i], 
                             type="h", col="purple")
                }
                        text(spec$mZ[m$idx[LABEL.abc]][i], spec$intensity[m$idx[LABEL.abc]][i],
                            m$label[LABEL.abc][i], pos=3, col="purple")
                        
                }
                
                if (length(m$idx[LABEL.xyz]) > 0) {
                for (i in 1:length(spec$mZ[m$idx[LABEL.xyz]])) {
                        lines(spec$mZ[m$idx[LABEL.xyz]][i], spec$intensity[m$idx[LABEL.xyz]][i], 
                              type="h", col="blue")
                        text(spec$mZ[m$idx[LABEL.xyz]][i], spec$intensity[m$idx[LABEL.xyz]][i],
                             m$label[LABEL.xyz][i], pos=3, col="blue")
                        
                }
                }
                
                if (ion.axes) {
                        if (length(m$idx[LABEL.abc]) > 0) {
                                axis(1, spec$mZ[m$idx[LABEL.abc]], m$label[LABEL.abc], 
                                     las = 2)
                        }
                        axis(2)
                        if (length(m$idx[LABEL.xyz]) > 0) {
                                axis(3, spec$mZ[m$idx[LABEL.xyz]], m$label[LABEL.xyz], 
                                     col.axis = "blue", las = 2)
                        }       
                       
                }
               
                par(op)
                return(m)
        }

################################################################################
#function: findUnoccupiedAndGlcNAc
#  find the gPSMs that are unoccupied and that contain GlcNAc
#  gPSMvalidator throws out unoccupied 
#  gPSMvalidator often throws out GlcNAc only because don;t always contain the 
# Y1 and Y0 ions.
################################################################################
findUnoccupiedAndGlcNAc <- function(IDPdb) {
        
        
        ###################  make a data.frame of IDs containing GlcNAc only (+203) ####
        
        GlcNAc.glycoforms <- subset(IDPdb, GlycanMass == 203)
       
        ### Make a data.frame of IDs that do not contain a glycan modification ########
        IDPdb.unoccupied <- subset(IDPdb, GlycanMass == 0)
        
        
 
        # rbind IDPdb.unoccupied with GlcNAc.glycoforms
        
        unoccupiedAndGlcNAc <- rbind(IDPdb.unoccupied, GlcNAc.glycoforms)
        
        
        ###################  #######################################
        unoccupiedAndGlcNAc$query <- NA
        unoccupiedAndGlcNAc$search <- "unoccupiedAndGlcNAc"
        unoccupiedAndGlcNAc$GlycanMass[is.na(unoccupiedAndGlcNAc$GlycanMass)] <- 0
        
        
        return(unoccupiedAndGlcNAc)
        
        
}

######################################################################################
# function: glycosylationProfileTables
######################################################################################
#problems: Table 3 has hardcoded: Asn, site, sequon

glycosylationProfileTables <- function(dat) {
        
        t2S <- unique(dat$table2Sequence)
        dat$glyco <- as.integer(dat$GlycanMass)
        
        ##  How many times was that peptide observed?  ###########################
        count <- sapply(1:length(t2S), function(i) {
                length(dat$table2Sequence[dat$table2Sequence == t2S[i]])
        })
        
        ##  When did it start eluting? ##############################
        min.rt <- sapply(1:length(t2S), function(j) {
                min(dat$scan.time[dat$table2Sequence == t2S[j]])
        })        
        
        ##  When did it stop eluting?  ##############################
        max.rt <- sapply(1:length(t2S), function(k) {
                max(dat$scan.time[dat$table2Sequence == t2S[k]])
        })        
        
        ####  What is the peptide sequence associated with each peptide mass?  ############################
        peptideMass <- sapply(1:length(t2S), function(l) {
                unique(dat$MonoisotopicPeptideMass[dat$table2Sequence == t2S[l]])
        }) 
        
        n.sequence <- sapply(1:length(t2S), function(l) {
                unique(dat$n.sequence[dat$table2Sequence == t2S[l]])
        }) 
        
        ####  What is the Asn number?  ############################
        Asn <- sapply(1:length(t2S), function(n) {
                unique(dat$Asn[dat$table2Sequence == t2S[n]])
        }) 
        
        
        ##  what is the smallest glycan? ##############################
        glycan.S <- sapply(1:length(t2S), function(m) {
                min(dat$glyco[dat$table2Sequence == t2S[m]])
        }) 
        
        ##  what is the largest glycan?  ##############################
        glycan.L <- sapply(1:length(t2S), function(k) {
                max(dat$glyco[dat$table2Sequence == t2S[k]])
        })
        
        ####
        
        Table2 <- data.frame(Asn, peptideMass, n.sequence, t2S, count, min.rt, max.rt, glycan.L,
                             glycan.S, stringsAsFactors=F)
        
        oo <- order(Table2$Asn)
        Table2 <- Table2[oo,]
        
        Table2 <- data.frame( Table2$Asn, Table2$peptideMass,
                              Table2$n.sequence, Table2$t2S, 
                              Table2$count, Table2$min.rt, Table2$max.rt, 
                              Table2$glycan.S, Table2$glycan.L)
        
        names(Table2) <- c("Asn", "peptideMass", "n.sequence", 
                           "table2Sequence", "count", "min.rt",
                           "max.rt", "glycan.S", "glycan.L")
        
        
        ##########################################################
        # Table3: Heterogeneity of glycosylation
        ##########################################################
        
        dat$Asn <- as.factor(dat$Asn)
        
        site <- 1:13
        Asn <- c(8,28,60,114,127,144,156,185,188,211,256,267,298)
        sequon <- c("NQS", "NNS", "NNT", "NIT", "NVS", "NAT", "NLT.AD", "NFSNTS", "NTS","NST", "NLS", "NLT.AW", "NCS")
        
        ####  How many different glycoforms were observed for each glycosylation site?  ############################
        glycoforms <- sapply(1:length(Asn), function(m) {
                length(unique(dat$GlycanMass[dat$Asn == Asn[m]]))
        }) 
        
        
        ####  How many spectra were observed for that site?  ###########
        count <- as.data.frame(table(dat$Asn))
        names(count) <- c("Asn", "count")
        
        
        ################################
        
        Table3 <- data.frame(site, Asn, glycoforms, sequon)
        
        Table3 <- merge(Table3, count)
        
        Table3 <- data.frame(Table3$site, Table3$Asn, Table3$sequon, Table3$glycoforms,
                             Table3$count)
        names(Table3) <- c("site", "Asn", "sequon", "glycoforms", "count")
        
        dat$exact.precursor.mz <- ((dat$exact.mass + (dat$charge* 1.007276))/dat$charge)
        
        return(list(Table2=Table2, Table3=Table3, dat=dat))
        
}

##############################################################################
#   Function 19. Read.RQ. Read in quantitation information 
# and associate it with the SIC peaks file
###############################################################################

Read.RQ <- function (input="RQ/RQ_Chym1_ELUTE.csv",
                     dir="RQ/") {
        
        
        RQ <- read.csv(file=input, header=TRUE, stringsAsFactors=FALSE)
        RQ$exact.precursor.mz <- round(RQ$exact.precursor.mz, 4)
        
        # Make a data frame 'SIC.all' with columns 'title' and 'exact.precursor.mz'; 
        # it is a list of the SIC files and the exact mz
        
        title <- list.files(dir)
        peaks <- grep("peaks", title)
        
        data <- title[peaks]
        title <- title[peaks]
        
        # change: chr "ljz_20131022_MR_Chym1_ELUTE.mzML.sic.1096.0156.peaks" 
        #to: num  1096.0156
        
        data <- strsplit(data, split="sic.")
        data <- sapply(data, function(x) x[2])
        data <- strsplit(data, split=".peaks")
        data <- sapply(data, function(x) x[1])
        exact.precursor.mz <- as.numeric(data)
        
        # this is the title and scan number. 
        #data will be filled in after merging with IDPdb
        SIC.all <- data.frame(exact.precursor.mz, title) 
        RQ <- merge(  x = SIC.all, y = RQ,  by = "exact.precursor.mz", all.y = TRUE)
        RQ$title <- as.character(RQ$title)
        return(RQ)
}

##############################################################################
#   Function. 20. Read.SICData   Read in SIC data 
##############################################################################


Read.SICData <- function (dir="quantitation_Chym1_ELUTE", Asn=211) {
        
        read.dat <- function(file="ljz_20131022_MR_Chym1_ELUTE.mzML.sic.
                             519.5880.peaks.csv") 
        {                                
                
                dat <- read.csv(paste(dir, "/", file, sep=""), skip=1)
                dat <- dat[c("rt", "sumIntensity","peakMZ")]
                dat$rt.min <- dat$rt/60
                
                tt <- dat$peakMZ
                tt <- as.character(tt)
                tt <- strsplit(tt, split=".", fixed=T)
                tag <- sapply(tt, function(x) x[1])
                dat$tag <- as.numeric(tag)
                
                return(dat)
        }
        
        # Make a list 'MS2Data' containing all MS2 binary data
        # these are all of the titles matching ids in RQ        
        title.identified <- RQ$title[RQ$Asn==Asn]							
        
        SICData <- vector( mode="list", length=length(title.identified))
        SICData <- lapply( title.identified, read.dat )
        names(SICData) <- title.identified
        
        
        return(SICData)
}

##############################################################################
#   Function. 21. RQ of glycoforms
##############################################################################

glycoRQ <- function(RQ, Table2, dir, rt.min.minus, rt.min.plus) {
        
        rr <- numeric()
        
        ii <- 0
        
        for (i in 1:length(unique(RQ$Asn))) {
                
                ii = ii + 1
                
                Asn <- unique(RQ$Asn)[i]
                
                #Read in SIC data from the RQ directory for the glycosylation site of interest
                SICData <- Read.SICData(dir=dir, Asn=Asn)
                
                RQ.Asn <- RQ[RQ$Asn==Asn,]
                
                # subset data calculate stuff and write a data frame
                #Asn.min.rt <- Table2$min.rt[Table2$table2Sequence==unique(RQ.Asn$table2Sequence)]
                #Asn.max.rt <- Table2$max.rt[Table2$table2Sequence==unique(RQ.Asn$table2Sequence)]
                
               Asn.min.rt <- 
                        ((Table2$min.rt[Table2$table2Sequence==unique(RQ.Asn$table2Sequence)]) - rt.min.minus)
                Asn.max.rt <- 
                       ((Table2$min.rt[Table2$table2Sequence==unique(RQ.Asn$table2Sequence)]) + rt.min.plus)
             
                
                subset_SICData <- lapply (SICData, function(x) 
                        subset(x, x$rt.min >= Asn.min.rt & x$rt.min <= Asn.max.rt))
                
                exact.precursor.mz <- RQ.Asn$exact.precursor.mz
                mean_peakMZ <- sapply(subset_SICData, function(x) mean(x$peakMZ))
                sd_peakMZ <- sapply(subset_SICData, function(x) sd(x$peakMZ))
                sum_sumIntensity <-  sapply(subset_SICData, function(x) sum(x$sumIntensity))
                
                
                results <- data.frame(exact.precursor.mz, mean_peakMZ,sd_peakMZ,
                                      sum_sumIntensity)
                RQ.Asn <- merge(x=RQ.Asn, y=results, by="exact.precursor.mz")
                RQ.Asn$relativeRatio <- RQ.Asn$sum_sumIntensity/sum(RQ.Asn$sum_sumIntensity)
                
                
                rr <- rbind(rr, RQ.Asn)
                
        }
        return(as.data.frame(rr))
}

################################################################################
# pGlycoFilter
################################################################################
pGlycoFilter <- function(structure, data=NULL) {
        
        varname <- as.character(substitute(structure))
        if(!is.null(data)) {
                structure = data[[varname]]
        }
        #structure<-eval(substitute(structure),data, parent.frame())
        
        #part 1 rename structures to be consistent
        select <- grepl("(Man)3(GlcNAc)2", structure, fixed=T)
        structure[select] <- gsub("(HexNAc)1","(HexNAc)3",structure[select], fixed=T)
        structure[select] <- gsub("(HexNAc)2","(HexNAc)4",structure[select], fixed=T)
        structure[select] <- gsub("(Hex)6","(Hex)9",structure[select], fixed=T)
        structure[select] <- gsub("(Hex)5","(Hex)8",structure[select], fixed=T)
        structure[select] <- gsub("(Hex)4","(Hex)7",structure[select], fixed=T)
        structure[select] <- gsub("(Hex)3","(Hex)6",structure[select], fixed=T)
        structure[select] <- gsub("(Hex)2","(Hex)5",structure[select], fixed=T)
        structure[select] <- gsub("(Hex)1","(Hex)4",structure[select], fixed=T)
        
        select <- grepl("(Man)3(GlcNAc)2", structure, fixed=T) & !grepl("(Hex)", structure, fixed=T)
        structure[select] <- gsub("(Man)","(Hex)",structure[select], fixed=T)
        
        select <- grepl("(GlcNAc)", structure) & !grepl("(HexNAc)", structure, fixed=T)
        structure[select] <- gsub("(GlcNAc)","(HexNAc)",structure[select], fixed=T)
        
        select <- grepl("+", structure, fixed=T)
        structure[select] <- gsub("(Man)3(GlcNAc)2"," ",structure[select], fixed=T)
        structure[select] <- gsub("+"," ",structure[select], fixed=T)
        structure[select] <- gsub("(GlcNAc)2"," ",structure[select], fixed=T)
        structure[select] <- gsub("(Man)3"," ",structure[select], fixed=T)
        
        
        #part 2: filter structures
        getGroupCount <- function(structure, group) {
                if(grepl(sprintf("(%s)",group), structure,fixed=T)) {
                        pattern <- sprintf("^.*\\(%s\\)([1-9]?).*$",group)
                        res <- sub(pattern, "\\1", structure, perl=T)
                        if(res == "") {
                                return(1)
                        } else {
                                return(as.numeric(res))
                        }
                } else {
                        return(0)
                }
        }
        
        allowed_conditions <- data.frame(matrix(ncol=4, nrow=6))
        colnames(allowed_conditions) <- c("HexNAc", "Hex", "Deoxyhexose", "Pent")
        allowed_conditions[["HexNAc"]] <- c(list(1), list(2), list(2), list(2), list(3), list(4))
        allowed_conditions[["Hex"]] <- c(list(0), list(0), list(1:5), list(6:9), list(2:6), list(3:5))
        allowed_conditions[["Deoxyhexose"]] <- c(list(0:1), list(0:1), list(0:1), list(0), list(0:2), list(0:3))
        allowed_conditions[["Pent"]] <- c(list(0), list(0), list(0:1), list(0), list(0:1), list(0:1))
        struct_mat <- t(sapply(structure, function(s) {
                sapply(colnames(allowed_conditions), function(gr) {
                        getGroupCount(s, gr)
                })
        }))
        
        rallow_list <- rep(FALSE, length(structure))
        for(i in 1:nrow(struct_mat)) {
                for(r in 1:nrow(allowed_conditions)) {
                        rallow <- all(sapply(colnames(allowed_conditions), function(gr) {
                                range = allowed_conditions[r,gr][[1]]
                                struct_mat[i,gr] %in% range
                        }))
                        if(rallow) {
                                rallow_list[i] <- TRUE
                                break
                        }
                }
        }
        
        structure = gsub(" +"," ",structure, perl=T)
        structure = gsub(" $","",structure, perl=T)
        structure = gsub("^ ","",structure, perl=T)
        
        #should the different groups be sorted?
        
        if(is.null(data)) {
                return(structure[rallow_list])
        } else {
                data[[varname]] <- structure
                data <- data[rallow_list,]
                return(data)
        }
}

##############################################################################
#   Function 15. Read in data from glycomod search (IDPdb, in silico digest, 
#       list of MS2 data names
##############################################################################

Read.GlycoMod <- function (input=pGlycoFilter_output, 
                           ChainSaw,
                           spectrum.table=spectrum.table, 
                           dir="MS2Data") {
        
        spectrum.table <- read.csv(spectrum.table, skip=1)
        
        names(input) <- c("glycoform.mass", "mass.error.ppm", "structure", "type", 
                           "peptide.mass",
                           "sequence", "Exact.Mass", "mod", "links" )
        
        input$Exact.Mass <- as.numeric(input$Exact.Mass)
        input$mass.error.ppm <- as.numeric(input$mass.error.ppm)
        
        input$Observed.Mass <- input$Exact.Mass + (input$mass.error.ppm/ 1e+06*input$Exact.Mass)
        input$round.Observed.Mass <- round(input$Observed.Mass, digits=0)
        input$round.Observed.Mass.minus1 <- (input$round.Observed.Mass-1)
        input$round.Observed.Mass.plus1 <- (input$round.Observed.Mass+1)
        spectrum.table$Observed.Mass <- ((spectrum.table$precursorMZ *spectrum.table$charge)-
                                            (spectrum.table$charge * 1.007276))
        spectrum.table$round.Observed.Mass <- round(spectrum.table$Observed.Mass, digits=0)
        
        
        input_round <- merge (x=input, y=spectrum.table, by="round.Observed.Mass", all.x=TRUE)
        
        
        input_round_minus <- merge (x=input, y=spectrum.table, 
                                     by.x="round.Observed.Mass.minus1", by.y="round.Observed.Mass",
                                     all.x=TRUE)
        
        input_round_plus <- merge (x=input, y=spectrum.table, 
                                    by.x="round.Observed.Mass.plus1", by.y="round.Observed.Mass",
                                    all.x=TRUE)
        
        
        glycoModIDs <- rbind (input_round, input_round_minus, input_round_plus)
        
        glycoModIDs$analysis <- "GlycoMod"
        glycoModIDs$Q.Value <- 0
        glycoModIDs$mass.error <- glycoModIDs$Exact.Mass - glycoModIDs$Observed.Mass.y
        glycoModIDs$Scan.Time <- glycoModIDs$rt/60
        
        
        IDPdb <- glycoModIDs
        
        IDPdb <- IDPdb[c("id", "Scan.Time", "Observed.Mass.y", "precursorMZ",
                         "Exact.Mass", "mass.error","analysis", "charge", "Q.Value", 
                         "sequence",
                         "glycoform.mass", "peptide.mass", "structure")]
        
        names(IDPdb) <- c("scans", "Scan.Time", "Observed.Mass","precursorMZ", "Exact.Mass", 
                          "mass.error", 
                          "analysis", "charge", "Q.Value", "sequence", "glycoform.mass", 
                          "peptide.mass", "structure")        
        
        
        
        # make a ppm.mass.error column
        #IDPdb$mass.error/IDPdb$Exact.Mass *10^6
        
        ppm.mass.error <- (IDPdb$mass.error/IDPdb$Exact.Mass *10^6)
        IDPdb$ppm.mass.error <- ppm.mass.error
        IDPdb <- subset(IDPdb, IDPdb$ppm.mass.error <= 15 & IDPdb$ppm.mass.error >= -15)
        
        # change: factor 0.1.4586  to: numeric 4586
        sc <- IDPdb$scans
        sc <- as.character(sc)
        sc <- strsplit(sc, split=".", fixed=T)
        scans <- sapply(sc, function(x) x[3])
        scans <- as.numeric(scans)
        IDPdb$scans <- scans
        
        # Make a data frame 'MS2.all' with columns 'title' and 'scans'; 
        # it contains this info for all MS2 binary files
        
        title <- list.files(dir)
        
        # change: chr "ljz_20131022_MR_Chym1_HILIC_MS2.mzML.binary.sn1830.txt" 
        #to: num  10011
        
        sc <- strsplit(title, split="sn")
        scans <- sapply(sc, function(x) x[2])
        scans <- gsub(pattern=".txt", replacement="", x=scans)
        scans <- as.numeric(scans)
        
        #  title and scan number. data will be filled in after merging with IDPdb
        MS2.all <- data.frame(scans, title) 
        IDPdb <- merge(  x = MS2.all,y = IDPdb,  by = "scans", all.y = TRUE)
        IDPdb$title <- as.character(IDPdb$title)        
        
        
        ###add methionine oxidation as a variable modification
        M <- ChainSaw[grep("M", ChainSaw$sequence), ]
        mox <- gsub( "M", "m", M$sequence ,fixed=TRUE)
        mox <- gsub( "[A-Z]", "0", mox)
        mox <- strsplit(mox , split="")
        mox <- lapply(mox , function(x){gsub("m", "15.994915", x)})
        mox <- lapply(mox , function(x){as.numeric(x)})
        mox <- sapply(mox , function(x) {sum(x)})      
        M$mass <- mapply(sum, M$mass , mox)
        
        # change: chr "QHGFTMMNVYNSTK" to: chr "QHGFTMM16NVYNSTK"
        moxSeq <- M$sequence  
        moxSeq <- as.character(moxSeq)
        moxSeq <- gsub( "M", "M16",  moxSeq ,fixed=TRUE)
        M$sequence <- moxSeq
        
        ChainSaw <- rbind(ChainSaw, M)
        
        
        ChainSaw$peptide.mass <- round(ChainSaw$mass, 3)
        
        
        IDPdb <- merge( x = ChainSaw, y = IDPdb, by = "peptide.mass" )
        
        #Modify the column contents of IDPdb
        
        
        # change: chr "QHGFTMM[16]NVYN[1170]STK" to: chr "0000000200010000"
        modification <- IDPdb$sequence.x
        modification <- gsub( "M16", "2",  modification,fixed=TRUE)
        modification<- gsub( "N[ABCDEFGHIJKLMNOQRSTUVWXYZ]T", "100", modification)
        modification<- gsub( "N[ABCDEFGHIJKLMNOQRSTUVWXYZ]S", "100",modification)
        modification<- gsub( "NF", "10",  modification,fixed=TRUE)
        modification <- gsub( "[A-Z]", "0",modification)
        modification[1:length(modification)] <- 
                paste("0",modification[1:length(modification)],"0", sep="")
        IDPdb$modification <- modification
        
        IDPdb$charge <- as.numeric(IDPdb$charge)
        
        peptideSequence <- IDPdb$sequence.x
        peptideSequence <- gsub("M16", "M", peptideSequence)

        IDPdb$peptideSequence <- peptideSequence
        
        IDPdb$table2Sequence <- IDPdb$sequence.x
        
        # change: chr "Q[-17]HGFTMM[16]NVYN[1170]STK" to: num 1170
        
        gg <- IDPdb$glycoform.mass
        gg <- as.character(gg)
        gg <- strsplit(gg, split=".", fixed=T)
        GlycanMass <- sapply(gg, function(x) x[1])  
        
        IDPdb$GlycanMass <- as.numeric(GlycanMass)
        
        # Calculate Y1 values
        
        modMass <- IDPdb$sequence
        
        IDPdb$MonoisotopicY1mass <- mapply(sum, 203.079373 , IDPdb$mass)
        
        MonoisotopicY1mass <- IDPdb$MonoisotopicY1mass
        MonoisotopicY1mass[is.na(MonoisotopicY1mass)] <- 0
        IDPdb$MonoisotopicY1mass <- MonoisotopicY1mass
        
        IDPdb$MonoisotopicPeptideMass <- IDPdb$mass
        
        IDPdb <- IDPdb[!duplicated(IDPdb),]
        
        IDPdb$Group.Source.Spectrum <- IDPdb$scans
        IDPdb$Sequence <- IDPdb$sequence.x
        IDPdb$Precursor.m.z <- IDPdb$precursorMZ
        
        
        return(IDPdb) 
        
}

##############################################################################
# ## psm
##############################################################################


psm <- 
        function (sequence, spec, FUN = defaultIon, plot = FALSE, 
                  fi = fragmentIon(sequence, 
                                   FUN = FUN)[[1]], fragmentIonError = 0.02) 
        {
                n <- nchar(sequence)
                pim <- fi$y[nrow(fi)]
                by.mZ <- numeric()
                by.label <- character()
                fi.names <- names(fi)
                for (i in 1:ncol(fi)) {
                        by.mZ <- c(by.mZ, fi[, i])
                        by.label <- c(by.label, paste(fi.names[i], 1:n, sep = ""))
                }
                out <- .C("findNN_", nbyion = as.integer(length(by.mZ)), 
                          nmZ = as.integer(length(spec$mZ)), byion = as.double(by.mZ), 
                          mZ = as.double(spec$mZ), NN = as.integer(rep(-1, length(by.mZ))))
                mZ.error <- spec$mZ[out$NN + 1] - by.mZ
                if (plot == TRUE) {
                        plot(mZ.error[mZ.error.idx <- order(mZ.error)], 
                             ylim = c(-5 * fragmentIonError, 5 * fragmentIonError), 
                             pch = "o", cex = 0.5 )
                        
                        hits = (abs(mZ.error) < fragmentIonError)
                        nHits <- sum(hits)
                        sumMZerror = round(sum(abs(mZ.error[hits])), 2)
                        avgMZerror = round(sumMZerror/nHits, 2)
                        cover = round(nHits/(nrow(fi) * ncol(fi)), 2)
                        legend("topleft", paste(c("nHits", "sumMZerror", "avgMZerror", 
                                                  "cover"), as.character(c(nHits, sumMZerror, avgMZerror, 
                                                                           cover)), sep = "="))
                }
                return(list(mZ.Da.error = mZ.error, mZ.ppm.error = 1e+06 * 
                                    mZ.error/by.mZ, idx = out$NN + 1, label = by.label, score = -1, 
                            sequence = sequence, fragmentIon = fi))
        }

##############################################################################
#   Function 15. Restrict the IDPdb from glycoMod search to the retention time
# determined and recorded in table 2
##############################################################################
rt.restrict <- function(IDPdb, Table2, rt.min.minus, rt.min.plus) {
        
        rr <- numeric()
        
        ii <- 0
        
        for (i in 1:length(unique(IDPdb$table2Sequence))) {
                
                ii = ii + 1
                
table2Sequence <- unique(IDPdb$table2Sequence)[i]
                
                
IDPdb.table2Sequence <- IDPdb[IDPdb$table2Sequence==table2Sequence,]
                
# subset data calculate stuff and write a data frame
table2Sequence.min.rt <- 
        ((Table2$min.rt[Table2$table2Sequence==unique(IDPdb.table2Sequence$table2Sequence)]) - rt.min.minus)
table2Sequence.max.rt <- 
        ((Table2$min.rt[Table2$table2Sequence==unique(IDPdb.table2Sequence$table2Sequence)]) + rt.min.plus)
                
                
IDPdb.subset<-  subset(IDPdb.table2Sequence, IDPdb.table2Sequence$Scan.Time >= table2Sequence.min.rt 
                                       & IDPdb.table2Sequence$Scan.Time <= table2Sequence.max.rt)
                
                
                
                
                rr <- rbind(rr, IDPdb.subset)
                
        }
        return(as.data.frame(rr))
}
##### The End ###########################################################