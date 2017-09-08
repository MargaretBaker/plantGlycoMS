#' A function to calculate monoisotopic precursor mass
#'
#' this function calculates the monoisotopic precursor mass from then precursor mz
#' @param input a csv file with one column named precursorMZ
#' @keywords calculate
#' @export
#' @examples
#' spectrumTable()

spectrumTable <- function (input) {
        
spectrumTable <- read.csv(file=input, header=TRUE, skip=1)
      
spectrumTable$MonoPrecursorMass <- (spectrumTable$precursorMZ * 
                                            spectrumTable$charge - 
                                            (spectrumTable$charge * 1.007276 )) 
        
        return(spectrumTable)
        
}

#' A function to modify a chainsaw in silico digest
#'
#' this function adds the mass for carbamidomethylation and methionine oxidation and filters for glycopeptides
#' @param digest.N a chainsaw in silico digest with all glycosylation sites having a capital n (N)
#' @param digest.n a chainsaw in silico digest with all glycosylation sites having a lower case n (n)
#' @param carbamidomethylation 57 Da is added to the peptide mass for every cysteine present, default=TRUE
#' @param methionineOxidation 16 Da is added to the peptide mass for every methionine present, and a new row is added to the table so that every methionine containing peptide has a value for oxidized and not oxidized, default=TRUE
#' @param glycoOnly return a dataframe containing only data for peptide with glycosylation sites, default=TRUE
#' @keywords digest
#' @export
#' @examples
#' glycoChainSaw()

glycoChainSaw <- function(digest.N, digest.n, carbamidomethylation=TRUE, methionineOxidation=TRUE, glycoOnly =TRUE) {
  
  
  names(digest.n) <- c("n.sequence", "protein", "mass",  "missedCleavages",  "specificity",     
                       "nTerminusIsSpecific" ,"cTerminusIsSpecific" )
  
  digest.n$sequence <- gsub("n", "N", digest.n$n.sequence)
  
  digest.Nn <- merge(digest.N, digest.n, by="sequence")
  
  digest.Nn <- digest.Nn[c("sequence", "protein.x","mass.x","missedCleavages.x",   
                           "specificity.x","nTerminusIsSpecific.x","cTerminusIsSpecific.x","n.sequence")]         
  
  names(digest.Nn) <- c("sequence","protein","mass","missedCleavages",
                        "specificity","nTerminusIsSpecific","cTerminusIsSpecific" ,"n.sequence"  )
  
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
  
  if (methionineOxidation ==TRUE ) {
    
    MoxPeptides <- grep ("M", digest.Nn$sequence)
    MoxPeptideTable <- digest.Nn[MoxPeptides,]
    
    sequence <- MoxPeptideTable$sequence        
    sequence <- gsub( "M", "m", sequence ,fixed=TRUE)
    MoxPeptideTable$sequence <- sequence
    
    modMass <- MoxPeptideTable$sequence
    modMass <- gsub( "[A-Z]", "0", modMass)
    modMass <- strsplit(modMass , split="")
    modMass <- lapply(modMass , function(x){gsub("m", "15.994915", x)})
    modMass <- lapply(modMass , function(x){as.numeric(x)})
    modMass <- sapply(modMass , function(x) {sum(x)})
    MoxPeptideTable$mass <- mapply(sum, modMass , MoxPeptideTable$mass)
    
    digest.Nn <- rbind(MoxPeptideTable, digest.Nn)
    
  }
  
  if (glycoOnly==TRUE){
    glycoOnly <- grep("n", digest.Nn$n.sequence)
    digest.Nn <- digest.Nn[glycoOnly,]
    
  }
  
  return(digest.Nn)
}

#' A function to add the glycosylation motif for a MyriMatch database search config
#'
#' this function adds N!{P}[ST] * in front of a mass value
#' @param input a vector of dynamic modification masses
#' @keywords import data
#' @export
PreferredDeltaMasses <- function (input)
{
        
        deltaMass<- paste("N!{P}[ST] * ", input$glycoform)
        deltaMass <- unique(deltaMass)
        #cat(deltaMass, file="output/PreferredDeltaMasses.txt")

}

#' A function to import IDPicker results
#'
#' This function imports IDPicker results and calculates relevant values.
#' @param IDPdb a data.frame containing peptide spectrum matches
#' @param ChainSaw a data.frame of insilico digest results (see glycoChainSaw)
#' @param dir a directory containing MS2 binary data
#' @keywords import data
#' @export
#' @examples
#' Read.IDPdb()

Read.IDPdb <- function (IDPdb, ChainSaw,dir) {
        

        IDPdb <- IDPdb[, -c(2:6)]
       IDPdb <- IDPdb[-c(1,2),]
       
      # x <-  data.frame (IDPdb$Group.Source.Spectrum, IDPdb$Sequence)
     #  IDPdb <- IDPdb[!duplicated(x),]
        
        ppm.mass.error <- IDPdb$Mass.Error/IDPdb$Exact.Mass *10^6
        IDPdb$ppm.mass.error <- ppm.mass.error     
        
        exact.precursor.mz <- (IDPdb$Exact.Mass+(1.00727647*IDPdb$Charge))/IDPdb$Charge
        IDPdb$exact.precursor.mz <- exact.precursor.mz
        
        
        # change: factor 0.1.4586  to: numeric 4586 
        sc <- IDPdb$Group.Source.Spectrum
        sc <- strsplit(sc, split=".", fixed=T)
        scans <- sapply(sc, function(x) x[3])
        scans <- as.numeric(scans)
        IDPdb$Group.Source.Spectrum <- scans
        
        # Make a data frame 'MS2.all' with columns 'title' and 'scans'; 
        # it contains this info for all MS2 binary files
        
        title <- list.files(dir)
        
        #chr "MS2Data.mzML.binary.sn1830.txt" to: num  1830
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
        peptideSequence <- gsub(pattern="[]", replacement="", peptideSequence,fixed = TRUE)
     
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

#' A function to import MS2 binary data generated with msaccess (ProteoWizard).
#'
#' this function imports MS2 binary data
#' @param IDPdb a data.frame containing peptide spectrum matches
#' @param dir name of the directory containing the MS2 binary data
#' @keywords import data
#' @export
#' @examples
#' Read.MS2Data()

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

#' A function to calculate Y-ions
#'
#' this function calculates Y-ions
#' @param IDPdb a data.frame containing MonoisotopicPeptideMasses and corresponding MonoisotopicY1mass, and the charge state
#' @keywords calculate
#' @export
#' @examples MonoisotopicPeptideMass <- c(1000, 2000, 3000, 4000)
#' MonoisotopicY1mass <- c(1203.079, 2203.079, 3203.079, 4203.079)
#' charge <- c(2,5,3,4)
#' IDPdb <- data.frame(MonoisotopicY1mass,MonoisotopicPeptideMass, charge)
#' Ions <- calculateIons(IDPdb)

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

#' A function to calculate the Y1-ions
#'
#' this function calculates the m/z value of a Y1-ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y1(1000,4)

find_Y1 <- function(MonoisotopicY1mass, charge ) {
        
        Y1 <- ((MonoisotopicY1mass + ( charge * 1.007276 ) ) / charge)
        
        return(Y1)
        
}

#' A function to calculate the neutral loss of methane sulfenic acid from a Y-1 ion
#'
#' this function calculates the m/z value of the neutral loss of methane sulfenic acid from a Y-1 ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y1Mox(1000,4)

find_Y1Mox <- function(MonoisotopicY1mass, charge ) {
        
        Y1Mox <- (((MonoisotopicY1mass-64.10686) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y1Mox)
        
}

#' A function to calculate the neutral loss of a carbamidomethyl group from a Y-1 ion
#'
#' this function calculates the m/z value of the neutral loss of a carbamidomethyl group from a Y-1 ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y1Ca(1000,4)

find_Y1Ca <- function(MonoisotopicY1mass, charge ) {
        
        Y1Ca <- (((MonoisotopicY1mass-57.021464) + ( charge * 1.007276 ) ) 
                  / charge)
        
        return(Y1Ca)
        
}

#' A function to calculate the Y0-ions
#'
#' this function calculates the m/z value of a Y0-ion
#' @param MonoisotopicPeptideMass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y0(1000,4)

find_PepPlus <- function(MonoisotopicPeptideMass, charge ) {
        
        PP <- ((MonoisotopicPeptideMass + ( charge * 1.007276 ) ) / charge)
        
        return(PP)
        
}

#' A function to calculate the neutral loss of ammonia from a Y0- ion
#'
#' this function calculates the m/z value of the neutral loss of ammonia from a Y0- ion
#' @param MonoisotopicPeptide mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples
#' find_Y0NH3(1000,4)

find_Y0NH3 <- function(MonoisotopicPeptideMass, charge ) {
        
        Y0NH3 <- (((MonoisotopicPeptideMass- 17.026549) + ( charge * 1.007276 ) ) 
                  / charge)
        
        return(Y0NH3)
        
}

#' A function to calculate the Y1- plus fucose ion
#'
#' this function calculates the m/z value of the Y1 plus fucose ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y1F(1000,4)

find_Y1F <- function(MonoisotopicY1mass, charge ) {
        
        Y1F <- (((MonoisotopicY1mass+146.057909) + ( charge * 1.007276 ) ) 
               / charge)
        
        return(Y1F)
        
}

#' A function to calculate the Y2- plus fucose ion
#'
#' this function calculates the m/z value of the Y2 plus fucose ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y2F(1000,4)

find_Y2F <- function(MonoisotopicY1mass, charge ) {
        
        Y2F <- (((MonoisotopicY1mass+146.057909+203.079373) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y2F)
        
}

#' A function to calculate the Y3- plus fucose ion
#'
#' this function calculates the m/z value of the Y3 plus fucose ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y3F(1000,4)

find_Y3F <- function(MonoisotopicY1mass, charge ) {
        
        Y3F <- (((MonoisotopicY1mass+146.057909+203.079373+162.052824) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y3F)
        
}

#' A function to calculate the Y3- plus fucose and xylose ion
#'
#' this function calculates the m/z value of the Y1 plus fucose and xylose ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y3FX(1000,4)

find_Y3FX <- function(MonoisotopicY1mass, charge ) {
        
        Y3FX <- (((MonoisotopicY1mass+146.057909+203.079373+162.052824+132.042259) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y3FX)
        
}

#' A function to calculate the Y2-ion
#'
#' this function calculates the m/z value of the Y2-ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export
#' @examples 
#' find_Y2(1000,4)

find_Y2 <- function(MonoisotopicY1mass, charge ) {
        
        Y2 <- (((MonoisotopicY1mass+203.079373 ) + ( charge * 1.007276 ) ) 
                / charge)
        
        return(Y2)
        
}

#' A function to calculate the Y3-ion
#'
#' this function calculates the m/z value of the Y3-ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export 
#' @examples 
#' find_Y3(1000,4)

find_Y3 <- function(MonoisotopicY1mass, charge ) {
        
        Y3 <- (((MonoisotopicY1mass+ 203.079373 + 162.052824  ) + ( charge * 1.007276 ) ) 
               / charge)
        
        return(Y3)
        
}

#' A function to calculate the Y3 plus xylose-ion
#'
#' this function calculates the m/z value of the Y3 plus xylose-ion
#' @param MonoisotopicY1mass numeric value
#' @param charge numeric value
#' @keywords calculate
#' @export 
#' @examples 
#' find_Y3X(1000,4)

find_Y3X <- function(MonoisotopicY1mass, charge ) {
        
        Y3X <- (((MonoisotopicY1mass+ 203.079373 + 162.052824 +  132.042259 ) + ( charge * 1.007276 ) ) 
               / charge)
        
        return(Y3X)
        
}

#' A function to compile data for the 'data' argument in plantGlycoMS
#'
#' this function compiles data for the 'data' argument in plantGlycoMS
#' @param analysis a character string with the name of the analysis, defaults to "analysis"
#' @param sampleName a character string with the name of the sample, defaults to "results"
#' @param IDPdb a data.frame containing the peptide spectrum matches, defaults to IDPdb
#' @param MS2Data a list of MS2 binary text files, defaults to MS2Data
#' @param Ions, a list of ion m/z values for each peptide spectrum match, defaults to Ions
#' @keywords compile
#' @export
#' @examples 
#' compileData()

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
                                     Scan.Time=IDPdb$Scan.Time[i], 
                                     
                                     sequence = IDPdb$Sequence[i], 
                                     peptideSequence=IDPdb$peptideSequence[i], 
                                     n.sequence=IDPdb$n.sequence[i],
                                     table2Sequence=IDPdb$table2Sequence[i],
                                     
                                     exact.mass=IDPdb$Exact.Mass[i],
                                     obs.mass=IDPdb$Observed.Mass[i], 
                                     exact.precursor.mz=IDPdb$exact.precursor.mz[i],
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

#' A function to validate and annotate plant glycopeptide-spectrum matches
#'
#' this function validates and annotates plant glycopeptide-spectrum matches
#' @param data a list compiled with compileData function
#' @param modification vector containing the modification mass
#' @param modificationName a character vector containing the name of the modification
#' @param mZmarkerIons a vector containing the mass values of the oxonium ions
#' @param minMarkerIons the minimum number of oxonium ions required for a match to be valid, default is 2
#' @param itol_ppm the MS ion tolerance in ppm, default is 15
#' @param minMarkerIntensityRatio the minimum intensity ratio of marker ions, default is 2
#' @param peakplot02 if true, plots will be generated, default is TRUE
#' @param validate if true, gPSMs will be validated, if false, all will be plotted, default is FALSE
#' @keywords validation
#' @export
#' @examples 
#' gPSMvalidator()


gPSMvalidator <-
        
        function (data, modification, modificationName, mZmarkerIons, 
                  minMarkerIons = 2, itol_ppm = 15, minMarkerIntensityRatio = 2, 
                  peakplot02=TRUE, validate=FALSE) 
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
                        
                        idx <- protViz::findNN(mZmarkerIons, data[[i]]$mZ)
                        ppm.error <- 1e-06 * itol_ppm * data[[i]]$mZ[idx]
                        b <- (abs(mZmarkerIons - data[[i]]$mZ[idx]) < ppm.error)
                        sum.mZmarkerIons.intensity <- sum(data[[i]]$intensity[idx[b]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZmarkerIons.intensity
                        percent.mZmarkerIons <- round(100 * (sum.mZmarkerIons.intensity/
                                                                     (sum.mZmarkerIons.intensity + sum.intensity)), 1)
                        idx.ppm.error <- (mZmarkerIons[b]-data[[i]]$mZ[idx][b])/mZmarkerIons[b]*1e06
                        idx.mZ <- data[[i]]$mZ[idx][b]
                        
                        idY1 <- protViz::findNN(data[[i]]$Y1 , data[[i]]$mZ)
                        ppm.errorY1 <- 1e-06 * itol_ppm * data[[i]]$mZ[idY1]
                        c <- (abs(data[[i]]$Y1  - data[[i]]$mZ[idY1]) < ppm.errorY1)
                        sum.mZY1Ions.intensity <- sum(data[[i]]$intensity[idY1[c]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY1Ions.intensity
                        percent.mZY1Ions <- round(100 *(sum.mZY1Ions.intensity/
                                                                (sum.mZY1Ions.intensity +  sum.intensity)), 1)
                        idY1.ppm.error <- (data[[i]]$Y1[c]-data[[i]]$mZ[idY1][c])/data[[i]]$Y1[c]*1e06
                        idY1.mZ <- data[[i]]$mZ[idY1][c]
                        
                        idY1Mox <- protViz::findNN(data[[i]]$Y1Mox , data[[i]]$mZ)
                        ppm.errorY1Mox <- 1e-06 * itol_ppm * data[[i]]$mZ[idY1Mox]
                        bb <- (abs(data[[i]]$Y1Mox  - data[[i]]$mZ[idY1Mox]) < ppm.errorY1Mox)
                        sum.mZY1MoxIons.intensity <- sum(data[[i]]$intensity[idY1Mox[bb]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY1MoxIons.intensity
                        percent.mZY1MoxIons <- round(100 *(sum.mZY1MoxIons.intensity/
                                                                   (sum.mZY1MoxIons.intensity +  sum.intensity)), 1)
                        idY1Mox.ppm.error <- (data[[i]]$Y1Mox[bb]-data[[i]]$mZ[idY1Mox][bb])/data[[i]]$Y1Mox[bb]*1e06
                        idY1Mox.mZ <- data[[i]]$mZ[idY1Mox][bb]
                        
                        
                        idY1Ca <- protViz::findNN(data[[i]]$Y1Ca , data[[i]]$mZ)
                        ppm.errorY1Ca <- 1e-06 * itol_ppm * data[[i]]$mZ[idY1Ca]
                        ee <- (abs(data[[i]]$Y1Ca  - data[[i]]$mZ[idY1Ca]) < ppm.errorY1Ca)
                        sum.mZY1CaIons.intensity <- sum(data[[i]]$intensity[idY1Ca[ee]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY1CaIons.intensity
                        percent.mZY1CaIons <- round(100 *(sum.mZY1CaIons.intensity/
                                                                  (sum.mZY1CaIons.intensity +  sum.intensity)), 1)
                        idY1Ca.ppm.error <- (data[[i]]$Y1Ca[ee]-data[[i]]$mZ[idY1Ca][ee])/data[[i]]$Y1Ca[ee]*1e06
                        idY1Ca.mZ <- data[[i]]$mZ[idY1Ca][ee]
                        
                        idY2 <- protViz::findNN(data[[i]]$Y2 , data[[i]]$mZ)
                        ppm.errorY2 <- 1e-06 * itol_ppm * data[[i]]$mZ[idY2]
                        dd <- (abs(data[[i]]$Y2  - data[[i]]$mZ[idY2]) < ppm.errorY2)
                        sum.mZY2Ions.intensity <- sum(data[[i]]$intensity[idY2[dd]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY2Ions.intensity
                        percent.mZY2Ions <- round(100 *(sum.mZY2Ions.intensity/
                                                                (sum.mZY2Ions.intensity +  sum.intensity)), 1)
                        idY2.ppm.error <- (data[[i]]$Y2[dd]-data[[i]]$mZ[idY2][dd])/data[[i]]$Y2[dd]*1e06
                        idY2.mZ <- data[[i]]$mZ[idY2][dd]
                        
                        idY1F <- protViz::findNN(data[[i]]$Y1F , data[[i]]$mZ)
                        ppm.errorY1F <- 1e-06 * itol_ppm * data[[i]]$mZ[idY1F]
                        ff <- (abs(data[[i]]$Y1F  - data[[i]]$mZ[idY1F]) < ppm.errorY1F)
                        sum.mZY1FIons.intensity <- sum(data[[i]]$intensity[idY1F[ff]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY1FIons.intensity
                        percent.mZY1FIons <- round(100 *(sum.mZY1FIons.intensity/
                                                                 (sum.mZY1FIons.intensity +  sum.intensity)), 1)
                        idY1F.ppm.error <- (data[[i]]$Y1F[ff]-data[[i]]$mZ[idY1F][ff])/data[[i]]$Y1F[ff]*1e06
                        idY1F.mZ <- data[[i]]$mZ[idY1F][ff]
                        
                        idY2F <- protViz::findNN(data[[i]]$Y2F , data[[i]]$mZ)
                        ppm.errorY2F <- 1e-06 * itol_ppm * data[[i]]$mZ[idY2F]
                        gg <- (abs(data[[i]]$Y2F  - data[[i]]$mZ[idY2F]) < ppm.errorY2F)
                        sum.mZY2FIons.intensity <- sum(data[[i]]$intensity[idY2F[gg]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY2FIons.intensity
                        percent.mZY2FIons <- round(100 *(sum.mZY2FIons.intensity/
                                                                 (sum.mZY2FIons.intensity +  sum.intensity)), 1)
                        idY2F.ppm.error <- (data[[i]]$Y2F[gg]-data[[i]]$mZ[idY2F][gg])/data[[i]]$Y2F[gg]*1e06
                        idY2F.mZ <- data[[i]]$mZ[idY2F][gg]
                        
                        idY3F <- protViz::findNN(data[[i]]$Y3F , data[[i]]$mZ)
                        ppm.errorY3F <- 1e-06 * itol_ppm * data[[i]]$mZ[idY3F]
                        kk <- (abs(data[[i]]$Y3F  - data[[i]]$mZ[idY3F]) < ppm.errorY3F)
                        sum.mZY3FIons.intensity <- sum(data[[i]]$intensity[idY3F[kk]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY3FIons.intensity
                        percent.mZY3FIons <- round(100 *(sum.mZY3FIons.intensity/
                                                                 (sum.mZY3FIons.intensity +  sum.intensity)), 1)
                        idY3F.ppm.error <- (data[[i]]$Y3F[kk]-data[[i]]$mZ[idY3F][kk])/data[[i]]$Y3F[kk]*1e06
                        idY3F.mZ <- data[[i]]$mZ[idY3F][kk]
                        
                        idY3FX <- protViz::findNN(data[[i]]$Y3FX , data[[i]]$mZ)
                        ppm.errorY3FX <- 1e-06 * itol_ppm * data[[i]]$mZ[idY3FX]
                        ll <- (abs(data[[i]]$Y3FX  - data[[i]]$mZ[idY3FX]) < ppm.errorY3FX)
                        sum.mZY3FXIons.intensity <- sum(data[[i]]$intensity[idY3FX[ll]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY3FXIons.intensity
                        percent.mZY3FXIons <- round(100 *(sum.mZY3FXIons.intensity/
                                                                  (sum.mZY3FXIons.intensity +  sum.intensity)), 1)
                        idY3FX.ppm.error <- (data[[i]]$Y3FX[ll]-data[[i]]$mZ[idY3FX][ll])/data[[i]]$Y3FX[ll]*1e06
                        idY3FX.mZ <- data[[i]]$mZ[idY3FX][ll]
                        
                        idY3 <- protViz::findNN(data[[i]]$Y3 , data[[i]]$mZ)
                        ppm.errorY3 <- 1e-06 * itol_ppm * data[[i]]$mZ[idY3]
                        cc <- (abs(data[[i]]$Y3  - data[[i]]$mZ[idY3]) < ppm.errorY3)
                        sum.mZY3Ions.intensity <- sum(data[[i]]$intensity[idY3[cc]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY3Ions.intensity
                        percent.mZY3Ions <- round(100 *(sum.mZY3Ions.intensity/
                                                                (sum.mZY3Ions.intensity +  sum.intensity)), 1)
                        idY3.ppm.error <- (data[[i]]$Y3[cc]-data[[i]]$mZ[idY3][cc])/data[[i]]$Y3[cc]*1e06
                        idY3.mZ <- data[[i]]$mZ[idY3][cc]
                        
                        idY3X <- protViz::findNN(data[[i]]$Y3X , data[[i]]$mZ)
                        ppm.errorY3X <- 1e-06 * itol_ppm * data[[i]]$mZ[idY3X]
                        hh <- (abs(data[[i]]$Y3X  - data[[i]]$mZ[idY3X]) < ppm.errorY3X)
                        sum.mZY3XIons.intensity <- sum(data[[i]]$intensity[idY3X[hh]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY3XIons.intensity
                        percent.mZY3XIons <- round(100 *(sum.mZY3XIons.intensity/
                                                                 (sum.mZY3XIons.intensity +  sum.intensity)), 1)
                        idY3X.ppm.error <- (data[[i]]$Y3X[hh]-data[[i]]$mZ[idY3X][hh])/data[[i]]$Y3X[hh]*1e06
                        idY3X.mZ <- data[[i]]$mZ[idY3X][hh]
                        
                        idPepPlus <- protViz::findNN(data[[i]]$PepPlus , data[[i]]$mZ)
                        ppm.errorPepPlus <- 1e-06 * itol_ppm * data[[i]]$mZ[idPepPlus]
                        f <- (abs(data[[i]]$PepPlus  - data[[i]]$mZ[idPepPlus]) < 
                                      ppm.errorPepPlus)
                        sum.mZPepPlusIons.intensity <- sum(data[[i]]$intensity[idPepPlus[f]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZPepPlusIons.intensity
                        percent.mZPepPlusIons <- round(100 * (sum.mZPepPlusIons.intensity/
                                                                      (sum.mZPepPlusIons.intensity +  sum.intensity)), 1)
                        idPepPlus.ppm.error <- (data[[i]]$PepPlus[f]-data[[i]]$mZ[idPepPlus][f])/data[[i]]$PepPlus[f]*1e06
                        idPepPlus.mZ <- data[[i]]$mZ[idPepPlus][f]
                        
                        idY0NH3 <- protViz::findNN(data[[i]]$Y0NH3 , data[[i]]$mZ)
                        ppm.errorY0NH3 <- 1e-06 * itol_ppm * data[[i]]$mZ[idY0NH3]
                        mm <- (abs(data[[i]]$Y0NH3  - data[[i]]$mZ[idY0NH3]) < ppm.errorY0NH3)
                        sum.mZY0NH3Ions.intensity <- sum(data[[i]]$intensity[idY0NH3[mm]])
                        sum.intensity <- sum(data[[i]]$intensity) - sum.mZY0NH3Ions.intensity
                        percent.mZY0NH3Ions <- round(100 *(sum.mZY0NH3Ions.intensity/
                                                                   (sum.mZY0NH3Ions.intensity +  sum.intensity)), 1)
                        idY0NH3.ppm.error <- (data[[i]]$Y0NH3[mm]-data[[i]]$mZ[idY0NH3][mm])/data[[i]]$Y0NH3[mm]*1e06
                        idY0NH3.mZ <- data[[i]]$mZ[idY0NH3][mm]
                        
                        ido <- protViz::findNN(otherOxonium, data[[i]]$mZ)
                        ppm.errorOO <- 1e-06 * itol_ppm * data[[i]]$mZ[ido]
                        e <- (abs(otherOxonium - data[[i]]$mZ[ido]) < ppm.errorOO)
                        sum.otherOxonium.intensity <- sum(data[[i]]$intensity[ido[e]])
                        ido.ppm.error <- (otherOxonium[e]-data[[i]]$mZ[ido][e])/otherOxonium[e]*1e06
                        ido.mZ <-   data[[i]]$mZ[ido][e]     
                        
                        if(validate)  {   
                                
                                if ((data[[i]]$GlycanMass == 203) &
                                    (length(data[[i]]$mZ[idx[b]]) >= minMarkerIons) & 
                                    percent.mZmarkerIons > minMarkerIntensityRatio |
                                    
                                    (data[[i]]$GlycanMass == 0) |
                                        (length(data[[i]]$mZ[idx[b]]) >= minMarkerIons) & 
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
                                                   Scan.Time = data[[i]]$Scan.Time,
                                                   structure = data[[i]]$structure,
                                                   precursorMZ = data[[i]]$precursorMZ,
                                                   exact.mass = data[[i]]$exact.mass,
                                                   exact.precursor.mz = data[[i]]$exact.precursor.mz,
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
                                                
                                                
                                                
                                                fi <- protViz::fragmentIon(sequence = data[[i]]$peptideSequence, 
                                                                  FUN = defaultIons, 
                                                                  modified = substr(data[[i]]$modification, 2, 
                                                                                    nchar(data[[i]]$modification) - 1), 
                                                                  modification = modification)
                                                
                                                fi.by <- as.data.frame(cbind(b = fi[[1]]$b, y = fi[[1]]$y))
                                                
                                                check<- peakplot02(peptideSequence=data[[i]]$peptideSequence, spec = data[[i]], 
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
                                             ylim = c(-2* itol_ppm, 2* itol_ppm)
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
                                           Scan.Time = data[[i]]$Scan.Time,
                                           structure = data[[i]]$structure,
                                           exact.precursor.mz = data[[i]]$exact.precursor.mz,
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
                                        
                                        
                                        
                                        fi <- protViz::fragmentIon(sequence = data[[i]]$peptideSequence, 
                                                          FUN = defaultIons, 
                                                          modified = substr(data[[i]]$modification, 2, 
                                                                            nchar(data[[i]]$modification) - 1), 
                                                          modification = modification)
                                        
                                        fi.by <- as.data.frame(cbind(b = fi[[1]]$b, y = fi[[1]]$y))
                                        
                                        check<- peakplot02(data[[i]]$peptideSequence, spec = data[[i]], 
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
                                               idY2F.ppm.error, idY3F.ppm.error,
                                               idY3FX.ppm.error,
                                               idY3.ppm.error, idY3X.ppm.error,
                                               idPepPlus.ppm.error, 
                                               idY0NH3.ppm.error, ido.ppm.error,
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
                                     ylim = c(-2* itol_ppm, 2* itol_ppm))
                                
                                
                                
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

#' A function to define default ion values
#'
#' this function defines default ion values
#' @param b mass of a b-ion
#' @param y mass of a y-ion
#' @keywords protViz
#' @export
#' @examples 
#' defaultIons()

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

#' A function to import RQ data
#'
#' this function is used to import RQ data and associate it with SIC data
#' @param input a csv file containing exact precursor mass values
#' @param dir a directory containing SIC data
#' @keywords relative quantitation
#' @keywords import
#' @export
#' @examples 
#' Read.RQ()

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

#' A function to import SIC data
#'
#' this function is used to import with SIC data generated with msaccess (ProteoWizard)
#' @param dir a directory containing SIC data
#' @keywords relative quantitation
#' @keywords import
#' @export
#' @examples 
#' Read.SICData()

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

#' A function to plot MS data
#'
#' this function plots MS data
#' @param peptideSequence a peptide sequence
#' @param spec a spectrum
#' @param FUN a function to calculate ions, default is defaultIons
#' @param fi a fragment ion function, default is protViz::fragmentIon(peptideSequence, FUN = FUN)[[1]] 
#' @param main the main title of a plot, default is NULL
#' @param sub the subtitle of a plot, default is NULL
#' @param xlim the limits of the x axis, default is range(spec$mZ, na.rm = TRUE)
#' @param the limits of the y axis, default is ylim = range(spec$intensity, na.rm = TRUE)
#' @param itol the fragment ion mass tolerance in delta mz, default is 0.02
#' @param pattern.abc the abc pattern, default is "[abc].*"
#' @param pattern.xyz the xyz pattern, default is "[xyz].*", ion.axes = TRUE)
#' @keywords protViz
#' @export
#' @examples 
#' peakplot02()

peakplot02 <- 
        function (peptideSequence, spec, FUN = defaultIons, 
                  fi = protViz::fragmentIon(peptideSequence, FUN = FUN)[[1]], 
                  main = NULL, 
                  sub = NULL, 
                  xlim = range(spec$mZ, na.rm = TRUE), 
                  ylim = range(spec$intensity, na.rm = TRUE), 
                  itol = 0.02, pattern.abc = "[abc].*", 
                  pattern.xyz = "[xyz].*", ion.axes = TRUE) 
        {
                n <- nchar(peptideSequence)
                m <- psm02(peptideSequence, spec, FUN, fi = fi, plot = FALSE)
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

#' A function to import SIC data
#'
#' this function is used to import with SIC data generated with msaccess (ProteoWizard)
#' @param dir the directory with SIC data
#' @param Asn PLACEHOLDER
#' @keywords relative quantitation
#' @keywords import
#' @export
#' @examples 
#' Read.SICData()

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

#' A function to estimate relative abunance of glycoforms
#'
#' this function estimates relative abunance of glycoforms
#' @param RQ a data.frame
#' @param rtTable a data.frame with retention time information
#' @param dir a directory with SIC data
#' @param rt.min.minus retention time minus
#' @param rt.min.plus retention time plus
#' @keywords relative quantitation
#' @export
#' @examples 
#' glycoRQ()

glycoRQ <- function(RQ, rtTable, dir, rt.min.minus, rt.min.plus) {
        
        rr <- numeric()
        
        ii <- 0
        
        for (i in 1:length(unique(RQ$Asn))) {
                
                ii = ii + 1
                
                Asn <- unique(RQ$Asn)[i]
                
                #Read in SIC data from the RQ directory for the glycosylation site of interest
                SICData <- Read.SICData(dir=dir, Asn=Asn)
                
                RQ.Asn <- RQ[RQ$Asn==Asn,]
                
                
               Asn.min.rt <- 
                        ((rtTable$min.rt[rtTable$table2Sequence==unique(RQ.Asn$table2Sequence)]) - rt.min.minus)
                Asn.max.rt <- 
                       ((rtTable$min.rt[rtTable$table2Sequence==unique(RQ.Asn$table2Sequence)]) + rt.min.plus)
             
                
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

#' A function to validate a glycan composition prediction
#'
#' this function uses rules based on the plant N-glycan biosynthesis pathway to test if 
#' a predicted glycan composition is biologically possible
#' @param structure the glycan composition
#' @param data a csv file from the GlycoMod search
#' @keywords glycan
#' @export
#' @examples 
#' pGlycoFilter()

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

#' A function to import GlycoMod data
#'
#' this function is for importing GlycoMod data
#' @param input glycopeptide identifications from GlycoMod
#' @param ChainSaw an in silico digest
#' @param spectrum.table a summary table of MS data
#' @param dir a directory with MS2 binary data
#' @keywords import
#' @export
#' @examples 
#' Read.GlycoMod()

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
        
        
        exact.precursor.mz <- (IDPdb$Exact.Mass+(1.00727647*IDPdb$charge))/IDPdb$charge
        IDPdb$exact.precursor.mz <- exact.precursor.mz
        
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
        #modification<- gsub( "NF", "10",  modification,fixed=TRUE)
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
        
        oo <- order(IDPdb$scans)
        IDPdb <- IDPdb[oo,]
        
        
        return(IDPdb) 
        
}

#' A function to match peptides to spectra
#'
#' this function plots MS data
#' @param sequence the peptide sequence
#' @param spec a spectrum
#' @param FUN a function for calculating ions, default is defaultIon
#' @param plot should the psm be plotted, default is FALSE
#' @param fi how are fragment ions calculated, default is protViz::fragmentIon(sequence, FUN = FUN)[[1]]
#' @param fragmentIonError the mass error for fragment ions in mz, default is 0.02
#' @keywords protViz
#' @export
#' @examples 
#' psm02()


psm02 <- 
        function (sequence, spec, FUN = defaultIon, plot = FALSE, 
                  fi = protViz::fragmentIon(sequence, 
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
                out <- .C("__findNN_", nbyion = as.integer(length(by.mZ)), 
                          nmZ = as.integer(length(spec$mZ)), byion = as.double(by.mZ), 
                          mZ = as.double(spec$mZ), NN = as.integer(rep(-1, length(by.mZ))),
                          PACKAGE="protViz")
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

#' A function to filter gPSMs by retention time
#'
#' this function filters gPSMs by retention time
#' @param IDPdb a table with peptide spectrum match information generated with IDPicker
#' @param rtTable a table with information about identified glycopeptides and their retention times
#' @param rt.min.minus amount of minutes to subtract from the retention time minimum
#' @param rt.min.plus amount of minutes to add to the retention time maximum
#' @keywords retentionTime
#' @export
#' @examples 
#' rt.restrict()
		
rt.restrict <- function(IDPdb, rtTable, rt.min.minus, rt.min.plus) {
  
  rr <- numeric()
  
  ii <- 0
  
  for (i in 1:length(unique(IDPdb$table2Sequence))) {
    
    ii = ii + 1
    
    table2Sequence <- unique(IDPdb$table2Sequence)[i]
    
    
    IDPdb.table2Sequence <- IDPdb[IDPdb$table2Sequence==table2Sequence,]
    
    # subset data calculate stuff and write a data frame
    table2Sequence.min.rt <- 
      ((rtTable$min.rt[rtTable$table2Sequence==unique(IDPdb.table2Sequence$table2Sequence)]) - rt.min.minus)
    table2Sequence.max.rt <- 
      ((rtTable$min.rt[rtTable$table2Sequence==unique(IDPdb.table2Sequence$table2Sequence)]) + rt.min.plus)
    
    
    IDPdb.subset<-  subset(IDPdb.table2Sequence, IDPdb.table2Sequence$Scan.Time >= table2Sequence.min.rt 
                           & IDPdb.table2Sequence$Scan.Time <= table2Sequence.max.rt)
    
    
    
    
    rr <- rbind(rr, IDPdb.subset)
    
  }
  return(as.data.frame(rr))
}

#' A function to make a retention time table
#'
#' this function filters gPSMs by retention time
#' @param dat a data.frame with gPSM data
#' @keywords table
#' @export
#' @examples 
#' retentionTimeTable()

retentionTimeTable <- function(dat) {
  
  t2S <- unique(dat$table2Sequence)
  dat$glyco <- as.integer(dat$GlycanMass)
  
  ##  How many times was that peptide observed?  ###########################
  count <- sapply(1:length(t2S), function(i) {
    length(dat$table2Sequence[dat$table2Sequence == t2S[i]])
  })
  
  ##  When did it start eluting? ##############################
  min.rt <- sapply(1:length(t2S), function(j) {
    min(dat$Scan.Time[dat$table2Sequence == t2S[j]])
  })        
  
  ##  When did it stop eluting?  ##############################
  max.rt <- sapply(1:length(t2S), function(k) {
    max(dat$Scan.Time[dat$table2Sequence == t2S[k]])
  })        
  
  ####  What is the peptide sequence associated with each peptide mass?  ############################
  peptideMass <- sapply(1:length(t2S), function(l) {
    unique(dat$MonoisotopicPeptideMass[dat$table2Sequence == t2S[l]])
  }) 
  
  n.sequence <- sapply(1:length(t2S), function(l) {
    unique(dat$n.sequence[dat$table2Sequence == t2S[l]])
  }) 
  
  ####  What is the Asn number?  ############################
  #Asn <- sapply(1:length(t2S), function(n) {
  #  unique(dat$Asn[dat$table2Sequence == t2S[n]])
  # }) 
  
  
  ##  what is the smallest glycan? ##############################
  glycan.S <- sapply(1:length(t2S), function(m) {
    min(dat$glyco[dat$table2Sequence == t2S[m]])
  }) 
  
  ##  what is the largest glycan?  ##############################
  glycan.L <- sapply(1:length(t2S), function(k) {
    max(dat$glyco[dat$table2Sequence == t2S[k]])
  })
  
  ####
  
  Table2 <- data.frame(
    # Asn, 
    peptideMass, n.sequence, t2S, 
    count, min.rt, max.rt, glycan.L,
    glycan.S, stringsAsFactors=F)
  
  #  oo <- order(Table2$Asn)
  #  Table2 <- Table2[oo,]
  
  Table2 <- data.frame( #Table2$Asn, 
    Table2$peptideMass,
    Table2$n.sequence, Table2$t2S, 
    Table2$count, Table2$min.rt, Table2$max.rt, 
    Table2$glycan.S, Table2$glycan.L)
  
  names(Table2) <- c(#"Asn", 
    "peptideMass", "n.sequence", 
    "table2Sequence", "count", "min.rt",
    "max.rt", "glycan.S", "glycan.L")
  
  return(as.data.frame(Table2))
  # return(list(Table2=Table2, dat=dat))
  
}

#' A function to make the RQ input table
#'
#' this function makes the RQ input table
#' @param gPSMs.ALL a data.frame with all gPSM data
#' @keywords import
#' @export
#' @examples 
#' Read.RQinput()

Read.RQinput <- function (gPSMs.ALL) {
        
        
        
        dat <- subset(gPSMs.ALL, select = c(charge, table2Sequence, 
                                            GlycanMass, exact.precursor.mz))
        RQ_input <- unique(dat)
        
        Asn <- subset(RQ_input, select=c(table2Sequence, charge))
        uniqueAsn <- unique(Asn)
        uniqueAsn$Asn <- 1:nrow(uniqueAsn)
        RQ_input <- merge (uniqueAsn, RQ_input)
        
        
        
        return(RQ_input)
        
}

#' A function to change N to n for predicted glycosylation sites
#'
#' this function changes N to n for predicted glycosylation sites
#' @param fasta a character string containing the amino acid sequence of a protein. A fasta file can be imported with the R package 'seqinr'. 
#' @keywords digest
#' @export
#' @examples 
#' proteinSeq <- "DLQIGFYNQSCPSAESLVQQAVAAAFANNSGIAPGLIRMHFHDCFVRGCDASVLLDSTANNTAEKDAAPNNPSLRGFEVIAAAKSAVEAACPKTVSCADILAFAARDSAALAGNITYQVPSGRRDGNVSLASEALTNIPAPTFNATQLINSFAGKNLTADEMVTLSGAHSIGVSHCFSFLNRIYNFSNTSQVDPTLSSSYADLLRTKCPSNSTRFTPITVSLDIITPTVLDNRYYTGVQLTLGLLTSDQALVTEANLSAAVKNNADNLTAWVAKFAQAIVKMGQIQVLTGTQGEIRTNCSVVNSAS"
#' glycoChange(fasta=proteinSeq)
glycoChange <- function(fasta)
{
        
        
        glycoFasta <- lapply(fasta, function(x){gsub("NAS",
                                                     "nAS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NCS",
                                                          "nCS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NDS",
                                                          "nDS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NES",
                                                          "nES",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NFS",
                                                          "nFS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NGS",
                                                          "nGS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NHS",
                                                          "nHS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NIS",
                                                          "nIS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NKS",
                                                          "nKS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NLS",
                                                          "nLS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NMS",
                                                          "nMS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NNS",
                                                          "nNS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NQS",
                                                          "nQS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NRS",
                                                          "nRS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NSS",
                                                          "nSS",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NTS",
                                                          "nTS",x,fixed=F)})	
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NVS",
                                                          "nVS",x,fixed=F)})	
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NYS",
                                                          "nYS",x,fixed=F)})	
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NWS",
                                                          "nWS",x,fixed=F)})	
        
        
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NAT",
                                                          "nAT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NCT",
                                                          "nCT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NDT",
                                                          "nDT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NET",
                                                          "nET",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NFT",
                                                          "nFT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NGT",
                                                          "nGT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NHT",
                                                          "nHT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NIT",
                                                          "nIT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NKT",
                                                          "nKT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NLT",
                                                          "nLT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NMT",
                                                          "nMT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NNT",
                                                          "nNT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NQT",
                                                          "nQT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NRT",
                                                          "nRT",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NST",
                                                          "nST",x,fixed=F)})
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NTT",
                                                          "nTT",x,fixed=F)})	
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NVT",
                                                          "nVT",x,fixed=F)})	
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NYT",
                                                          "nYT",x,fixed=F)})	
        glycoFasta <- lapply(glycoFasta, function(x){gsub("NWT",
                                                          "nWT",x,fixed=F)})	
        return(glycoFasta)
}
