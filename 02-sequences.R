#+ setup, include=FALSE
library("knitr")
opts_chunk$set(cache=TRUE, cache.path=".cache/sequences/", fig.width=12)

#' # Read .fasta files
library("Biostrings")

proteom <- readAAStringSet("data/sequences/Human_uniprot_ref_can_2015_05_25.fasta.gz")
insulin <- proteom[[grep("INS_HUMAN", names(proteom))]]

#' # In-silico Digestion

library(cleaver)
insulin

## single sequence;
## enyzm == trypsin, split after K/R if not followed by P
cleavedInsulin <- cleave(insulin, enzym="trypsin")
cleavedInsulin

## different protease;
## enzym == pepsin, split after AA (not D/E) if followed by A/F/I/L/M/V
cleave(insulin, enzym="thermolysin")

## whole proteom
cleave(proteom)

#' Fragment calculation

library("MSnbase")
knitr::kable(calculateFragments(as.character(cleavedInsulin[3])))

#' Fragment information for spectra
## find path to a mzXML file
quantFile <- dir(system.file(package="MSnbase", dir="extdata"),
                 full.name=TRUE, pattern="mzXML$")
## find path to a mzIdentML file
identFile <- dir(system.file(package="MSnbase", dir="extdata"),
                 full.name=TRUE, pattern="dummyiTRAQ.mzid")
## create basic MSnExp
msexp <- readMSData(quantFile, verbose=FALSE, centroided.=TRUE)
## add ID information
msexp <- addIdentificationData(msexp, id=identFile)

plot(msexp[[2]], fData(msexp)$sequence[2])

#' Predicting Isotopes (toy example)
library("BRAIN")
useBRAIN(list(C=10, H=86))

#' Predicting Isotopes Insulin
library("RColorBrewer")
col <- paste0(brewer.pal(7, "Set1"))

## create empty plot area
plot(NA, xlim=c(170, 4300), ylim=c(0, 1),
     xlab="mass", ylab="relative intensity",
     main="tryptic digested insulin - isotopic distribution")

## loop through peptides
for (i in seq(along=cleavedInsulin)) {
  ## count C, H, N, O, S atoms in current peptide
  atoms <- BRAIN::getAtomsFromSeq(cleavedInsulin[[i]])
  ## calculate isotopic distribution
  d <- useBRAIN(atoms)
  ## draw peaks
  lines(d$masses, d$isoDistr, type="h", col=col[i])
}

#+ plotinsulinzoom, echo=FALSE
#' mz = 170:280
## create empty plot area
plot(NA, xlim=c(170, 280), ylim=c(0, 1),
     xlab="mass", ylab="relative intensity",
     main="tryptic digested insulin - isotopic distribution")

## loop through peptides
for (i in seq(along=cleavedInsulin)) {
  atoms <- BRAIN::getAtomsFromSeq(cleavedInsulin[[i]])
  d <- useBRAIN(atoms)
  lines(d$masses, d$isoDistr, type="h", col=col[i])
}
