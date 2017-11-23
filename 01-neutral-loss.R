#+ setup, include=FALSE
library("knitr")
opts_chunk$set(cache=TRUE, cache.path=".cache/")

#' # MSnbase example
#' Setup
library("MSnbase")

msdir <- file.path("data", "neutral-loss")
msfile <- file.path(
    msdir,
    "LUM2_01470_KS_L1-5-2_EC17-123-150In2_55oC_NLtrig.mzML.gz"
)

#ms <- readMSData(msfile, mode="onDisk")
#saveRDS(ms, file=file.path(msdir, "mnltrig.rds"))
ms <- readRDS(file.path(msdir, "mnltrig.rds"))

ms1 <- filterMsLevel(ms, 1)
ms2 <- filterMsLevel(ms, 2)
ms3 <- filterMsLevel(ms, 3)

#' # TIC plot
plot(rtime(ms1) / 60, tic(ms1),
     type="l", col="red",
     xlim=c(40, 120),
     main="TIC", xlab="rt", ylab="intenisty")

#' # BPI plot
plot(rtime(ms1) / 60, fData(ms1)$basePeakIntensity,
     type="l", col="blue",
     xlim=c(40, 120),
     main="BPI", xlab="rt", ylab="intensity")

#' # MS/MS per minute
ms2pm <- round(rtime(ms2) / 60)

hist(ms2pm)
plot(table(ms2pm), type="l",
     main="MS2 per min", xlab="rt", ylab="number of MS2 spectra")

#' # NL trigger

#' own top5 function
.topMz <- function(x, n=5) {
    mz(x)[order(intensity(x), decreasing=TRUE)[1:n]]
}

#' top5 for each spectrum
top5mz <- spectrapply(ms2, .topMz)
# convert to matrix
top5mz <- do.call(rbind, top5mz)
#' get precursor mz
precmz <- precursorMz(ms2)
#' &Delta; mz
deltamz <- precmz - top5mz
ph <- c(phosphate=80, phosphoNL=97.9763)

#' How many phosphates?
ph80 <- as.logical(rowSums(abs(deltamz - ph["phosphate"]) < 0.5))
sum(ph80)

#' How many phosphoNL?
ph97 <- as.logical(rowSums(abs(deltamz - ph["phosphoNL"]) < 0.5))
sum(ph97)
