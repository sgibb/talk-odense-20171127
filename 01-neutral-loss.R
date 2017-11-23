#+ setup, include=FALSE
library("knitr")
opts_chunk$set(cache=TRUE, cache.path=".cache/", fig.width=12)

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

#' # Ion Injection Time plot
hist(fData(ms1)$injectionTime, main="Ion Injection Time MS1", xlim=c(0, 0.055))
hist(fData(ms2)$injectionTime, main="Ion Injection Time MS2", xlim=c(0, 0.055))
hist(fData(ms3)$injectionTime, main="Ion Injection Time MS3", xlim=c(0, 0.055))

plot(rtime(ms1) / 60, fData(ms1)$injectionTime,
     type="l",
     xlim=c(40, 120),
     main="Ion Injection Time MS1", xlab="rt", ylab="intensity")
plot(rtime(ms2) / 60, fData(ms2)$injectionTime,
     xlim=c(40, 120),
     main="Ion Injection Time MS2", xlab="rt", ylab="intensity")
plot(rtime(ms3) / 60, fData(ms3)$injectionTime,
     xlim=c(40, 120),
     main="Ion Injection Time MS3", xlab="rt", ylab="intensity")

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

#' MS level vs rtime
library("RColorBrewer")
col <- paste0(brewer.pal(5, "Set1"), "80")
plot(rtime(ms) / 60, jitter(msLevel(ms), factor=2),
     col=col[msLevel(ms)], pch=20,
     ylim=c(5, 0.5),
     main="MS level vs rtime", xlab="rt", yaxt="n", ylab="")
points(rtime(ms2)[ph80] / 60, jitter(rep(4, sum(ph80))), col=col[4], pch=20)
points(rtime(ms2)[ph97] / 60, jitter(rep(5, sum(ph97))), col=col[5], pch=20)
legend("top", legend=c(paste0("MS", 1:3), "Ph80", "Ph97"),
       col=col, pch=20, bty="n", horiz=TRUE)
axis(2, at=-(1:5), labels=c(paste0("MS", 1:3), "Ph80", "Ph97"))
