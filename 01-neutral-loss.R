#+ setup, include=FALSE
library("knitr")
opts_chunk$set(cache=TRUE, cache.path=".cache/neutral-loss/", fig.width=12)

#' # MSnbase example
#' Setup
library("MSnbase")

# create some colors
library("RColorBrewer")
col <- paste0(brewer.pal(5, "Set1"))
names(col) <- c("MS1", "MS2", "MS3", "NL80", "NL97")

msdir <- file.path("data", "neutral-loss")
msfile <- file.path(
    msdir,
    "LUM2_01470_KS_L1-5-2_EC17-123-150In2_55oC_NLtrig.mzML.gz"
)

#+ readmsdata, eval=FALSE, echo=1
ms <- readMSData(msfile, mode="onDisk")
saveRDS(ms, file=file.path(msdir, "mnltrig.rds"))

#+ readrds, eval=TRUE, echo=FALSE
ms <- readRDS(file.path(msdir, "mnltrig.rds"))

ms1 <- filterMsLevel(ms, 1)
ms2 <- filterMsLevel(ms, 2)
ms3 <- filterMsLevel(ms, 3)

#' # TIC plot
plot(rtime(ms1) / 60, tic(ms1),
     type="l", col="red",
     xlim=c(40, 120),
     main="TIC", xlab="rt [min]", ylab="intenisty")

#' # BPI plot
plot(rtime(ms1) / 60, fData(ms1)$basePeakIntensity,
     type="l", col="blue",
     xlim=c(40, 120),
     main="BPI", xlab="rt [min]", ylab="intensity")

#' # Ion Injection Time plot
hist(fData(ms1)$injectionTime,
     main="Ion Injection Time MS1", xlim=c(0, 0.055), col=col["MS1"])
hist(fData(ms2)$injectionTime,
     main="Ion Injection Time MS2", xlim=c(0, 0.055), col=col["MS2"])
hist(fData(ms3)$injectionTime,
     main="Ion Injection Time MS3", xlim=c(0, 0.055), col=col["MS3"])

plot(rtime(ms1) / 60, fData(ms1)$injectionTime,
     type="l", col=col["MS1"],
     xlim=c(40, 120),
     main="Ion Injection Time MS1", xlab="rt [min]", ylab="intensity")
plot(rtime(ms2) / 60, fData(ms2)$injectionTime,
     xlim=c(40, 120), col=col["MS2"],
     main="Ion Injection Time MS2", xlab="rt [min]", ylab="intensity")
plot(rtime(ms3) / 60, fData(ms3)$injectionTime,
     xlim=c(40, 120), col=col["MS3"],
     main="Ion Injection Time MS3", xlab="rt [min]", ylab="intensity")

#' # MS1 per minute
ms1pm <- round(rtime(ms1) / 60)

hist(ms1pm, col=col["MS1"])
plot(table(ms1pm), type="l", col=col["MS1"],
     main="MS1 per min", xlab="rt [min]", ylab="number of MS1 spectra")

#' # MS2 per minute
ms2pm <- round(rtime(ms2) / 60)

hist(ms2pm, col=col["MS2"])
plot(table(ms2pm), type="l", col=col["MS2"],
     main="MS2 per min", xlab="rt [min]", ylab="number of MS2 spectra")

#' # MS3 per minute
ms3pm <- round(rtime(ms3) / 60)

hist(ms3pm, col=col["MS3"])
plot(table(ms3pm), type="l", col=col["MS3"],
     main="MS3 per min", xlab="rt [min]", ylab="number of MS3 spectra")

#' # NL trigger

#' own top5 function
.topN <- function(x, n=5) {
    order(intensity(x), decreasing=TRUE)[1:n]
}
.topMz <- function(x, n=5) {
    mz(x)[.topN(x, n)]
}

#' top5 for each spectrum
top5mz <- spectrapply(ms2, .topMz)
# convert to matrix
top5mz <- do.call(rbind, top5mz)
#' get precursor mz
precmz <- precursorMz(ms2)
#' &Delta; mz
deltamz <- precmz - top5mz
nl <- c(phosphate=80, phosphoNL=97.9763)

#' How many phosphates?
nl80 <- as.logical(rowSums(abs(deltamz - nl["phosphate"]) < 0.5))
sum(nl80)

#' How many phosphoNL?
nl97 <- as.logical(rowSums(abs(deltamz - nl["phosphoNL"]) < 0.5))
sum(nl97)

#' MS level vs rtime
#+ mslevelvsrt, fig.width=12, fig.height=12
## add alpha channel
cola <- paste0(col, 80)
plot(rtime(ms) / 60, jitter(msLevel(ms), factor=2),
     col=cola[msLevel(ms)], pch=20,
     ylim=c(5, 0.5),
     main="MS level vs rtime", xlab="rt [min]", yaxt="n", ylab="")
points(rtime(ms2)[nl80] / 60, jitter(rep(4, sum(nl80))), col=cola[4], pch=20)
points(rtime(ms2)[nl97] / 60, jitter(rep(5, sum(nl97))), col=cola[5], pch=20)
legend("top", legend=names(col), col=cola, pch=20, bty="n", horiz=TRUE)
axis(2, at=-(1:5), labels=names(col))

#' spectra vs rtime
#+ msnvsrt, fig.width=12
plot(NA, col=col["MS1"],
     ylim=c(0.1, 1200), xlim=c(0, 120),
     main="# spectra vs rtime", xlab="rt [min]", ylab="# of spectra")
l <- list(ms1pm, ms2pm, ms3pm,
          round(rtime(ms2)[nl80] / 60), round(rtime(ms2)[nl97] / 60))
for (i in seq(along=l)) {
    lines(table(l[[i]]), col=col[i], type="l")
}
legend("top", legend=names(col), col=col, pch=20, bty="n", horiz=TRUE)

#' # Compare NL trigger
#' ms2 precursor with CID energy from MS3 filter string
ms2df <- data.frame(scanId=scanIndex(ms2),
                    pcCid=gsub(".*ms2 *([0-9.]+@cid[0-9.]+) .*", "\\1",
                               fData(ms2)$filterString),
                    stringsAsFactors=FALSE)
ms3df <- data.frame(scanId=scanIndex(ms3),
                    pcCid=gsub(".*ms3 *([0-9.]+@cid[0-9.]+) .*", "\\1",
                               fData(ms3)$filterString),
                    stringsAsFactors=FALSE)
nldf <- merge(ms3df, ms2df,
              all.x=TRUE, all.y=FALSE,
              sort=FALSE, by="pcCid", suffixes=c(".ms3", ".ms2"))
knitr::kable(head(nldf))

# calc scan difference
nldf$delta <- nldf$scanId.ms3 - nldf$scanId.ms2
knitr::kable(head(nldf))

#' just keep entries where MS3 was acquired after MS2
nldf <- nldf[nldf$delta > 0,]

#' sort by scanId.ms3 followed by delta and remove duplicated
nldf <- nldf[order(nldf$scanId.ms3, nldf$delta),]
knitr::kable(head(nldf))

nldf <- nldf[!duplicated(nldf$scanId.ms3),]
knitr::kable(head(nldf))

#' compare nl trigger
ms2idNL <- list(own=sort(scanIndex(ms2)[nl80 | nl97]),
                thermo=nldf$scanId.ms2)

#' plot sets
library("UpSetR")
upset(fromList(ms2idNL), order.by="freq")

#' find unique ones
com <- intersect(ms2idNL$own, ms2idNL$thermo)
uOwn <- setdiff(ms2idNL$own, ms2idNL$thermo)
uThermo <- setdiff(ms2idNL$thermo, ms2idNL$own)


#+ plotms2, include=FALSE
plotMs2 <- function(s, xlim=range(mz(s)), ylim=range(intensity(s)), tol=0.5, topn=5) {
    plot(NA, xlim=xlim, ylim=ylim, xlab="m/z", ylab="intensity")
    lines(mz(s), intensity(s), type="h", col=col["MS2"])
    n <- .topN(s, topn)
    top5mz <- mz(s)[n]
    sel <- unique(c(which(abs(precursorMz(s) - (top5mz + nl["phosphate"])) < tol),
                    which(abs(precursorMz(s) - (top5mz + nl["phosphoNL"])) < tol)))
    nt <- n[sel]

    if (length(nt)) {
        x <- mz(s)[nt]
        y <- intensity(s)[nt]
        nlcol <- col[ifelse(precursorMz(s) - x > nl["phosphate"], "NL97", "NL80")]
        lines(x, y, col=nlcol, lwd=1.5, type="h")
        points(x, y, col=nlcol, pch=20)
        text(x, y, paste0(round(x, 4), "(top ", sel, ")"), pos=3, col=nlcol)
        lg <- paste(c("precursorMz", "NL mz"),
                    round(c(precursorMz(s), mz(s)[nt]), 4), sep=": ")
    } else {
        mi <- min(intensity(s)[n])
        abline(v=precursorMz(s) - nl, col=col[c("NL80", "NL97")], lty=2)
        text(precursorMz(s) - nl, max(intensity(s)), c("NL80", "NL97"),
             pos=2, col=col[c("NL80", "NL97")])
        abline(h=mi, col="#808080", lty=4)
        text(xlim[1], mi, paste("top", topn), col="#808080", pos=3)
        lg <- paste("precursorMz", round(precursorMz(s), 4), sep=": ")
    }
    legend("topleft", legend=lg, bty="n")
}

#' plot some common ones
plotMs2(ms[[com[1]]])
plotMs2(ms[[com[2]]])

#' plot some unique ones
plotMs2(ms[[uOwn[1]]])
plotMs2(ms[[uOwn[2]]])
plotMs2(ms[[uThermo[1]]])
plotMs2(ms[[uThermo[2]]])

#' minimal differences
i <- match(com, scanIndex(ms2))
minDiffCom <- apply(deltamz[i,], 1, function(r)min(abs(r - rep(nl, each=length(r)))))
i <- match(uOwn, scanIndex(ms2))
minDiffOwn <- apply(deltamz[i,], 1, function(r)min(abs(r - rep(nl, each=length(r)))))
i <- match(uThermo, scanIndex(ms2))
minDiffThermo <- apply(deltamz[i,], 1, function(r)min(abs(r - rep(nl, each=length(r)))))

brks <- seq(0, 2.1, by=0.1)

hist(minDiffCom, breaks=brks,
     main="Minimal Differences in precursorMz / NL triggered in both")
hist(minDiffOwn, breaks=brks,
     main="Minimal Differences in precursorMz / NL that should be triggered")
hist(minDiffThermo, breaks=brks,
     main="Minimal Differences in precursorMz / NL triggered just in Thermo")

#' # Increase NL threshold
threshold <- 1

#' How many phosphates?
nl80 <- as.logical(rowSums(abs(deltamz - nl["phosphate"]) < threshold))
sum(nl80)

#' How many phosphoNL?
nl97 <- as.logical(rowSums(abs(deltamz - nl["phosphoNL"]) < threshold))
sum(nl97)

#' compare nl trigger
ms2idNL <- list(own=sort(scanIndex(ms2)[nl80 | nl97]),
                thermo=nldf$scanId.ms2)

#' plot sets
upset(fromList(ms2idNL), order.by="freq")
