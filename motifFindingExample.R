library(trena)
library(trenaViz)
library(org.Mm.eg.db)
library(MotifDb)
#trena <- Trena("mm10")

PORT.RANGE <- 8000:8020
if(!exists("tv")) {
   tv <- trenaViz(PORT.RANGE, quiet=FALSE)
   setGenome(tv, "mm10")
   }

showGenomicRegion(tv, "TCF7")
showGenomicRegion(tv, "chr11:52282407-52283222")

# adjust as you wish, then
current.region <- parseChromLocString(getGenomicRegion(tv))
tbl.regions <- data.frame(chrom=current.region$chrom,
                          start=current.region$start,
                          end=current.region$end,
                          stringsAsFactors=FALSE)

print(tbl.regions)

# creating just this minimal data.frame:
#   chrom    start      end
# 1 chr11 52282407 52283222

# configure the MotifMatcher, a class within the trena package
pfms.mouse <- query(query(MotifDb, "jaspar2016"), "mmus")
pfms.human <- query(query(MotifDb, "jaspar2016"), "hsap")
pfms <- as.list(c(pfms.human, pfms.mouse))

motifMatcher <- MotifMatcher(genomeName="mm10", pfms)

tbl.motifs <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions, pwmMatchMinimumAsPercentage=98)

tbl.toDisplay <- tbl.motifs[, c("chrom", "motifStart", "motifEnd", "motifName", "motifScore")]

addBedTrackFromDataFrame(tv, trackName="test", tbl.toDisplay, color="blue")
