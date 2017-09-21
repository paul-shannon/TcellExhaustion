library(trena)
library(trenaViz)
library(colorspace)
library(annotate)
library(org.Mm.eg.db)
library(MotifDb)
#--------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#--------------------------------------------------------------------------------
stopifnot(packageVersion("trena")    >= "0.99.149")
stopifnot(packageVersion("trenaViz") >= "0.99.11")

if(!exists("trena"))
   trena <- Trena("mm10")

PORT.RANGE <- 8000:8020
if(!exists("tv")) {
   tv <- trenaViz(PORT.RANGE, quiet=FALSE)
   setGenome(tv, "mm10")
   }
#----------------------------------------------------------------------------------------------------
readHamidsTCellNetwork <- function()
{
   tbl <- read.table("cytoscapeEdgeAnnot_11July2017.txt", sep="\t", header=TRUE, as.is=TRUE)
   tbl.xtab <- as.data.frame(table(c(tbl$Source, tbl$Target)))
   tbl.xtab <- tbl.xtab[order(tbl.xtab$Freq, decreasing=TRUE),]
   head(tbl.xtab)
      # DNMT3A is the most connected gene, a methyl transferase.
      # extended location: chr12:3,751,728-3,970,655  (218kb)

   tbl

} # readHamidsTCellNetwork
#----------------------------------------------------------------------------------------------------
# gene <- "Dnmt3a" chrom <- "chr12" loc.start <- 3751728 loc.end   <- 3970655
# ifng: Chr10:118441047-118445892 bp, + strand gene <- "Ifng" chrom <- "chr10" loc.start <- 118441047 loc.end   <- 118445892
#----------------------------------------------------------------------------------------------------
#    The are 3 replicates for each of expression and ATAC-seq peaks for 'naive' T cells:
#
#     GSM2365761        RNA_N1
#     GSM2365762        RNA_N2
#     GSM2365763        RNA_N3
#
#     GSM2365799        ATAC_N1
#     GSM2365800        ATAC_N2
#     GSM2365801        ATAC_N3
#
#   These define the starting state of the T cells (point of
#   reference). At this point the cells are mature but "inexperienced".
#
#   Then there are expression and peak data sets with identifier like
#   RNA_L5_1, etc. The L stands for Liver tumor and the number after
#   the L is how many days after the naive cells (starting state) were
#   stimulated with antigen. What I am doing is building a gene
#   regulatory network of interactions that start in state naive and go
#   through states L5, L7, ... till L28. The data includes later days
#   (35 and 60), in these late stages a second differentiation process
#   starts (memory cell formation) and for now we can skip that stage.
#
#   So you have 4 developmental stages (L5, 7, 14, 21) with RNA and
#   ATAC-seq. The genes that I think are the most likely players are
#   listed in the attached Cytoscape annot file.
#
#----------------------------------------------------------------------------------------------------
calculateATACregions <- function(chrom, loc.start, loc.end, display)
{
   atac.files <- grep("_ATAC_", list.files("./data/"), v=TRUE)
   short.names <- unlist(lapply(strsplit(atac.files, "_"),
                          function(tokens)
                             if(grepl("^N", tokens[3]))
                                return(tokens[3])
                             else
                               sprintf("%s-%s", tokens[3], tokens[4])))

   names(atac.files) <- short.names
   showGenomicRegion(tv, sprintf("%s:%d-%d", chrom, loc.start, loc.end))

   atac.file.count <- length(atac.files)
   colors <- rainbow_hcl(atac.file.count)
   i <- 0
   tbls.regions <- list()

   for(id in names(atac.files)[1:atac.file.count]){
      i <- i + 1
      filename <- file.path("./data/", atac.files[[id]])
      tbl.tmp <- read.table(filename, header=TRUE, sep="\t", as.is=TRUE)
        # Tcf7 is on the minus strand
      track.start <- loc.start
      track.end   <- loc.end
      tbl.gene <- subset(tbl.tmp, chr==chrom & start >= track.start & end <= track.end)
      if(nrow(tbl.gene) == 0){
         printf("no atac regions found for %s:%d-%d in %s", chrom, track.start, track.end, filename)
         next;
         }
      tbl.gene$sample = id
      colnames(tbl.gene)[7] <- "score"
      #printf("colors[%d]: %s", i, colors[i])
      tbls.regions[[id]] <- tbl.gene
      # if(display)
        # addBedGraphTrackFromDataFrame(tv, id, tbl.gene[, c(1,2,3,7,5)], color=colors[i], minValue=0, maxValue=750)
      }

   tbl.regions <- do.call(rbind, tbls.regions)[, c(1,2,3,7,8)]
   tbl.regions <- tbl.regions[order(tbl.regions$start, decreasing=FALSE),]
   rownames(tbl.regions) <- NULL
   colnames(tbl.regions)[1] <- "chrom"   # better than "chr"
   tbl.regions

} # calculateATACregions
#----------------------------------------------------------------------------------------------------
# plot ataq-seq scores for each open region
# a possible trend:
#          first three samples are high:    N1-3
#     next seven are about half as much:    L5, L7
#     final six are at about 10% of max:    L14, L21
#----------------------------------------------------------------------------------------------------
plotAtaqSeqScores <- function(tbl.regions)
{
   plot(tbl.regions[1:16, "score"], type="b", col="blue")
   plot(tbl.regions[17:32, "score"], type="b", col="darkgreen")
   plot(tbl.regions[33:48, "score"], type="b", col="red")

} # plotAtaSeqScores
#----------------------------------------------------------------------------------------------------
# what motifs are found in these open chromatin regions?
#----------------------------------------------------------------------------------------------------
findAndDisplayMotifs <- function(tbl.regions, pwmMatchMinimumAsPercentage, source, trackName)
{
   stopifnot(source %in% c("MotifDb", "TFClass"))
   pfms.mouse <- query(query(MotifDb, "jaspar2016"), "mmus")
   pfms.human <- query(query(MotifDb, "jaspar2016"), "hsap")
   pfms <- as.list(c(pfms.human, pfms.mouse))
   mm <- MotifMatcher(genomeName="mm10", pfms)

   tbl.regions.uniq <- unique(tbl.regions[, 1:3])
   print(tbl.regions.uniq)
   tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.regions.uniq, pwmMatchMinimumAsPercentage=pwmMatchMinimumAsPercentage)

   shortMotifs <- unlist(lapply(strsplit(tbl.motifs$motifName, "-"), function(tokens) tokens[length(tokens)]))
   tbl.motifs$shortMotif <- shortMotifs
   tbl.motifs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source=source, expand.rows=TRUE)

   motifs.without.tfs <- which(is.na(tbl.motifs$geneSymbol))
   if(length(motifs.without.tfs) > 0){
      printf("%d/%d motifs had no TF/geneSymbol", length(motifs.without.tfs), nrow(tbl.motifs))
      tbl.motifs <- tbl.motifs[-motifs.without.tfs,]
      }


   motif.tfs <- sort(unique(tbl.motifs$geneSymbol))
   addBedTrackFromDataFrame(tv, trackName,
                            tbl.motifs[, c("chrom", "motifStart", "motifEnd", "motifName", "motifRelativeScore")],
                            color="red")

   invisible(tbl.motifs)

} # findAndDisplayMotifs
#----------------------------------------------------------------------------------------------------
# read in the early-stage RNA, N1-3, L5_1-3, in which the ataq-seq data suggests the chromatin
# is open
#----------------------------------------------------------------------------------------------------
readExpressionFiles <- function()
{
   rna.files <- grep("RNA", list.files("./data/"), v=TRUE)
   if(length(rna.files) == 0){
      stop(sprintf("cannot find rna files from current working directory, %s", getwd()))
      }
   printf ("read %d RNA files", length(rna.files))

   sample.name.tokens <- strsplit(rna.files, "_")
   sample.names <- unlist(lapply(sample.name.tokens, function(tokens) sprintf("%s.%s", tokens[3], tokens[4])))
   names(rna.files) <- sample.names
   rna.files <- as.list(rna.files)

   tbls.rna <- list()
   for(i in seq_len(length(rna.files))){
      file.name <- file.path("./data", rna.files[[i]])
      sample.name <- names(rna.files)[i]
      tbl.tmp <- read.table(file.name, sep="\t", header=FALSE)
      colnames(tbl.tmp) <- c("geneID", sample.names[i])
      rownames(tbl.tmp) <- tbl.tmp$geneID
      tbl.tmp <- tbl.tmp[, 2, drop=FALSE]   # drop the geneID column
      tbls.rna[[i]] <- tbl.tmp
      }

   #----------------------------------------------------------------------------------------------------
   # ensure that all geneIDs are the same, in the same order
   # then combine column-wise
   #----------------------------------------------------------------------------------------------------
   for(i in seq_len(length(tbls.rna)))
      stopifnot(all(rownames(tbls.rna[[1]]) == rownames(tbls.rna[[i]])))

   tbl.rna <- do.call(cbind, tbls.rna)
   gene.symbols <- select(org.Mm.eg.db, keys=rownames(tbl.rna), keytype="ENTREZID", columns=c("SYMBOL"))
   tbl.rna$gene <- gene.symbols$SYMBOL
   dups <- which(duplicated(tbl.rna$gene))
   if(length(dups) > 0)
      tbl.rna <- tbl.rna[-dups, , drop=FALSE]
   nas <- which(is.na(tbl.rna$gene))
   if(length(nas) > 0)
      tbl.rna <- tbl.rna[-nas, , drop=FALSE]

   rownames(tbl.rna) <- tbl.rna$gene
   gene.column <- grep("gene", colnames(tbl.rna))
   mtx.rna <- as.matrix(tbl.rna[, -gene.column])
   rownames(mtx.rna) <- toupper(rownames(mtx.rna))
   zeros <- which(rowSums(mtx.rna) == 0)
   if(length(zeros) > 0)
      mtx.rna <- mtx.rna[-zeros, , drop=FALSE]

   printf("--- mtx.rna, %d x %d", nrow(mtx.rna), ncol(mtx.rna))
   printf("mtx.rna before asinh transform: ")
   print(fivenum(mtx.rna))
   mtx.rna <- asinh(mtx.rna)
   printf("mtx.rna after asinh transform: ")
   print(fivenum(mtx.rna))
   invisible(mtx.rna)

} # readExpressionFiles
#----------------------------------------------------------------------------------------------------
makeModel <- function(trena, targetGene, targetGene.tss, tbl.motifs, mtx.rna, pcaMaxThreshold=1.0)
{
   printf("possible motifs for which we have expression: %d/%d",
          length(intersect(tbl.motifs$geneSymbol, rownames(mtx.rna))), length(unique(tbl.motifs$geneSymbol)))

  chromosome <- "chr11"

  #solver.names <- c("lasso", "lassopv", "pearson", "randomForest", "ridge", "spearman")
  solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")

  tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifs, mtx.rna)
  tbl.geneModel.strong <- subset(tbl.geneModel, pcaMax >= pcaMaxThreshold)

  tbl.regulatoryRegions.strong <- subset(tbl.motifs, geneSymbol %in% tbl.geneModel.strong$gene)   # 71 rows

    # add two additional column to describe the regulatory regions, given that the TCF7 tss is 52283013
    #  motifName chrom motifStart motifEnd strand
    #   MA0516.1 chr11   52282407 52282421      -
    #
    #    1) the distance (+ or -) to the targetGene's TSS
    #    2) an id, e.g., "TCF7.fp.-.000606.MA0516.1"
    #----------------------------------------------------------------------------------------------------

  distance <- tbl.regulatoryRegions.strong$motifStart - targetGene.tss
  direction <- rep("upstream", length(distance))
  direction[which(distance < 0)] <- "downstream"
  tbl.regulatoryRegions.strong$distance.from.tss <- distance
  tbl.regulatoryRegions.strong$id <- sprintf("%s.fp.%s.%06d.%s", targetGene, direction, abs(distance), tbl.regulatoryRegions.strong$motifName)

  #save(tbl.geneModel, tbl.geneModel.strong, tbl.regulatoryRegions.strong, file="tbl.geneModel.8sep2017.9am.RData")

  tbl.tfBindingFreq <- as.data.frame(table(tbl.regulatoryRegions.strong$geneSymbol))
  tbl.tfBindingFreq <- tbl.tfBindingFreq[order(tbl.tfBindingFreq$Freq, decreasing=TRUE),]
  colnames(tbl.tfBindingFreq) <- c("gene", "binding.sites")
  tbl.geneModel.strong <- merge(tbl.geneModel.strong, tbl.tfBindingFreq, by="gene")
  tbl.geneModel.strong <- tbl.geneModel.strong[order(tbl.geneModel.strong$pcaMax, decreasing=TRUE),]
  list(model=tbl.geneModel.strong, regions=tbl.regulatoryRegions.strong)

} # makeModel
#----------------------------------------------------------------------------------------------------
simple.demo <- function()
{
   mtx.rna <- readExpressionFiles()
   targetGene <- "TCF7"
   targetGene.tss <- 52283013

   showGenomicRegion(tv, "chr11:52,282,097-52,283,776")   # about 1300 kb straddling a TSS

     # get that region in a list, a reusable data structure
     # this also allows you to change the current region of interest interactively
   current.region <- parseChromLocString(getGenomicRegion(tv))
   removeTracksByName(tv, getTrackNames(tv)[-1])   # remove all but the first track,  Gencode

   tbl.regions.atac <- with(current.region,
                            calculateATACregions(chrom=chrom, loc.start=start,  loc.end=end, display=TRUE))

      # pshannon added next four lines, emulating the ones which follow

   rowsN <- tbl.regions.atac[1:3, ]
   newN <- rowsN[1, ]
   newN$score <- mean(rowsN$score)
   newN$sample <- "N"

   rows5 <- tbl.regions.atac[4:7, ]
   new5 <- rows5[1, ]
   new5$score <- mean(rows5$score)
   new5$sample <- "TD5"

   rows7 <- tbl.regions.atac[8:10, ]
   new7 <- rows7[1, ]
   new7$score <- mean(rows7$score)
   new7$sample <- "TD7"

   rows14 <- tbl.regions.atac[11:13, ]
   new14 <- rows14[1, ]
   new14$score <- mean(rows14$score)
   new14$sample <- "TD14"

   rows21 <- tbl.regions.atac[14:16, ]
   new21 <- rows21[1, ]
   new21$score <- mean(rows21$score)
   new21$sample <- "TD21"

   HB.tbl.regions.atac <- rbind(newN, new5, new7, new14, new21)
   HB.tbl.regions.atac.toDisplay <- HB.tbl.regions.atac
   colnames(HB.tbl.regions.atac.toDisplay)[1] <- "chr"
   HB.tbl.regions.atac.toDisplay$score <- round(HB.tbl.regions.atac.toDisplay$score)

   colors <- rainbow_hcl(nrow(HB.tbl.regions.atac.toDisplay))
     # clear off all tracks but Gencode, the first track
   removeTracksByName(tv, getTrackNames(tv)[-1])

   for (i in 1:nrow(HB.tbl.regions.atac)) {
     addBedGraphTrackFromDataFrame(tv,
                                   HB.tbl.regions.atac.toDisplay$sample[i],
                                   HB.tbl.regions.atac.toDisplay[i,],       # pshannon added [i,]
                                   displayMode = "COLLAPSED",
                                   color=colors[i],
                                   minValue=0,
                                   maxValue=max(HB.tbl.regions.atac.toDisplay$score))
    }

     # prepare to identify motifs in these regions: collapse the many tracks into the two
     # genomic regions they cover
   tbl.regions.collapsed <- unique(HB.tbl.regions.atac[, c("chrom", "start", "end")])

      # the called function collects and provides motif pfms from MotifDb,
      # in this case, all jaspar2016 human and mouse.
      # other choices are certanly reasonable
      # this step is slow because each motif in the pfms list is matched against
      # the sequence in the two regions
   tbl.motifs <- findAndDisplayMotifs(tbl.regions.collapsed,
                                        pwmMatchMinimumAsPercentage=100,
                                        source="MotifDb", # pshannon capitialized leading 'M'
                                        trackName=sprintf("%s.%d%%", targetGene, 100))
   # now make a trena model.
       # the pcaMaxThreshold is generous: several low-significance gnees are therefore included
       # but because the region is small, let's leave them in for now.
   model <- makeModel(trena, targetGene, targetGene.tss, tbl.motifs, mtx.rna, pcaMaxThreshold=0.5)

       # one option is to now examine model$model, keeping only those (for instance) with a
       # random forest score > 10 or ....
   model$model <- subset(model$model, rf.score > 10)
       # eliminate all motifs whose transcription factors are not named in the gene model
   model$regions <- subset(model$regions, geneSymbol %in% model$model$gene)
       # just one model, but we can make more, and view them on demand in cytoscape.js
   models <- list(simple=model)

   g <- buildMultiModelGraph(tv, targetGene="TCF7", models) # comment out for speed
   g.lo <- addGeneModelLayout(tv, g, xPos.span=1500)
   setGraph(tv, g.lo, names(models))
   setStyle(tv, "style.js")
   fit(tv)

   list(tv=tv, model=model)

} # simple.demo
#----------------------------------------------------------------------------------------------------
hamids.new.mouse.stringent.motifs <- function(match.threshold)
{
   current.region <- parseChromLocString(getGenomicRegion(tv))
   tbl.regions <- data.frame(chrom=current.region$chrom,
                             start=current.region$start,
                             end=current.region$end,
                             stringsAsFactors=FALSE)

   #organism.rows <- grep ('Mmusculus', values(MotifDb)$organism, ignore.case=TRUE)
   #source.rows <- grep ('JASPAR', values(MotifDb)$dataSource, ignore.case=TRUE)
   #mouse.jaspar.rows <- intersect (organism.rows, source.rows)	# 278 PWMs
   #mouse.jaspar.pfms <- MotifDb [mouse.jaspar.rows]

   mouse.jaspar.pfms <- query(query(MotifDb, "JASPAR"), "mmus")  # pshannon simplification

   pfms <- as.list(mouse.jaspar.pfms)
   motifMatcher <- MotifMatcher(genomeName="mm10", pfms)

   tbl.motifs <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions,
						pwmMatchMinimumAsPercentage=match.threshold) # pshannon loosened constraint

   tbl.toDisplay <- tbl.motifs[, c("chrom", "motifStart", "motifEnd",
                                   "motifName", "motifScore")]
   tmp <- strsplit(tbl.toDisplay$motifName, "-")
   tbl.toDisplay$motifName <- unlist(lapply(tmp, function(x) x[3]))

   trackName <- sprintf("mouse.%d", match.threshold)
   addBedTrackFromDataFrame(tv, trackName=trackName, tbl.toDisplay, color="blue")

} # hamids.new.mouse.stringent.motifs
#----------------------------------------------------------------------------------------------------

# >  source(explore.R)  # the latest version, which you will have just pulled down
# >  x <- simple.demo()
# >  names(x) # “tv” “model”
# >  names(x$model) # “model” “regions”   model$model is the gene model data.frame.  regions is all the motif info
# >  dim(x$model$regions) # 49 17
