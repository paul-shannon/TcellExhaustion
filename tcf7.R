library(trena)
library(trenaViz)
library(colorspace)
library(annotate)
library(org.Mm.eg.db)
library(MotifDb)
library(biomaRt)
library(RUnit)
#--------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#--------------------------------------------------------------------------------
stopifnot(packageVersion("trena")    >= "0.99.192")
stopifnot(packageVersion("trenaViz") >= "0.99.23")
stopifnot(packageVersion("MotifDb")  >= "1.19.18")

if(!exists("trena"))
   trena <- Trena("mm10", quiet=TRUE)


PORT.RANGE <- 8000:8020
#tv <- FALSE   # no browser visualization needed

if(!exists("tv")){
   tv <-  trenaViz(PORT.RANGE, quiet=TRUE)
   setGenome(tv, "mm10")
   }

if(!exists("mm10.mart")){
   mm10.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
   }

if(!exists("pfms")){
   pfms.jaspar <- query(query(MotifDb, "jaspar"), "mmus")
   pfms.jolma  <- query(query(MotifDb, "jolma2013"), "mmus")
   pfms <- as.list(c(pfms.jaspar, pfms.jolma))   # 411
   }


#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test.transcriptAwareProximalPromoter()
   test.buildModel.allTranscripts()

} # runTests
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

   printf("--- mtx.rna, %d x %d", nrow(mtx.rna), ncol(mtx.rna))
   printf("mtx.rna before asinh transform: ")
   print(fivenum(mtx.rna))
   mtx.rna <- asinh(mtx.rna)

   variance <- apply(mtx.rna, 1, var)
   low.variance.genes <- which(variance < 0.25)
   if(length(low.variance.genes) > 0)
      mtx.rna <- mtx.rna[-low.variance.genes, , drop=FALSE]

   printf("mtx.rna, %d x %d range of values after asinh transform: ", nrow(mtx.rna), ncol(mtx.rna))
   print(fivenum(mtx.rna))
   printf("mtx.rna, range of variance after asinh transform: ")
   print(fivenum(apply(mtx.rna, 1, var)))

   invisible(mtx.rna)

} # readExpressionFiles
#----------------------------------------------------------------------------------------------------
if(!exists("mtx.rna"))
   mtx.rna <- readExpressionFiles()
#----------------------------------------------------------------------------------------------------
readHamidsTCellNetwork <- function()
{
   tbl <- read.table("data/cytoscapeEdgeAnnot_11July2017.txt", sep="\t", header=TRUE, as.is=TRUE)
   tbl.xtab <- as.data.frame(table(c(tbl$Source, tbl$Target)))
   tbl.xtab <- tbl.xtab[order(tbl.xtab$Freq, decreasing=TRUE),]
   head(tbl.xtab)
      # DNMT3A is the most connected gene, a methyl transferase.
      # extended location: chr12:3,751,728-3,970,655  (218kb)

   tbl

} # readHamidsTCellNetwork
#----------------------------------------------------------------------------------------------------
calculateATACregions <- function(chrom, loc.start, loc.end, display=FALSE)
{
   atac.files <- grep("_ATAC_", list.files("./data/"), v=TRUE)
   short.names <- unlist(lapply(strsplit(atac.files, "_"),
                          function(tokens)
                             if(grepl("^N", tokens[3]))
                                return(tokens[3])
                             else
                               sprintf("%s-%s", tokens[3], tokens[4])))

   names(atac.files) <- short.names
   if(display)
      showGenomicRegion(tv, sprintf("%s:%d-%d", chrom, loc.start, loc.end))

   atac.file.count <- length(atac.files)
   colors <- rainbow_hcl(atac.file.count)
   i <- 0
   tbls.regions <- list()

   for(id in names(atac.files)[1:atac.file.count]){
      i <- i + 1
      filename <- file.path("./data/", atac.files[[id]])
      tbl.tmp <- read.table(filename, header=TRUE, sep="\t", as.is=TRUE)
      track.start <- loc.start
      track.end   <- loc.end
      tbl.gene <- subset(tbl.tmp, chr==chrom & start >= track.start & end <= track.end)
      if(nrow(tbl.gene) == 0){
         printf("no atac regions found for %s:%d-%d in %s", chrom, track.start, track.end, filename)
         next;
         }
      tbl.gene$sample = id
      colnames(tbl.gene)[7] <- "score"
      tbls.regions[[id]] <- tbl.gene
      if(display)
         addBedGraphTrackFromDataFrame(tv, id, tbl.gene[, c(1,2,3,7,5)], color=colors[i], minValue=0, maxValue=750)
      }
   tbl.regions <- do.call(rbind, tbls.regions)
   if(is.null(tbl.regions))
      return(data.frame())

   tbl.regions <- tbl.regions[, c(1,2,3,7,8)]
   tbl.regions <- tbl.regions[order(tbl.regions$start, decreasing=FALSE),]
   rownames(tbl.regions) <- NULL
   colnames(tbl.regions)[1] <- "chrom"   # better than "chr"
   tbl.regions

} # calculateATACregions
#----------------------------------------------------------------------------------------------------
findMotifs <- function(pfms, tbl.regions, pwmMatchMinimumAsPercentage, source, display=FALSE, trackName=NA_character_)
{
   mm <- MotifMatcher(genomeName="mm10", pfms)

   tbl.regions.uniq <- unique(tbl.regions[, 1:3])
   # print(tbl.regions.uniq)
   tbl.motifs <- findMatchesByChromosomalRegion(mm, tbl.regions.uniq, pwmMatchMinimumAsPercentage=pwmMatchMinimumAsPercentage)
   tbl.motifs$seq <- unlist(lapply(tbl.motifs$seq,
                                   function(seq){
                                      if(nchar(seq) > 16){
                                         seq <- sprintf("%s...", substring(seq, 1, 13))
                                      };
                                      return(seq)
                                   }))

   if(nrow(tbl.motifs) == 0){
      printf("--- no match of pfms in supplied regions at %d%%", pwmMatchMinimumAsPercentage)
      return(data.frame())
      }

   shortMotifs <- unlist(lapply(strsplit(tbl.motifs$motifName, "-"), function(tokens) tokens[length(tokens)]))
   tbl.motifs$shortMotif <- shortMotifs
   tbl.motifs <- associateTranscriptionFactors(MotifDb, tbl.motifs, source="MotifDb", expand.rows=TRUE)

   motifs.without.tfs <- which(is.na(tbl.motifs$geneSymbol))

   if(length(motifs.without.tfs) > 0){
      printf("%d/%d motifs had no TF/geneSymbol, removing", length(motifs.without.tfs), nrow(tbl.motifs))
      tbl.motifs <- tbl.motifs[-motifs.without.tfs,]
      }


   motif.tfs <- sort(unique(tbl.motifs$geneSymbol))
   if(display){
      addBedTrackFromDataFrame(tv, trackName,
                               tbl.motifs[, c("chrom", "motifStart", "motifEnd", "motifName", "motifRelativeScore")],
                               color="red")
      }

   invisible(tbl.motifs)

} # findMotifs
#----------------------------------------------------------------------------------------------------
transcriptAwareProximalPromoter <- function(targetGene, upstream, downstream)
{
     # selected from tbl.mart[-grep("homolog", tbl.mart$name),].  see also p_value, transcript_gencode_basic

   column.names <- c("ensembl_transcript_id", "chromosome_name", "transcription_start_site",
                     "transcript_start", "transcript_end",
                     "mgi_symbol", "strand", "transcript_tsl", "transcript_length")

   tbl.geneInfo <- getBM(attributes=column.names, filters="mgi_symbol", value=targetGene, mart=mm10.mart)
   if(nrow(tbl.geneInfo) == 0)
      return(NA)

   tbl.geneInfo$transcript_tsl <- sub("tsl", "", tbl.geneInfo$transcript_tsl)
   index <- grep("_tsl", colnames(tbl.geneInfo))
   colnames(tbl.geneInfo)[index] <- "support"

      # make sure all transcripts are on the same strand
   strand <- unique(tbl.geneInfo$strand)
   stopifnot(length(strand) == 1)

   chromosome <- sprintf("chr%s", unique(tbl.geneInfo$chromosome_name[1]))
   stopifnot(length(chromosome) == 1)
   tbl.geneInfo$chromosome_name <- chromosome

      # assume + strand
   tss <- tbl.geneInfo$transcription_start_site
   start <- tss - upstream
   end   <- tss + downstream

   if(strand == -1){
      tss <- tbl.geneInfo$transcription_start_site
      start <- tss - downstream
      end   <- tss + upstream
      }

   tbl.geneInfo$promoter_start <- start
   tbl.geneInfo$promoter_end <- end

   tbl.geneInfo

} # transcriptAwareProximalPromoter
#----------------------------------------------------------------------------------------------------
test.transcriptAwareProximalPromoter <- function()
{
   printf("--- test.transcriptAwareProximalPromoter")
   upstream <- 2000
   downstream <- 200

   tbl.pp <- transcriptAwareProximalPromoter("Prdm1", upstream, downstream)
   checkTrue(all((tbl.pp$promoter_end - tbl.pp$promoter_start) == upstream + downstream))
   checkTrue(nrow(tbl.pp) >= 3)
   checkTrue(ncol(tbl.pp) >= 10)

} # test.transcriptAwareProximalPromoter
#----------------------------------------------------------------------------------------------------
buildModel <- function(targetGene, transcript=NA, chrom, start, end, pfms, motifMatchThreshold, display=FALSE)
{
   if(display == TRUE && !exists("tv")){
      tv <<-  trenaViz(PORT.RANGE, quiet=TRUE)
      setGenome(tv, "mm10")
      }

      # increase the atac search region in order to catch those which exceed our possibly small proximal promoter

   tbl.regions.atac <- calculateATACregions(chrom=chrom, loc.start=start-5000,  loc.end=end+5000, display)

   if(nrow(tbl.regions.atac) == 0)
      return(list(model=data.frame(), regulatoryRegions=data.frame()))

      # an opportunity here to filter on sample and/or score, i.e., keep only thos
      # regions which have an atac sscore > 100.
      # right now we keep them all, ignoring score, ignoring sample, uniqifying on chrom:start-end

   tbl.candidateRegions <- unique(tbl.regions.atac[, 1:3])
   tbl.atacInPromoter <- as.data.frame(GenomicRanges::intersect(GRanges(tbl.candidateRegions),
                                                                GRanges(seqnames=chrom, IRanges(start=start, end=end))))
   if(nrow(tbl.atacInPromoter) == 0)
      return(list(model=data.frame(), regulatoryRegions=data.frame()))

   colnames(tbl.atacInPromoter) <- c("chrom", "start", "end", "width", "strand")

   if(display){
      featureName <- targetGene
      if(!is.na(transcript))
         featureName <- sprintf("%s.%s", featureName, transcript)
      atacTrackName <- sprintf("atacInPromoter.%s", featureName)
      promoterTrackName <- sprintf("promoter.%s", featureName)
      addBedTrackFromDataFrame(tv, promoterTrackName, tbl.promoter, color="#C29DDE")
      addBedTrackFromDataFrame(tv, atacTrackName, tbl.atacInPromoter, color="#C29DDE")
      } # display
   tbl.motifsInRegulatoryRegions  <- findMotifs(pfms, tbl.atacInPromoter, motifMatchThreshold)
   printf("%d motifs in atac seq region about threshold %d", nrow(tbl.motifsInRegulatoryRegions), motifMatchThreshold)
   if(display){
      addBedTrackFromDataFrame(tv, "motifs",
                       tbl.motifsInRegulatoryRegions[, c("chrom", "motifStart", "motifEnd", "motifName", "motifRelativeScore")],
                       color="red")
      } # display

   solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")
   printf("--- creating gene model")
   suppressWarnings(
      tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifsInRegulatoryRegions, mtx.rna)
      )
   if(nrow(tbl.geneModel) > 0)
      tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rfScore, decreasing=TRUE),]

   regRegions.rows.to.keep <- sort(unlist(lapply(tbl.geneModel$gene,
                                                 function(geneSymbol){
                                                    grep(geneSymbol, toupper(tbl.motifsInRegulatoryRegions$geneSymbol))})))

   tbl.regulatoryRegions <- tbl.motifsInRegulatoryRegions[regRegions.rows.to.keep,]
   list(model=tbl.geneModel, regulatoryRegions=tbl.regulatoryRegions)

} # buildModel
#----------------------------------------------------------------------------------------------------
buildModel.allTranscripts <- function(targetGene, tbl.promotersByTranscript, pfms, motifMatchThreshold)
{
   coi <- c("chromosome_name", "promoter_start", "promoter_end", "ensembl_transcript_id")
   chrom <- tbl.promotersByTranscript$chromosome_name[1]
   promoters.start <- min(tbl.promotersByTranscript$promoter_start)
   promoters.end   <- max(tbl.promotersByTranscript$promoter_end)
   shoulder <- 5000
   tbl.regions.atac <- calculateATACregions(chrom=chrom, promoters.start - shoulder,
                                            promoters.end + shoulder)

   tbl.regions.atac.unique <- unique(tbl.regions.atac[, 1:3])

   tbl.atacInPromoter <- as.data.frame(intersect(GRanges(tbl.promotersByTranscript[, coi]),
                                                 GRanges(tbl.regions.atac.unique)))

   printf("%s: found %d intersection/s between transcript promoters and %d atac-seq regions",
          targetGene, nrow(tbl.atacInPromoter), nrow(tbl.regions.atac.unique))

   if(nrow(tbl.atacInPromoter) == 0)
      return(list(model=data.frame(), regulatoryRegions=data.frame()))

   colnames(tbl.atacInPromoter) <- c("chrom", "start", "end", "width", "strand")

   tbl.motifsInRegulatoryRegions  <- findMotifs(pfms, tbl.atacInPromoter, motifMatchThreshold)
   tbl.motifsInRegulatoryRegions$geneSymbol  <- toupper(tbl.motifsInRegulatoryRegions$geneSymbol)
   printf("%d motifs in atac seq region above threshold %d", nrow(tbl.motifsInRegulatoryRegions), motifMatchThreshold)
   printf("tf candidates with expression data: %d/%d",
          length(intersect(rownames(mtx.rna), unique(tbl.motifsInRegulatoryRegions$geneSymbol))),
          length(unique(tbl.motifsInRegulatoryRegions$geneSymbol)))


   printf("--- creating gene model for %s, using promoters from %d transcripts, motif match >= %d",
          targetGene, nrow(tbl.promotersByTranscript), motifMatchThreshold)

   solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")
   suppressWarnings(
      tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifsInRegulatoryRegions, mtx.rna)
      )

   if(nrow(tbl.geneModel) > 0)
      tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rfScore, decreasing=TRUE),]

   regRegions.rows.to.keep <- sort(unlist(lapply(tbl.geneModel$gene,
                                                 function(geneSymbol){
                                                    grep(geneSymbol, toupper(tbl.motifsInRegulatoryRegions$geneSymbol))})))

   tbl.regulatoryRegions <- unique(tbl.motifsInRegulatoryRegions[regRegions.rows.to.keep,])
   printf("returning regulatory model of %d tfs for target gene %s, based on %d unique motif/location/strand locations",
          nrow(tbl.geneModel), targetGene, nrow(tbl.regulatoryRegions))

   list(model=tbl.geneModel, regulatoryRegions=tbl.regulatoryRegions)

} # buildModel.allTranscripts
#----------------------------------------------------------------------------------------------------
tcf7.buildModel.allTranscripts <- function()
{
   targetGene <- "TCF7"
   upstream <- 2000
   downstream <- 200
   motifMatchThreshold <- 90

   tbl.promotersByTranscript <- transcriptAwareProximalPromoter(targetGene, upstream, downstream)

   if(display){
      if(!exists("tv")){
         tv <<-  trenaViz(PORT.RANGE, quiet=TRUE)
         setGenome(tv, "mm10")
         }
      removeTracksByName(tv, getTrackNames(tv)[-1])
      region.comprehensive <- sprintf("%s:%d-%d",
                                      unique(tbl.promotersByTranscript$chromosome_name),
                                         min(tbl.promotersByTranscript$transcript_start) - 100000,
                                         max(tbl.promotersByTranscript$transcript_end)   + 100000)
      showGenomicRegion(tv, region.comprehensive)
      columns.of.interest <- c("chromosome_name", "promoter_start", "promoter_end", "ensembl_transcript_id")
      tbl.promoterRegions <- tbl.promotersByTranscript[, columns.of.interest]
      addBedTrackFromDataFrame(tv, "promoters", tbl.promoterRegions, color="blue", displayMode="EXPANDED")
      chromLoc <- parseChromLocString(getGenomicRegion(tv))
      tbl.regions.atac <- calculateATACregions(chrom=chromLoc$chrom, loc.start=chromLoc$start,
                                               loc.end=chromLoc$end, display=FALSE)
      tbl.regions.atac <- unique(tbl.regions.atac[, 1:3])
      addBedTrackFromDataFrame(tv, "ATAC", tbl.regions.atac, color="darkGreen", displayMode="EXPANDED")
      tbl.atacInPromoter <- as.data.frame(GenomicRanges::intersect(GRanges(tbl.regions.atac), GRanges(tbl.promoterRegions)))
      addBedTrackFromDataFrame(tv, "ATAC.promoter", tbl.atacInPromoter, color="orange", displayMode="EXPANDED")
      } # if display

   x <- buildModel.allTranscripts(targetGene, tbl.promotersByTranscript, pfms, 90)
   tbl.bed <- x$regulatoryRegions[, c("chrom", "motifStart", "motifEnd", "motifName")]
   removeTracksByName(tv, "motifs")
   addBedTrackFromDataFrame(tv, "motifs", tbl.bed, displayMode="EXPANDED", color="magenta")

} # tcf7.buildModel.allTranscripts
#----------------------------------------------------------------------------------------------------
test.buildModel.allTranscripts <- function(display=FALSE)
{
   printf("--- test.buildModel.allTranscripts")
   targetGene <- "BACH2"
   upstream <- 2000
   downstream <- 200
   motifMatchThreshold <- 90

   tbl.promotersByTranscript <- transcriptAwareProximalPromoter(target.gene, upstream, downstream)

   if(display){
      if(!exists("tv")){
         tv <<-  trenaViz(PORT.RANGE, quiet=TRUE)
         setGenome(tv, "mm10")
         }
      removeTracksByName(tv, getTrackNames(tv)[-1])
      region.comprehensive <- sprintf("%s:%d-%d",
                                      unique(tbl.promotersByTranscript$chromosome_name),
                                         min(tbl.promotersByTranscript$transcript_start) - 100000,
                                         max(tbl.promotersByTranscript$transcript_end)   + 100000)
      showGenomicRegion(tv, region.comprehensive)
      columns.of.interest <- c("chromosome_name", "promoter_start", "promoter_end", "ensembl_transcript_id")
      tbl.promoterRegions <- tbl.promoters[, columns.of.interest]
      addBedTrackFromDataFrame(tv, "promoters", tbl.promoterRegions, color="blue", displayMode="EXPANDED")
      chromLoc <- parseChromLocString(getGenomicRegion(tv))
      tbl.regions.atac <- calculateATACregions(chrom=chromLoc$chrom, loc.start=chromLoc$start,
                                               loc.end=chromLoc$end, display=FALSE)
      tbl.regions.atac <- unique(tbl.regions.atac[, 1:3])
      addBedTrackFromDataFrame(tv, "ATAC", tbl.regions.atac, color="darkGreen", displayMode="EXPANDED")
      tbl.atacInPromoter <- as.data.frame(GenomicRanges::intersect(GRanges(tbl.regions.atac), GRanges(tbl.promoterRegions)))
      addBedTrackFromDataFrame(tv, "ATAC.promoter", tbl.atacInPromoter, color="orange", displayMode="EXPANDED")
      } # if display


   model.bach2 <- buildModel.allTranscripts(targetGene, tbl.promotersByTranscript, pfms, motifMatchThreshold)

   target.gene.2 <- "PRDM1"
   tbl.promotersByTranscript.prdm1 <- transcriptAwareProximalPromoter(target.gene.2, upstream, downstream)
   model.prdm1.90 <- buildModel.allTranscripts(target.gene.2, tbl.promotersByTranscript.prdm1, pfms, 90)
   model.prdm1.85 <- buildModel.allTranscripts(target.gene.2, tbl.promotersByTranscript.prdm1, pfms, 85)

   tbl.promotersByTranscript.lef1 <- transcriptAwareProximalPromoter("LEF1", upstream, downstream)
   model.lef1 <- buildModel.allTranscripts("LEF1", tbl.promotersByTranscript.lef1, pfms, 90)


} # test.buildModel.allTranscripts
#----------------------------------------------------------------------------------------------------
buildManyModels <- function (genes.of.interest, upstream, downstream, motifMatchThreshold, display=FALSE)
{
   failures <- c()   # misnamed, lacking atac-seq regions in proximal promoters, no tfs for motifs, ....

   model.tbls <- list()
   region.tbls <- list()

   for(targetGene in genes.of.interest){
      printf("%s", " ")
      printf("--- building model for %s (%d, %d, %d)", targetGene, upstream, downstream,
             motifMatchThreshold)
      tbl.promotersByTranscript <- transcriptAwareProximalPromoter(targetGene, upstream, downstream)
      result <- buildModel.allTranscripts(targetGene, tbl.promotersByTranscript, pfms, motifMatchThreshold)
      if(nrow(result$model) == 0){
         printf("no regulators found for %s", targetGene);
         failures <- c(failures, targetGene)
         next;
         }
      result$model$targetGene <- targetGene
      result$regulatoryRegions$targetGene <- targetGene
      model.tbls[[targetGene]] <- result$model[order(result$model$rfScore, decreasing=TRUE),]
      region.tbls[[targetGene]] <- result$regulatoryRegions
      printf("--- model for %s", targetGene)
      print(result$model)
      } # for targetGene

   tbl.models <- do.call(rbind, model.tbls)
   rownames(tbl.models) <- NULL
   tbl.regions <- do.call(rbind, region.tbls)
   rownames(tbl.regions) <- NULL
   list(models=tbl.models, regions=tbl.regions, failures=failures)

} # buildManyModels
#----------------------------------------------------------------------------------------------------
run <- function()
{
   upstream=2000
   downstream=200
   motifMatchThreshold=90

   target.genes <- unique(tbl.network$Target)
   models.all <- buildManyModels(target.genes, upstream, downstream, motifMatchThreshold, display=FALSE)
   invisible(models.all)

} # run
#----------------------------------------------------------------------------------------------------
modelVersion <- function(pfms, tbl.regulatoryRegions, motifMatchThreshold, mtx)
{
   tbl.motifsInRegulatoryRegions  <- findMotifs(pfms, tbl.regulatoryRegions, motifMatchThreshold)
   tbl.motifs <- associateTranscriptionFactors(MotifDb, tbl.motifsInRegulatoryRegions, source="MotifDb", expand.rows=TRUE)
   tbl.motifs$geneSymbol <- toupper(tbl.motifs$geneSymbol)

   #solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman", "sqrtlasso", "lassopv")
   solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")
   candidates <- unique(intersect(tbl.motifs$geneSymbol, rownames(mtx)))
   save(targetGene, mtx, candidates, solver.names, file=sprintf("solver.bug.%s.RData", gsub(" ", ".", Sys.time(), fixed=TRUE)))
   #suppressWarnings(
      tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifs, mtx)
   #   )

   if(nrow(tbl.geneModel) > 0)
      tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rfScore, decreasing=TRUE),]

   targetGene.tss <- with(tbl.promoters, ifelse(strand[1] == "+",
                                                max(transcription_start_site),
                                                min(transcription_start_site)))

   tbl.motifs.strong <- subset(tbl.motifs, geneSymbol %in% tbl.geneModel$gene)
   distance <- tbl.motifs.strong$motifStart - targetGene.tss
   direction <- rep("upstream", length(distance))
   direction[which(distance < 0)] <- "downstream"
   tbl.motifs.strong$distance.from.tss <- distance
   tbl.motifs.strong$id <- sprintf("%s.fp.%s.%06d.%s", targetGene, direction, abs(distance), tbl.motifs.strong$motifName)

   list(model=tbl.geneModel, regions=tbl.motifs.strong)

} # modelVersion
#------------------------------------------------------------------------------------------------------------------------
explore.one.gene <- function(targetGene)
{
   targetGene <- "TCF7"

   upstream=2000
   downstream=200
   motifMatchThreshold=90

   removeTracksByName(tv, getTrackNames(tv)[-1])

   tbl.promoters <- transcriptAwareProximalPromoter(targetGene, upstream, downstream)

      # our standard approach is to look for binding motifs in ATAC-seq regions
      # in proximal promoters.
      # a useful (sometimes humbling) contrast is to ignore open chromatin,
      # and consider all binding sites in ALL dna in the gene's TSS +/- upstream & downstream
      # calculate the region first

   strand <- tbl.promoters$strand[1]
   targetGene.tss <- with(tbl.promoters, ifelse(strand == "+",
                                                max(transcription_start_site),
                                                min(transcription_start_site)))
   targetGene.upstream <- ifelse(strand == "+", targetGene.tss - upstream, targetGene.tss + upstream)
   targetGene.downstream <- ifelse(strand == "+", targetGene.tss + downstream, targetGene.tss - downstream)

   tbl.fullRegionAllDNA <- data.frame(chrom=tbl.promoters$chromosome_name[1],
                                      start=min(c(targetGene.upstream, targetGene.downstream)),
                                      end=max(c(targetGene.upstream, targetGene.downstream)),
                                      stringsAsFactors=FALSE)


      # show the whole gene
   showGenomicRegion(tv, targetGene)
   chromLoc <- parseChromLocString(getGenomicRegion(tv))

      # expand view to gene +/- 5kb
   showGenomicRegion(tv, sprintf("%s:%d-%d", chromLoc$chrom, chromLoc$start-5000, chromLoc$end+10000))

      # show promoters in this broad view
   tbl.promoters.bed <- tbl.promoters[, c("chromosome_name", "promoter_start", "promoter_end", "ensembl_transcript_id")]
   addBedTrackFromDataFrame(tv, "promoters", tbl.promoters.bed, color="blue", displayMode="EXPANDED")

      # survey the ATAC-seq regions across theis broad view.
   chromLoc <- parseChromLocString(getGenomicRegion(tv))
   tbl.regions.atac <- calculateATACregions(chrom=chromLoc$chrom, loc.start=chromLoc$start,
                                            loc.end=chromLoc$end, display=FALSE)
   tbl.regions.atac <- unique(tbl.regions.atac[, 1:3])
   addBedTrackFromDataFrame(tv, "ATAC", tbl.regions.atac, color="darkGreen", displayMode="EXPANDED")

      # now zoom in on the proximal promoter region: this is where we will build the model
   with(tbl.promoters.bed, showGenomicRegion(tv, sprintf("%s:%d-%d", chromosome_name[1],
                                                         min(promoter_start) - 500,
                                                         max(promoter_end) + 500)))

   tbl.atacInPromoter <- as.data.frame(GenomicRanges::intersect(GRanges(tbl.regions.atac), GRanges(tbl.promoters.bed)))
   colnames(tbl.atacInPromoter) <- c("chrom", "start", "end", "width", "strand")
   addBedTrackFromDataFrame(tv, "ATAC/promoter", tbl.atacInPromoter, color="magenta", displayMode="EXPANDED")

   pfms.jaspar <- query(query(MotifDb, "jaspar"),     "mmus")
   pfms.jolma  <- query(query(MotifDb, "jolma2013"),  "mmus")
   pfms.j <- as.list(c(pfms.jaspar, pfms.jolma))      # 411

   pfms.mus <- as.list(query(MotifDb, "mmus"))                 # 1251



   mtx.early <- mtx.rna[, 1:19]
   variance.early <- apply(mtx.early, 1, var)
   early.genes.low.variance <- which(variance.early < 0.1)
   if(length(early.genes.low.variance) > 0)
      mtx.early <- mtx.early[-early.genes.low.variance,,drop=FALSE]
   mtx.late  <- mtx.rna[, 20:38]
   variance.late <- apply(mtx.late, 1, var)
   late.genes.low.variance <- which(variance.late < 0.1)
   if(length(late.genes.low.variance) > 0)
      mtx.late <- mtx.late[-late.genes.low.variance,,drop=FALSE]

   models <- list()
   models[["allMousePFMs97"]] <- modelVersion(pfms.mus, tbl.atacInPromoter, 97, mtx.rna)
   models[["allMousePFMs90"]] <- modelVersion(pfms.mus, tbl.atacInPromoter, 90, mtx.rna)
   models[["JMousePFMs90"]] <- modelVersion(pfms.j, tbl.atacInPromoter, 90, mtx.rna)
   models[["allMousePFMs-90-early"]] <- modelVersion(pfms.mus, tbl.atacInPromoter, 90, mtx.early)
   models[["allMousePFMs-90-late"]] <- modelVersion(pfms.mus, tbl.atacInPromoter, 90, mtx.late)
   models[["allDNAmousePFMs90"]] <- modelVersion(pfms.mus, tbl.fullRegionAllDNA, 90, mtx.rna)

   g <- buildMultiModelGraph(tv, targetGene, models)

   xCoordinate.span <- 1500
   g.lo <- addGeneModelLayout(tv, g, xPos.span=xCoordinate.span)
   setGraph(tv, g.lo, names(models))
   setStyle(tv, "style.js")
   fit(tv)

} # explore.one.gene
#----------------------------------------------------------------------------------------------------
if(!exists("tbl.network"))
   tbl.network <- readHamidsTCellNetwork()
