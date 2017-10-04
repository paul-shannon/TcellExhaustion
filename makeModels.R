library(trena)
library(trenaViz)
library(colorspace)
library(annotate)
library(org.Mm.eg.db)
library(MotifDb)
library(biomaRt)
#--------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#--------------------------------------------------------------------------------
stopifnot(packageVersion("trena")    >= "0.99.182")
stopifnot(packageVersion("trenaViz") >= "0.99.21")
stopifnot(packageVersion("MotifDb")  >= "1.19.13")

if(!exists("trena"))
   trena <- Trena("mm10", quiet=TRUE)

PORT.RANGE <- 8000:8020
tv <- FALSE   # no browser visualization needed

if(!exists("tv")) {
   tv <- trenaViz(PORT.RANGE, quiet=TRUE)
   setGenome(tv, "mm10")
   }

if(!exists("mm10.mart")){
   mm10.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
   }


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
proximalPromoter <- function(targetGene, upstream, downstream)
{
   tbl.geneInfo <- getBM(attributes=c("chromosome_name", "transcription_start_site", "mgi_symbol", "strand"),
                         filters="mgi_symbol", value=targetGene, mart=mm10.mart)

   if(nrow(tbl.geneInfo) == 0)
      return(NA)

      # make sure all transcripts are on the same strand
   strand <- unique(tbl.geneInfo$strand)
   stopifnot(length(strand) == 1)

   chromosome <- sprintf("chr%s", unique(tbl.geneInfo$chromosome_name[1]))
   stopifnot(length(chromosome) == 1)

      # assume + strand
   tss <- min(tbl.geneInfo$transcription_start_site)
   start <- tss - upstream
   end   <- tss + downstream

   if(strand == -1){
      tss <- max(tbl.geneInfo$transcription_start_site)
      start <- tss - downstream
      end   <- tss + upstream
      }

   data.frame(chrom=chromosome, start=start, end=end, stringsAsFactors=FALSE)

} # proximalPromoter
#----------------------------------------------------------------------------------------------------
buildModel <- function(targetGene, upstream, downstream, pfms, motifMatchThreshold, display=FALSE)
{
   tbl.promoter <- proximalPromoter(targetGene, upstream, downstream)

   if(all(is.na(tbl.promoter))){
      printf("no promoter found for gene '%s', skipping", targetGene)
      return(list(model=data.frame(), regulatoryRegions=data.frame()))
      }

      # increase the atac search region in order to catch those which
      # exceed our possibly small proximal promoter
   tbl.regions.atac <- with(tbl.promoter,
                            calculateATACregions(chrom=chrom, loc.start=start-5000,  loc.end=end+5000, display=FALSE))
   if(nrow(tbl.regions.atac) == 0)
     return(list(model=data.frame(), regulatoryRegions=data.frame()))

      # an opportunity here to filter on sample and/or score, i.e., keep only thos
      # regions which have an atac sscore > 100.
      # right now we keep them all, ignoring score, ignoring sample, uniqifying on chrom:start-end

   tbl.candidateRegions <- unique(tbl.regions.atac[, 1:3])
   if(display){
      removeTracksByName(tv, getTrackNames(tv)[-1])
      showGenomicRegion(tv, targetGene)
      #with(tbl.promoter, showGenomicRegion(tv, sprintf("%s:%d-%d", chrom, start-300, end+300)))
      addBedTrackFromDataFrame(tv, "atac", tbl.candidateRegions, color="#C29DDE")
      } # display

   tbl.motifsInRegulatoryRegions  <- findMotifs(pfms, tbl.candidateRegions, motifMatchThreshold)
   printf("%d motifs before subset", nrow(tbl.motifsInRegulatoryRegions))
   tbl.motifs <- subset(tbl.motifsInRegulatoryRegions, motifStart >= tbl.promoter$start &
                                                       motifEnd   <= tbl.promoter$end)
   printf("%d motifs after subset", nrow(tbl.motifs))
   if(display){
      addBedTrackFromDataFrame(tv, "motifs",
                       tbl.motifs[, c("chrom", "motifStart", "motifEnd", "motifName", "motifRelativeScore")],
                       color="red")
      } # display

   solver.names <- c("lasso", "pearson", "randomForest", "ridge", "spearman")
   tbl.geneModel <- createGeneModel(trena, targetGene, solver.names, tbl.motifsInRegulatoryRegions, mtx.rna)
   if(nrow(tbl.geneModel) > 0)
      tbl.geneModel <- tbl.geneModel[order(tbl.geneModel$rf.score, decreasing=TRUE),]

   list(model=tbl.geneModel, regulatoryRegions=tbl.motifsInRegulatoryRegions)

} # buildModel
#----------------------------------------------------------------------------------------------------
buildManyModels <- function (genes.of.interest, upstream, downstream, motifMatchThreshold, display=FALSE)
{
   failures <- c()   # misnamed, lacking atac-seq regions in proximal promoters, no tfs for motifs, ....

   pfms.jaspar <- query(query(MotifDb, "jaspar"), "mmus")
   pfms.jolma  <- query(query(MotifDb, "jolma2013"), "mmus")
   pfms <- as.list(c(pfms.jaspar, pfms.jolma))   # 411

   model.tbls <- list()
   region.tbls <- list()

   for(targetGene in genes.of.interest){
      printf("--- building model for %s (%d, %d, %d)", targetGene, upstream, downstream,
             motifMatchThreshold)
      result <- buildModel(targetGene, upstream, downstream, pfms, motifMatchThreshold, display)
         # two fields in result:
         #   a data.frame for the regulatory model
         #   a data.frame for the regulatory regions
         # add a "targetGene" column to both
      if(nrow(result$model) == 0){
         printf("no regulators found for %s", targetGene);
         failures <- c(failures, targetGene)
         next;
         }
      result$model$targetGene <- targetGene
      result$regulatoryRegions$targetGene <- targetGene
      model.tbls[[targetGene]] <- result$model[order(result$model$rf.score, decreasing=TRUE),]
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
if(!exists("tbl.network"))
   tbl.network <- readHamidsTCellNetwork()


