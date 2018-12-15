#------------------------------------------------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#------------------------------------------------------------------------------------------------------------------------
library(RCyjs)
stopifnot(packageVersion("RCyjs") >= "1.7.15")
#------------------------------------------------------------------------------------------------------------------------

setwd("/Volumes/hbolouri/TCE/CyJS/July2017")
# setwd("/Volumes/bolouri_h/user/hbolouri/TCE/CyJS/July2017")

# nodeAnnotationFile <- "nodes.txt"
# networkFile <- "edges.txt"
# nodeAnnotationFile <- "CytoscapeNodeAnnot_11July2017.txt"
# networkFile <- "cytoscapeEdgeAnnot_11July2017.txt"
# nodeAnnotationFile <- "cancerNodeAnnot_17July2017.txt"
# networkFile <- "cancerEdgeAnnot_17July2017.txt"

# tbl.edges <- read.table(networkFile, sep="\t", header=TRUE, as.is=TRUE)
# # tbl.edges <- tbl.edges[ ,-(5:10)]

# tbl.nodes <- read.table(nodeAnnotationFile, sep="\t", header=TRUE, as.is=TRUE)


tbl.edges <- get(load("validatedEdgesTbl.RData"))
tbl.nodes <- get(load("validatedNodesTbl.RData"))

colnames(tbl.nodes)[1] <- "Symbol" #### mis-labeled in cancer data! ####


nodes.in.edges <- sort(unique(c(tbl.edges$Source, tbl.edges$Target)))
nodes.from.nodesFile <- as.character(tbl.nodes$Symbol)
all.nodes <- sort(unique(c(nodes.from.nodesFile, nodes.in.edges)))
g <- graphNEL(all.nodes, edgemode="directed")
nodeDataDefaults(g, attr = "label") <- "default node label"
nodeDataDefaults(g, attr = "nodeType") <- "default node type"
nodeDataDefaults(g, attr = "expression") <- 0.0
edgeDataDefaults(g, attr = "edgeType") <- "undefined"
edgeDataDefaults(g, attr = "score") <- 0.0

source.nodes <- tbl.edges$Source
target.nodes <- tbl.edges$Target
edgeTypes <- tbl.edges$edgeType

g <- addEdge(source.nodes, target.nodes, g)
edgeData(g, source.nodes, target.nodes, "edgeType") <- edgeTypes

tbl.node.types <- data.frame(node=unique(c(source.nodes, target.nodes)),
							type="protein", stringsAsFactors=FALSE)

nodeData(g, tbl.node.types$node, "nodeType") <- tbl.node.types$type

exNodes <- get(load("exhaustionNodes.RData"))
nodeDataDefaults(g, attr="type") <- "undefined"
nodeData(g, exNodes, attr="type") <- "ex"

rcy <- RCyjs(10000:10020, title="Celgene", graph=g)
restoreLayout(rcy, "layout.RData")
fit(rcy)
# httpSetStyle(rcy, "style.js")
httpSetStyle(rcy, "styleCancer_21July2017.js")

savePNG(rcy, "noExpDataNet.png")


expression.conditions <- grep("Exp", colnames(tbl.nodes), value=TRUE)
score.conditions <- grep("Scores", colnames(tbl.edges), value=TRUE)

setNodeAttributes(rcy, "expression", tbl.nodes$Symbol, tbl.nodes[, expression.conditions[1]])

if (length(expression.conditions) == length(score.conditions)) { 
for (i in 1:length(expression.conditions)) {
     expression.condition <- expression.conditions[i]
     score.condition <- score.conditions[i]      
     printf("displaying node condition %s, edge condition %s", expression.condition, score.condition)

     nodeNames <- tbl.nodes$Symbol
     values    <- tbl.nodes[, expression.condition]
     setNodeAttributes(rcy, "expression", nodeNames, values)

     sourceNodes <- tbl.edges$Source
     targetNodes <- tbl.edges$Target
     edgeTypes   <- tbl.edges$edgeType
     values      <- as.numeric(tbl.edges[, score.condition])
     setEdgeAttributes(rcy, "score", sourceNodes, targetNodes, edgeTypes, values)
     redraw(rcy)
	 # httpSetStyle(rcy, paste0("style", i, ".js")) # May need to re-draw after this!
	 savePNG(rcy, paste0(expression.condition, ".png"))
     Sys.sleep(2)
     } # for i
  } # if lengths match



saveLayout(rcy, "layout.RData")


