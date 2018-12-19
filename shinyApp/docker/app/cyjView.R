library(cyjShiny)
library(later)
#----------------------------------------------------------------------------------------------------
tbl.nodes <- get(load("validatedNodesTbl.RData"))
tbl.edges <- get(load("validatedEdgesTbl.RData"))
colnames(tbl.edges)[1:4] <- c("source", "target", "reference", "interaction")
tbl.edges$score <- 0

expression.conditionNames <- grep("Exp", colnames(tbl.nodes), value=TRUE)
score.conditionNames <- grep("Scores", colnames(tbl.edges), value=TRUE)
#----------------------------------------------------------------------------------------------------
node.conditions <- c("avExp.aPDL1.v.TIL",
                     "avExp.PDL1plusPDL1minus",
                     "avExp.Tim3PlusTim3Minus",
                     "avExp.TILN2CMV",
                     "avExp.TILN2naive",
                     "avExp.TILN2EBV",
                     "avExp.TILN2tumorPBMC",
                     "avExp.acuteD5.2N.89307",
                     "avExp.acuteD7.2N.89307",
                     "avExp.T.D5.2N.89307",
                     "avExp.T.D7.2N.89307",
                     "avExp.T.D14.2N.89307",
                     "avExp.T.D21.2N.89307")

edge.conditions = c("aPDL1TILScores",
                    "PDL1plusPDL1minusScores",
                    "Tim3Plus2Tim3MinusScores",
                    "TILN2CMVScores",
                    "TILN2naiveScores",
                    "TILN2EBVScores",
                    "TILN2tumorPBMCScores",
                    "acuteD5.2NScores",
                    "acuteD7.2NScores",
                    "T.D5.2NScores",
                    "T.D7.2NScores",
                    "T.D14.2NScores",
                    "T.D21.2NScores")

condition.names <- c("aPDL1 vs TIL",
                     "PD1+ vs PD1-",
                     "TIM3+vs TIM3-",
                     "TILN vs CMV",
                     "TILN vs naive",
                     "TILN vs EBV",
                     "TILN vs tumorPBMC",
                     "Acute Day 5 vs naive",
                     "Acute Day 7 vs naive",
                     "TIL Day 5 vs naive",
                     "TIL Day 7 vs naive",
                     "TIL Day 14 vs naive",
                     "TIL Day 21 vs naive")

tbl.conditions <- data.frame(node=node.conditions, edge=edge.conditions, stringsAsFactors=FALSE)
rownames(tbl.conditions) <- condition.names
stopifnot(all(tbl.conditions$node %in% colnames(tbl.nodes)))
stopifnot(all(tbl.conditions$edge %in% colnames(tbl.edges)))

# readHamidsTCellNetwork <- function()
# {
#    tbl <- read.table("cytoscapeEdgeAnnot_11July2017.txt", sep="\t", header=TRUE, as.is=TRUE)
#    tbl.xtab <- as.data.frame(table(c(tbl$Source, tbl$Target)))
#    tbl.xtab <- tbl.xtab[order(tbl.xtab$Freq, decreasing=TRUE),]
#    head(tbl.xtab)
#       # DNMT3A is the most connected gene, a methyl transferase.
#       # extended location: chr12:3,751,728-3,970,655  (218kb)
#
#    colnames(tbl)[1:2] <- c("source", "target")
#    tbl
#
# } # readHamidsTCellNetwork
#----------------------------------------------------------------------------------------------------
styleList <- c("",
               "cancer 21jul2017" = "styleCancer_21July2017.js",
               "generic style"    = "genericStyle.js")

graph <- dataFramesToJSON(tbl.edges, tbl.nodes)
all.nodes <- sort(unique(c(tbl.edges$source, tbl.edges$target)))
tbl.hamidsCuratedLayout  <- get(load("hamidsCuratedLayout.RData"))
#----------------------------------------------------------------------------------------------------
ui = shinyUI(fluidPage(

  sidebarLayout(
      sidebarPanel(
          selectInput("loadStyleFile", "Select Style: ", choices=styleList),
          selectInput("doLayout", "Select Layout:",
                      choices=c(" ",
                                "Hamid's curated layout",
                                "cose",
                                "cola",
                                "circle",
                                "concentric",
                                "breadthfirst",
                                "grid",
                                "random",
                                "dagre",
                                "cose-bilkent")),

          selectInput("setNodeAndEdgeAttributes", "Select Condition:",
                      choices=c("", rownames(tbl.conditions))),
          selectInput("selectName", "Select Node by ID:", choices = c("", all.nodes)),
          actionButton("sfn", "Select First Neighbor"),
          actionButton("fit", "Fit Graph"),
          actionButton("fitSelected", "Fit Selected"),
          actionButton("clearSelection", "Deselect Nodes"),
          HTML("<br>"),
          #actionButton("loopConditions", "Loop Conditions"),
          #HTML("<br>"),
          actionButton("getSelectedNodes", "Get Selected Nodes"),
          HTML("<br><br>"),
          htmlOutput("selectedNodesDisplay"),
          width=2
      ),
      mainPanel(cyjShinyOutput('cyjShiny', height=800), width=10
      )
  ) # sidebarLayout
))
#----------------------------------------------------------------------------------------------------
server = function(input, output, session)
{
    observeEvent(input$fit, ignoreInit=TRUE, {
       fit(session, 80)
       })

    observeEvent(input$setNodeAndEdgeAttributes, ignoreInit=TRUE, {
       condition.name <- input$setNodeAndEdgeAttributes
       if(condition.name != ""){
          node.condition.name <- tbl.conditions[condition.name, "node"]
          edge.condition.name <- tbl.conditions[condition.name, "edge"]
          node.vector <- tbl.nodes[, node.condition.name]
          edge.vector <- tbl.edges[, edge.condition.name]
          attribute <- "expression"
          node.ids <- tbl.nodes$Source
          setNodeAttributes(session, attributeName=attribute, nodes=node.ids, values=as.numeric(node.vector))
          setEdgeAttributes(session,
                            attributeName="score",
                            sourceNodes=tbl.edges$source,
                            targetNodes=tbl.edges$target,
                            interactions=tbl.edges$interaction,
                            values=as.numeric(edge.vector))
          }
       })

    observeEvent(input$loadStyleFile,  ignoreInit=TRUE, {
      if(input$loadStyleFile != "")
         tryCatch({
            loadStyleFile(input$loadStyleFile)
            }, error=function(e){
                 printf("loadStyleFile warning: %s", e$message)
                 showModal(modalDialog(title="cyjShiny error", e$message))
                 }
               )
       })

    observeEvent(input$doLayout,  ignoreInit=TRUE,{
       strategy <- input$doLayout
       if(strategy != " "){
          if(strategy == "Hamid's curated layout"){
             setNodePositions(session, tbl.hamidsCuratedLayout)
          } else {
             doLayout(session, strategy)
	          }
          later(function() {updateSelectInput(session, inputId="doLayout", selected=" ")}, 1);
          }
       })

    observeEvent(input$selectName,  ignoreInit=TRUE,{
       session$sendCustomMessage(type="selectNodes", message=list(input$selectName))
       })

    observeEvent(input$sfn,  ignoreInit=TRUE,{
       session$sendCustomMessage(type="sfn", message=list())
       })

    observeEvent(input$fitSelected,  ignoreInit=TRUE,{
       fitSelected(session, 100)
       })

    observeEvent(input$getSelectedNodes, ignoreInit=TRUE, {
       output$selectedNodesDisplay <- renderText({" "})
       getSelectedNodes(session)
       })

    observeEvent(input$clearSelection,  ignoreInit=TRUE, {
       session$sendCustomMessage(type="clearSelection", message=list())
       })

    observeEvent(input$loopConditions, ignoreInit=TRUE, {
       condition.names <- c("gal1RGexp", "gal4RGexp", "gal80Rexp")
       for(condition.name in condition.names){
          expression.vector <- tbl.mrna[, condition.name]
          setNodeAttributes(session, attributeName="lfc", nodes=yeastGalactodeNodeIDs, values=expression.vector)
          Sys.sleep(1)
          } # for condition.name
       updateSelectInput(session, "setNodeAttributes", selected="gal1RGexp")
       })

    observeEvent(input$selectedNodes, {
       newNodes <- input$selectedNodes;
       output$selectedNodesDisplay <- renderText({
          paste(newNodes)
          })
       })

    output$value <- renderPrint({input$action})

    output$cyjShiny <- renderCyjShiny(
       cyjShiny(graph, "cola", "yeastGalactoseStyle.js")
       )

} # server
#----------------------------------------------------------------------------------------------------
app <- shinyApp(ui = ui, server = server)
