library(cyjShiny)
#----------------------------------------------------------------------------------------------------
readHamidsTCellNetwork <- function()
{
   tbl <- read.table("cytoscapeEdgeAnnot_11July2017.txt", sep="\t", header=TRUE, as.is=TRUE)
   tbl.xtab <- as.data.frame(table(c(tbl$Source, tbl$Target)))
   tbl.xtab <- tbl.xtab[order(tbl.xtab$Freq, decreasing=TRUE),]
   head(tbl.xtab)
      # DNMT3A is the most connected gene, a methyl transferase.
      # extended location: chr12:3,751,728-3,970,655  (218kb)

   colnames(tbl)[1:2] <- c("source", "target")
   tbl

} # readHamidsTCellNetwork
#----------------------------------------------------------------------------------------------------
styleList <- c("",
               "cancer 21jul2017" = "received.15dec2018/styleCancer_21July2017.js",
               "generic style" = "style.js")

condition <- c("", "gal1RGexp", "gal4RGexp", "gal80Rexp")
tbl.edges <- readHamidsTCellNetwork()
tbl.dumb <- data.frame(source="A", target="B", interaction="simple", stringsAsFactors=FALSE)
tbl.dumb <- tbl.edges[1:3, c("source", "target", "interaction")]
graph <- dataFramesToJSON(tbl.edges)
all.nodes <- sort(unique(c(tbl.edges$source, tbl.edges$target)))
load("received.15dec2018/layout.RData")
tbl.hamidsCuratedLayout <- tbl.layout
#----------------------------------------------------------------------------------------------------
ui = shinyUI(fluidPage(

  #tags$head(
  #        tags$link(rel = "stylesheet", type = "text/css",
  #                  href = "http://maxcdn.bootstrapcdn.com/font-awesome/4.2.0/css/font-awesome.min.css")),

  sidebarLayout(
      sidebarPanel(
          selectInput("loadStyleFile", "Select Style: ", choices=styleList),
          selectInput("doLayout", "Select Layout:",
                      choices=c("",
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

          selectInput("setNodeAttributes", "Select Condition:", choices=condition),
          selectInput("selectName", "Select Node by ID:", choices = c("", all.nodes)),
          actionButton("sfn", "Select First Neighbor"),
          actionButton("fit", "Fit Graph"),
          actionButton("fitSelected", "Fit Selected"),
          actionButton("clearSelection", "Deselect Nodes"),
          HTML("<br>"),
          actionButton("loopConditions", "Loop Conditions"),
          HTML("<br>"),
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

    observeEvent(input$setNodeAttributes, ignoreInit=TRUE, {
       attribute <- "lfc"
       expression.vector <- switch(input$setNodeAttributes,
                                   "gal1RGexp" = tbl.mrna$gal1RGexp,
                                   "gal4RGexp" = tbl.mrna$gal4RGexp,
                                   "gal80Rexp" = tbl.mrna$gal80Rexp)
       setNodeAttributes(session, attributeName=attribute, nodes=yeastGalactodeNodeIDs, values=expression.vector)
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
       if(strategy == "Hamid's curated layout"){
          setNodePositions(session, tbl.hamidsCuratedLayout)
       } else {
          doLayout(session, strategy)
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
