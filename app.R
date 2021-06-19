library(shiny)
library(ape)
library(phytools)


trimHS.maxC <- function (N, HS, n, check.unique= FALSE) {
    trim.int <- function (x, HS, n) {
        HS.LUT <- which(HS == 1, arr.in= TRUE)
        HS.LUT <- cbind(HS.LUT, 1:nrow(HS.LUT))
        df <- as.data.frame(HS.LUT)
        hs.lut <- subset(df[sample(nrow(df)), ],
                         !duplicated(row) & !duplicated(col))
        if (nrow(hs.lut) < n) hs <- NULL else {
            hs.lut <- hs.lut[sample(nrow(hs.lut), n), ]
            hs <- diag(nrow(hs.lut))
            rownames(hs) <- rownames(HS[hs.lut[ ,1], ])
            colnames(hs) <- colnames(HS[ ,hs.lut[ ,2]])
            return(hs)
        }
    }
    trim.HS <- lapply(1:N, trim.int, HS= HS, n= n )
    if (check.unique == TRUE) trim.HS <- unique(trim.HS)
    if (length(trim.HS) < N) warning("No. of trimmed H-S assoc. matrices < No. of runs")
    return(trim.HS)
}

One2one.f <- function (hs, reps= 1e+4) {
    HS.LUT <- which(hs ==1, arr.in=TRUE)
    HS.LUT <- cbind(HS.LUT,1:nrow(HS.LUT))
    df <- as.data.frame(HS.LUT)
    V <- rep(NA,reps)
    for(i in 1:reps){
        hs.lut <- subset(df[sample(nrow(df)),],
                         !duplicated(row) & !duplicated(col))
        n <- sum(HS)
        while (n >0) {
            n <- n-1;
            if (nrow(hs.lut) == n) break
        }
        V[i]<- n
    }
    V <- min(V)
    return(V)
}

ui <- fluidPage(
    titlePanel("Congruence analisis using Random TaPas"),
    sidebarLayout(
        sidebarPanel(
            radioButtons("type", label = "Select the format of the phylogenies:",
                         choices = c("Newick", "Nexus")),
            fileInput("phylo1", label = "Select a file with the hosts phylogeny:"),
            fileInput("mphylo1", label = "Select a file with the hosts posterior probability trees:"),
            #selectInput("outgroup1", label = "outgroup of hosts phylogeny:", choices = ""),
            fileInput("phylo2", label = "Select another file with the symbionts phylogeny:"),
            fileInput("mphylo2", label = "Select a file with the symbionts posterior probability trees:"),
            #selectInput("outgroup2", label = "outgroup of symbionts phylogeny:", choices = ""),
            fileInput("relmat", label = "Select a file with an association matrix:"),
            radioButtons("globalfit", label = "Global fit statistics", choices = c("PACo", "PACo and GeoD"), choiceValues = c("PACo","PAGD"),
                         selected = "PACo"),
            sliderInput("reps", label = "Number of runs in Random Tapas:", min = 1e+3, max = 1e+4, step = 1e+3, value = 1e+3),
            radioButtons("dupli", label = "Remove duplicate associations:",
                         choiceNames = c("Yes", "No"), choiceValues = c(TRUE, FALSE), selected = TRUE),
            sliderInput("n", label = "Number of unique associations:", min = 1, max = 20, step = 1, value = 5),
            radioButtons("paco.ss", label = "Procrustes superposition:",
                         choiceNames = c("Symmetric", "Asymmetric"), choiceValues = c(TRUE, FALSE), selected = FALSE),
            radioButtons("paco.warn", label = "show PACo's warnings:",
                         choiceNames = c("Yes", "No"), choiceValues = c(TRUE, FALSE), selected = FALSE),
            selectInput("paco.ei.corr", label = "Eigenvalues correction:", 
                        choices = list("None"="none", "Lingoes"="lingoes", "Cailliez"="cailliez", "Element‐wise square‐root of the phylogenetic distances"="sqrt.D"), selected = "sqrt.D"),
            numericInput("percentile", label = "Percentile used to compute the frequency of each h-s association:",value = 0.01,
                        min = 0.0, max = 1, step = 0.01),
            radioButtons("link.res.fq", label = "Compute the residual frequency values for each h-s association:",
                         choiceNames = c("Yes", "No"), choiceValues = c(TRUE, FALSE), selected = FALSE),
            radioButtons("below.p", label = "Compute the frequency for the h-s associations above the percentile threshold:",
                         choiceNames = c("Yes", "No"), choiceValues = c(TRUE, FALSE), selected = FALSE),
            HTML('<p style="font-size:20px"><b>GRAPHICAL OPTIONS<b><p>'),
            radioButtons("colscale", label = "Color scale:", choiceNames = c("Diverging", "Sequential"), 
                         choiceValues = c("diverging", "sequential"), selected = "sequential"),
            sliderInput("n.breaks", label = "Number of shades of color:", min = 3, max = 100, step = 1, value = 50),
            radioButtons("node.tag", label = "Plot fast maximum likelihood estimators of ancestral states of each node:", choiceNames = c("Yes", "No"), 
                         choiceValues = c(TRUE, FALSE), selected = TRUE),
            HTML('<p style="font-size:20px"><b>Generate a report in HTML format, (It may take a while)<b><p>'),
            downloadButton("report.download", label = "Generate report")
            
        ),
        mainPanel(plotOutput("n_plot"),
                  textOutput("max_n"),
                  plotOutput("tangle", height = "auto"))
        
)
)
    
server <- function(input, output, session) {
    
    
    tree1<-reactive({
        file1<-input$phylo1
        file1<-file1$datapath
        if (!is.null(file1)) {
            if(input$type=="Newick"){
                if(class(try(read.tree(file = file1)))!="try-error"){
                    read.tree(file = file1)}else{NULL}
            }else{
                if(class(try(read.nexus(file = file1)))!="try-error"){
                    read.nexus(file = file1)}else{NULL}
            }
            
        }
    })
    tree1label<-reactive({
        tree<-tree1()
        tree$tip.label
    })
    #observe({
    #    updateSelectInput(session, "outgroup1", label = "outgroup of phylogeny 1:", choices = tree1label())
    #})
    tree2<-reactive({
        file2<-input$phylo2
        file2<-file2$datapath
        if (!is.null(file2)) {
            if(input$type=="Newick"){
                if(class(try(read.tree(file = file2)))!="try-error"){
                read.tree(file = file2)}else{NULL}
            }else{
                if(class(try(read.nexus(file = file2)))!="try-error"){
                    read.nexus(file = file2)}else{NULL}
            }
            
        }
    })
    
    HS<-reactive({
        path2HS<-input$relmat
        path2HS<-path2HS$datapath
        if(!is.null(path2HS)){
            HS<-read.table(path2HS,row.names = 1, header = T)
            as.matrix(HS)
        }else{NULL}
        
    })   
    tree2label<-reactive({
        tree<-tree2()
        tree$tip.label
    })
    
    #observe({
    #    updateSelectInput(session, "outgroup2", label = "outgroup of phylogeny 2:", choices = tree2label())
    #})
    
    output$max_n<-renderText({
        if(!is.null(HS())){
            n_max<-sum(HS())
            sprintf("The number of associations in A is %i", n_max)
        }else{NULL}
    })
    
    output$n_plot<-renderPlot({
        if(!is.null(HS())){
            hs<-HS()
            X <- 1:12
            Y <- rep(NA, length(X))
            for(i in 1:length(X)) {
                THS <- trimHS.maxC(input$reps, HS = hs, n=X[i], check.unique=TRUE)
                Y[i] <- length(THS)  
        }
            plot(X,Y, type="b", xlab="Number of unique H-S associations",
                 ylab="Number of runs accomplished")
            axis(side = 1, at = X,labels = T)
        }
    })
    
    output$tangle<-renderPlot({
        TreeH<-tree1()
        TreeS<-tree2()
        HS<-HS()
        HS.lut <- which(HS ==1, arr.ind=TRUE)
        linkhs <- cbind(rownames(HS)[HS.lut[,1]], colnames(HS)[HS.lut[,2]])
        obj <- cophylo(TreeH,TreeS, linkhs)
        par()
        plot.cophylo(obj, link.type="curved")
    }, height = function() {
        session$clientData$output_tangle_width
    })
    
    output$report.download <- downloadHandler(
           filename = "report.html",
           content = function(file) {
               phylo1<-input$phylo1
               phylo1<-phylo1$datapath
               phylo2<-input$phylo2
               phylo2<-phylo2$datapath
               mphylo1<-input$mphylo1
               mphylo1<-mphylo1$datapath
               mphylo2<-input$mphylo2
               mphylo2<-mphylo2$datapath
               relmat<-input$relmat
               relmat<-relmat$datapath
               params = list(type = input$type, phylo1 = phylo1,
                             mphylo1 = mphylo1, phylo2 = phylo2,
                             mphylo2 = mphylo2, relmat = relmat,
                             reps = input$reps, dupli = input$dupli, n = input$n,
                             paco.ss = input$paco.ss, paco.warn = input$paco.warn,
                             paco.ei.corr = input$paco.ei.corr, percentile = input$percentile,
                             link.res.fq = input$link.res.fq, below.p = input$below.p,
                             colscale = input$colscale, n.breaks = input$n.breaks,
                             node.tag = input$node.tag
                             )
               if(input$globalfit=="PACo"){
             rmarkdown::render(input = "ReportPACo.Rmd", params = params,
                               envir = new.env(parent = globalenv()), 
                               output_file = file)}else{
                                   rmarkdown::render(input = "ReportPACoGD.Rmd", params = params,
                                                     envir = new.env(parent = globalenv()), 
                                                     output_file = file)                       
                               }
           }
         )
    
}

shinyApp(ui, server)

