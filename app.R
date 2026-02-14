

Sys.setenv(RENV_CONFIG_ENABLED = "FALSE")
options(shiny.minified = TRUE)

#Libraries----
library(shiny)
library(shinyjs)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(shinydashboard)
library(shinyWidgets)

# shared app logger (terminal + browser sessions)----
app_log_callbacks <- list()
app_log_callback_id <- 0L

register_app_log_callback <- function(callback) {
    app_log_callback_id <<- app_log_callback_id + 1L
    id <- as.character(app_log_callback_id)
    app_log_callbacks[[id]] <<- callback
    id
}

unregister_app_log_callback <- function(id) {
    app_log_callbacks[[id]] <<- NULL
}

emit_app_log <- function(text) {
    line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text)
    message(line)
    if (length(app_log_callbacks) > 0) {
        for (cb in app_log_callbacks) {
            try(cb(line), silent = TRUE)
        }
    }
    invisible(line)
}

#l ightweight startup data
Gene_desp <- read.csv("DATA/gene_annotations/gene_annotation_sorted.csv")
t3 <- read.csv("DATA/sc_meta.csv")

# load Seurat object at startup
emit_app_log("Loading Islets object...")
.islets_load_time <- system.time({
    islets_obj <- readRDS("DATA/Islets4-slimmed.Rds")
})
emit_app_log(sprintf("Islets loaded in %.2fs", .islets_load_time[["elapsed"]]))
islets_cache <- function() islets_obj

# cache expensive per-gene computations
gene_fetch_cache <- new.env(parent = emptyenv())
gene_dotplot_cache <- new.env(parent = emptyenv())

get_gene_fetch <- function(gene) {
    key <- toupper(gene)
    if (!exists(key, envir = gene_fetch_cache, inherits = FALSE)) {
        emit_app_log(sprintf("FetchData cache miss for %s, fetching...", key))
        elapsed <- system.time({
            islets <- islets_cache()
            dat <- FetchData(islets, vars = key)
            dat$celltypes <- islets$CellTypes
            dat$age <- islets$Age
            assign(key, dat, envir = gene_fetch_cache)
        })
        emit_app_log(sprintf("FetchData for %s done in %.2fs", key, elapsed[["elapsed"]]))
    }
    get(key, envir = gene_fetch_cache, inherits = FALSE)
}

get_gene_dotplot_data <- function(gene) {
    key <- toupper(gene)
    if (!exists(key, envir = gene_dotplot_cache, inherits = FALSE)) {
        emit_app_log(sprintf("DotPlot cache miss for %s, computing...", key))
        elapsed <- system.time({
            islets <- islets_cache()
            dp <- DotPlot(islets, feature = key)
            assign(key, dp$data, envir = gene_dotplot_cache)
        })
        emit_app_log(sprintf("DotPlot for %s done in %.2fs", key, elapsed[["elapsed"]]))
    }
    get(key, envir = gene_dotplot_cache, inherits = FALSE)
}

# precompute static UMAP
static_umap_celltype_plot <- local({
    p <- NULL
    function() {
        if (is.null(p)) {
            islets <- islets_cache()
            p <<- DimPlot(islets, reduction = "umap", label = TRUE, pt.size = 1,
                          cols = c('Alpha' = '#F8766D', 'Beta' = '#39B600', 'Delta' = '#D89000',
                                   'Gamma' = '#A3A500', 'Epsilon' = '#00BF7D', 'Acinar' = '#00BFC4',
                                   'Ductal' = '#00B0F6', 'Endothelial' = '#9590FF',
                                   'Stellate' = '#E76BF3', 'Immune' = '#FF62BC')) +
                NoLegend() +
                theme(text = element_text(size = 18, face = "bold"),
                      axis.text = element_text(size = 18, face = "bold"),
                      axis.line.x = element_line(color = "black", size = 0.8),
                      axis.line.y = element_line(color = "black", size = 0.8))
        }
        p
    }
})

#Single-cell gene expression atlas of human pancreatic islets

#1.Header----
header <- dashboardHeader(titleWidth = "100%",
                          # Set height of dashboardHeader
                          tags$li(class = "dropdown",
                                  tags$style(".main-header {max-height: 200px}"),
                                  tags$style(".main-header .logo {height: 100px}")
                          )
)


#webpage links to the images
anchor <- tags$header(
    tags$a(href='https://www.powersbrissovaresearch.org/',
           tags$img(src='logo-4.png', width='200',style="float:left; margin:0 70px 10px 20px;" )),
    tags$a(href='https://cds.vanderbilt.edu',
           tags$img(src='CDS-logo-600x85.png', width='200',style="float:right; margin-left:70px; margin-top:15px; height:auto;" )),
    'Single cell gene expression atlas of human pancreatic islets',
    style = "color: #2b6cb3;
                           float:left;
                           /*font-family: Avenir Light;*/
                           font-size: 25px;
                           padding:20px;
                           font-weight: bold"
)


header$children[[2]]$children <- tags$div(
    tags$head(tags$style(HTML(".name { background-color: white } Gene-label { font-size:80%;} "))),
    anchor,
    class = 'name')


#2.User Interface----

#*  Dashboard header----
ui<-dashboardPage(header,
                  title = "Single cell gene expression atlas of human pancreatic islets - Powers & Brissova Research Group",
                  skin = "black",
                  #* Dashboard sidebar ----
                  dashboardSidebar(width = 300,

                                   sidebarMenu(
                                       id = "tabs",
                                       menuItem("Home", tabName = "home", selected = T),
                                       menuItem(
                                           startExpanded = TRUE,
                                           selectizeInput(inputId = "Gene",
                                                          label = "Enter Official Gene Symbol",
                                                          choices=NULL
                                           )
                                       ),

                                       menuItem("Violinplot", tabName = "vlnplot", icon = icon("vp")),
                                       menuItem("Umap", tabName = "umap", icon = icon("ump")),
                                        menuItem("Dotplot", tabName = "dotplot", icon = icon("dp")),
                                        menuItem("Expression values", tabName = "cellno", icon = icon("cellno")),
                                        menuItem("Manuscript", icon = icon("Manuscript"), href ="https://insight.jci.org/articles/view/151621"),
                                        menuItem("Experimental Summary", tabName = "expsum", icon = icon("ES"))
                                    )
                  ),

                  #*  Dashboardbody----
                  dashboardBody(
                       useShinyjs(),

                       tags$head(
                           tags$script(HTML("
                               Shiny.addCustomMessageHandler('rlog', function(msg) {
                                 console.log('[R]', msg);
                               });
                           ")),
                           tags$title("Single cell gene expression atlas of human pancreatic islets"),
                           tags$style(HTML('
                                            /* body */
                                            body, .wrapper, .skin-black .wrapper,
                                            .skin-black .content-wrapper,
                                            .skin-black .main-footer,
                                            .skin-black .wrapper {
                                            background-color: #FFFFFF !important;
                                            }
                                            html { background-color: #FFFFFF !important; }
                                            .content-wrapper, .right-side {
                                            background-color: #FFFFFF !important;
                                            max-width: 1400px;
                                            }

                                            /* main sidebar */
                                            .skin-blue .main-sidebar { font-size: 20px;
                                                            background-color: #F0F8FF;
                                            }
                                            .main-sidebar { font-size: 20px; }


                                            .left-side, .main-sidebar {
                                            	padding-top: 110px;
                                            }
                                            section.sidebar .shiny-input-container {
                                            	padding-left: 0 !important;
                                            }

                                            /* Gene-label { font-size:70%;} */

                                            /* image {max-width: 60%; width: 60%; height: auto; } */

                                            /* fix for spinner showing up in right of plots in large monitors
                                            not elegant, but quick fix */
                                            .loading-spinner { left:25% !important;}

                                            header { padding-top:20px 0 0 0 !important;}

                                            '))),



                      #For violinplot tab
                      tabItems(
                          tabItem(tabName = "home",
                                  fluidPage(
                                      verticalLayout(tags$h2("Welcome!"),
                                                     hr(),
                                                     tags$h4("This app provides interactive access to our single cell RNA-Seq data that is reported in:"),
                                                     tags$div(
                                                         HTML("<p style = 'font-size:20px;'><u><b><a href='https://insight.jci.org/articles/view/151621'>Combinatorial transcription factor profiles predict mature and functional human islet α and β cells.</a></u></b><br></p> Shristi Shrestha*, Diane C. Saunders*, John T. Walker*, Joan Camunas-Soler, Xiao-Qing Dai, Rachana Haliyur, Radhika Aramandla, Greg Poffenberger, Nripesh Prasad, Rita Bottino, Roland Stein, Jean-Philippe Cartailler, Stephen C. J. Parker, Patrick E. MacDonald, Shawn E. Levy, Alvin C. Powers, Marcela Brissova, <b><a href='https://insight.jci.org/articles/view/151621'>JCI Insight. 2021 doi: 10.1172/jci.insight.151621</a></b><br> *first co-authors <blockquote style='font-size:15px'> Abstract <br> Islet-enriched transcription factors (TFs) exert broad control over cellular processes in pancreatic α and β cells and changes in their expression are associated with developmental state and diabetes. However, the implications of heterogeneity in TF expression across islet cell populations are not well understood. To define this TF heterogeneity and its consequences for cellular function, we profiled >40,000 cells from normal human islets by scRNA-seq and stratified α and β cells based on combinatorial TF expression. Subpopulations of islet cells co-expressing ARX/MAFB (α cells) and MAFA/MAFB (β cells) exhibited greater expression of key genes related to glucose sensing and hormone secretion relative to subpopulations expressing only one or neither TF. Moreover, all subpopulations were identified in native pancreatic tissue from multiple donors. By Patch-seq, MAFA/MAFB co-expressing β cells showed enhanced electrophysiological activity. Thus, these results indicate combinatorial TF expression in islet α and β cells predicts highly functional, mature subpopulations.</blockquote>
          	                                                      "),
                                                         tags$div(HTML("<b>Data Availability:</b>")),
                                                         tags$h5("Processed data as seurat object now available at",tags$a(href='https://zenodo.org/record/7626110',"zenodo")),
                                                         tags$h5("Raw data available at GEO",tags$a(href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183568',"GSE183568")),
                                                         tags$h5("All scripts to this interactive explorer now available at",tags$a(href='https://github.com/Powers-Brissova-Research-Group/single-cell-app-islets',"Github"))
                                                     ),
                                                     hr(),
                                                     tags$div(
                                                         tags$h4("What's New - February 14, 2026"),
                                                         tags$ul(
                                                             tags$li("Resizable plots on all pages with sidebar controls"),
                                                             tags$li("Faster violin plots - individual cell points off by default, toggle on via sidebar"),
                                                             tags$li("UMAP: square aspect ratio, adjustable point size, viridis color schemes"),
                                                             tags$li("Gene info moved to sidebar on all plot pages"),
                                                             tags$li("Expression values table with adjustable decimal precision"),
                                                             tags$li("Performance logging added"),
                                                             tags$li("Code refactoring and optimizations")
                                                         ),
                                                         tags$h5("More details at ",
                                                                  tags$a(href='https://github.com/Powers-Brissova-Research-Group/single-cell-app-islets#releases',"Github"))
                                                     )
                                      )
                                  )

                          ),
                          tabItem(tabName = "vlnplot",
                                  fluidRow(
                                      column(9,
                                          addSpinner(plotOutput("plot1", height = "auto"), spin = "dots", color = "#2b6cb3"),
                                          plotOutput("plot2", height = "auto")
                                      ),
                                      column(3,
                                          wellPanel(
                                              tags$h4("Plot Options"),
                                              sliderInput("VlnHeight", "Plot height (px):", min = 300, max = 800, value = 400, step = 50),
                                              sliderInput("VlnWidth", "Plot width (px):", min = 400, max = 1000, value = 800, step = 50),
                                              hr(),
                                              checkboxInput("ShowBoxplot", "Show boxplot", value = TRUE),
                                              checkboxInput("ShowDots", "Show individual cells (slower to render)", value = FALSE),
                                              conditionalPanel(
                                                  condition = "input.ShowDots",
                                                  sliderInput("Cellsize",
                                                              "Cell Size:",
                                                              min=0.1,
                                                              max=1,
                                                              value=0.1,
                                                              step = 0.1)
                                              )
                                          ),
                                          wellPanel(
                                              tags$h4("Gene Info"),
                                              tableOutput("geneinfo_vlnplot")
                                          )
                                      )
                                  )
                          ),

                          #For UMAP plot tab
                          tabItem(tabName = "umap",
                                  fluidRow(
                                      column(9,
                                          addSpinner(plotOutput("plot3", width = "auto", height = "auto"), spin = "dots", color = "#2b6cb3"),
                                          addSpinner(plotOutput("plot4", width = "auto", height = "auto"), spin = "dots", color = "#2b6cb3")
                                      ),
                                      column(3,
                                          wellPanel(
                                              tags$h4("Plot Options"),
                                              sliderInput("UmapSize", "Plot size (px):", min = 300, max = 1500, value = 800, step = 100),
                                              sliderInput("UmapPtSize", "Point size:", min = 0.1, max = 3, value = 1, step = 0.1),
                                              selectInput("UmapColorScheme", "Color scheme (feature plot):",
                                                          choices = c("Grey-Red" = "greyred",
                                                                      "Viridis" = "viridis",
                                                                      "Magma" = "magma",
                                                                      "Inferno" = "inferno",
                                                                      "Plasma" = "plasma",
                                                                      "Cividis" = "cividis"),
                                                          selected = "greyred")
                                          ),
                                          wellPanel(
                                              tags$h4("Gene Info"),
                                              tableOutput("geneinfo_umap")
                                          )
                                      )
                                  )
                          ),

                          #For dotplot tab
                          tabItem(tabName = "dotplot",
                                  fluidRow(
                                      column(9,
                                          addSpinner(plotOutput("plot5", height = "auto"), spin = "dots", color = "#2b6cb3")
                                      ),
                                      column(3,
                                          wellPanel(
                                              tags$h4("Plot Options"),
                                              sliderInput("DotSize", "Plot size (px):", min = 400, max = 1200, value = 800, step = 50)
                                          ),
                                          wellPanel(
                                              tags$h4("Gene Info"),
                                              tableOutput("geneinfo_dotplot")
                                          )
                                      )
                                  )
                          ),

                          #For table on expression counts
                          tabItem(tabName = "cellno",
                                  fluidRow(
                                      column(9,
                                          DT::dataTableOutput("table2")
                                      ),
                                      column(3,
                                          wellPanel(
                                              tags$h4("Table Options"),
                                              sliderInput("DecimalPoints", "Decimal places:", min = 0, max = 10, value = 2, step = 1)
                                          ),
                                          wellPanel(
                                              tags$h4("Gene Info"),
                                              tableOutput("geneinfo_cellno")
                                          )
                                      )
                                  )
                          ),

                          #Add Manuscript link
                          tabItem(tabName = "manuscript",
                                  h2(plotOutput("manuscript"))
                          ),

                          tabItem(tabName = "expsum",
                                  fluidPage(
                                      verticalLayout(br(),
                                                     tags$h3(HTML(paste0("<b>","Single Cell RNA-seq Metadata","</b>")) ),
                                                     tableOutput("table3"),
                                                     tags$h6("Metadata format standardized according to",tags$a(href='https://www.nature.com/articles/s41587-020-00744-z',"Fullgrabe et al.,2020"))

                                      )
                                  )
                          )



                      )))


#3.Server----
server<-function(input, output,session)

    #* table 1 -Gene Description----
{
    append_runtime_log <- function(line) {
        session$sendCustomMessage("rlog", line)
    }

    callback_id <- register_app_log_callback(append_runtime_log)
    session$onSessionEnded(function() {
        unregister_app_log_callback(callback_id)
    })

    initial_load <- reactiveVal(TRUE)

    observe({

        updateSelectizeInput(
            session,
            'Gene',
            choices = Gene_desp$hgnc_symbol,
            server = TRUE,
            selected = "INS"
        )


        # if user interacts with gene filter and is not on umap/vlnplot/dotplot page, then take
        # them to vlnplot, otherwise do not change tab
        observeEvent(input$Gene, {
            emit_app_log(sprintf("Gene selected: %s", toupper(input$Gene)))

            # skip redirect on initial default gene load so we stay on landding  page
            if (initial_load()) {
                initial_load(FALSE)
                return()
            }

            if (input$tabs == "home" || input$tabs == "cellno" || input$tabs == "expsum" || input$tabs == "cellno") {
                # it requires an ID of sidebarMenu (in this case)
                updateTabsetPanel(session, inputId="tabs", selected="vlnplot")
            }

        }, ignoreInit = TRUE) # ignore initial load, since nothing is actually clicked

    })

    gene_info_transposed <- reactive({
        req(input$Gene)
        row <- Gene_desp %>% filter(hgnc_symbol == input$Gene)
        if (nrow(row) == 0) return(data.frame(Field = character(), Value = character()))
        data.frame(Field = names(row), Value = unlist(row[1, ], use.names = FALSE))
    })

    output$geneinfo_vlnplot <- renderTable({ gene_info_transposed() }, colnames = FALSE)
    output$geneinfo_umap <- renderTable({ gene_info_transposed() }, colnames = FALSE)
    output$geneinfo_dotplot <- renderTable({ gene_info_transposed() }, colnames = FALSE)
    output$geneinfo_cellno <- renderTable({ gene_info_transposed() }, colnames = FALSE)

    #* plot 1 (Violinplot by Cell types)
    output$plot1 <- renderPlot({

        req(input$Gene)
        emit_app_log(sprintf("Rendering violin plot (cell types) for %s", toupper(input$Gene)))
        elapsed <- system.time({
        celltype_data <- get_gene_fetch(input$Gene)

        p <- ggplot(celltype_data, aes(x = celltypes, y = celltype_data[,1], fill = celltypes)) +
            geom_violin(trim = T, alpha = 1, width =1, size=1, scale="width") +
            theme_minimal()+
            labs(x="", y= "ln(UMI -per-10,000 +1)", title = toupper(input$Gene))+
            { if (input$ShowDots) geom_jitter(width = 0.2, size = input$Cellsize, alpha = 0.6, height = 0) } +
            { if (input$ShowBoxplot) geom_boxplot(width = 0.1,size=1, outlier.shape = NA, alpha = 0, na.rm = TRUE, position = position_dodge(width = 1), color = "#5A5A5A") } +
            scale_fill_manual(values = c('Alpha'='#F8766D','Beta'='#39B600','Delta'='#D89000','Gamma'='#A3A500','Epsilon'='#00BF7D','Acinar'='#00BFC4','Ductal'='#00B0F6','Endothelial'='#9590FF','Stellate'='#E76BF3','Immune'='#FF62BC'))+
            guides(fill = guide_legend(title = "Cell types"))+  # Add legend title
            theme(text = element_text(size = 18, face="bold"),
                  axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
                  axis.text = element_text(size = 18, face="bold"),
                  axis.line.x = element_line(color="black", size = 0.8),
                  axis.line.y = element_line(color="black", size = 0.8))
        })
        emit_app_log(sprintf("Violin plot (cell types) for %s done in %.2fs", toupper(input$Gene), elapsed[["elapsed"]]))
        p
    }, height = function() input$VlnHeight, width = function() input$VlnWidth)


    #* plot2 (Violin plot by Age)
    output$plot2 <- renderPlot({
        req(input$Gene)
        emit_app_log(sprintf("Rendering violin plot (age) for %s", toupper(input$Gene)))
        elapsed <- system.time({
        age_data <- get_gene_fetch(input$Gene)

        p <- ggplot(age_data, aes(x = age, y = age_data[,1], fill = age)) +
            geom_violin(trim = T, alpha = 0.7, width =1, size=1, scale="area") +
            theme_minimal()+
            labs(x="", y= "ln(UMI -per-10,000 +1)", title = toupper(input$Gene))+
            { if (input$ShowDots) geom_jitter(width = 0.2, size = input$Cellsize, alpha = 0.6, height = 0) } +
            { if (input$ShowBoxplot) geom_boxplot(width = 0.1,size=1, outlier.shape = NA, alpha = 0, na.rm = TRUE, position = position_dodge(width = 1), color = "#5A5A5A") } +
            scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"))+
            guides(fill = guide_legend(title = "Donor age"))+  # Add legend title
            theme(text = element_text(size = 18, face="bold"),
                  axis.text = element_text(size = 18, face="bold"),
                  axis.line.x = element_line(color="black", size = 0.8),
                  axis.line.y = element_line(color="black", size = 0.8))
        })
        emit_app_log(sprintf("Violin plot (age) for %s done in %.2fs", toupper(input$Gene), elapsed[["elapsed"]]))
        p
    }, height = function() input$VlnHeight, width = function() input$VlnWidth)


    #* plot3 (umap-celltypes)
    output$plot3 <- renderPlot({
        emit_app_log("Rendering UMAP (cell types)")
        elapsed <- system.time({
            islets <- islets_cache()
            p <- DimPlot(islets, reduction = "umap", label = TRUE, pt.size = input$UmapPtSize,
                          cols = c('Alpha' = '#F8766D', 'Beta' = '#39B600', 'Delta' = '#D89000',
                                   'Gamma' = '#A3A500', 'Epsilon' = '#00BF7D', 'Acinar' = '#00BFC4',
                                   'Ductal' = '#00B0F6', 'Endothelial' = '#9590FF',
                                   'Stellate' = '#E76BF3', 'Immune' = '#FF62BC')) +
                NoLegend() +
                coord_fixed(ratio = 1) +
                theme(text = element_text(size = 18, face = "bold"),
                      axis.text = element_text(size = 18, face = "bold"),
                      axis.line.x = element_line(color = "black", size = 0.8),
                      axis.line.y = element_line(color = "black", size = 0.8))
        })
        emit_app_log(sprintf("UMAP (cell types) done in %.2fs", elapsed[["elapsed"]]))
        p
    }, height = function() input$UmapSize, width = function() input$UmapSize)


    #* plot 4 (umap)
    output$plot4 <- renderPlot({
        req(input$Gene)
        emit_app_log(sprintf("Rendering UMAP (feature) for %s", toupper(input$Gene)))
        elapsed <- system.time({
        islets <- islets_cache()
        color_scheme <- input$UmapColorScheme
        p <- FeaturePlot(islets, features = toupper(input$Gene), pt.size = input$UmapPtSize,
                         cols = c("lightgrey", "red")) +
            coord_fixed(ratio = 1) +
            { if (color_scheme != "greyred") scale_color_viridis_c(option = color_scheme) } +
            theme(text = element_text(size = 18,face="bold"),
                  axis.text = element_text(size = 18,face="bold"),
                  axis.line.x = element_line(color="black", size = 0.8),
                  axis.line.y = element_line(color="black", size = 0.8))
        })
        emit_app_log(sprintf("UMAP (feature) for %s done in %.2fs", toupper(input$Gene), elapsed[["elapsed"]]))
        p
    }, height = function() input$UmapSize, width = function() input$UmapSize)



    #* plot5 (Dotplot) ----
    Known.markers <- c("GCG","INS","SST","PPY","GHRL", "PRSS1", "KRT19","PECAM1", "PDGFRB", "HLA-DRA")
    output$plot5 <- renderPlot({
        req(input$Gene)
        emit_app_log(sprintf("Rendering dotplot for %s", toupper(input$Gene)))
        elapsed <- system.time({
        islets <- islets_cache()
        selected_markers <- Known.markers[!(Known.markers %in% input$Gene)]
        p <- DotPlot(islets, feature=c(selected_markers,toupper(input$Gene)), dot.scale = 8)+
            coord_flip()+
            scale_color_gradient2(low = "blue", high = "red",mid = "white")+
            labs(y="", x="Genes")+RotatedAxis()+
            #geom_point(shape = 21) +
            scale_y_discrete(position = "right") +
            theme_light() +
            guides(x =  guide_axis(angle = 45))+
            guides(color = guide_colorbar(title = 'Mean expression z-score'),size = guide_legend("% Cells expressed"))+
            theme(text = element_text(size = 15,face="bold"),
                  axis.text = element_text(size = 15,face="bold"))
        })
        emit_app_log(sprintf("Dotplot for %s done in %.2fs", toupper(input$Gene), elapsed[["elapsed"]]))
        p
    }, height = function() round(input$DotSize * 0.8), width = function() input$DotSize)


    #* Table2 (Expression Values) ----
    output$table2 <- DT::renderDataTable({
        req(input$Gene)
        dp <- input$DecimalPoints
        emit_app_log(sprintf("Rendering expression values table for %s", toupper(input$Gene)))
        elapsed <- system.time({
        tble <- get_gene_dotplot_data(input$Gene)
        tble <- tble %>%rename(Gene = features.plot,
                               CellType = id,
                               'Mean expression z-score'= avg.exp.scaled,
                               '% Cells expressed'= pct.exp,
                               'Average expression'= avg.exp)
        tble <- tble[ , c(3, 4, 2,1,5)]
        tble <- tble %>% filter_all(any_vars(. %in% toupper(input$Gene)))
        tble <- tble %>% mutate(across(where(is.numeric), ~round(.x, dp)))
        })
        emit_app_log(sprintf("Expression values table for %s done in %.2fs", toupper(input$Gene), elapsed[["elapsed"]]))
        tble
    })


    #* Table3 (Experimental Summary) ----
    output$table3 <- renderTable({
        t3
    })
}



# shinyApp----
shinyApp(ui=ui,server=server)
