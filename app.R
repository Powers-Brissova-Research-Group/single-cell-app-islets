


#Libraries----
library(shiny)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(shinythemes)
library(shinydashboard)
library(shinyBS)
library(dittoSeq)

#Might need this to deploy----
library(BiocManager)
options(repos = BiocManager::repositories())
#options(repos = BiocInstaller::biocinstallRepos())
#options(repos = c("CRAN" = "https://cran.rstudio.com/", "BioCsoft" = "https://bioconductor.org/packages/3.8/bioc", "BioCann" = "https://bioconductor.org/packages/3.8/data/annotation"))
#getOption("repos")

#load data----
load("DATA/Islets2.Rda")

#Fetch gene id, gene description info from biomart
# ensembl = useMart(
#   "ensembl", 
#   host = "uswest.ensembl.org",
#   dataset = "hsapiens_gene_ensembl" )
# Gene_desp <- getBM( attributes = c("ensembl_gene_id",'hgnc_symbol', "description",'gene_biotype', 'chromosome_name', 'start_position', 'end_position','source'),
#                   filters = "hgnc_symbol",
#                   values = rownames(Islets),
#                   mart = ensembl,
#                   useCache = FALSE)
#write.csv(Gene_desp, "DATA/biomart_annotation.csv")

#saved the output to csv
Gene_desp<-read.csv("DATA/biomart_annotation.csv",header = TRUE, row.names = 1)


#Single-cell gene expression atlas of human pancreatic islets

#1.Header----
header <- dashboardHeader(titleWidth = "100%",
  # Set height of dashboardHeader
  tags$li(class = "dropdown",
          tags$style(".main-header {max-height: 100px}"),
          tags$style(".main-header .logo {height: 60px}")))


#add powerslab webpage link to the powerslab image
anchor <- tags$li(a(href='https://www.powersbrissovaresearch.org/',
                 tags$img(src='logo-4.png', height='60', width='200'),
                 'Single cell gene expression atlas of human pancreatic islets',
                 style = "color: #4292C6;
                          font-family: Avenir Light;
                          font-size: 40px;
                          padding-top:10px; padding-bottom:10px;
                          font-weight: bold"))


header$children[[2]]$children[[2]] <-tags$div( align="left",#height='50px',
                                               tags$head(tags$style(HTML(".name { background-color: white }"))),
                                               anchor,
                                               class = 'name')

#2.User Interface----

#*  Dashboard header----
ui<-dashboardPage(header,
                 skin = "black",
                 
                 #* Dashboard sidebar ----
                 dashboardSidebar(width = 250,
                                  sidebarMenu(
                                    selectizeInput(inputId = "Gene",
                                              label = "Enter Official Gene Symbol",
                                              options = list(placeholder='Enter Gene'),choices=NULL),
                                    menuItem("Violinplot", tabName = "vlnplot", icon = icon("vp")),
                                    menuItem("Umap", tabName = "umap", icon = icon("ump")),
                                    menuItem("Dotplot", tabName = "dotplot", icon = icon("dp")),
                                    menuItem("Expression values", tabName = "cellno", icon = icon("cellno")),
                                    menuItem("Manuscript", icon = icon("Manuscript"), href ="https://www.biorxiv.org/content/10.1101/2021.02.23.432522v1"),
                                    tipify(menuItem("Experimental Summary", tabName = "Con"),"Explore gene expression level across pancreatic islet cell types from total of 44,953 single cells.Data derived from droplet based 10x Chromium single cell libraries sequenced on NovaSeq 100bp PE resulting in ~146,000 reads per cell.",
                                           placement="bottom", trigger = "hover"),
                                    tipify(menuItem("Contact Us", tabName = "Con"),"shristi.shrestha@vumc.org",
                                           placement="bottom", trigger = "hover")
                  )
                ),
                
#*  Dashboardbody----         
      dashboardBody(tags$head(tags$style(HTML('
                                            /* body */
                                            .content-wrapper, .right-side {
                                            background-color: #FFFFFF;
                                            }
                                            /* main sidebar */
                                            .skin-blue .main-sidebar { font-size: 20px;
                                                            background-color: #F0F8FF;
                                            }

                                            .main-sidebar { font-size: 20px; }

                                            '))),
                    
                    #For gene description table
                    mainPanel(tableOutput("table1")), 
                    
                    #For violinplot tab
                    tabItems(
                      tabItem(tabName = "vlnplot", 
                              fluidPage(
                                verticalLayout(plotOutput("plot1"),
                                               plotOutput("plot2"),
                                               sidebarPanel(
                                                 sliderInput("Cellsize", 
                                                             "Increase Cell Size:",
                                                             min=-1,
                                                             max=1,
                                                             value=0.1,
                                                             step = 0.1,
                                                             animate=TRUE),
                                                 hr(),
                                                 helpText("set to -1 to remove dots(cells)")
                                                 )
                                               )
                                )
                              ),
                      
                      #For UMAP plot tab
                      tabItem(tabName = "umap",
                              fluidPage(
                                verticalLayout(plotOutput("plot3"),
                                               plotOutput("plot4")
                                               )
                                )
                              ),
                      
                      #For dotplot tab
                      tabItem(tabName = "dotplot",
                              fluidPage(
                                verticalLayout(tags$h6("Note: if you see an error, enter a gene that is not a cell type marker to avoid duplication in dotplot"),
                                               plotOutput("plot5")
                                               )
                                )
                              ),
                      
                      #For table on expression counts
                      tabItem(tabName = "cellno",
                              fluidPage(
                                verticalLayout(br(),
                                               tableOutput("table2")
                                               )
                                )
                              ),
                      
                      #Add Manuscript link
                      tabItem(tabName = "manuscript",
                              h2(plotOutput("manuscript"))
                              ),
                      
                      #Add Contact us email
                      tabItem(tabName = "contactus",
                              h2(plotOutput("contactus"))
                              )
    
    )))


#3.Server----
server<-function(input, output,session) 
  
#* table 1 -Gene Description----
{
  updateSelectizeInput(session, 'Gene', choices = rownames(Islets), server = TRUE)
  output$table1<-renderTable({
    req(input$Gene)
    #rownames(Gene_desp) <- make.names(Gene_desp[,1], unique = TRUE)
    Gene_desp[toupper(input$Gene),]
    })


#* plot 1 (Violinplot by Cell types)----
  output$plot1<- renderPlot({
    req(input$Gene)
    dittoPlot(Islets, toupper(input$Gene), group.by = "CellTypes",
              plots = c("jitter", "vlnplot", "boxplot"), # <- order matters
              color.panel=c('Alpha'='#F8766D','Beta'='#39B600','Delta'='#D89000','Gamma'='#A3A500','Epsilon'='#00BF7D','Acinar'='#00BFC4','Ductal'='#00B0F6','Endothelial'='#9590FF','Stellate'='#E76BF3','Immune'='#FF62BC'),
              # change the color and size of jitter points
              jitter.color = "black", jitter.size = input$Cellsize,

              # change the outline color and width, and remove the fill of boxplots
              boxplot.color = "#636363", boxplot.width = 0.1,
              boxplot.fill = FALSE,

              # change how the violin plot widths are normalized across groups
              vlnplot.scaling = "width"
    )+ylab("ln(UMI -per-10,000 +1)")+
          theme(text = element_text(size = 18,
                                    face="bold"),
                axis.text = element_text(size = 18,
                                         face="bold"),
                axis.line.x = element_line(color="black",
                                           size = 0.8),
                axis.line.y = element_line(color="black",
                                           size = 0.8))
      
    }, height = 400, width = 750)
  
  
#* plot2 (Violin plot by Age)----
  output$plot2<- renderPlot({
    req(input$Gene)
    dittoPlot(Islets, toupper(input$Gene), group.by = "Age",
              plots = c("jitter", "vlnplot", "boxplot"), # <- order matters

              # change the color and size of jitter points
              jitter.color = "black", jitter.size = input$Cellsize,

              # change the outline color and width, and remove the fill of boxplots
              boxplot.color = "#636363", boxplot.width = 0.1,
              boxplot.fill = FALSE,
              legend.title = "Age",
              vlnplot.scaling = "width" 
    )+ylab("ln(UMI -per-10,000 +1)")+
      theme(text = element_text(size = 18,
                                face="bold"),
            axis.text = element_text(size = 18,
                                     face="bold"),
            axis.line.x = element_line(color="black",
                                       size = 0.8),
            axis.line.y = element_line(color="black",
                                       size = 0.8))
    
  }, height = 400, width = 700)
  
 
#* plot3 (umap-celltypes)----
  output$plot3<- renderPlot({
    DimPlot(Islets , reduction = "umap", label = TRUE,pt.size = 1,cols = c('Alpha'='#F8766D','Beta'='#39B600','Delta'='#D89000','Gamma'='#A3A500','Epsilon'='#00BF7D','Acinar'='#00BFC4','Ductal'='#00B0F6','Endothelial'='#9590FF','Stellate'='#E76BF3','Immune'='#FF62BC') )+NoLegend()+
      theme(text = element_text(size = 18,face="bold"),
            axis.text = element_text(size = 18,face="bold"),
            axis.line.x = element_line(color="black", size = 0.8),
            axis.line.y = element_line(color="black", size = 0.8))
  }, height = 400, width = 600)

   
#* plot 4 (umap) ----
  output$plot4<- renderPlot({
    req(input$Gene)
    FeaturePlot(Islets, features = toupper(input$Gene), pt.size =1, cols = c("lightgrey","red"))+
      theme(text = element_text(size = 18,face="bold"),
            axis.text = element_text(size = 18,face="bold"),
            axis.line.x = element_line(color="black", size = 0.8),
            axis.line.y = element_line(color="black", size = 0.8))
  }, height = 400, width = 600)

  
  
#* plot5 (Dotplot) ---- 
  Known.markers<-c("GCG","INS","SST","PPY","GHRL", "PRSS1", "KRT19","PECAM1", "PDGFRB", "HLA-DRA")
  output$plot5<- renderPlot({
  req(input$Gene)
  DotPlot(Islets, feature=c(Known.markers,toupper(input$Gene)), dot.scale = 8)+ 
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
  
  }, height = 600, width = 800)
  

#* Table2 (Expression Values) ----
  output$table2<-renderTable({
    req(input$Gene)
    ttt<-DotPlot(Islets, feature=toupper(input$Gene))
    tble<-ttt$data
    tble<-tble %>%rename(Gene= features.plot, 
                         CellType=id, 
                         'Mean expression z-score'=avg.exp.scaled,
                         '% Cells expressed'= pct.exp,
                         'Average expression'= avg.exp)
    tble<-tble[ , c(3, 4, 2,1,5)]
    tble<-tble %>% filter_all(any_vars(. %in% toupper(input$Gene)))
    tble
  })
  
}



#4. shinyApp----
shinyApp(ui=ui,server=server)
