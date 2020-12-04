
#Might need this to deploy----
#library(BiocManager)
#options(repos = BiocManager::repositories())
#options(repos = BiocInstaller::biocinstallRepos())
#options(repos = c("CRAN" = "https://cran.rstudio.com/", "BioCsoft" = "https://bioconductor.org/packages/3.8/bioc", "BioCann" = "https://bioconductor.org/packages/3.8/data/annotation"))
#getOption("repos")

#Libraries----
library(shiny)
library(Seurat)
library(ggplot2)
library(Matrix)
library(dplyr)
library(shinythemes)
library(shinydashboard)
library(shinyBS)
library(data.table)
library(dittoSeq)

#load data----
load("DATA/Islets.Rda")
Gene_desp<-read.csv("DATA/annotation.csv",header = TRUE)


#Single-cell gene expression atlas of human pancreatic islets

header <- dashboardHeader(titleWidth = "100%",
  # Set height of dashboardHeader
  tags$li(class = "dropdown",
          tags$style(".main-header {max-height: 100px}"),
          tags$style(".main-header .logo {height: 60px}")
  ))
anchor <- tags$li(a(href='https://www.powersbrissovaresearch.org/',
                 tags$img(src='logo-4.png', height='60', width='200'),
                 'Single cell gene expression atlas of human pancreatic islets',
                 style = "color: #4292C6;
                          font-family: Avenir Light;
                          font-size: 40px;
                          padding-top:10px; padding-bottom:10px;
                          font-weight: bold"))

header$children[[2]]$children[[2]] <-tags$div( align="left", 
                                               #height='50px',
  tags$head(tags$style(HTML(".name { background-color: white }"))),
  anchor,
  class = 'name')


ui<-dashboardPage(header,
                 skin = "black",
      dashboardSidebar(
                  width = 250,
          
                  sidebarMenu(
                     textInput(inputId = "Gene",
                               label = "Enter Official Gene Symbol"),
                          menuItem("Violinplot", tabName = "vlnplot", icon = icon("vp")),
                          menuItem("Umap", tabName = "umap", icon = icon("ump")),
                          menuItem("Dotplot", tabName = "dotplot", icon = icon("dp")),
                          menuItem("Expression values", tabName = "cellno", icon = icon("cellno")),
                    tipify(menuItem("Experimental Summary", tabName = "Con"),"Explore gene expression level across pancreatic islet cell types from total of 44,953 single cells.Data derived from droplet based 10x Chromium single cell libraries sequenced on NovaSeq 100bp PE resulting in ~146,000 reads per cell. Expression values = natural log of transcript per 10,000 molecules. Average expression/Scaled expression = mean expression z-score, Percent Expressed= % cells expressing the gene of interest",
                           placement="bottom", trigger = "hover"),
                    tipify(menuItem("Manuscript", tabName = "man"),"Manucript Link",
                           placement="bottom", trigger = "hover"),
                    tipify(menuItem("Contact Us", tabName = "Con"),"shristi.shrestha@vumc.org",
                           placement="bottom", trigger = "hover")
                  )
                ),
                
                
      dashboardBody(
                  tags$head(tags$style(HTML('
                                            /* body */
                                            .content-wrapper, .right-side {
                                            background-color: #FFFFFF;
                                            }
                                            /* main sidebar */
                                            .skin-blue .main-sidebar { font-size: 20px;
                                                            background-color: #F0F8FF;
                                            }

                                            .main-sidebar { font-size: 20px; }

                                            '
                  ))),
                                            
                
                  
                  mainPanel(
                    tableOutput("table1")
                  ),     
                  
  
  tabItems(
    
    tabItem(tabName = "vlnplot", 
            fluidPage(
              verticalLayout(plotOutput("plot1"),
                          plotOutput("plot2"),
                          sidebarPanel(
                            sliderInput("Cellsize", "Increase Cell Size:",
                                        min=-1,
                                        max=1,
                                        value=0.1,
                                        step = 0.1,
                                        
                                        animate=TRUE
                            ),
                            hr(),
                            helpText("set to -1 to remove dots(cells)")
                                      )

                             )
                        )
  
              ),
  
    tabItem(tabName = "umap",
            fluidPage(
              verticalLayout(plotOutput("plot3"),
                             plotOutput("plot4")
                             )
            )
    ),
    
    tabItem(tabName = "dotplot",
            fluidPage(
              verticalLayout(tags$h6("Note: enter a gene that is not a cell type marker to avoid duplication in dotplot"),
                             plotOutput("plot5")
                             
              ))

    ),
    
    tabItem(tabName = "cellno",
            fluidPage(
              verticalLayout(br(),
                             tableOutput("table3")
                             )
              )
    ),
    
    tabItem(tabName = "manuscript",
            h2(plotOutput("manuscript"))
    ),
    tabItem(tabName = "contactus",
            h2(plotOutput("contactus"))
           
    )
    
    )))

server<-function(input, output) 
  
  #table 1
{
  output$table1<-renderTable({
    req(input$Gene)
    rownames(Gene_desp) <- make.names(Gene_desp[,1], unique = TRUE)
    Gene_desp[toupper(input$Gene),]
  })


#plot 1 (Violinplot by Cell types)
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
    )+
          theme(text = element_text(size = 18,
                                    face="bold"),
                axis.text = element_text(size = 18,
                                         face="bold"),
                axis.line.x = element_line(color="black",
                                           size = 0.8),
                axis.line.y = element_line(color="black",
                                           size = 0.8))
      
    }, height = 400, width = 750)
  
  
  #plot2 (Violin plot by Age)
  output$plot2<- renderPlot({
    req(input$Gene)
    dittoPlot(Islets, toupper(input$Gene), group.by = "Age",
              plots = c("jitter", "vlnplot", "boxplot"), # <- order matters

              # change the color and size of jitter points
              jitter.color = "black", jitter.size = input$Cellsize,

              # change the outline color and width, and remove the fill of boxplots
              boxplot.color = "#636363", boxplot.width = 0.1,
              boxplot.fill = FALSE,

              vlnplot.scaling = "width" 
    )+
      theme(text = element_text(size = 18,
                                face="bold"),
            axis.text = element_text(size = 18,
                                     face="bold"),
            axis.line.x = element_line(color="black",
                                       size = 0.8),
            axis.line.y = element_line(color="black",
                                       size = 0.8))
    
  }, height = 400, width = 700)
  
 
#plot3 (umap-celltypes)----
  output$plot3<- renderPlot({
    DimPlot(Islets , reduction = "umap", label = TRUE,pt.size = 1,cols = c('Alpha'='#F8766D','Beta'='#39B600','Delta'='#D89000','Gamma'='#A3A500','Epsilon'='#00BF7D','Acinar'='#00BFC4','Ductal'='#00B0F6','Endothelial'='#9590FF','Stellate'='#E76BF3','Immune'='#FF62BC') )+NoLegend()+
      theme(text = element_text(size = 18,face="bold"),
            axis.text = element_text(size = 18,face="bold"),
            axis.line.x = element_line(color="black", size = 0.8),
            axis.line.y = element_line(color="black", size = 0.8))
  }, height = 400, width = 600)

   
#plot 4 (umap) ----
  output$plot4<- renderPlot({
    req(input$Gene)
    FeaturePlot(Islets, features = toupper(input$Gene), pt.size =1, cols = c("lightgrey", "red","blue"), order=TRUE,min.cutoff = "q10", max.cutoff = "q90")+
      theme(text = element_text(size = 18,face="bold"),
            axis.text = element_text(size = 18,face="bold"),
            axis.line.x = element_line(color="black", size = 0.8),
            axis.line.y = element_line(color="black", size = 0.8))
  }, height = 400, width = 600)

  
  
#plot5 (Dotplot) ---- 
  Known.markers<-c("GCG","INS","SST","PPY","GHRL", "PRSS1", "KRT19","PECAM1", "PDGFRB", "HLA-DRA")
  output$plot5<- renderPlot({
  req(input$Gene)
  DotPlot(Islets, feature=c(Known.markers,toupper(input$Gene)), dot.scale = 8)+ 
    coord_flip()+ 
    scale_color_gradient2(low = "blue", high = "red",mid = "white")+
    labs(y=" ", x="Genes", color = "Expression z-score")+RotatedAxis()+
    theme(text = element_text(size = 15,face="bold"),axis.text = element_text(size = 15,face="bold"))
  
  }, height = 700, width = 900)
  

  #Table3 (Expression Values) ----
  dataset <- reactive({
    Geneofinterest<-subset(Islets, toupper(input$Gene)>0)
    
    return(Geneofinterest)
  })
  

  output$table3<-renderTable({
    req(input$Gene)
    ttt<-DotPlot(Islets, feature=toupper(input$Gene))
    tble<-ttt$data
    tble<-tble %>%rename(Gene= features.plot, 
                         CellType=id, 
                         AverageExpression_Scaled=avg.exp.scaled,
                         '% cells expressed'= pct.exp,
                         AverageExpression= avg.exp)
    tble<-tble[ , c(3, 4, 2,1,5)]
    tble<-tble %>% filter_all(any_vars(. %in% toupper(input$Gene)))
    tble
  })
  
}




shinyApp(ui=ui,server=server)
