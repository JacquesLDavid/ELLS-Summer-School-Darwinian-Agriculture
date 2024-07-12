# authors Jacques David, Nicolas Salas
# date : July 10th 2024
# Why : ELLS Summer school Darwinian Agriculture
# Still to be bone 
# Why is it so long ?
# Add the genetic advance in the summary


library(shiny)

# Load the necessary packages
if(!requireNamespace("tidyverse", quietly = TRUE)){ install.packages("tidyverse")}
if(!requireNamespace("MASS", quietly = TRUE)){ install.packages("MASS")}
if(!requireNamespace("ggplot2", quietly = TRUE)){ install.packages("ggplot2")}
if(!requireNamespace("sommer", quietly = TRUE)){ install.packages("sommer")}
if(!requireNamespace("data.table", quietly = TRUE)){ install.packages("data.table")}
if(!requireNamespace("shinycssloaders", quietly = TRUE)){ install.packages("shinycssloaders")}
if(!requireNamespace("schoolmath", quietly = TRUE)){ install.packages("schoolmath")}
if(!requireNamespace("corrplot", quietly = TRUE)){ install.packages("corrplot")}

library(MASS) # For mvrnorm
library(ggplot2) # For plots
library(sommer)
library(data.table)
library(tidyverse)
library(shinycssloaders)
library(corrplot)



# UI
ui <- fluidPage(
  tags$style(HTML("
            .table-container {
                display: flex;
                justify-content: center;
            }
            table {
                margin: auto;
            }
        ")),
  titlePanel("Population Evolution Simulation"),
  sidebarLayout(
    sidebarPanel(
		p("GENERAL PARAMETERS:"),
      checkboxInput("setSeed", "Set seed for reproducibility", value = FALSE),
      numericInput("seedValue", "Seed value", value = 12345),
      checkboxInput("Mean","Mean effects (checked) or Sum effects (unchecked) for IGE calculation", value=FALSE),
      checkboxInput("Asreml","Use of Asreml to make inference", value=FALSE),
		  checkboxInput("is.spherical","'Spherical' neighbordhood or 'flat' neighborhood", value=TRUE),
	  p("Make the product of N x Rep reasonably close to 1000 at the maximum"),
      sliderInput("N", "Number of genotypes (N)", min = 10, max = 500, value = 50),
      sliderInput("rep", "Number of rep per genotype  (rep)", min = 1, max = 100, value = 2),
      sliderInput("N_sim", "Number of simulations  (N_sim)", min = 1, max = 100, value = 1),
	  tags$hr(),
      p("GENETIC VARIANCES:"),
      numericInput("varG11", "Genetic variance DGE", value = 1),
      numericInput("varG22", "Genetic variance IGE", value = 0.125),
      sliderInput("r", "Genetic correlation DGE : IGE", min = -1, max = 1, value = 0, step = 0.1),
	  tags$hr(),
      p("ENVIRONMENTAL VARIANCES:"),
      numericInput("varE11", "Environmental variance DGE", value = 1),
      numericInput("varE22", "Environmental variance IGE", value = 0.125),
      div(style = "text-align: center;font-weight: bold; color: blue;",
          actionButton("goButton", label = div(style = "color: red; font-weight: bold;", "Step 1. Run Simulation"))
      ),
      tags$hr(),
      p("SELECTION"),
      sliderInput("p", "Selection pressure", min = 0.01, max = 1, value = 0.5, step = 0.1),          
      sliderInput("b_DGE", "Selection index : Weight on DGE (IGE =1-DGE)", min = -1, max = 1, value = 0.5, step = 0.1),
      div(style = "text-align: center;",
          actionButton("SelButton",
                       label = div(style = "color: blue; font-weight: bold;", "Step 2. Breed !"))
      )
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Population characteristics", 
                           tags$hr(),
                           p("Relation between True Direct Genetic Effect and the True Indirect Genetic Effect"),
                           plotOutput("plotTRUE_DGE_IGE"),
                           tags$hr(),
                           p("Quality of estimates of DGE"),
                           plotOutput("plotTRUE_DGE_PRED_DGE"),
                           tags$hr(),
                           p("Quality of estimates of iGE"),
                           plotOutput("plotTRUE_IGE_PRED_IGE"),
                           tags$hr(),
                           p("Relation between estimates of DGE and IGE"),
                           plotOutput("plotPRED_DGE_IGE"),
                  ),
                  
                  tabPanel("Selection", 
                           tags$hr(),
                           p("Individual Mass selection of phenotypic values"),
                           plotOutput("plotMass_differential"),
                           tags$hr(),
                           p("Genetic advance on phenotypic values under mass selection"),
                           plotOutput("plotMass_Selection"),
                           tags$hr(),
                           p("Genetic advance on phenotypic values under DGEvsIGE index selection "),
                           plotOutput("plotIndex_Selection"),
                  ),
                  tabPanel("Summary",
                           div(style = "text-align: center;",
                               h2("Before selection"),
                               h4("TRUE Mean and Variance parammeters"),
                               br(),
                               h4("Variance Table"),
                               div(class = "table-container",
                                   tableOutput("table_TrueOutput")
                               ),
                               br(),
                               h4("Correlation Values"),
                               plotOutput("corrplot"),
                               br(),
                               h2("After selection"),
                               br(),
                               h3("Mass Selection"),
                               div(class = "table-container",
                                   tableOutput("table_True_selOutput")
                               ),
                               br(),
                               h3("Index selection"),
                               div(class = "table-container",
                                   tableOutput("table_True_sel_IOutput")
                               ),
                               br(),
                               h3("Phenotypic gain"),
                               div(class = "table-container",
                                   tableOutput("TABLE_gain")
                               ),
                           ),
                           
                  )
                  
      )
    )
  )
)

# Server
server <- function(input, output) {
  
  voisin_fct=function(DATA,mat_grid){
    for (i in 1:nrow(mat_grid)){
      if(i==1){
        for (j in 1:ncol(mat_grid)){
          if(j==1){
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j+1]]+1
          }
          else if(j==ncol(mat_grid)){
            
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,1]]+1
          }
          else{
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[nrow(mat_grid),j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j+1]]+1
          }
        }
      }
      else if (i==nrow(mat_grid)){
        for (j in 1:(ncol(mat_grid))){
          if(j==1){
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j+1]]+1
          }
          else if(j==ncol(mat_grid)){
            
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,1]]+1
          }
          else{
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[1,j+1]]+1
          }
        }
      }
      
      else{
        for (j in 1:(ncol(mat_grid))){
          if(j==1){
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,ncol(mat_grid)]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,ncol(mat_grid)]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j+1]]+1
            
          }
          else if(j==ncol(mat_grid)){
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,1]]+1
          }
          else{
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j-1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j-1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i-1,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i,j+1]]+1
            DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j+1]]=DATA[which(DATA$Row==i&DATA$Column==j),mat_grid[i+1,j+1]]+1
          }
        }
      }
    }
    return(DATA)
  }
  voisin_fct_sans_recollage <- function(DATA, mat_grid) {
    nrows <- nrow(mat_grid)
    ncols <- ncol(mat_grid)
    
    for (i in 1:nrows) {
      for (j in 1:ncols) {
        voisins <- list(
          c(i-1, j-1), c(i-1, j), c(i-1, j+1),
          c(i, j-1),             c(i, j+1),
          c(i+1, j-1), c(i+1, j), c(i+1, j+1)
        )
        
        for (v in voisins) {
          vi <- v[1]
          vj <- v[2]
          
          if (vi >= 1 && vi <= nrows && vj >= 1 && vj <= ncols) {
            focal <- DATA[which(DATA$Row == i & DATA$Column == j), "Focal"]
            DATA[which(DATA$Row == i & DATA$Column == j), mat_grid[vi, vj]] <- 
              DATA[which(DATA$Row == i & DATA$Column == j), mat_grid[vi, vj]] + 1
          }
        }
      }
    }
    return(DATA)
  }
  assign("voisin_fct",voisin_fct,globalenv())
  assign("voisin_fct_sans_recollage",voisin_fct_sans_recollage,globalenv())
  
  observeEvent(input$goButton, {
    output$plotTRUE_DGE_IGE<-renderPlot(ggplot())
    output$plotTRUEvsPRED_DGE<-renderPlot(ggplot())
    output$plotTRUEvsPRED_IGE<-renderPlot(ggplot())
    output$plotPred_DGE_IGE<-renderPlot(ggplot())
    output$plotMass_differential<-renderPlot(ggplot())
    output$plotMass_Selection<-renderPlot(ggplot())
    output$plotIndex_Selection<-renderPlot(ggplot())
    is.spherical=input$is.spherical
    if(input$setSeed) {
      set.seed(input$seedValue)
    }
    
    Asreml=input$Asreml
    Mean=input$Mean
    # Asreml=FALSE
    # Mean=FALSE
    if(Asreml){
      library(asreml)
      asreml.options(workspace="2gb",maxit=5,ai.sing=TRUE)
      Asreml=TRUE
    }
    else{
      Asreml=FALSE
    }
    
    # V_env_DGE=0
    # V_env_IGE=0
    # V_geno=1
    # V_voisin=1/8*V_geno
    # cor_geno_voisin=0
    # cov_geno_voisin=cor_geno_voisin*sqrt(V_geno*V_voisin)
    # N_geno=100
    # N_rep=2
    # b_DGE=0.5
    # p=0.1
    
    V_env_DGE=input$varE11
    V_env_IGE=input$varE22
    V_geno=input$varG11
    V_voisin=input$varG22
    cor_geno_voisin=input$r
    cov_geno_voisin=cor_geno_voisin*sqrt(V_geno*V_voisin)
    N_geno=input$N
    N_rep=input$rep
    N_sim=input$N_sim
    N_row=N_col=0
    
    if((V_voisin==0)&(V_env_IGE==0)){
      V_env_IGE=10^-6
    }
    
    if((V_geno==0)&(V_env_DGE==0)){
      V_env_DGE=10^-6
    }
    
    
    c=0
    if (ceiling(sqrt(N_geno*N_rep))!=sqrt(N_geno*N_rep)){
      N_row=N_col=ceiling(sqrt(N_geno*N_rep))
    }else {
      N_row=N_col=sqrt(N_geno*N_rep)
    }
    
    assign("N_geno",N_geno,envir = globalenv())
    assign("N_rep",N_rep,envir = globalenv())
    assign("N_row",N_row,envir = globalenv())
    assign("N_col",N_col,envir = globalenv())
    assign("Mean",Mean,envir = globalenv())
    assign("Asreml",Asreml,envir=globalenv())
    assign("is.spherical",is.spherical,envir=globalenv())
    
    SIM <- c()
    calc_SIM <- c()
    SIM_DGE=c()
    SIM_IGE=c()
    SIM_cov=c()
    
    SIM_mean <- c()
    SIM_mean_DGE=c()
    SIM_mean_IGE=c()
    
    # Pré-allocation
    mu <- c(0, 0)
    G <- matrix(c(V_geno, cov_geno_voisin, cov_geno_voisin, V_voisin), ncol=2, nrow=2)
    E <- matrix(c(V_env_DGE, 0, 0, V_env_IGE), ncol=2, nrow=2)
    P <- G + E
    assign("mu",mu,envir = globalenv())
    assign("G",G,envir=globalenv())
    assign("E",E,envir = globalenv())
    assign("P",P,envir = globalenv())
    
    Focal <- sprintf("G%03d", 1:N_geno)
    rep <- sprintf("%03d", 1:N_rep)
    grid <- expand.grid(Row = 1:N_row, Column = 1:N_col)
    
    for (i in 1:N_sim) {
      print(i)
      
      df <- mvrnorm(n = N_geno, mu, G)
      colnames(df) <- c("DGE", "IGE")
      
      output$plotTRUE_DGE_IGE=renderPlot({ggplot(df,aes(DGE,IGE))+
        geom_point()+
        labs(x="TRUE DGE",y="TRUE IGE")+
        theme_minimal()
      })
      df_E <- mvrnorm(n = (N_geno * N_rep), mu, E)
      assign("df_E",df_E,globalenv())
      
      DGE <- df[, 1]
      assign("DGE",DGE,globalenv())
      IGE <- df[, 2]
      assign("IGE",IGE,globalenv())
      
      DATA=expand_grid("Focal"=as.character(paste0("G",sprintf("%03d", 1:N_geno))),"rep"=as.character(sprintf("%03d", 1:N_rep)))
      
      grid=expand_grid("Row"=factor(1:N_row),"Column"=factor(1:N_col))
      grid=grid[sample(1:(N_row*N_col),(N_geno*N_rep)),]
      
      mat_grid=matrix("vide",nrow=N_row,ncol=N_col)
      DATA=cbind(DATA,grid)
      
      for (i in 1:(nrow(mat_grid))){
        for (j in 1:(ncol(mat_grid))){
          if (!is_empty(DATA[(which(DATA$Row==i&DATA$Column==j)),"Focal"])){
            mat_grid[i,j]=DATA[(which(DATA$Row==i&DATA$Column==j)),"Focal"]
          }
        }
      }
      mat_grid[grep(mat_grid,pattern="vide")]=sample(DATA$Focal,length(mat_grid[grep(mat_grid,pattern="vide")]))
      
      matrice_voisin=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA$Focal)))
      DATA=cbind(matrice_voisin,DATA,data.frame("vide"=NA))
      
      if (is.spherical){
        DATA=voisin_fct(DATA,mat_grid)
      }else{
        DATA=voisin_fct_sans_recollage(DATA,mat_grid)
      }
      
      Zg <- model.matrix(~Focal - 1, DATA)
      dimnames(Zg)[[2]] <- paste0("G", sprintf("%03d", 1:N_geno))
      
      if (Mean) {
        Zv <- as.matrix(DATA[, 1:N_geno]) / rowSums(DATA[1:N_geno])
        DATA[, 1:N_geno] <- Zv
      } else {
        Zv <- as.matrix(DATA[, 1:N_geno])
      }
      
      Pheno <- Zg %*% DGE + Zv %*% IGE + df_E[, 1] + df_E[, 2]
      
      SIM <- c(SIM,var(Pheno))
      calc_SIM <- c(calc_SIM,round(var(DGE) + 8 * var(IGE) + 8 * mean(tcrossprod(IGE + DGE)) * (2 * cov(DGE, IGE) + 7 * var(IGE)) / (N_col * N_row), 3))
      SIM_cov=c(SIM_cov,cov(DGE,IGE))
      SIM_DGE=c(SIM_DGE,var(DGE))
      SIM_IGE=c(SIM_IGE,var(IGE))
      
      SIM_mean <- c(SIM_mean,mean(Pheno))
      SIM_mean_DGE=c(SIM_mean_DGE,mean(DGE))
      SIM_mean_IGE=c(SIM_mean_IGE,mean(IGE))
    }
    
    TABLE_TRUE <- data.frame(
      "Effect" = c("DGE", "IGE","Cov_DGE_IGE", "Pheno"),
      "Variance" = c(mean(SIM_DGE), mean(SIM_IGE), mean(SIM_cov),mean(SIM)),
      "Mean"=c(mean(SIM_mean_DGE),mean(SIM_mean_IGE),NA,mean(SIM_mean))
    )
    
    output$table_TrueOutput <- renderTable({
      TABLE_TRUE
    })
    # mass phenotypic selection
    DATA$Pheno <- as.vector(Pheno)
    DATA$Focal <- as.factor(DATA$Focal)
    
    # Estimation of BLUP 
    if (Asreml) {
      Mod <- try(asreml(fixed = Pheno ~ 1,
                        random = ~str(~Focal + grp(Voisin), ~us(2):id(Focal)),
                        group = list(Voisin = 1:N_geno),
                        residual = ~units,
                        data = DATA))
      
      tmp_DGE <- data.frame("DGE_pred" = summary(Mod, coef = TRUE)$coef.random[1:N_geno, 1],
                            "Focal" = str_split(names(summary(Mod, coef = TRUE)$coef.random[1:N_geno, 1]), pattern = "_", simplify = TRUE)[, 2])
      tmp_IGE <- data.frame("IGE_pred" = summary(Mod, coef = TRUE)$coef.random[(N_geno + 1):(2 * N_geno), 1],
                            "Focal" = str_split(names(summary(Mod, coef = TRUE)$coef.random[1:N_geno, 1]), pattern = "_", simplify = TRUE)[, 2])
      pred <- merge(tmp_DGE, tmp_IGE, by = "Focal")
    } else {
      Mod <- try(mmer(fixed = Pheno ~ 1,
                      random = ~vsr(Focal) + vsr(Zv),
                      rcov = ~units,
                      data = DATA, nIters = 4))
      
      if (length(Mod) < 3) {
        Mod <- NA
        pred <- NA
      } else {
        DGE_pred <- data.frame("Focal" = names(randef(Mod)$`u:Focal`$Pheno), "DGE_pred" = as.numeric(randef(Mod)$`u:Focal`[[1]]))
        IGE_pred <- data.frame("Focal" = names(randef(Mod)$`u:Focal`$Pheno), "IGE_pred" = as.numeric(randef(Mod)$`u:Zv`[[1]]))
        pred <- merge(DGE_pred, IGE_pred, by = "Focal")
      }
    }
    
    output$plotTRUE_DGE_PRED_DGE=renderPlot({ggplot(df,aes(x=DGE,y=pred$DGE_pred))+
      geom_point()+
      labs(x="TRUE DGE",y="PRED DGE")+
      theme_minimal()
    })
    
    output$plotTRUE_IGE_PRED_IGE=renderPlot({ggplot(df,aes(x=IGE,y=pred$IGE_pred))+
      geom_point()+
      labs(x="TRUE IGE",y="PRED IGE")+
      theme_minimal()
  })
    output$plotPRED_DGE_IGE=renderPlot({ggplot(pred,aes(x=IGE_pred,y=DGE_pred))+
      geom_point()+
      labs(x="PRED IGE",y="PRED DGE")+
      theme_minimal()
})
    assign("pred",pred,globalenv())
    assign("DATA", DATA, envir = globalenv())
    assign("N_geno", N_geno, envir = globalenv())
    assign("Pheno", Pheno, envir = globalenv())
    assign("Zg",Zg,globalenv())
    assign("Zv", Zv, envir = globalenv())
    
    output$corrplot=renderPlot({
      corrplot::corrplot(cor(data.frame("TRUE_DGE"=DGE,"TRUE_IGE"=IGE,"PRED_DGE"=pred$DGE_pred,"PRED_IGE"=pred$IGE_pred)),type = "lower",addCoef.col = "black")
    })
  })
  
  observeEvent(input$SelButton,{
    b_DGE=input$b_DGE
    p<-input$p
    
    #####Index selection
    
    pred$I <- b_DGE*pred$DGE_pred+(1-b_DGE)*pred$IGE_pred
    
    # list of selected genotypes (numbers will be different from mass phenotype selection)
    sel = which(pred$I>quantile(pred$I,1-p))
    length(sel)
    
    R_I<-c(mean(DGE[sel]),mean(IGE[sel]))
    
    # R=G%*%solve(P)%*%S
    
    mus=c(mu[1]+R_I[1],mu[2]+R_I[2])
    
    Focal_sel_I <- sprintf("G%03d", 1:N_geno)
    rep_sel_I <- sprintf("%03d", 1:N_rep)
    grid <- expand.grid(Row = 1:N_row, Column = 1:N_col)
    
    
    df_sel_I=mvrnorm(N_geno,mus,G)
    colnames(df_sel_I)=c("DGE","IGE")
    df_E_sel_I <- mvrnorm(n = (N_geno * N_rep), mus, E)
    
    DGE_sel_I <- df_sel_I[, 1]
    IGE_sel_I <- df_sel_I[, 2]
    
    DATA_sel_I=expand_grid("Focal"=as.character(paste0("G",sprintf("%03d", 1:N_geno))),"rep"=as.character(sprintf("%03d", 1:N_rep)))
    
    grid=expand_grid("Row"=factor(1:N_row),"Column"=factor(1:N_col))
    grid=grid[sample(1:(N_row*N_col),(N_geno*N_rep)),]
    
    mat_grid=matrix("vide",nrow=N_row,ncol=N_col)
    DATA_sel_I=cbind(DATA_sel_I,grid)
    
    for (i in 1:(nrow(mat_grid))){
      for (j in 1:(ncol(mat_grid))){
        if (!is_empty(DATA_sel_I[(which(DATA_sel_I$Row==i&DATA_sel_I$Column==j)),"Focal"])){
          mat_grid[i,j]=DATA_sel_I[(which(DATA_sel_I$Row==i&DATA_sel_I$Column==j)),"Focal"]
        }
      }
    }
    mat_grid[grep(mat_grid,pattern="vide")]=sample(DATA_sel_I$Focal,length(mat_grid[grep(mat_grid,pattern="vide")]))
    
    matrice_voisin=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA_sel_I$Focal)))
    DATA_sel_I=cbind(matrice_voisin,DATA_sel_I,data.frame("vide"=NA))
    
    DATA_sel_I=voisin_fct(DATA_sel_I,mat_grid)
    
    Zg_sel_I <- model.matrix(~Focal - 1, DATA_sel_I)
    dimnames(Zg_sel_I)[[2]] <- paste0("G", sprintf("%03d", 1:N_geno))
    
    if (Mean) {
      Zv_sel_I <- as.matrix(DATA_sel_I[, 1:N_geno]) / 8
      DATA_sel_I[, 1:N_geno] <- Zv
    } else {
      Zv_sel_I <- as.matrix(DATA_sel_I[, 1:N_geno])
    }
    
    Pheno_sel_I <- Zg_sel_I %*% DGE_sel_I + Zv_sel_I %*% IGE_sel_I + df_E_sel_I[, 1] + df_E_sel_I[, 2]
    
    # Combine the two vectors into a dataframe
    combinedData_I <- rbind(data.frame(Value = Pheno, Phase = "Before selection"),
                            data.frame(Value = Pheno_sel_I, Phase = "After selection"))
    
    mean_before_sel_I=round(mean(combinedData_I[combinedData_I$Phase=="Before selection","Value"]),2)
    mean_after_sel_I=round(mean(combinedData_I[combinedData_I$Phase=="After selection","Value"]),2)
    
    output$plotIndex_Selection <- renderPlot({
      ggplot(combinedData_I, aes(x = Value, fill = Phase)) +
        geom_density(alpha = 0.5, adjust = 1.5) +
        geom_vline(aes(xintercept=mean_before_sel_I),color="blue",size=1)+
        geom_vline(aes(xintercept=mean_after_sel_I),color="red",size=1)+# Adjustment for smoothing, 'adjust' controls the smoothing
        scale_fill_manual(values = c("Before selection" = "blue", "After selection" = "red")) +
        labs(x = "Phenotypic Value", y = "Density", title = "Index selection") +
        theme_minimal() +
        theme(legend.title = element_blank())+
        xlim(-5*var(Pheno),+5*var(Pheno))+
        annotate("text",x=mean_after_sel_I+1.5,y=0.3,label=paste0('Delta_mu_pheno = ', mean_after_sel_I - mean_before_sel_I))# Remove the legend title
    })
    
    ###Mass selection
    
    # Selected individuals
    N_sel <- round(nrow(DATA) * p)
    list_sel <- order(DATA$Pheno, decreasing = TRUE)[1:N_sel]
    
    # Définir le seuil
    threshold <- min(DATA$Pheno[list_sel])
    
    # graphDATA
    graphDATA<-c()
    graphDATA$Pheno <-DATA$Pheno
    threshold <- min(graphDATA$Pheno[list_sel])
    graphDATA$color <- ifelse(graphDATA$Pheno < threshold, "left", "right")
    
    # histogram appearing with S the selection differential
    output$plotMass_differential <- renderPlot({
      ggplot(as.data.frame(graphDATA), aes(x = Pheno, fill = color)) +
        geom_histogram(binwidth = 0.2, color = "black", boundary = threshold) +
        scale_fill_manual(values = c("left" = "red", "right" = "blue")) +
        geom_vline(xintercept = threshold, linetype = "dashed", color = "black", size = 1) +
        labs(title = "mass selection on phenotypic values",
             x = "Pheno",
             y = "Frequency") +
        theme_minimal()
    })
    
    # Selection differential
    S<-c()
    S[1] <- mean((Zg%*%DGE + df_E[,1])[list_sel])
    S[2] <- mean((Zg%*%IGE + df_E[,2])[list_sel])
    
    # Genetic gain
    R=G%*%solve(P)%*%S
    mus=c(mu[1]+R[1],mu[2]+R[2])
    
    Focal_sel <- sprintf("G%03d", 1:N_geno)
    rep_sel <- sprintf("%03d", 1:N_rep)
    grid <- expand.grid(Row = 1:N_row, Column = 1:N_col)

    
    df_sel=mvrnorm(N_geno,mus,G)
    colnames(df_sel)=c("DGE","IGE")
    df_E_sel <- mvrnorm(n = (N_geno * N_rep), mus, E)
    
    DGE_sel <- df_sel[, 1]
    IGE_sel <- df_sel[, 2]
    
    DATA_sel=expand_grid("Focal"=as.character(paste0("G",sprintf("%03d", 1:N_geno))),"rep"=as.character(sprintf("%03d", 1:N_rep)))
    
    grid=expand_grid("Row"=factor(1:N_row),"Column"=factor(1:N_col))
    grid=grid[sample(1:(N_row*N_col),(N_geno*N_rep)),]
    
    mat_grid=matrix("vide",nrow=N_row,ncol=N_col)
    DATA_sel=cbind(DATA_sel,grid)
    
    for (i in 1:(nrow(mat_grid))){
      for (j in 1:(ncol(mat_grid))){
        if (!is_empty(DATA_sel[(which(DATA_sel$Row==i&DATA_sel$Column==j)),"Focal"])){
          mat_grid[i,j]=DATA_sel[(which(DATA_sel$Row==i&DATA_sel$Column==j)),"Focal"]
        }
      }
    }
    mat_grid[grep(mat_grid,pattern="vide")]=sample(DATA_sel$Focal,length(mat_grid[grep(mat_grid,pattern="vide")]))
    
    matrice_voisin=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA_sel$Focal)))
    DATA_sel=cbind(matrice_voisin,DATA_sel,data.frame("vide"=NA))
    
    DATA_sel=voisin_fct(DATA_sel,mat_grid)
    
    Zg_sel <- model.matrix(~Focal - 1, DATA_sel)
    dimnames(Zg_sel)[[2]] <- paste0("G", sprintf("%03d", 1:N_geno))
    
    if (Mean) {
      Zv_sel <- as.matrix(DATA_sel[, 1:N_geno]) / 8
      DATA_sel[, 1:N_geno] <- Zv
    } else {
      Zv_sel <- as.matrix(DATA_sel[, 1:N_geno])
    }
    
    Pheno_sel <- Zg_sel %*% DGE_sel + Zv_sel %*% IGE_sel + df_E_sel[, 1] + df_E_sel[, 2]
    
   
    # Combine the two vectors into a dataframe
    combinedData <- rbind(data.frame(Value = Pheno, Phase = "Before selection"),
                          data.frame(Value = Pheno_sel, Phase = "After selection"))
    
    mean_before_sel=round(mean(combinedData[combinedData$Phase=="Before selection","Value"]),2)
    mean_after_sel=round(mean(combinedData[combinedData$Phase=="After selection","Value"]),2)
    
    output$plotMass_Selection <- renderPlot({
      ggplot(combinedData, aes(x = Value, fill = Phase)) +
        geom_density(alpha = 0.5, adjust = 1.5) +
        geom_vline(aes(xintercept=mean_before_sel),color="blue",size=1)+
        geom_vline(aes(xintercept=mean_after_sel),color="red",size=1)+# Adjustment for smoothing, 'adjust' controls the smoothing
        scale_fill_manual(values = c("Before selection" = "blue", "After selection" = "red")) +
        labs(x = "Phenotypic Value", y = "Density", title = "Phenotypic mass selection") +
        theme_minimal() +
        xlim(-5*var(Pheno),+5*var(Pheno))+
        theme(legend.title = element_blank())+# Remove the legend title
        annotate("text",x=mean_after_sel+1.5,y=0.3,label=paste0('Delta_mu_pheno = ', mean_after_sel - mean_before_sel))# Remove the legend title
      
    })
    
    
    TABLE_TRUE_sel_I <- data.frame(
      "Effect" = c("DGE", "IGE","Cov_DGE_IGE", "Pheno"),
      "Variance" = c(var(DGE_sel_I), var(IGE_sel_I), cov(DGE_sel_I,IGE_sel_I),var(Pheno_sel_I)),
      "Mean"=c(mean(DGE_sel_I),mean(IGE_sel_I),NA,mean(Pheno_sel_I))
    )
    
    output$table_True_sel_IOutput <- renderTable({
      TABLE_TRUE_sel_I
    })
    
    TABLE_TRUE_sel <- data.frame(
      "Effect" = c("DGE", "IGE","Cov_DGE_IGE", "Pheno"),
      "Variance" = c(var(DGE_sel), var(IGE_sel), cov(DGE_sel,IGE_sel),var(Pheno_sel)),
      "Mean"=c(mean(DGE_sel),mean(IGE_sel),NA,mean(Pheno_sel))
    )
    
    output$table_True_selOutput <- renderTable({
      TABLE_TRUE_sel
    })
    
    TABLE_gain=data.frame(
      "Selection_type"=c("Mass selection","Index selection"),
      "Gain"=c((mean(Pheno_sel)-mean(Pheno)),(mean(Pheno_sel_I)-mean(Pheno)))
    )
    
    output$TABLE_gain <- renderTable({
      TABLE_gain
    })
    
  })
  
}

# Run the Shiny application
shinyApp(ui = ui, server = server)

#Minor change to test git push to branch test_PB