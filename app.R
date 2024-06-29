library(shiny)

# Load the necessary packages
if(!requireNamespace("tidyverse", quietly = TRUE)){ install.packages("tidyverse")}
if(!requireNamespace("MASS", quietly = TRUE)){ install.packages("MASS")}
if(!requireNamespace("ggplot2", quietly = TRUE)){ install.packages("ggplot2")}
if(!requireNamespace("sommer", quietly = TRUE)){ install.packages("sommer")}
if(!requireNamespace("data.table", quietly = TRUE)){ install.packages("data.table")}
if(!requireNamespace("shinycssloaders", quietly = TRUE)){ install.packages("shinycssloaders")}
if(!requireNamespace("schoolmath", quietly = TRUE)){ install.packages("schoolmath")}

library(MASS) # For mvrnorm
library(ggplot2) # For plots
library(sommer)
library(data.table)
library(tidyverse)
library(shinycssloaders)



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
      checkboxInput("setSeed", "Set seed for reproducibility", value = FALSE),
      numericInput("seedValue", "Seed value", value = 12345),
      #checkboxInput("Mean","Mean effects (checked) or Sum effects (unchecked) for IGE calculation", value=TRUE),
      #checkboxInput("Asreml","Use of Asreml to make inference", value=FALSE),
      sliderInput("N", "Number of genotypes (N)", min = 10, max = 500, value = 450),
      sliderInput("rep", "Number of rep per genotype  (rep)", min = 1, max = 100, value = 2),
      numericInput("varG11", "Genetic variance DGE", value = 1),
      numericInput("varG22", "Genetic variance IGE", value = 0.125),
      sliderInput("r", "Genetic correlation DGE : IGE", min = -1, max = 1, value = 0, step = 0.1),
      numericInput("varE11", "Environmental variance DGE", value = 1),
      numericInput("varE22", "Environmental variance IGE", value = 0.125),
      div(style = "text-align: center;",
          actionButton("goButton", "Run Simulation")
      ),
      sliderInput("p", "Selection pressure", min = 0.01, max = 1, value = 0.1, step = 0.1),          
      sliderInput("b_DGE", "Index weight for DGE (weight for IGE = 1 - weight for DGE)", min = -1, max = 1, value = 0.5, step = 0.1),
      div(style = "text-align: center;",
          actionButton("SelButton","Make Selection")
      )
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Population characteristics", 
                           plotOutput("plotTRUE_DGE_IGE"),
                           plotOutput("plotTRUEvsPRED_DGE"),
                           plotOutput("plotTRUEvsPRED_IGE"),
                           plotOutput("plotPred_DGE_IGE"),
                  ),
                  tabPanel("Selection", 
                           plotOutput("plotMass_differential"),
                           plotOutput("plotMass_Selection"),
                           plotOutput("plotIndex_Selection"),
                  ),
                  tabPanel("Summary",
                           div(style = "text-align: center;",
                               h2("Before selection"),
                               h4("TRUE Mean and Variance parammeters"),
                               textOutput("TRUEVarOutputDGE"),
                               textOutput("TRUEMeanOutputDGE"),
                               textOutput("TRUEVarOutputIGE"),
                               textOutput("TRUEMeanOutputDGE"),
                               textOutput("TRUEVarOutputEnv_DGE"),
                               textOutput("TRUEMeanOutputDGE"),
                               textOutput("TRUEVarOutputEnv_IGE"),
                               textOutput("TRUEMeanOutputDGE"),
                               textOutput("summaryOutput"),
                               br(),
                               h4("Variance Table"),
                               div(class = "table-container",
                                   tableOutput("tableOutput")
                               ),
                               br(),
                               h4("Correlation Values"),
                               textOutput("TRUECorTRUE_DGE_IGE"),
                               textOutput("TRUECorTRUE_DGE_PRED_DGE"),
                               textOutput("TRUECorTRUE_IGE_PRED_IGE"),
                               textOutput("TRUECorPRED_DGE_PRED_IGE"),
                               br(),
                               h2("After selection"),
                               br(),
                               h3("Mass Selection"),
                               br(),
                               h3("Index selection")
                           ),
                           
                  )
                  
      )
    )
  )
)

# Server
server <- function(input, output) {
  observeEvent(input$goButton, {
    output$plotTRUE_DGE_IGE<-renderPlot(ggplot())
    output$plotTRUEvsPRED_DGE<-renderPlot(ggplot())
    output$plotTRUEvsPRED_IGE<-renderPlot(ggplot())
    output$plotPred_DGE_IGE<-renderPlot(ggplot())
    output$plotMass_differential<-renderPlot(ggplot())
    output$plotMass_Selection<-renderPlot(ggplot())
    output$plotIndex_Selection<-renderPlot(ggplot())
    
    if(input$setSeed) {
      set.seed(input$seedValue)
    }
    #Asreml=input$Asreml
    #Mean=input$Mean
    Asreml=FALSE
    Mean=FALSE
    if(Asreml){
      library(asreml)
      asreml.options(workspace="2gb",maxit=5,ai.sing=TRUE)
      Asreml=TRUE
    }
    else{
      Asreml=FALSE
    }
    
    # V_env_DGE=0.1
    # V_env_IGE=0
    # V_geno=0.75
    # V_voisin=V_geno
    # cor_geno_voisin=-0.6
    # cov_geno_voisin=cor_geno_voisin*sqrt(V_geno*V_voisin)
    # N_geno=25
    # N_rep=75
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
    N_row=N_col=0
    
    c=0
    if (ceiling(sqrt(N_geno*N_rep))!=sqrt(N_geno*N_rep)){
      N_row=N_col=ceiling(sqrt(N_geno*N_rep))
    }
    else {
      N_row=N_col=sqrt(N_geno*N_rep)
    }
    
    assign("N_geno",N_geno,envir = globalenv())
    assign("N_rep",N_rep,envir = globalenv())
    assign("N_row",N_row,envir = globalenv())
    assign("N_col",N_col,envir = globalenv())
    assign("Mean",Mean,envir = globalenv())
    assign("Asreml",Asreml,envir=globalenv())
    
    mu=c(0,0)
    assign("mu",mu,envir = globalenv())
    
    G=matrix(c(V_geno,cov_geno_voisin,
               cov_geno_voisin,V_voisin),
             ncol=2,nrow=2)
    
    assign("G",G,envir = globalenv())
    
    E=matrix(c(V_env_DGE,0,0,V_env_IGE), ncol=2,nrow=2)
    
    assign("E",E,envir = globalenv())
    
    P=G+E
    
    assign("P",P,envir = globalenv())
    
    df=mvrnorm(n=N_geno,mu,G)
    colnames(df)=c("DGE","IGE")
    df_E=mvrnorm(n=(N_geno*N_rep),mu,E)
    
    DGE=as.matrix(df[,1])
    assign("DGE",DGE,envir = globalenv())
    # var(DGE)
    
    IGE=as.matrix(df[,2])
    assign("IGE",IGE,envir = globalenv())
    # var(IGE)
    output$TRUEVarOutputDGE <- renderText({
      paste0("V(DGE) = ",round(var(DGE),3))
    })
    
    output$TRUEVarOutputIGE <- renderText({
      paste0("V(IGE) = ",round(var(IGE),3))
    })
    
    output$TRUEVarOutputEnv_DGE <- renderText({
      paste0("V(Env_DGE) = ",round(var(df_E[,1]),3))
    })
    
    output$TRUEVarOutputEnv_IGE <- renderText({
      paste0("V(Env_IGE) = ",round(var(df_E[,2]),3))
    })
    
    #output$TRUEMeanOutputDGE <- renderText({
    #  paste0("Mean(DGE) = ",round(Mean(DGE),3))
    #})
    
    #output$TRUEMeanOutputIGE <- renderText({
    #  paste0("Mean(IGE) = ",round(Mean(IGE),3))
    #})
    
    #output$TRUEMeanOutputEnv_DGE <- renderText({
    #  paste0("Mean(Env_DGE) = ",round(Mean(df_E[,1]),3))
    #})
    
    #output$TRUEMeanOutputEnv_IGE <- renderText({
    #  paste0("Mean(Env_IGE) = ",round(Mean(df_E[,2]),3))
    #})
    
    DATA=expand_grid("Focal"=as.character(paste0("G",sprintf("%03d", 1:N_geno))),"rep"=as.character(sprintf("%03d", 1:N_rep)))
    
    grid=expand_grid("Row"=factor(1:N_row),"Column"=factor(1:N_col))
    grid=grid[sample(1:(N_row*N_col),(N_geno*N_rep)),]
    
    mat_grid=matrix("vide",nrow=N_row+2,ncol=N_col+2)
    DATA=cbind(DATA,grid)
    
    for (i in 1:(nrow(mat_grid)-2)){
      for (j in 1:(ncol(mat_grid)-2)){
        if (!is_empty(DATA[(DATA$Row==i&DATA$Column==j),"Focal"])){
          mat_grid[i+1,j+1]=DATA[(DATA$Row==i&DATA$Column==j),"Focal"]
        }
      }
    }
    mat_grid[grep(mat_grid,pattern="vide")]=sample(DATA$Focal,length(mat_grid[grep(mat_grid,pattern="vide")]))
    
    matrice_voisin=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA$Focal)))
    DATA=cbind(matrice_voisin,DATA,data.frame("vide"=NA))
    
    for (i in 2:(nrow(mat_grid)-1)){
      for (j in 2:(ncol(mat_grid)-1)){
        if (!is_empty(DATA[(DATA$Row==(i-1)&DATA$Column==(j-1)),"Focal"])){
          DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i-1,j-1]]=DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i-1,j-1]]+1
          DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i,j-1]]=DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i,j-1]]+1
          DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i+1,j-1]]=DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i+1,j-1]]+1
          DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i-1,j]]=DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i-1,j]]+1
          DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i+1,j]]=DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i+1,j]]+1
          DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i-1,j+1]]=DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i-1,j+1]]+1
          DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i,j+1]]=DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i,j+1]]+1
          DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i+1,j+1]]=DATA[DATA$Row==(i-1)&DATA$Column==(j-1),mat_grid[i+1,j+1]]+1
        }
      }
    }
    
    Zg=model.matrix(~Focal-1,DATA)
    dimnames(Zg)[[2]]=paste0("G",sprintf("%03d", 1:N_geno))
    
    if(Mean==TRUE){
      Zv=as.matrix(DATA[,1:N_geno])/8
      DATA[,1:N_geno]=Zv
    }
    else{
      Zv=as.matrix(DATA[,1:N_geno])
    }
    
    assign("Zg",Zg,envir=globalenv())
    assign("Zv",Zv,envir=globalenv())
    assign("df_E",df_E,envir=globalenv())
    
    Pheno=Zg%*%DGE+Zv%*%IGE+df_E[,1]+df_E[,2]
    
    DATA$Pheno=as.vector(Pheno)
    DATA$Focal=as.factor(DATA$Focal)
    
    
    # mass phenotypic selection
    DATA$Pheno=as.vector(Pheno)
    DATA$Focal=as.factor(DATA$Focal)
    
    assign("DATA",DATA,envir = globalenv())
    assign("N_geno",N_geno,envir = globalenv())
    assign("Pheno",Pheno,envir = globalenv())
    # Estimation of BLUP 
    if (Asreml){
      Mod=asreml(fixed = Pheno~1,
                    random = ~str(~Focal+grp(Voisin),~us(2):id(Focal)),
                    group=list(Voisin=1:N_geno),
                    residual = ~units,
                    data=DATA)
      
      tmp_DGE=data.frame("DGE_pred"=summary(Mod,coef=TRUE)$coef.random[1:N_geno,1],
                         "Focal"=str_split(names(summary(Mod,coef=TRUE)$coef.random[1:N_geno,1]),pattern = "_",simplify = TRUE)[,2])
      tmp_IGE=data.frame("IGE_pred"=summary(Mod,coef=TRUE)$coef.random[(N_geno+1):(2*N_geno),1],
                         "Focal"=str_split(names(summary(Mod,coef=TRUE)$coef.random[1:N_geno,1]),pattern = "_",simplify = TRUE)[,2])
      pred=merge(tmp_DGE,tmp_IGE,by="Focal")
    }
    else{
      Mod=mmer(fixed = Pheno~1,
                  random= ~vsr(Focal)+vsr(Zv),
                  rcov = ~units,
                  data=DATA,nIters = 4 )
      
      DGE_pred=data.frame("Focal"=names(randef(Mod)$`u:Focal`$Pheno),"DGE_pred"=as.numeric(randef(Mod)$`u:Focal`[[1]]))
      IGE_pred=data.frame("Focal"=names(randef(Mod)$`u:Focal`$Pheno),"IGE_pred"=as.numeric(randef(Mod)$`u:Zv`[[1]]))
      pred=merge(DGE_pred,IGE_pred,by="Focal")
    }
    assign("pred",pred,envir = globalenv())
    assign("Mod",Mod,envir = globalenv())
    # summary(Mod)
    
    #plot(randef(Mod)$`u:Focal`$Pheno,randef(Mod)$`u:Zv`$Pheno)
    
    # plot(pred$DGE_pred,pred$IGE_pred)
    # plot(pred$DGE_pred,DGE)
    # plot(pred$IGE_pred,IGE)
    output$TRUECorTRUE_DGE_IGE <- renderText({
      paste0("r(TRUE DGE vs. IGE) = ",round(cor(DGE,IGE),2))
    })
    
    output$plotTRUEvsPRED_DGE <- renderPlot({
      plot(DGE,pred$DGE_pred, main = "TRUE vs PRED DGE")
    })
    
    output$TRUECorTRUE_DGE_PRED_DGE <- renderText({
      paste0("r(TRUE DGE vs. PRED DGE) = ",round(cor(DGE,pred$DGE_pred),2))
    })
    
    output$plotTRUEvsPRED_IGE <- renderPlot({
      plot(IGE,pred$IGE_pred, main = "TRUE vs PRED IGE")
    })
    
    output$TRUECorTRUE_IGE_PRED_IGE <- renderText({
      paste0("r(TRUE IGE vs. PRED IGE) = ",round(cor(IGE,pred$IGE_pred),2))
    })
    
    output$TRUECorPRED_DGE_PRED_IGE <- renderText({
      paste0("r(PRED DGE vs. PRED IGE) = ",round(cor(pred$DGE_pred,pred$IGE_pred),2))
    })
    
    
    output$plotTRUE_DGE_IGE <- renderPlot({
      ggplot() +
        geom_point(data = as.data.frame(df), aes(x = DGE, y = IGE), color = "gray") +  # General points
        xlab("DGE") +
        ylab("IGE") +
        ggtitle("True DGE vs. IGE")+
        theme_bw()
    })
    output$plotPred_DGE_IGE <- renderPlot({
      ggplot() +
        geom_point(data = as.data.frame(pred), aes(x = DGE_pred, y = IGE_pred), color = "gray") +  # General points
        xlab("DGE_pred") +
        ylab("IGE_pred") +
        ggtitle("DGE vs. IGE")+
        theme_bw()
    })
    
    table=summary(Mod)$varcomp
    table$Effect=c("DGE","IGE","Env")
    table=table[,c(5,1:4)]
    
    # Display the mean of phenotypes as an example of summary output
    output$tableOutput <- renderTable({
      table
    })
  })
  observeEvent(input$SelButton,{
    b_DGE=input$b_DGE
    p<-input$p
    
    #####Index selection
    
    pred$I <- b_DGE*pred$DGE_pred+(1-b_DGE)*pred$IGE_pred
    
    cor(pred$DGE_pred, pred$I)
    cor(pred$IGE_pred, pred$I)
    
    # list of selected genotypes (numbers will be different from mass phenotype selection)
    sel = which(pred$I>quantile(pred$I,1-p))
    length(sel)
    
    R<-c(mean(DGE[sel]),mean(IGE[sel]))
    
    # R=G%*%solve(P)%*%S
    
    mus=c(mu[1]+R[1],mu[2]+R[2])
    
    df_sel_I=mvrnorm(N_geno,mus,G)
    colnames(df_sel_I)=c("DGE","IGE")
    
    DGE_sel_I=as.matrix(df_sel_I[,1])
    IGE_sel_I=as.matrix(df_sel_I[,2])
    df_E_sel_I=mvrnorm(n=(N_geno*N_rep),mu,E)
    
    DATA_sel_I=expand_grid("Focal_sel_I"=as.character(paste0("G_sel_I",sprintf("%03d", 1:N_geno))),"rep"=as.character(sprintf("%03d", 1:N_rep)))
    
    grid_sel_I=expand_grid("Row"=factor(1:N_row),"Column"=factor(1:N_col))
    grid_sel_I=grid_sel_I[sample(1:(N_row*N_col),(N_geno*N_rep)),]
    
    mat_grid_sel_I=matrix("vide",nrow=N_row+2,ncol=N_col+2)
    
    DATA_sel_I=cbind(DATA_sel_I,grid_sel_I)
    
    for (i in 1:(nrow(mat_grid_sel_I)-2)){
      for (j in 1:(ncol(mat_grid_sel_I)-2)){
        if (!is_empty(DATA_sel_I[(DATA_sel_I$Row==i&DATA_sel_I$Column==j),"Focal_sel_I"])){
          mat_grid_sel_I[i+1,j+1]=DATA_sel_I[(DATA_sel_I$Row==i&DATA_sel_I$Column==j),"Focal_sel_I"]
        }
      }
    }
    mat_grid_sel_I[grep(mat_grid_sel_I,pattern="vide")]=sample(DATA_sel_I$Focal,length(mat_grid_sel_I[grep(mat_grid_sel_I,pattern="vide")]))
    
    matrice_voisin=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA_sel_I$Focal)))
    DATA_sel_I=cbind(matrice_voisin,DATA_sel_I,data.frame("vide"=NA))
    
    for (i in 2:(nrow(mat_grid_sel_I)-1)){
      for (j in 2:(ncol(mat_grid_sel_I)-1)){
        if (!is_empty(DATA_sel_I[(DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1)),"Focal_sel_I"])){
          DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i-1,j-1]]=DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i-1,j-1]]+1
          DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i,j-1]]=DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i,j-1]]+1
          DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i+1,j-1]]=DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i+1,j-1]]+1
          DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i-1,j]]=DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i-1,j]]+1
          DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i+1,j]]=DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i+1,j]]+1
          DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i-1,j+1]]=DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i-1,j+1]]+1
          DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i,j+1]]=DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i,j+1]]+1
          DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i+1,j+1]]=DATA_sel_I[DATA_sel_I$Row==(i-1)&DATA_sel_I$Column==(j-1),mat_grid_sel_I[i+1,j+1]]+1
        }
      }
    }
    
    Zg_sel_I=model.matrix(~Focal_sel_I-1,DATA_sel_I)
    dimnames(Zg_sel_I)[[2]]=paste0("G_sel_I",sprintf("%03d", 1:N_geno))
    
    if(Mean==TRUE){
      Zv_sel_I=as.matrix(DATA_sel_I[,1:N_geno])/8
    }
    else{
      Zv_sel_I=as.matrix(DATA_sel_I[,1:N_geno])
    }
    
    Pheno_sel_I=Zg_sel_I%*%DGE_sel_I+Zv_sel_I%*%IGE_sel_I+df_E_sel_I[,1]+df_E_sel_I[,2]
    
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
    
    
    # find the corresponding genotype
    geno_sel<- as.numeric(sapply(as.character(DATA$Focal[list_sel]), function(x) substring(x, 2)) )
    sizeBubble<-table(geno_sel)
    
    # plot(df[,1], df[,2], col="grey")
    # points(df[geno_sel,1], df[geno_sel,2], col="red", cex=sizeBubble)
    
    # Selection differential
    S<-c()
    S[1] <- mean((Zg%*%DGE + df_E[,1])[list_sel])
    S[2] <- mean((Zg%*%IGE + df_E[,2])[list_sel])
    
    # Genetic gain
    R=G%*%solve(P)%*%S
    mus=c(mu[1]+R[1],mu[2]+R[2])
    
    df_sel=mvrnorm(N_geno,mus,G)
    colnames(df_sel)=c("DGE","IGE")
    
    DGE_sel=as.matrix(df_sel[,1])
    IGE_sel=as.matrix(df_sel[,2])
    df_E_sel=mvrnorm(n=(N_geno*N_rep),mu,E)
    
    DATA_sel=expand_grid("Focal_sel"=as.character(paste0("G_sel",sprintf("%03d", 1:N_geno))),"rep"=as.character(sprintf("%03d", 1:N_rep)))
    
    grid_sel=expand_grid("Row"=factor(1:N_row),"Column"=factor(1:N_col))
    grid_sel=grid_sel[sample(1:(N_row*N_col),(N_geno*N_rep)),]
    
    mat_grid_sel=matrix("vide",nrow=N_row+2,ncol=N_col+2)
    
    DATA_sel=cbind(DATA_sel,grid_sel)
    
    for (i in 1:(nrow(mat_grid_sel)-2)){
      for (j in 1:(ncol(mat_grid_sel)-2)){
        if (!is_empty(DATA_sel[(DATA_sel$Row==i&DATA_sel$Column==j),"Focal_sel"])){
          mat_grid_sel[i+1,j+1]=DATA_sel[(DATA_sel$Row==i&DATA_sel$Column==j),"Focal_sel"]
        }
      }
    }
    mat_grid_sel[grep(mat_grid_sel,pattern="vide")]=sample(DATA_sel$Focal,length(mat_grid_sel[grep(mat_grid_sel,pattern="vide")]))
    
    matrice_voisin=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA_sel$Focal)))
    DATA_sel=cbind(matrice_voisin,DATA_sel,data.frame("vide"=NA))
    
    for (i in 2:(nrow(mat_grid_sel)-1)){
      for (j in 2:(ncol(mat_grid_sel)-1)){
        if (!is_empty(DATA_sel[(DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1)),"Focal_sel"])){
          DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i-1,j-1]]=DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i-1,j-1]]+1
          DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i,j-1]]=DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i,j-1]]+1
          DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i+1,j-1]]=DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i+1,j-1]]+1
          DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i-1,j]]=DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i-1,j]]+1
          DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i+1,j]]=DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i+1,j]]+1
          DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i-1,j+1]]=DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i-1,j+1]]+1
          DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i,j+1]]=DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i,j+1]]+1
          DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i+1,j+1]]=DATA_sel[DATA_sel$Row==(i-1)&DATA_sel$Column==(j-1),mat_grid_sel[i+1,j+1]]+1
        }
      }
    }
    
    Zg_sel=model.matrix(~Focal_sel-1,DATA_sel)
    dimnames(Zg_sel)[[2]]=paste0("G_sel",sprintf("%03d", 1:N_geno))
    
    if(Mean==TRUE){
      Zv_sel=as.matrix(DATA_sel[,1:N_geno])/8
    }
    else{
      Zv_sel=as.matrix(DATA_sel[,1:N_geno])
    }
    
    Pheno_sel=Zg_sel%*%DGE_sel+Zv_sel%*%IGE_sel+df_E_sel[,1]+df_E_sel[,2]
    
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
        theme(legend.title = element_blank())+# Remove the legend title
        annotate("text",x=mean_after_sel+1.5,y=0.3,label=paste0('Delta_mu_pheno = ', mean_after_sel - mean_before_sel))# Remove the legend title
      
    })
    
    output$summaryOutput <- renderText({
      paste("Summary of calculations...")
    })
    
  })
    
}

# Run the Shiny application
shinyApp(ui = ui, server = server)
