library(shiny)

# Load the necessary packages
if(!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if(!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
if(!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if(!requireNamespace("sommer", quietly = TRUE)) install.packages("sommer")
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if(!requireNamespace("shinycssloaders", quietly = TRUE)) install.packages("shinycssloaders")


library(MASS) # For mvrnorm
library(ggplot2) # For plots
library(sommer)
library(data.table)
library(tidyverse)
library(shinycssloaders)


# UI
ui <- fluidPage(
  titlePanel("Population Evolution Simulation"),
  sidebarLayout(
    sidebarPanel(
      actionButton("goButton", "Run Simulation"),
      sliderInput("N", "Number of genotypes (N)", min = 10, max = 500, value = 450),
      sliderInput("rep", "Number of rep per genotype  (rep)", min = 1, max = 20, value = 2),
      numericInput("varG11", "Genetic variance DGE", value = 0.75),
      numericInput("varG22", "Genetic variance IGE", value = 0.225),
      sliderInput("r", "Genetic correlation DGE : IGE", min = -1, max = 1, value = -0.6, step = 0.1),
      numericInput("varE11", "Environmental variance DGE", value = 1),
      numericInput("varE22", "Environmental variance IGE", value = 0.07),
      sliderInput("p", "Selection pressure", min = 0.01, max = 1, value = 0.1, step = 0.1),          
      sliderInput("b_DGE", "Index weight for DGE (weight for IGE = 1 - weight for DGE)", min = 0, max = 1, value = 0.5, step = 0.1),
      checkboxInput("setSeed", "Set seed for reproducibility", value = FALSE),
      numericInput("seedValue", "Seed value", value = 12345)
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Graph", 
                           plotOutput("plotTRUE_DGE_IGE"),
                           plotOutput("plotTRUEvsPRED_DGE"),
                           plotOutput("plotTRUEvsPRED_IGE"),
                           plotOutput("plotMass_differential"),
                           plotOutput("plotMass_Selection"),
                           plotOutput("plotIndex_Selection"),
    ),
                  
                  tabPanel("Summary", verbatimTextOutput("summaryOutput"))
      )
    )
  )
)

# Server
server <- function(input, output) {
  observeEvent(input$goButton, {
    if(input$setSeed) {
      set.seed(input$seedValue)
    }

    # V_env_DGE=0.1
    # V_env_IGE=0
    # V_geno=0.75
    # V_voisin=V_geno
    # cor_geno_voisin=-0.6
    # cov_geno_voisin=cor_geno_voisin*sqrt(V_geno*V_voisin)
    # N_geno=450
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
    b_DGE=input$b_DGE
    p<-input$p
    
    mu=c(0,0)
    
    G=matrix(c(V_geno,cov_geno_voisin,
                   cov_geno_voisin,V_voisin),
                   ncol=2,nrow=2)

    E=matrix(c(V_env_DGE,0,0,V_env_IGE), ncol=2,nrow=2)

    P=G+E
    
    df=mvrnorm(n=N_geno,mu,G)
    colnames(df)=c("DGE","IGE")

    Pheno=mvrnorm(n=(N_geno*N_rep),mu,E)
    
    DGE=as.matrix(df[,1])
    
    # var(DGE)

    IGE=as.matrix(df[,2])
    # var(IGE)
    
    DATA=expand_grid("Focal"=as.character(paste0("G",1:(N_geno))),"rep"=as.character(1:N_rep))

    grid=expand_grid("Row"=factor(1:sqrt(N_geno*N_rep)),"Column"=factor(1:sqrt(N_geno*N_rep)))
    grid=grid[sample(1:(N_geno*N_rep),(N_geno*N_rep)),]
    
    mat_grid=matrix("vide",nrow=sqrt(N_geno*N_rep)+2,ncol=sqrt(N_geno*N_rep)+2)
    DATA=cbind(DATA,grid)
    
    for (i in 1:(nrow(mat_grid)-2)){
      for (j in 1:(ncol(mat_grid)-2)){
        mat_grid[i+1,j+1]=DATA[(DATA$Row==i&DATA$Column==j),"Focal"]
      }
    }

    mat_grid[1,]=sample(DATA$Focal,length(mat_grid[1,]))
    mat_grid[sqrt(N_geno*N_rep)+2,]=sample(DATA$Focal,length(mat_grid[sqrt(N_geno*N_rep)+2,]))
    mat_grid[,1]=sample(DATA$Focal,length(mat_grid[,1]))
    mat_grid[,sqrt(N_geno*N_rep)+2]=sample(DATA$Focal,length(mat_grid[,sqrt(N_geno*N_rep)+2]))
    
    matrice_voisin=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA$Focal)))
    DATA=cbind(matrice_voisin,DATA,data.frame("vide"=NA))
    
    for (i in 2:(nrow(mat_grid)-1)){
      for (j in 2:(nrow(mat_grid)-1)){
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
    
    Zg=model.matrix(~Focal-1,DATA)
    Zv=as.matrix(DATA[,1:N_geno])
    
    Pheno=Zg%*%DGE+Zv%*%IGE+df_E[,1]+df_E[,2]
    
    DATA$Pheno=as.vector(Pheno)
    DATA$Focal=as.factor(DATA$Focal)
    
    Zv=Zv[,levels(DATA$Focal)]
 
    
    # mass phenotypic selection
    DATA$Pheno=as.vector(Pheno)
    DATA$Focal=as.factor(DATA$Focal)
    
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
    S[2] <- mean((Zv%*%IGE + df_E[,2])[list_sel])
    
    # Genetic gain  
    R=G%*%solve(P)%*%S
    mus=c(mu[1]+R[1],mu[2]+R[2])
        
    df_sel=mvrnorm(N_geno,mus,G)
    colnames(df_sel)=c("DGE","IGE")
    
    DGE_sel=as.matrix(df_sel[,1])
    IGE_sel=as.matrix(df_sel[,2])
    df_E_sel=mvrnorm(n=(N_geno*N_rep),mu,E)
    
    DATA_sel=expand_grid("Focal_sel"=as.character(paste0("G_sel",1:(N_geno))),"rep"=as.character(1:N_rep))
    
    grid_sel=expand_grid("Row"=factor(1:sqrt(N_geno*N_rep)),"Column"=factor(1:sqrt(N_geno*N_rep)))
    grid_sel=grid_sel[sample(1:(N_geno*N_rep),(N_geno*N_rep)),]
    
    mat_grid_sel=matrix("vide",nrow=sqrt(N_geno*N_rep)+2,ncol=sqrt(N_geno*N_rep)+2)
    
    DATA_sel=cbind(DATA_sel,grid_sel)
    
    for (i in 1:(nrow(mat_grid_sel)-2)){
      for (j in 1:(ncol(mat_grid_sel)-2)){
        mat_grid_sel[i+1,j+1]=DATA_sel[(DATA_sel$Row==i&DATA_sel$Column==j),"Focal_sel"]
      }
    }
    #mat_grid[1,]=sample(DATA$Focal,length(mat_grid[1,]))
    #mat_grid[sqrt(N_geno*N_rep)+2,]=sample(DATA$Focal,length(mat_grid[sqrt(N_geno*N_rep)+2,]))
    #mat_grid[,1]=sample(DATA$Focal,length(mat_grid[,1]))
    #mat_grid[,sqrt(N_geno*N_rep)+2]=sample(DATA$Focal,length(mat_grid[,sqrt(N_geno*N_rep)+2]))
    
    matrice_voisin_sel=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA_sel$Focal_sel)))
    DATA_sel=cbind(matrice_voisin_sel,DATA_sel,data.frame("vide"=NA))
    
    for (i in 2:(nrow(mat_grid_sel)-1)){
      for (j in 2:(nrow(mat_grid_sel)-1)){
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
    
    Zg_sel=model.matrix(~Focal_sel-1,DATA_sel)
    Zv_sel=as.matrix(DATA_sel[,1:N_geno])
    Pheno_sel=Zg_sel%*%DGE_sel+Zv_sel%*%IGE_sel+df_E_sel[,1]+df_E_sel[,2]
    
    # Combine the two vectors into a dataframe
    combinedData <- rbind(data.frame(Value = Pheno, Phase = "Before selection"),
                          data.frame(Value = Pheno_sel, Phase = "After selection"))
    
    output$plotMass_Selection <- renderPlot({
      ggplot(combinedData, aes(x = Value, fill = Phase)) +
        geom_density(alpha = 0.5, adjust = 1.5) + # Adjustment for smoothing, 'adjust' controls the smoothing
        scale_fill_manual(values = c("Before selection" = "blue", "After selection" = "red")) +
        labs(x = "Phenotypic Value", y = "Density", title = "Phenotypic mass selection") +
        theme_minimal() +
        theme(legend.title = element_blank())  # Remove the legend title
    })
    
  
    
    
    ################# 
    # Estimation of BLUP 
    
    modèle_sommer_no_cov=mmer(fixed = Pheno~1,
                              random= ~vsr(Focal)+vsr(as.matrix(DATA[,1:N_geno])),
                              rcov = ~units,
                              data=DATA,nIters = 4 )
    
    # summary(modèle_sommer_no_cov)
    
    # plot(randef(modèle_sommer_no_cov)$`u:Focal`$Pheno,randef(modèle_sommer_no_cov)$`u:DATA:N_geno`$Pheno)
    
    DGE_pred=data.frame("Focal"=names(randef(modèle_sommer_no_cov)$`u:Focal`$Pheno),"DGE_pred"=as.numeric(randef(modèle_sommer_no_cov)$`u:Focal`[[1]]))
    IGE_pred=data.frame("Focal"=names(randef(modèle_sommer_no_cov)$`u:Focal`$Pheno),"IGE_pred"=as.numeric(randef(modèle_sommer_no_cov)$`u:DATA:N_geno`[[1]]))
    
    pred=merge(DGE_pred,IGE_pred,by="Focal")
    
    # plot(pred$DGE_pred,pred$IGE_pred)
    # plot(pred$DGE_pred,DGE)
    # plot(pred$IGE_pred,IGE)
    
    output$plotTRUEvsPRED_DGE <- renderPlot({
      plot(DGE,pred$DGE_pred, main = "TRUE vs PRED DGE")
    })
    
    output$plotTRUEvsPRED_IGE <- renderPlot({
      plot(IGE,pred$IGE_pred, main = "TRUE vs PRED IGE")
    })
    
    output$plotTRUE_DGE_IGE <- renderPlot({
      ggplot() +
        geom_point(data = as.data.frame(df), aes(x = DGE, y = IGE), color = "gray") +  # General points
        xlab("DGE") +
        ylab("IGE") +
        ggtitle("True DGE vs. IGE")+
        theme_bw()
    })


    
    # Index selection index
    
    b_DGE<-1
    pred$I <- b_DGE*pred$DGE_pred+(1-b_DGE)*pred$IGE_pred

    cor(pred$DGE_pred, pred$I)
    cor(pred$IGE_pred, pred$I)
    cor(pred$IGE_pred, pred$DGE_pred)
    
    hist(pred$I)
    
    # list of selected genotypes (numbers will be different from mass phenotype selection)
    sel = which(pred$I>quantile(pred$I,1-p))
    length(sel)
    
    R[1]<-mean(DGE[sel])
    R[2]<-mean(IGE[sel])
  
    # R=G%*%solve(P)%*%S
    
    mus=c(mu[1]+R[1],mu[2]+R[2])
    
    df_sel=mvrnorm(N_geno,mus,G)
    colnames(df_sel)=c("DGE","IGE")
    
    DGE_sel=as.matrix(df_sel[,1])
    IGE_sel=as.matrix(df_sel[,2])
    df_E_sel=mvrnorm(n=(N_geno*N_rep),mu,E)
    
    DATA_sel=expand_grid("Focal_sel"=as.character(paste0("G_sel",1:(N_geno))),"rep"=as.character(1:N_rep))
    
    grid_sel=expand_grid("Row"=factor(1:sqrt(N_geno*N_rep)),"Column"=factor(1:sqrt(N_geno*N_rep)))
    grid_sel=grid_sel[sample(1:(N_geno*N_rep),(N_geno*N_rep)),]
    
    mat_grid_sel=matrix("vide",nrow=sqrt(N_geno*N_rep)+2,ncol=sqrt(N_geno*N_rep)+2)
    
    DATA_sel=cbind(DATA_sel,grid_sel)
    
    for (i in 1:(nrow(mat_grid_sel)-2)){
      for (j in 1:(ncol(mat_grid_sel)-2)){
        mat_grid_sel[i+1,j+1]=DATA_sel[(DATA_sel$Row==i&DATA_sel$Column==j),"Focal_sel"]
      }
    }
    #mat_grid[1,]=sample(DATA$Focal,length(mat_grid[1,]))
    #mat_grid[sqrt(N_geno*N_rep)+2,]=sample(DATA$Focal,length(mat_grid[sqrt(N_geno*N_rep)+2,]))
    #mat_grid[,1]=sample(DATA$Focal,length(mat_grid[,1]))
    #mat_grid[,sqrt(N_geno*N_rep)+2]=sample(DATA$Focal,length(mat_grid[,sqrt(N_geno*N_rep)+2]))
    
    matrice_voisin_sel=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(DATA_sel$Focal_sel)))
    DATA_sel=cbind(matrice_voisin_sel,DATA_sel,data.frame("vide"=NA))
    
    for (i in 2:(nrow(mat_grid_sel)-1)){
      for (j in 2:(nrow(mat_grid_sel)-1)){
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
    
    Zg_sel=model.matrix(~Focal_sel-1,DATA_sel)
    Zv_sel=as.matrix(DATA_sel[,1:N_geno])
    Pheno_sel=Zg_sel%*%DGE_sel+Zv_sel%*%IGE_sel+df_E_sel[,1]+df_E_sel[,2]
  
    output$plotIndex_Selection <- renderPlot({
      ggplot(combinedData, aes(x = Value, fill = Phase)) +
        geom_density(alpha = 0.5, adjust = 1.5) + # Adjustment for smoothing, 'adjust' controls the smoothing
        scale_fill_manual(values = c("Before selection" = "blue", "After selection" = "red")) +
        labs(x = "Phenotypic Value", y = "Density", title = "Index selection") +
        theme_minimal() +
        theme(legend.title = element_blank())  # Remove the legend title
    })
    

    
  output$plotGenotype <- renderPlot({
      ggplot() +
        geom_point(data = as.data.frame(DATA), aes(x = DGE, y = IGE), color = "gray") +  # General points
        geom_point(data = as.data.frame(df_sel), aes(x = DGE, y = IGE), color = "red") +  # Special points
        xlab("DGE") +
        ylab("IGE") +
        ggtitle("Genotype Distribution")+
        theme_bw()
    })
    
    
    
    
    output$plotGenotype <- renderPlot({
      ggplot() +
        geom_point(data = as.data.frame(df), aes(x = DGE, y = IGE), color = "gray") +  # General points
        geom_point(data = as.data.frame(df_sel), aes(x = DGE, y = IGE), color = "red") +  # Special points
        xlab("DGE") +
        ylab("IGE") +
        ggtitle("Genotype Distribution")+
        theme_bw()
    })
    
    output$summaryOutput <- renderText({
      paste("Summary of calculations...")
    })
    
    # Combine the two vectors into a dataframe
    combinedData <- rbind(data.frame(Value = Pheno, Phase = "Before selection"),
                          data.frame(Value = Pheno_sel, Phase = "After selection"))
    
    output$plotPhenotype <- renderPlot({
      ggplot(combinedData, aes(x = Value, fill = Phase)) +
        geom_density(alpha = 0.5, adjust = 1.5) + # Adjustment for smoothing, 'adjust' controls the smoothing
        scale_fill_manual(values = c("Before selection" = "blue", "After selection" = "red")) +
        labs(x = "Phenotypic Value", y = "Density", title = "Distribution of Phenotypic Values Before and After Selection") +
        theme_minimal() +
        theme(legend.title = element_blank())  # Remove the legend title
    })
    
    # Display the mean of phenotypes as an example of summary output
    output$summaryOutput <- renderText({
      paste("Mean phenotypes:", mean(Pheno))
    })
  })
}

# Run the Shiny application
shinyApp(ui = ui, server = server)
