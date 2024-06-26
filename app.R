library(shiny)
library(MASS) # For mvrnorm
library(ggplot2) # For plots
library(sommer)
library(data.table)
library(tidyverse)

# UI
ui <- fluidPage(
  titlePanel("Population Evolution Simulation"),
  sidebarLayout(
    sidebarPanel(
      actionButton("goButton", "Run Simulation"),
      sliderInput("N", "Number of genotypes (N)", min = 10, max = 100, value = 20),
      sliderInput("rep", "Number of rep per genotype  (rep)", min = 1, max = 20, value = 5),
      numericInput("varG11", "Genetic variance DGE", value = 0.75),
      numericInput("varG22", "Genetic variance IGE", value = 0.07),
      sliderInput("r", "Genetic correlation DGE : IGE", min = -1, max = 1, value = 0, step = 0.1),
      numericInput("varE11", "Environmental variance DGE", value = 0.7),
      numericInput("varE22", "Environmental variance IGE", value = 0.07),
      sliderInput("p", "Selection pressure", min = 0.01, max = 1, value = 0.5, step = 0.1),          
      sliderInput("b_DGE", "Index weight for DGE (weight for IGE = 1 - weight for DGE)", min = 0, max = 1, value = 0.5, step = 0.1),
      checkboxInput("setSeed", "Set seed for reproducibility", value = FALSE),
      numericInput("seedValue", "Seed value", value = 12345)
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Graph", plotOutput("plotGenotype"), plotOutput("plotPhenotype")),
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

    V_env_DGE=1
    V_env_IGE=0
    V_geno=0.75
    V_voisin=0.3*V_geno
    cor_geno_voisin=-0.6
    cov_geno_voisin=cor_geno_voisin*sqrt(V_geno*V_voisin)
    N_geno=450
    N_rep=2
    b_DGE=0.5
    p=0.1
    
    V_env_DGE=input$varE11
    V_env_IGE=input$varE22
    V_geno=input$varG11
    V_voisin=input$varG22
    cor_geno_voisin=input$r
    cov_geno_voisin=cor_geno_voisin*sqrt(V_geno*V_voisin)
    N_geno=input$N
    N_rep=input$rep
    b_DGE=input$b_DGE
    
    mu=c(0,0)
    G=matrix(c(V_geno,cov_geno_voisin,
                   cov_geno_voisin,V_voisin),
                 ncol=2,nrow=2)
    E=matrix(c(V_env_DGE,0,
               0,V_env_IGE),
             ncol=2,nrow=2)
    P=G+E
    
    df=mvrnorm(n=N_geno,mu,G)
    colnames(df)=c("DGE","IGE")
    df_E=mvrnorm(n=(N_geno*N_rep),mu,E)
    
    DGE=as.matrix(df[,1])
    var(DGE)
    IGE=as.matrix(df[,2])
    var(IGE)
    
    données=expand_grid("Focal"=as.character(paste0("G",1:(N_geno))),"rep"=as.character(1:N_rep))
    grid=expand_grid("Row"=factor(1:sqrt(N_geno*N_rep)),"Column"=factor(1:sqrt(N_geno*N_rep)))
    grid=grid[sample(1:(N_geno*N_rep),(N_geno*N_rep)),]
    mat_grid=matrix("vide",nrow=sqrt(N_geno*N_rep)+2,ncol=sqrt(N_geno*N_rep)+2)
    données=cbind(données,grid)
    
    for (i in 1:(nrow(mat_grid)-2)){
      for (j in 1:(ncol(mat_grid)-2)){
        mat_grid[i+1,j+1]=données[(données$Row==i&données$Column==j),"Focal"]
      }
    }
    mat_grid[1,]=sample(données$Focal,length(mat_grid[1,]))
    mat_grid[sqrt(N_geno*N_rep)+2,]=sample(données$Focal,length(mat_grid[sqrt(N_geno*N_rep)+2,]))
    mat_grid[,1]=sample(données$Focal,length(mat_grid[,1]))
    mat_grid[,sqrt(N_geno*N_rep)+2]=sample(données$Focal,length(mat_grid[,sqrt(N_geno*N_rep)+2]))
    
    matrice_voisin=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(données$Focal)))
    données=cbind(matrice_voisin,données,data.frame("vide"=NA))
    
    for (i in 2:(nrow(mat_grid)-1)){
      for (j in 2:(nrow(mat_grid)-1)){
        données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i-1,j-1]]=données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i-1,j-1]]+1
        données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i,j-1]]=données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i,j-1]]+1
        données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i+1,j-1]]=données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i+1,j-1]]+1
        données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i-1,j]]=données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i-1,j]]+1
        données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i+1,j]]=données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i+1,j]]+1
        données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i-1,j+1]]=données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i-1,j+1]]+1
        données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i,j+1]]=données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i,j+1]]+1
        données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i+1,j+1]]=données[données$Row==(i-1)&données$Column==(j-1),mat_grid[i+1,j+1]]+1
      }
    }
    
    Zg=model.matrix(~Focal-1,données)
    Zv=as.matrix(données[,1:N_geno])
    
    Pheno=Zg%*%DGE+Zv%*%IGE+df_E[,1]+df_E[,2]
    données$Pheno=as.vector(Pheno)
    données$Focal=as.factor(données$Focal)
    Zv=Zv[,levels(données$Focal)]
  
    
    modèle_sommer_no_cov=mmer(fixed = Pheno~1,
                              random= ~vsr(Focal)+vsr(as.matrix(données[,1:N_geno])),
                              rcov = ~units,
                              data=données,nIters = 4 )
    
    summary(modèle_sommer_no_cov)
    
    plot(randef(modèle_sommer_no_cov)$`u:Focal`$Pheno,randef(modèle_sommer_no_cov)$`u:données:N_geno`$Pheno)
    
    DGE_pred=data.frame("Focal"=names(randef(modèle_sommer_no_cov)$`u:Focal`$Pheno),"DGE_pred"=as.numeric(randef(modèle_sommer_no_cov)$`u:Focal`[[1]]))
    IGE_pred=data.frame("Focal"=names(randef(modèle_sommer_no_cov)$`u:Focal`$Pheno),"IGE_pred"=as.numeric(randef(modèle_sommer_no_cov)$`u:données:N_geno`[[1]]))
    pred=merge(DGE_pred,IGE_pred,by="Focal")
    plot(pred$DGE_pred,pred$IGE_pred)
    plot(pred$DGE_pred,DGE)
    plot(pred$IGE_pred,IGE)
    
    pred$I=b_DGE*pred$DGE_pred+(1-b_DGE)*pred$IGE_pred
    sel = pred[pred$I>quantile(pred$I,1-p),]
    
    S=c(mean(sel$DGE_pred),mean(sel$IGE_pred))
    
    R=G%*%solve(P)%*%S
    
    mus=c(mu[1]+R[1],mu[2]+R[2])
    
    df_sel=mvrnorm(N_geno,mus,G)
    colnames(df_sel)=c("DGE","IGE")
    
    DGE_sel=as.matrix(df_sel[,1])
    var(DGE_sel)
    IGE_sel=as.matrix(df_sel[,2])
    var(IGE_sel)
    df_E_sel=mvrnorm(n=(N_geno*N_rep),mu,E)
    
    données_sel=expand_grid("Focal_sel"=as.character(paste0("G_sel",1:(N_geno))),"rep"=as.character(1:N_rep))
    grid_sel=expand_grid("Row"=factor(1:sqrt(N_geno*N_rep)),"Column"=factor(1:sqrt(N_geno*N_rep)))
    grid_sel=grid_sel[sample(1:(N_geno*N_rep),(N_geno*N_rep)),]
    mat_grid_sel=matrix("vide",nrow=sqrt(N_geno*N_rep)+2,ncol=sqrt(N_geno*N_rep)+2)
    données_sel=cbind(données_sel,grid_sel)
    
    for (i in 1:(nrow(mat_grid_sel)-2)){
      for (j in 1:(ncol(mat_grid_sel)-2)){
        mat_grid_sel[i+1,j+1]=données_sel[(données_sel$Row==i&données_sel$Column==j),"Focal_sel"]
      }
    }
    #mat_grid[1,]=sample(données$Focal,length(mat_grid[1,]))
    #mat_grid[sqrt(N_geno*N_rep)+2,]=sample(données$Focal,length(mat_grid[sqrt(N_geno*N_rep)+2,]))
    #mat_grid[,1]=sample(données$Focal,length(mat_grid[,1]))
    #mat_grid[,sqrt(N_geno*N_rep)+2]=sample(données$Focal,length(mat_grid[,sqrt(N_geno*N_rep)+2]))
    
    matrice_voisin_sel=matrix(0,nrow = N_geno*N_rep,ncol=N_geno,dimnames = list(1:(N_geno*N_rep),unique(données_sel$Focal_sel)))
    données_sel=cbind(matrice_voisin_sel,données_sel,data.frame("vide"=NA))
    
    for (i in 2:(nrow(mat_grid_sel)-1)){
      for (j in 2:(nrow(mat_grid_sel)-1)){
        données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i-1,j-1]]=données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i-1,j-1]]+1
        données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i,j-1]]=données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i,j-1]]+1
        données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i+1,j-1]]=données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i+1,j-1]]+1
        données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i-1,j]]=données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i-1,j]]+1
        données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i+1,j]]=données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i+1,j]]+1
        données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i-1,j+1]]=données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i-1,j+1]]+1
        données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i,j+1]]=données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i,j+1]]+1
        données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i+1,j+1]]=données_sel[données_sel$Row==(i-1)&données_sel$Column==(j-1),mat_grid_sel[i+1,j+1]]+1
      }
    }
    
    Zg_sel=model.matrix(~Focal_sel-1,données_sel)
    Zv_sel=as.matrix(données_sel[,1:N_geno])
    Pheno_sel=Zg_sel%*%DGE_sel+Zv_sel%*%IGE_sel+df_E_sel[,1]+df_E_sel[,2]
    
    output$plotGenotype <- renderPlot({
      ggplot() +
        geom_density_2d(data = as.data.frame(df), aes(x = DGE, y = IGE), color = "gray") +  # General points
        geom_density_2d(data = as.data.frame(df_sel), aes(x = DGE, y = IGE), color = "red") +  # Special points
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
