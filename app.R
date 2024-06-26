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
      sliderInput("N", "Number of individuals (N)", min = 10, max = 100, value = 50),
      sliderInput("rep", "Number of repetitions per genotypes (rep)", min = 2, max = 50, value = 50),
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
    N_geno=50
    N_rep=50
    b_DGE=0.5
    
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
    
    données=expand_grid("Focal"=as.character(paste0("G",1:(N_geno))),"rep"=as.character(paste0("G",1:N_rep)))
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
    
    modèle_sommer_no_cov=mmer(fixed = Pheno~1,
                              random= ~vsr(Focal)+vsr(Zv),
                              rcov = ~units,
                              data=données,nIters = 4 )
    
    summary(modèle_sommer_no_cov)
    
    modèle_sommer_cov=mmer(fixed = Pheno~1,
                                random= ~gvsr(Focal,Zv),
                                rcov = ~units,
                                data=données,nIters = 4)
    
    summary(modèle_sommer_cov)
    
    DGE_pred=data.frame("Focal"=names(randef(modèle_sommer_no_cov)$`u:Focal`$Pheno),"DGE_pred"=as.numeric(randef(modèle_sommer_no_cov)$`u:Focal`[[1]]))
    IGE_pred=data.frame("Focal"=names(randef(modèle_sommer_no_cov)$`u:Zv`$Pheno),"IGE_pred"=as.numeric(randef(modèle_sommer_no_cov)$`u:Zv`[[1]]))
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
    
    ggplot() +
      geom_density_2d(data = as.data.frame(df), aes(x = DGE, y = IGE), color = "gray") +  # General points
      geom_density_2d(data = as.data.frame(df_sel), aes(x = DGE, y = IGE), color = "red") +  # Special points
      xlab("DGE") +
      ylab("IGE") +
      ggtitle("Genotype Distribution")+
      theme_bw()
    
    
     # Prepare G and E matrices
    
    G <- matrix(c(input$varG11, covdgE_IGE, covdgE_IGE, input$varG22), nrow = 2, byrow = TRUE)
    E <- matrix(c(input$varE11, 0, 0, input$varE22), nrow = 2, byrow = TRUE)
    
    P=G+E
    
    N <- input$N
    sde <- input$sde
    
    P <- G + E
    mu <- c(0, 0)
    
    # Generate genotypes and neighborhoods
    genotype <- mvrnorm(n = N^2, mu = mu, Sigma = G)
    enviro <- mvrnorm(n = N^2, mu = mu, Sigma = E)
    voisinage <- matrix(sample(1:N^2, 4*N^2, replace = TRUE), nrow = 2*N, ncol = 2*N)
    
    # Load the necessary package for qnorm, dnorm and pnorm functions
    if(!requireNamespace("stats", quietly = TRUE)) install.packages("stats")
    
    # Function to calculate selection intensity (i) from the selection proportion (p)
    calculate_selection_intensity <- function(p) {
      # Z-score corresponding to the non-selected proportion p
      z = qnorm(1 - p)
      
      # Calculate selection intensity i using the correct formula
      i = dnorm(z) / (1 - pnorm(z))
      
      return(i)
    }
    
    # Example of use to calculate i when the top 10% are selected
    
    i = calculate_selection_intensity(input$p)
    
    # Display selection intensity
    i
    
    # Valeur_G function
    Valeur_G <- function(x, y, genotype, voisinage) {
      geno <- c()
      for (i in (x-1):(x+1)) {
        for (j in (y-1):(y+1)) {
          if (!(i == x && j == y)) {
            geno <- c(geno, voisinage[i, j])
          }
        }
      }
      # addition of genetic and environmental effects
      DGE <- genotype[voisinage[x, y], 1] + enviro[voisinage[x, y], 1]
      IGE <- sum(genotype[geno, 2] + enviro[geno, 2])
      
      return(DGE + IGE)
    }
    
    # Calculations and display of results (to be adapted as needed)
    Pheno <- c()
    tri_geno <- c()
    cpt <- 0
    
    for (x in c(5:(N+5))) {
      for (y in c(5:(N+5))) {
        cpt <- cpt + 1
        Pheno[cpt] <- Valeur_G(x, y, genotype, voisinage)
        tri_geno[cpt] <- voisinage[x, y]
      }
    }
    
    # Before selection
    
    # mean(Pheno)
    # var(Pheno)
    
    # hist(Pheno)
    
    # Selection of 50%: phenotypic selection differential
    # mean(Pheno[order(Pheno, decreasing = TRUE)][1:50])
    
    # Selected individuals
    N_sel <- round(cpt * input$p)
    
    tri_geno[order(Pheno, decreasing = TRUE)][1:N_sel]
    geno_sel <- genotype[tri_geno[order(Pheno, decreasing = TRUE)][1:N_sel],]
    
    # Selection differential on DGE and IGE
    S <- c()
    S[1] <- mean(genotype[tri_geno[order(Pheno, decreasing = TRUE)][1:N_sel], 1])
    S[2] <- mean(genotype[tri_geno[order(Pheno, decreasing = TRUE)][1:N_sel], 2])
    
    # Response to selection on DGE and IGE
    R <- G %*% solve(P) %*% S
    # R
    
    # DGE and IGE of the new population
    mus <- mu
    mus[1] <- mu[1] + R[1]
    mus[2] <- mu[2] + R[2]
    
    # New population after selection
    # Calculate the new performance in mixture
    
    # Generate a cloud of these data with phenotypic data
    new_geno <- mvrnorm(n = N^2, mus, G)
    
    # Distribute them on an excess square grid
    new_voisinage <- matrix(sample(1:N^2, 4*N^2, replace = TRUE), nrow = 2*N, ncol = 2*N)
    
    new_Pheno <- c()
    new_tri_geno <- c()
    cpt <- 0
    
    # Extract an internal square in the large matrix
    for (x in c(5:(N+5))) {
      for (y in c(5:(N+5))) {
        cpt <- cpt + 1
        new_Pheno[cpt] <- Valeur_G(x, y, new_geno, new_voisinage)
        new_tri_geno[cpt] <- new_voisinage[x, y]
      }
    }
    
    # Example: Calculation and display of the average G values and phenotypic values
    # These are examples and should be adapted based on the specific calculations of your script
    
    output$summaryOutput <- renderText({
      paste("Summary of calculations...")
    })
    
    
    output$plotGenotype <- renderPlot({
      ggplot() +
        geom_point(data = as.data.frame(genotype), aes(x = V1, y = V2), color = "gray") +  # General points
        geom_point(data = as.data.frame(geno_sel), aes(x = V1, y = V2), color = "red") +  # Special points
        xlab("DGE") +
        ylab("IGE") +
        ggtitle("Genotype Distribution")
    })
    
    # Combine the two vectors into a dataframe
    combinedData <- rbind(data.frame(Value = Pheno, Phase = "Before selection"),
                          data.frame(Value = new_Pheno, Phase = "After selection"))
    
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
