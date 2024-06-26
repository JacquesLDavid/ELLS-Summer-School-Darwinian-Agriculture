library(shiny)
library(MASS) # For mvrnorm
library(ggplot2) # For plots

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
    
    
    # G <- matrix(c(0.7, -0.15, -0.15, 0.3), nrow = 2, byrow = TRUE)
    # E <- matrix(c(0.7, 0,0, 0.07), nrow = 2, byrow = TRUE)
    # N<-100
    # rep <- 2
    # r <- -0.5
    # covdgE_IGE <-r* (G[1,1] * G[2,2] )^0.5
    # G <- matrix(c(G[1,1], covdgE_IGE, covdgE_IGE, G[2,2]), nrow = 2, byrow = TRUE)

    
    
    # Prepare G and E matrices
    covdgE_IGE <- input$r * (input$varG11 * input$varG22)^0.5
    
    G <- matrix(c(input$varG11, covdgE_IGE, covdgE_IGE, input$varG22), nrow = 2, byrow = TRUE)
    E <- matrix(c(input$varE11, 0, 0, input$varE22), nrow = 2, byrow = TRUE)
    
    N <- input$N

    P <- G + E
    mu <- c(0, 0)
    
    # size of the square 
    size <- rep*N + 5
    
    Enviro <- mvrnorm(n = size^2 , mu = mu, Sigma = E)
    enviro_DGE <- matrix(Enviro[,1], nrow = size, ncol = size)
    enviro_IGE <- matrix(Enviro[,2], nrow = size, ncol = size)
    
    # the realized number of rep per genotype is not compulsarily rep 
    alea<-sample(1:N, size^2, replace = TRUE)
    voisinage <- matrix(alea, nrow = size, ncol = size)
    
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
    
    for (x in c(3:(N*rep+3))) {
      for (y in c(3:(N*rep+3))) {
        cpt <- cpt + 1
        Pheno[cpt] <- Valeur_G(x, y, genotype, voisinage)
        tri_geno[cpt] <- voisinage[x, y]
      }
    }
    
    length(Pheno)
    
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
