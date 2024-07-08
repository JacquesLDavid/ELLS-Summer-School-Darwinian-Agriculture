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
      checkboxInput("Mean","Mean effects (checked) or Sum effects (unchecked) for IGE calculation", value=FALSE),
      checkboxInput("Asreml","Use of Asreml to make inference", value=FALSE),
      sliderInput("N", "Number of genotypes (N)", min = 10, max = 500, value = 450),
      sliderInput("rep", "Number of rep per genotype  (rep)", min = 1, max = 100, value = 2),
      sliderInput("N_sim", "Number of simulations  (N_sim)", min = 1, max = 100, value = 1),
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
                           plotOutput("plotTRUE_DGE_PRED_DGE"),
                           plotOutput("plotTRUE_IGE_PRED_IGE"),
                           plotOutput("plotPRED_DGE_IGE"),
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
                               textOutput("TRUEMeanOutputIGE"),
                               textOutput("TRUEVarOutputEnv_DGE"),
                               textOutput("TRUEMeanOutputEnv_DGE"),
                               textOutput("TRUEVarOutputEnv_IGE"),
                               textOutput("TRUEMeanOutputEnv_IGE"),
                               textOutput("MeanPheno"),
                               textOutput("VarPheno"),
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
                               textOutput("MeanPhenoSel"),
                               textOutput("VarPhenoSel"),
                               textOutput("SelectionDifferential_DGE"),
                               textOutput("SelectionDifferential_IGE"),
                               textOutput("SelectionResponse_DGE"),
                               textOutput("SelectionResponse_IGE"),
                               br(),
                               h3("Index selection"),
                               textOutput("MeanPhenoSel_I"),
                               textOutput("VarPhenoSel_I"),
                               textOutput("SelectionResponse_DGE_I"),
                               textOutput("SelectionResponse_IGE_I")
                               
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
    
    count_neighbors <- function(mat) {
      n_row <- nrow(mat)
      n_col <- ncol(mat)
      offsets <- expand.grid(di = c(-1, 0, 1), dj = c(-1, 0, 1))
      offsets <- offsets[!(offsets$di == 0 & offsets$dj == 0), ]
      
      unique_genotypes <- unique(mat)
      neighbor_counts <- array(0, dim = c(n_row, n_col, length(unique_genotypes)))
      dimnames(neighbor_counts)[[3]] <- unique_genotypes
      
      for (k in 1:nrow(offsets)) {
        ni <- (row(mat) + offsets$di[k] - 1) %% n_row + 1
        nj <- (col(mat) + offsets$dj[k] - 1) %% n_col + 1
        for (i in 1:n_row) {
          for (j in 1:n_col) {
            current_genotype <- mat[ni[i, j], nj[i, j]]
            if (!is.na(match(current_genotype, unique_genotypes))) {
              neighbor_counts[i, j, match(current_genotype, unique_genotypes)] <- neighbor_counts[i, j, match(current_genotype, unique_genotypes)] + 1
            }
          }
        }
      }
      
      neighbor_counts
    }
    
    assign("count_neighbors",count_neighbors,globalenv())
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
    N_sim=input$N_sim
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
    
    SIM <- numeric(N_sim)
    calc_SIM <- numeric(N_sim)
    
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
      
      DATA <- expand.grid(Focal = Focal, rep = rep)
      grid_indices <- sample(N_row * N_col, N_geno * N_rep, replace = FALSE)
      grid_sample <- grid[grid_indices, ]
      DATA <- cbind(grid_sample, DATA)
      
      mat_grid <- matrix("vide", nrow = N_row, ncol = N_col)
      mat_grid[cbind(grid_sample$Row, grid_sample$Column)] <- as.character(DATA$Focal)
      
      empty_positions <- which(mat_grid == "vide", arr.ind = TRUE)
      mat_grid[empty_positions] <- sample(DATA$Focal, nrow(empty_positions), replace = TRUE)
      
      neighbor_counts <- count_neighbors(mat_grid)
      
      DATA <- cbind(matrix(0, nrow = nrow(DATA), ncol = length(Focal)), DATA)
      colnames(DATA)[1:N_geno] <- Focal
      
      for (k in 1:nrow(DATA)) {
        i <- DATA$Row[k]
        j <- DATA$Column[k]
        for (geno in Focal) {
          if (geno %in% dimnames(neighbor_counts)[[3]]) {
            DATA[k, geno] <- neighbor_counts[i, j, geno]
          }
        }
      }
      
      Zg <- model.matrix(~Focal - 1, DATA)
      dimnames(Zg)[[2]] <- paste0("G", sprintf("%03d", 1:N_geno))
      
      if (Mean) {
        Zv <- as.matrix(DATA[, 1:N_geno]) / 8
        DATA[, 1:N_geno] <- Zv
      } else {
        Zv <- as.matrix(DATA[, 1:N_geno])
      }
      
      Pheno <- Zg %*% DGE + Zv %*% IGE + df_E[, 1] + df_E[, 2]
      
      SIM[i] <- var(Pheno)
      calc_SIM[i] <- round(var(DGE) + 8 * var(IGE) + 8 * mean(tcrossprod(IGE + DGE)) * (2 * cov(DGE, IGE) + 7 * var(IGE)) / (N_col * N_row), 3)
    }
    
    TABLE_TRUE <- data.frame(
      "Effect" = c("DGE", "IGE", "Pheno", "calc_SIM"),
      "Variance" = c(V_geno, V_voisin, mean(SIM), mean(calc_SIM))
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
    
    DATA_sel_I <- expand.grid(Focal = Focal_sel_I, rep = rep_sel_I)
    grid_indices <- sample(N_row * N_col, N_geno * N_rep, replace = FALSE)
    grid_sample <- grid[grid_indices, ]
    DATA_sel_I <- cbind(grid_sample, DATA_sel_I)
    
    mat_grid <- matrix("vide", nrow = N_row, ncol = N_col)
    mat_grid[cbind(grid_sample$Row, grid_sample$Column)] <- as.character(DATA_sel_I$Focal)
    
    empty_positions <- which(mat_grid == "vide", arr.ind = TRUE)
    mat_grid[empty_positions] <- sample(DATA_sel_I$Focal, nrow(empty_positions), replace = TRUE)
    
    neighbor_counts <- count_neighbors(mat_grid)
    
    DATA_sel_I <- cbind(matrix(0, nrow = nrow(DATA_sel_I), ncol = length(Focal_sel_I)), DATA_sel_I)
    colnames(DATA_sel_I)[1:N_geno] <- Focal_sel_I
    
    for (k in 1:nrow(DATA_sel_I)) {
      i <- DATA_sel_I$Row[k]
      j <- DATA_sel_I$Column[k]
      for (geno in Focal_sel_I) {
        if (geno %in% dimnames(neighbor_counts)[[3]]) {
          DATA_sel_I[k, geno] <- neighbor_counts[i, j, geno]
        }
      }
    }
    
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
    
    DATA_sel <- expand.grid(Focal = Focal_sel, rep = rep_sel)
    grid_indices <- sample(N_row * N_col, N_geno * N_rep, replace = FALSE)
    grid_sample <- grid[grid_indices, ]
    DATA_sel <- cbind(grid_sample, DATA_sel)
    
    mat_grid <- matrix("vide", nrow = N_row, ncol = N_col)
    mat_grid[cbind(grid_sample$Row, grid_sample$Column)] <- as.character(DATA_sel$Focal)
    
    empty_positions <- which(mat_grid == "vide", arr.ind = TRUE)
    mat_grid[empty_positions] <- sample(DATA_sel$Focal, nrow(empty_positions), replace = TRUE)
    
    neighbor_counts <- count_neighbors(mat_grid)
    
    DATA_sel <- cbind(matrix(0, nrow = nrow(DATA_sel), ncol = length(Focal_sel)), DATA_sel)
    colnames(DATA_sel)[1:N_geno] <- Focal_sel
    
    for (k in 1:nrow(DATA_sel)) {
      i <- DATA_sel$Row[k]
      j <- DATA_sel$Column[k]
      for (geno in Focal_sel) {
        if (geno %in% dimnames(neighbor_counts)[[3]]) {
          DATA_sel[k, geno] <- neighbor_counts[i, j, geno]
        }
      }
    }
    
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
        theme(legend.title = element_blank())+# Remove the legend title
        annotate("text",x=mean_after_sel+1.5,y=0.3,label=paste0('Delta_mu_pheno = ', mean_after_sel - mean_before_sel))# Remove the legend title
      
    })
    
  })
  
}

# Run the Shiny application
shinyApp(ui = ui, server = server)
