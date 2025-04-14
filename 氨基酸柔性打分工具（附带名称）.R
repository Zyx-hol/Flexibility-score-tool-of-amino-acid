library(ggplot2)
library(shiny)

# BNORM 系数表（与论文一致）/BNORM coefficient table (consistent with the paper)
bnorm0 <- c(1.093, 1.169, 1.142, 1.055, 1.033, 1.094, 1.165, 1.073, 1.117, 1.038,
            1.041, 0.967, 0.982, 0.982, 0.961, 1.002, 0.930, 0.960, 0.925, 0.947)
bnorm1 <- c(1.082, 1.048, 1.042, 1.085, 1.089, 1.036, 1.028, 1.051, 1.006, 1.028,
            0.946, 0.961, 0.952, 0.927, 0.930, 0.892, 0.912, 0.878, 0.917, 0.862)
bnorm2 <- c(1.057, 0.923, 0.923, 0.932, 0.932, 0.933, 0.885, 0.934, 0.930, 0.901,
            0.892, 0.921, 0.894, 0.913, 0.837, 0.872, 0.914, 0.925, 0.803, 0.804)

amino_acids <- c('K', 'S', 'G', 'P', 'D', 'E', 'Q', 'T', 'N', 'R',
                 'A', 'L', 'H', 'V', 'Y', 'I', 'F', 'C', 'W', 'M')

# 固定权重/Fixed weight
fixed_weights <- c(0.25, 0.50, 0.75, 1.00, 0.75, 0.50, 0.25)

# 计算柔性得分函数/Calculate the flexibility score function
calculate_flexibility <- function(seq) {
  seq <- toupper(seq)
  n <- nchar(seq)
  scores <- rep(NA, n)
  
  for (i in 1:n) {
    aa <- substr(seq, i, i)
    if (!(aa %in% amino_acids)) next
    
    # 统计相邻刚性氨基酸数目/Counting the number of adjacent rigid amino acids
    left <- if (i > 1) substr(seq, i-1, i-1) else NA
    right <- if (i < n) substr(seq, i+1, i+1) else NA
    rigid_neighbors <- sum(c(left, right) %in% c('A','L','H','V','Y','I','F','C','W','M'), na.rm = TRUE)
    bnorm_values <- switch(as.character(rigid_neighbors),
                           "0" = bnorm0, "1" = bnorm1, "2" = bnorm2)
    
    # 滑动窗口/sliding window
    start <- max(1, i - 3)
    end <- min(n, i + 3)
    window <- strsplit(substr(seq, start, end), "")[[1]]
    window_length <- length(window)
    
    # weight
    weight_start <- if (i - start >= 3) 1 else (4 - (i - start))
    weights <- fixed_weights[weight_start:(weight_start + window_length - 1)]
    
    # 计算加权和/calculate a weighted sum
    weighted_sum <- 0
    for (j in 1:window_length) {
      aa_window <- window[j]
      idx <- match(aa_window, amino_acids)
      if (!is.na(idx)) {
        weighted_sum <- weighted_sum + bnorm_values[idx] * weights[j]
      }
    }
    
    scores[i] <- weighted_sum / 4.0
  }
  
  return(scores)
}

# UI 界面/ UI screen
ui <- fluidPage(
  titlePanel("Flexibility score tool of amino acid"),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      textInput("sequence", "Enter the amino acid sequence (single-letter code):", value = ""),
      actionButton("calculate", "Calculate the flexibility score"),
      br(), br(),
      strong("Mean flexibility score (6 amino acids before and after exclusion):"),
      textOutput("avgFlexibility"),
      br(),
      helpText("note：Flexibility score > 1 flexibility，< 1 rigid ")
    ),
    mainPanel(
      plotOutput("flexibilityPlot", height = "400px"),
      br(),
      plotOutput("trimmedFlexibilityPlot", height = "400px")  
    )
  )
)

# Server 逻辑/ Server logic
server <- function(input, output) {
  observeEvent(input$calculate, {
    seq <- gsub("[^A-Za-z]", "", input$sequence)
    if (nchar(seq) < 1) {
      showModal(modalDialog(title = "error", "Please enter valid amino acid sequence!", easyClose = TRUE))
      return()
    }
    
    scores <- calculate_flexibility(seq)
    n <- length(scores)
    
    # 平均值计算（排除前后3个氨基酸）/Calculation of the average (3 amino acids before and after exclusion)
    if (n >= 7) {
      start_avg <- 4
      end_avg <- n - 3
      avg_scores <- scores[start_avg:end_avg]
      avg_value <- mean(avg_scores, na.rm = TRUE)
    } else {
      avg_value <- NA
    }
    
    # 全图输出/output figures
    output$flexibilityPlot <- renderPlot({
      df <- data.frame(Position = 1:n, Flexibility = scores)
      ggplot(df, aes(x = Position, y = Flexibility)) +
        geom_line(color = "steelblue", linewidth = 1) +
        geom_point(color = "darkred", size = 2) +
        geom_hline(yintercept = 1.0, color = "red", linetype = "dashed", linewidth = 1) +
        labs(
          title = "Amino acid Flexibility score curve for each position",
          x = "Residue No. ",
          y = "Flexibility score"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 16, vjust = -0.5),
          axis.title.y = element_text(size = 16, vjust = 2.5),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14)
        )
    })
    
    # 截去首尾的图 + 添加氨基酸标签/Cut the beginning and end of the figure + add amino acid tags
    output$trimmedFlexibilityPlot <- renderPlot({
      if (n >= 7) {
        trimmed_positions <- 4:(n - 3)
        trimmed_scores <- scores[trimmed_positions]
        trimmed_aas <- unlist(strsplit(seq, ""))[trimmed_positions]
        
        trimmed_df <- data.frame(
          Position = trimmed_positions,
          Flexibility = trimmed_scores,
          AminoAcid = trimmed_aas
        )
        
        ggplot(trimmed_df, aes(x = Position, y = Flexibility, label = AminoAcid)) +
          geom_line(color = "seagreen", linewidth = 1) +
          geom_point(color = "red", size = 2) +
          geom_text(vjust = -1, size = 5, color = "black", fontface = "bold") + 
          geom_hline(yintercept = 1.0, color = "red", linetype = "dashed", linewidth = 1) +
          labs(
            title = expression(italic("SPECIES(enter the name)")),
            x = "Residue No. ",
            y = "Flexibility score"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
            axis.title.x = element_text(size = 18, vjust = -0.5),
            axis.title.y = element_text(size = 18, vjust = 2.5 ),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            plot.margin = margin(10, 20, 30, 50)
          )
      }
    })
    
    output$avgFlexibility <- renderText({
      if (is.na(avg_value)) {
        "Unable to calculate effective average (insufficient sequence length or too many invalid amino acids) ."
      } else {
        round(avg_value, 5)
      }
    })
  })
}

# 启动应用/Start the tool
shinyApp(ui = ui, server = server)
