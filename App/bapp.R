# Funktion zur Berechnung des HDIs
hdi <- function(sample_vector, credMass = 0.95) {
  sorted_vector <- sort(sample_vector)
  ciIdxInc <- floor(credMass * length(sorted_vector))
  nCIs <- length(sorted_vector) - ciIdxInc
  ciWidth <- rep(0, nCIs)
  for (i in 1:nCIs) {
    ciWidth[i] <- sorted_vector[i + ciIdxInc] - sorted_vector[i]
  }
  hdi_min <- sorted_vector[which.min(ciWidth)]
  hdi_max <- sorted_vector[which.min(ciWidth) + ciIdxInc]
  return(c(hdi_min, hdi_max))
}

library(shiny)
library(ggplot2)
library(bayestestR)

ui <- fluidPage(
  titlePanel("Beta Distribution"),
  
  sidebarLayout(
    sidebarPanel(
      sliderInput("a",
                  "Parameter a:",
                  min = 1,
                  max = 100,
                  value = 5),
      
      sliderInput("b",
                  "Parameter b:",
                  min = 1,
                  max = 100,
                  value = 5),
      
      numericInput("N",
                   "Number of observations (N):",
                   value = 50, 
                   min = 2),
      
      numericInput("z",
                   "Number of successes (z):",
                   value = 25, 
                   min = 0)
    ),
    
    mainPanel(
      plotOutput("distPlot"),
      htmlOutput("formula"),
      htmlOutput("ratio"),
      htmlOutput("hdi"),
      htmlOutput("data")
    )
  )
)

server <- function(input, output) {
  
  output$distPlot <- renderPlot({
    x <- seq(0, 1, length.out = 1000)
    
    y_prior <- dbeta(x, input$a, input$b)
    y_likelihood <- x^input$z * (1-x)^(input$N - input$z)  # Bernoulli-Verteilung
    scaling_factor <- max(y_prior) / max(y_likelihood)
    y_likelihood <- y_likelihood * scaling_factor
    y_posterior <- dbeta(x, input$z + input$a, input$N - input$z + input$b)
    
    df <- data.frame(x = x, y_prior = y_prior, y_likelihood = y_likelihood, y_posterior = y_posterior)
    
    p <- ggplot(df, aes(x = x)) +
      geom_line(aes(y = y_prior), color = "blue") +
      geom_line(aes(y = y_likelihood), color = "green") +
      geom_line(aes(y = y_posterior), color = "red") +
      labs(x = expression(theta), y = "Density", 
           color = "Verteilung") +
      theme_minimal()
    
    p <- p + annotate("text", x = 0.85, y = max(y_prior, na.rm = TRUE), label = "Prior (Beta)", hjust = 1, color = "blue") +
      annotate("text", x = 0.85, y = max(y_likelihood, na.rm = TRUE) - 0.3, label = "Likelihood (Bernoulli)", hjust = 1, color = "green") +
      annotate("text", x = 0.85, y = max(y_posterior, na.rm = TRUE), label = "Posterior (Beta)", hjust = 1, color = "red")
    
    # HDI und Modus hinzufügen
    hdi_prior <- hdi(rbeta(10000, input$a, input$b))
    hdi_posterior <- hdi(rbeta(10000, input$z + input$a, input$N - input$z + input$b))
    mode <- x[which.max(y_posterior)]
    p <- p + annotate("segment", x = hdi_prior[1], xend = hdi_prior[2], y = -0.1, yend = -0.1, colour = "blue", size = 1) +
      annotate("text", x = mean(hdi_prior), y = -0.25, label = "95% HDI (Prior)", colour = "blue") +
      annotate("segment", x = hdi_posterior[1], xend = hdi_posterior[2], y = -0.5, yend = -0.5, colour = "red", size = 1) +
      annotate("text", x = mean(hdi_posterior), y = -0.65, label = "95% HDI (Posterior)", colour = "red") +
      geom_vline(aes(xintercept = mode), linetype = "dashed", color = "black") +
      annotate("text", x = mode, y = max(y_posterior) + 0.2, label = "Modus", hjust = ifelse(mode < 0.5, 1, 0))
    return(p)
  })
  
  
  
  output$formula <- renderUI({
    withMathJax(paste0("$$ P(\\theta|z,N) = Beta(\\theta|z + ", input$a, ",N - z + ", input$b, ") $$"))
  })
  
  output$ratio <- renderUI({
    ratio <- (input$a + input$b) / (input$a + input$b + input$N)
    withMathJax(paste0("$$ \\text{Ratio (a + b / a + b + N)}: ", round(ratio, 2), " $$"))
  })
  
  output$hdi <- renderUI({
    hdi <- hdi(rbeta(10000, input$z + input$a, input$N - input$z + input$b))
    withMathJax(paste0("$$ 95\\% \\text{HDI}: [", round(hdi[1], 2), ", ", round(hdi[2], 2), "] $$"))
  })
  
  output$data <- renderUI({
    withMathJax(paste0("$$ \\text{Data : z} = ", input$z, ", N = ", input$N, " $$"))
  })
}

shinyApp(ui = ui, server = server)

