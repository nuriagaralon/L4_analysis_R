#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# The application is, to start, to deal with the survival data of SydLab one.

# Load necessary libraries
library(shiny)
library(tidyverse)
library(survival)
library(survminer)

# Define the UI of the app
ui <- fluidPage(
  
  # Title
  titlePanel("Survival Analysis and Plotting"),
  
  # Sidebar Layout
  sidebarLayout(
    sidebarPanel(
      # File input to upload the dataset
      fileInput("datafile", "Choose CSV File", accept = ".csv"),
      hr(),
      # Choose condition
      checkboxGroupInput("checkGroup", label = "Conditions", 
                         choices = NULL),
      hr(),
      # Button to generate plot
      actionButton("generate", "Generate Survival Plot"),
      hr(),
      # Show downloaded plot image
      downloadButton("download_plot", "Download Plot")
    ),
    
    mainPanel(
      # Output the survival plot
      plotOutput("survPlot")
    )
  )
)

# Define the server logic
server <- function(input, output) {
  
  # Reactive expression to load the data file when uploaded
  data_surv <- reactive({
    req(input$datafile)  # Ensure that file is uploaded
    # Read the CSV file
    path_surv <- input$datafile$datapath
    
    # Load headers and data
    headers <- read.csv(path_surv, nrows = 2, header = FALSE)
    headers <- sapply(headers, paste, collapse = "_")
    table_surv <- read.csv(file = path_surv, skip = 2, header = FALSE)
    names(table_surv) <- substr(headers, 11, 100)
    
    # Transform the data
    table_surv |>
      pivot_longer(
        cols = everything(),
        names_to = c("condition", "type"),
        names_pattern = "(.*)_(deaths|event)",
        values_to = "value"
      ) |>
      pivot_wider(
        names_from = type,
        values_from = value
      ) |>
      unnest(c(deaths, event)) |>
      drop_na()
  })
  
  # Generate the survival plot when the button is pressed
  observeEvent(input$generate, {
    req(data_surv())  # Ensure the data is loaded
    
    # Fit survival model
    km_fit <- survfit(Surv(deaths, event) ~ condition, data = data_surv())
    
    # Create the survival plot
    surv_plot <- ggsurvplot(
      km_fit,
      data = data_surv(),
      pval = TRUE,
      legend.title = "",
      legend.labs = levels(factor(data_surv()$condition))
    )
    
    # Store the plot in a reactive variable
    output$survPlot <- renderPlot({
      surv_plot
    })
    
    # Save the plot as an image
    output$download_plot <- downloadHandler(
      filename = function() { "survival_plot.png" },
      content = function(file) {
        ggsave(file, plot = surv_plot, width = 17, height = 15, dpi = 1000, units = "cm")
      }
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
