
# Load necessary libraries
library(shiny)
library(shinyjs)  # For dynamic notifications
library(isoviz)
library(dplyr)
library(shiny)
library(shinyjs)

# Increase max request size to 100MB
options(shiny.maxRequestSize = 100 * 1024^2)
library(shiny)
library(shinyjs)
library(shiny)
library(shinyjs)
library(shiny)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  
  titlePanel("Intron-Exon Visualization and Guide Table using isoviz"),
  includeCSS("PacBio-themed CSS.css"),  
  sidebarLayout(
    sidebarPanel(
      fileInput("psl_file", "Upload Genome .psl File", accept = c(".psl")),
      fileInput("gene_trans_file", "Upload Gene-Transcript Conversion .txt File", accept = c(".txt")),
      fileInput("junction_file", "Upload Junction File (.junc.txt)", accept = c(".txt")),
      fileInput("intron_annot_file", "Upload Intron Annotations .rda File", accept = c(".rda")),  # Corrected .rda file input
      
      textInput("gene_name", "Enter Gene Name:", value = "RBFOX2"),
      textInput("ensembl_id", "Enter Gene Ensembl ID:", value = "ENSG00000100320"),
      
      fileInput("junction_list_file", "Upload Junctions List (.txt)", accept = c(".txt")),  # Corrected file input for junction list
      
      textInput("cell_type", "Enter Cell Type:", value = "Custom"),  # Text input for cell type
      
      numericInput("junc_usage", "Minimum Junction Usage:", value = 5, min = 1),
      
      actionButton("plot", "Generate Plot"),
      actionButton("generate_table", "Generate Guide Table"),
      
      downloadButton("download_plot", "Download Plot"),
      downloadButton("download_table", "Download Guide Table")
    ),
    
    mainPanel(
      # Add tabset for plot and guide table
      tabsetPanel(
        tabPanel("Plot",
                 plotOutput("iso_plot"),
                 h4("Your exon-intron visualization will appear here.")
        ),
        tabPanel("Guide Table",
                 tableOutput("guide_table"),
                 h4("Your guide table will appear here.")
        )
      )
    )
  )
)

library(shiny)
library(shinyjs)


library(shiny)
library(isoviz)
library(dplyr)

library(shiny)
library(isoviz)
library(dplyr)

server <- function(input, output, session) {
  
  # Helper function to show notifications and return the ID
  showCustomNotification <- function(message, type = "message") {
    id <- showNotification(message, type = type, duration = NULL)
    return(id)
  }
  
  # Observe for plot generation
  observeEvent(input$plot, {
    # Display start notification
    plotNotificationId <- showCustomNotification("Plot generation started...", type = "message")
    
    # Ensure files are uploaded
    req(input$psl_file)
    req(input$gene_trans_file)
    req(input$junction_file)
    req(input$intron_annot_file)
    
    tryCatch({
      # Read the uploaded files
      file_path <- input$psl_file$datapath
      gene_trans <- input$gene_trans_file$datapath
      junctions_file <- input$junction_file$datapath
      intron_annotations_file <- input$intron_annot_file$datapath
      
      # Load the .rda file to get the gencode_intron_all_data object
      load(intron_annotations_file)
      
      # Check if gencode_intron_all_data is loaded
      if (!exists("gencode_intron_all_data")) {
        stop("gencode_intron_all_data object not found in .rda file!")
      }
      
      # Process uploaded PSL and gene transcript files
      all_coordinates <- isoviz_coords(file_path, gene_trans)
      exon_coords <- all_coordinates[[1]]
      intron_coords <- all_coordinates[[2]]
      
      # Run minicutter to get clusters from the uploaded junction file
      intron_clusts <- isoviz_minicutter(junctions_file)
      
      # Filter data based on the user-provided gene input
      gene_exons <- filter(exon_coords, gene_name == input$gene_name)
      gene_introns <- filter(intron_coords, gene_name == input$gene_name)
      
      # Check if a junction file was uploaded
      junctions_to_include <- NULL
      if (!is.null(input$junction_list_file)) {
        # Read user-uploaded junctions file
        junctions_to_include <- readLines(input$junction_list_file$datapath)
      } else {
        # Default junctions if no file is uploaded
        junctions_to_include <- c("junc178147", "junc178149", "junc178135", "junc178136", "junc178145", "junc178146")
      }
      
      # Map junctions to the gene
      mapped_junctions <- isoviz_map_junctions(
        cell_type = input$cell_type,  # Now user-provided cell type
        gene_introns, intron_clusts, gencode_intron_all_data
      )
      
      # Render the plot
      output$iso_plot <- renderPlot({
        isoviz_plot_juncs_to_iso(
          mapped_junctions, gene_exons, gene_introns,
          cell_type = input$cell_type,
          junc_usage = input$junc_usage, 
          intron_scale = "no"
        )
      })
      
      # Download plot as a PNG
      output$download_plot <- downloadHandler(
        filename = function() { paste("intron_exon_plot", Sys.Date(), ".png", sep = "") },
        content = function(file) {
          png(file)
          isoviz_plot_juncs_to_iso(
            mapped_junctions, gene_exons, gene_introns,
            cell_type = input$cell_type,
            junc_usage = input$junc_usage, 
            intron_scale = "no"
          )
          dev.off()
        }
      )
      
      # Remove plot generation notification and show completion
      removeNotification(plotNotificationId)
      showCustomNotification("Plot generation complete!", type = "message")
    }, error = function(e) {
      removeNotification(plotNotificationId)
      showCustomNotification(paste("Error: ", e$message), type = "error")
    })
  })
  
  # Observe for guide table generation
  observeEvent(input$generate_table, {
    # Display start notification
    tableNotificationId <- showCustomNotification("Guide table generation started...", type = "message")
    
    # Ensure files are uploaded
    req(input$junction_file)
    
    tryCatch({
      junctions_file <- input$junction_file$datapath
      intron_clusts <- isoviz_minicutter(junctions_file)
      
      # Check if a junction file was uploaded
      junctions_to_include <- NULL
      if (!is.null(input$junction_list_file)) {
        # Read user-uploaded junctions file
        junctions_to_include <- readLines(input$junction_list_file$datapath)
      } else {
        # Default junctions if no file is uploaded
        junctions_to_include <- c("junc178147", "junc178149", "junc178135", "junc178136", "junc178145", "junc178146")
      }
      
      # Generate the guide table using the Ensembl ID input
      guide_table <- isoviz_get_guide_predictions(
        gene = input$ensembl_id, 
        leafcutter_input = intron_clusts,
        guides_per_junction = 5,
        include_specific_junctions = junctions_to_include,
        output_format = "dataframe"
      )
      
      # Render the guide table in the UI
      output$guide_table <- renderTable({
        guide_table
      })
      
      # Allow downloading the guide table as a CSV
      output$download_table <- downloadHandler(
        filename = function() { paste("guide_table", Sys.Date(), ".csv", sep = "") },
        content = function(file) {
          write.csv(guide_table, file, row.names = FALSE)
        }
      )
      
      # Remove guide table generation notification and show completion
      removeNotification(tableNotificationId)
      showCustomNotification("Guide table generation complete!", type = "message")
    }, error = function(e) {
      removeNotification(tableNotificationId)
      showCustomNotification(paste("Error: ", e$message), type = "error")
    })
  })
}


shinyApp(ui=ui,server = server)

