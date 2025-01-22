# Load necessary libraries
library(shiny)
library(shinyjs)  # For dynamic notifications
library(isoviz)
library(dplyr)
library(markdown)

# Increase max request size to 100MB
options(shiny.maxRequestSize = 100 * 1024^2)

# User Interface
ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  
  titlePanel("Intron-Exon Visualization and Guide Table using isoviz"),
  includeCSS("PacBio-themed CSS.css"),  
  sidebarLayout(
    sidebarPanel(
      fileInput("psl_file", "Upload Genome .psl File", accept = c(".psl")),
      fileInput("gene_trans_file", "Upload Gene-Transcript Conversion .txt File", accept = c(".txt")),
      fileInput("junction_file", "Upload Junction File (.junc.txt)", accept = c(".txt")),
      fileInput("intron_annot_file", "Upload Intron Annotations .rda File", accept = c(".rda")),
      
      textInput("gene_name", "Enter Gene Name:", value = "RBFOX2"),
      textInput("ensembl_id", "Enter Gene Ensembl ID:", value = "ENSG00000100320"),
      
      fileInput("junction_list_file", "Upload Junctions List (.txt)", accept = c(".txt")),
      textInput("cell_type", "Enter Cell Type:", value = "Custom"),
      numericInput("junc_usage", "Minimum Junction Usage:", value = 5, min = 1),
      
      actionButton("plot", "Generate Plot"),
      actionButton("generate_table", "Generate Guide Table"),
      
      downloadButton("download_plot", "Download Plot"),
      downloadButton("download_table", "Download Guide Table")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot",
                 plotOutput("iso_plot"),
                 h4("Your exon-intron visualization will appear here.")
        ),
        tabPanel("Guide Table",
                 tableOutput("guide_table"),
                 h4("Your guide table will appear here.")
        ),
        tabPanel("Help",  # New Help Tab
                 includeMarkdown("README.md")  # Display README.md content
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Helper function to show notifications and return the ID
  showCustomNotification <- function(message, type = "message") {
    id <- showNotification(message, type = type, duration = NULL)
    return(id)
  }
  
  # Observe for plot generation
  observeEvent(input$plot, {
    plotNotificationId <- showCustomNotification("Plot generation started...", type = "message")
    withProgress(message = "Generating plot...", value = 0, {
      req(input$psl_file, input$gene_trans_file, input$junction_file, input$intron_annot_file)
      
      tryCatch({
        incProgress(0.1, detail = "Loading input files...")
        
        file_path <- input$psl_file$datapath
        gene_trans <- input$gene_trans_file$datapath
        junctions_file <- input$junction_file$datapath
        intron_annotations_file <- input$intron_annot_file$datapath
        
        load(intron_annotations_file)
        
        if (!exists("gencode_intron_all_data")) {
          stop("gencode_intron_all_data object not found in .rda file!")
        }
        
        incProgress(0.2, detail = "Processing coordinates...")
        
        all_coordinates <- isoviz_coords(file_path, gene_trans)
        exon_coords <- all_coordinates[[1]]
        intron_coords <- all_coordinates[[2]]
        
        intron_clusts <- isoviz_minicutter(junctions_file)
        incProgress(0.4, detail = "Filtering data by gene name...")
        
        gene_exons <- filter(exon_coords, gene_name == input$gene_name)
        gene_introns <- filter(intron_coords, gene_name == input$gene_name)
        
        junctions_to_include <- NULL
        if (!is.null(input$junction_list_file)) {
          junctions_to_include <- readLines(input$junction_list_file$datapath)
        } else {
          junctions_to_include <- c("junc178147", "junc178149", "junc178135", "junc178136", "junc178145", "junc178146")
        }
        
        incProgress(0.6, detail = "Mapping junctions...")
        
        mapped_junctions <- isoviz_map_junctions(
          cell_type = input$cell_type,
          gene_introns, intron_clusts, gencode_intron_all_data
        )
        
        output$iso_plot <- renderPlot({
          isoviz_plot_juncs_to_iso(
            mapped_junctions, gene_exons, gene_introns,
            cell_type = input$cell_type,
            junc_usage = input$junc_usage, 
            intron_scale = "no"
          )
        })
        
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
        
        incProgress(1, detail = "Plot generation complete!")
        removeNotification(plotNotificationId)
        showCustomNotification("Plot generation complete!", type = "message")
        
      }, error = function(e) {
        removeNotification(plotNotificationId)
        showCustomNotification(paste("Error: ", e$message), type = "error")
      })
    })
  })
  
  # Observe for guide table generation
  observeEvent(input$generate_table, {
    tableNotificationId <- showCustomNotification("Guide table generation started...", type = "message")
    req(input$junction_file)
    
    tryCatch({
      junctions_file <- input$junction_file$datapath
      intron_clusts <- isoviz_minicutter(junctions_file)
      
      junctions_to_include <- NULL
      if (!is.null(input$junction_list_file)) {
        junctions_to_include <- readLines(input$junction_list_file$datapath)
      } else {
        junctions_to_include <- c("junc178147", "junc178149", "junc178135", "junc178136", "junc178145", "junc178146")
      }
      
      guide_table <- isoviz_get_guide_predictions(
        gene = input$ensembl_id, 
        leafcutter_input = intron_clusts,
        guides_per_junction = 5,
        include_specific_junctions = junctions_to_include,
        output_format = "dataframe"
      )
      
      output$guide_table <- renderTable({
        guide_table
      })
      
      output$download_table <- downloadHandler(
        filename = function() { paste("guide_table", Sys.Date(), ".csv", sep = "") },
        content = function(file) {
          write.csv(guide_table, file, row.names = FALSE)
        }
      )
      
      removeNotification(tableNotificationId)
      showCustomNotification("Guide table generation complete!", type = "message")
    }, error = function(e) {
      removeNotification(tableNotificationId)
      showCustomNotification(paste("Error: ", e$message), type = "error")
    })
  })
}

# Launch the App
shinyApp(ui = ui, server = server)
