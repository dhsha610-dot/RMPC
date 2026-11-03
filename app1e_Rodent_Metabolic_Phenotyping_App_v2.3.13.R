# ====================================================================
# AUTOMATIC PACKAGE INSTALLATION AND LOADING
# ====================================================================
# This section automatically checks for required packages and installs
# any that are missing, making it easier for new users to run the app

cat("🔍 Checking for required packages...\n")

# List of required packages
required_packages <- c(
    "shiny",        # Web application framework
    "shinyFiles",   # File system access for Shiny
    "tidyverse",    # Data manipulation (includes dplyr, ggplot2, tidyr, etc.)
    "ggpubr",       # Publication-ready plots
    "RColorBrewer", # Color palettes
    "shinyjs",      # JavaScript operations in Shiny
    "readxl",       # Reading Excel files
    "DT"            # Interactive data tables
)

# Function to check and install packages
install_if_missing <- function(package) {
    if (!require(package, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("📦 Installing missing package: ", package, "\n"))
        install.packages(package, dependencies = TRUE, repos = "https://cran.rstudio.com/")
        library(package, character.only = TRUE)
        cat(paste0("✅ Successfully installed and loaded: ", package, "\n"))
    } else {
        cat(paste0("✅ Package already installed: ", package, "\n"))
    }
}

# Check and install all required packages
for (pkg in required_packages) {
    install_if_missing(pkg)
}

cat("\n🎉 All required packages are installed and loaded!\n")
cat("========================================\n\n")

# Load necessary libraries (this is now guaranteed to work)
library(shiny)
library(shinyFiles)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(shinyjs)
library(readxl)
library(DT)

# r2.2, 3/14/2026
library(pracma)
#r2.2.4, 3/26/2026, using GAM model
install_if_missing("mgcv")
library(mgcv)
#r2.2.5, 3/26/2026, add calibration of baseline
#r2.2.6, 3/27/2026, add Analysis selection and parapeters for Gam
# Pairwise difference smooths with gratia
library(gratia)
packageVersion("gratia")


# added Two-period ANCOVA comparison
#r2.2.7, 3/29/2026, add bodymass to Gam
#r2.2.8, 3/29/2026, add bodymass to two-period comparison - not working yet
library(emmeans)

#r2.2.9, 3/31/2026, added Gam pairwise comparison
#r2.2.10, 4/1/2026, separating hard cutoff from Hampel filter
#r2.2.11, 4/2/2026, two-period time cut both before and after 
library(lubridate)
#r2.2.12  added mean line and band for before/after two periods
#r2.2.13  added start time/injection time selection
#r2.2.14  added before/after change ANCOVA
#r2.2.15  added before/after GLM with interaction of Group*bodyMass for both/bothc
#r2.2.15.1  added before/after w/different length of windows
#r2.2.15.2  added before/after w/different length of windows, dynamically conditionally on method 
#r2.2.15.3  adding  before glm, after glm
#r2.2.16  added before glm, after glm (ANCOVA/Interaction)
#r2.2.18  added zero body mass
#r2.2.19  simplified the statistics analysis code
#r2.3.0  added two plot areas
#r2.3.1  fixed display for multiple int_p_val ( of two more condition group)
#r2.3.2  moved statistics out of redernplots
#r2.3.3  moved gam statistics out of redernplots
#r2.3.4  added ancova statistics results text render place but not for gam model(not showup)
#r2.3.5  added body mass effect
#r2.3.6  added Data Table
#r2.3.7  Update PDF dowlaod
#r2.3.8  
#r2.3.9  export post-analysis data fitted
#r2.3.10 exploring the model fitness AIC BIC R2
#r2.3.11 included baseline into fitted data (csv) download
# 6/18/2026, fixed sometime the calR file doesn't includethe column: Date_Time_1
# 2.3.12 7/9/2026, add parameters to the name of fitted data export csv file
# 2.3.13 9/17/2026, add condition to update the plots/analysis
# 2.3.13 10/8/2026, change data input from a folder to from multiple csv files

# install.packages(c("rmarkdown", "knitr", "tinytex"))
# tinytex::install_tinytex()
 



# r2.1, 3/14/2026

# Ensure metric_<cage> columns exist for a set of metric prefixes and cage numbers
ensure_metric_cage_cols <- function(df, prefixes, cage_numbers, fill = -1) {
  # sanitize inputs
  prefixes <- unique(trimws(as.character(prefixes)))
  prefixes <- prefixes[nzchar(prefixes)]
  cage_numbers <- sort(unique(suppressWarnings(as.integer(cage_numbers))))
  cage_numbers <- cage_numbers[!is.na(cage_numbers)]
  if (length(prefixes) == 0 || length(cage_numbers) == 0) return(df)
  
  # Build full set of expected names: e.g., si13c_1, vo2_1, ..., si13c_12, ...
  expected <- as.vector(outer(prefixes, cage_numbers, function(p, k) paste0(p, "_", k)))
  missing  <- setdiff(expected, names(df))
  
  if (length(missing) > 0) {
    for (col in missing) df[[col]] <- fill
  }
  df
}


# Set maximum file upload size to 50GB globally
options(shiny.maxRequestSize = 50 * 1024^3)  # 50GB in bytes

# Define the user interface (UI)
ui_app1e <- fluidPage(
    useShinyjs(),  # Enable shinyjs
    # titlePanel("🐭 Rodent Metabolic Phenotyping Analysis Platform "),
    
    
    titlePanel(
      div(
        # style = "display:flex; align-items:center;",
        style = "text-align:center;",
        span("🐭 Rodent Metabolic Phenotyping Analysis Platform "),
        tags$img(
          src = "pennmed.png",
          height = "45px",
          style = "margin-right:12px;"
        )
        
      )
    ),
    
    
    sidebarLayout(
      sidebarPanel(
        width = 4,
          
            h4("1. Select Data Type"),
            radioButtons("data_type", "Choose Analysis Type:",
                        choices = c("CLAMS Metabolic Cage Data" = "clams",
                                  "DSI Telemetry Data" = "dsi"),
                        selected = "clams"),
            
            # hr(),
            
            h4("2. Experiment Setup"),
            textInput("experiment_name", "Experiment Name", "LRAP_Cohort1"),
        
            # r2.2.13
            selectInput(
              inputId = "start_time_col",
              label   = "Select start time (Minute=0):",
              choices = c(
                "Start Time"        = "start_time",
                "Injection Time 1"  = "injection_time1",
                "Injection Time 2"  = "injection_time2"
              ),
              selected = "start_time"
            ),
        
            
            conditionalPanel(
                condition = "input.data_type == 'clams'",
                numericInput("duration_of_experiment", "Experiment Duration (minutes)", 150),
                numericInput("duration_before_experiment", "Baseline Duration (minutes)", 60)
            ),
        
            
            conditionalPanel(
                condition = "input.data_type == 'dsi'",
                numericInput("dsi_duration_of_experiment", "Experiment Duration (minutes)", 150),
                numericInput("dsi_duration_before_experiment", "Baseline Duration (minutes)", 60),
                numericInput("start_row", "Excel Start Row:", value = 3, min = 1, max = 20)
            ),
            
            h4("3. Data Input"),
            
            # CLAMS Data Input
            # conditionalPanel(
            #     condition = "input.data_type == 'clams'",
            #     shinyDirButton("data_dir", "Select CLAMS Data Folder", "Please select a folder"),
            #     verbatimTextOutput("dir_path_display")
            # ),
        # 10/8/2026
        conditionalPanel(
          condition = "input.data_type == 'clams'",
          fileInput(
            "clams_files",
            "Upload CLAMS CSV files",
            multiple = TRUE,
            accept = c(".csv")
          ),
          verbatimTextOutput("clams_files_display")
        ),
        
        

            # DSI Data Input  
            conditionalPanel(
                condition = "input.data_type == 'dsi'",
                fileInput("excel_file", "Upload DSI Excel File", accept = c(".xlsx", ".xls")),
                uiOutput("sheet_selector"),
                numericInput("start_row", "Data Start Row (exclude headers)", 3, min = 1),
                checkboxInput("use_metadata", "Use Metadata File", value = FALSE)
            ),
            
            # hr(),
            
            # Metadata section
            conditionalPanel(
                condition = "input.data_type == 'clams' || input.use_metadata == true",
                downloadButton("download_metadata_template", "Download Metadata Template"),
                fileInput("metadata_file", "Upload Metadata File (CSV)", multiple = FALSE, accept = ".csv")
            ),
            
            h4("4. Plot Options"),
            
            # CLAMS-specific plot options
            conditionalPanel(
                condition = "input.data_type == 'clams'",
                selectInput("clams_y_var", "Y-axis Variable:",
                           choices = c("¹³C in exhaled CO₂" = "si13c",
                                     "VO₂ (Oxygen Consumption)" = "vo2",
                                     "VCO₂ (CO₂ Production)" = "vco2",
                                     "Respiratory Quotient (RQ)" = "rq",
                                     "Energy Expenditure (kcal/hr)" = "kcal_hr",
                                     "Water Vapor (VH₂O)" = "vh2o",
                                     "Food Intake (up)" = "foodupa",
                                     "Water Intake A" = "waterupa",
                                     "Water Intake B" = "waterupb",
                                     "Body Mass" = "bodymass",
                                     "X-axis Activity" = "xbreak",
                                     "Y-axis Activity" = "ybreak",
                                     "Z-axis Activity" = "zbreak",
                                     "Pedometer Counts" = "pedmeters",
                                     "Total Activity" = "allmeters",
                                     "Body Temperature" = "bodytemp"),
                           selected = "si13c"),
                
                radioButtons("plot_color_by", "Color Plot By:",
                            choices = c("Samples", "Groups" = "Conditions"),
                            selected = "Conditions"),
                
                checkboxInput("use_faceting", "Enable Faceted Plots", value = FALSE),
                
                conditionalPanel(
                    condition = "input.use_faceting == true",
                    radioButtons("facet_by", "Facet (Split) Plots By:",
                               choices = c("Conditions" = "Conditions", 
                                         "Sample Dataset" = "Sample_Dataset",
                                         "Individual Samples" = "Samples"),
                               selected = "Conditions"),
                    numericInput("facet_cols", "Number of Columns:", value = 2, min = 1, max = 4)
                )
            ),
            
            # DSI-specific plot options  
            conditionalPanel(
                condition = "input.data_type == 'dsi'",
                selectInput("dsi_y_var", "Y-axis Variable:", choices = c(), selected = NULL),
                radioButtons("dsi_plot_group_by", "Group Lines By:",
                           choices = c("Mouse ID" = "mouse", "Condition" = "condition"),
                           selected = "mouse", inline = TRUE),
                checkboxInput("dsi_plot_facet", "Facet by Group", value = FALSE)
            ),
            
            #r2.2.5
            h4("5. Data Clean"),
            
            # r2.2.7
            # --- Body Mass selection (CLAMS only) ---

            
        conditionalPanel(
          condition = "input.data_type == 'clams'",
          selectInput(
            "bodymass_source",
            "5.0 Body Mass",

            choices = c("No Mass (0)" = "zeromass", 
                        "Use original CLAMS bodymass column" = "original",
                        "Lean"   = "lean",
                        "Fat"    = "fat",
                        "Total"  = "total"
            ),
            selected = "lean"
            
          ),
          # helpText(
          #   "When Lean/Fat/Total is selected, the 'bodymass' column in the combined CLAMS data ",
          #   "will be replaced per sample using the metadata columns ",
          #   "‘BodyMass_Lean’, ‘BodyMass_Fat’, and ‘BodyMass_Total’ (grams)."
          # )
        ),
       
        
        
            
            conditionalPanel(
              condition = "input.clams_y_var == 'si13c'",
              
              radioButtons(
                inputId = "baseline_mode",
                label   = "5.1 Baseline Calibration",
                choices = c(
                  "None" = "none",
                  "Baseline-adjusted" = "adjust"
                ),
                selected = "none"
              ),
              
              conditionalPanel(
                condition = "input.baseline_mode == 'adjust'",
                
                numericInput(
                  inputId = "baseline_st",
                  label   = "Baseline start (Minutes)",
                  value   = -60
                )
              )
            ),
        
            #add hard cutoff, r2.2.10
            
            conditionalPanel(
              condition = "input.clams_y_var == 'si13c'",
              
              selectInput(
                "boundary_cutoff",
                "5.2 Boundary cutoff",
                choices = c(
                  "No cutoff" = "none",
                  "Yes" = "cutoff"
                ),
                selected = "none"
              ),
              
              conditionalPanel(
                condition = "input.boundary_cutoff == 'cutoff'",
                
                # tags$hr(),
                strong("Bondary Settings"),
                
                #add hard cutoff
                numericInput("lower_cutoff", "5.2.1 Lower cutoff", value = -50,min=-250,step=-50),
                numericInput("upper_cutoff", "5.2.2 Upper cutoff", value = 1200,min=500,step=50)
              )
            ),
        
        
            
            #r2.2, DS, 3/22/2026
            
            #Add filter selector in the UI
            
            conditionalPanel(
              condition = "input.clams_y_var == 'si13c'",
              
              selectInput(
                "signal_filter",
                "5.3 Outlier removal",
                choices = c(
                  "No filter" = "none",
                  "Hampel filter (robust outlier removal)" = "hampel"
                ),
                selected = "none"
              ),
            
            conditionalPanel(
              condition = "input.signal_filter == 'hampel'",
              
              # tags$hr(),
              strong("Hampel Filter Settings"),
              
              numericInput("t_split", "5.3.1 Filter split time (minutes)", -30),
              numericInput("t_start_pre", "5.3.2 Baseline filter start (minutes)", value = -60, step = 1),
              numericInput("hampel_k", "5.3.3 Hampel window size (k)", value = 5, min = 1, step = 1),
              numericInput("hampel_t0", "5.3.4 Hampel threshold (t0, number of MADs)", value = 3, min = 0.1, step = 0.1)
            )
            ),
            
            # r2.2.6
            h4("6. Statistical Analysis"),
            
            radioButtons(
              inputId = "analysis_model",
              label   = "Select Analysis Method",
              choices = c(
                "None" = "none",
                "GAM (smooth time curves)" = "gam",
                "Two-period Mean comparison (ANCOVA)" = "twoperiod",
                "Two-period Mean comparison (GLM)" = "twoperiodglm"
              ),
              selected = "twoperiodglm"
            ),
            
            # --- GAM options ---
            conditionalPanel(
              condition = "input.analysis_model == 'gam'",
              
              numericInput(
                inputId = "gam_k",
                label   = "Spline flexibility (k)",
                value   = 10,
                min     = 3,
                max     = 50,
                step    = 1
              )
            ),
            
            # --- Two-period analysis options ---
            conditionalPanel(
              condition = "['twoperiod', 'twoperiodglm'].includes(input.analysis_model)",
              
              # numericInput(
              #   inputId = "cut_time",
              #   label   = "Cut point (Minutes)",
              #   value   = 100
              # ),
              
              sliderInput("cut_time", "Cut point (Minutes)", 0, 150, 150),
        
              
              selectInput(
                "two_periods",
                "Before or After",
                choices = c(
                  "Both before & after (change)" = "bothc",
                  "Both before & after" = "both",
                  "Before only" = "before",
                  "After only" = "after"
                ),
                selected = "before"
              ),
              
              
              
              
              # numericInput(
              #   inputId = "time_windowb4",
              #   label   = "Window length (Minutes/before)",
              #   value   = 30,
              #   min     = 1
              # ),
              # 
              # numericInput(
              #   inputId = "time_window",
              #   label   = "Window length (Minutes/after)",
              #   value   = 20,
              #   min     = 1
              # ),
              
              conditionalPanel(
                condition = "['both', 'bothc'].includes(input.two_periods)",

                numericInput(
                  inputId = "time_windowb4",
                  label   = "Window length (Minutes/before)",
                  value   = 30,
                  min     = 1
                ),

                numericInput(
                  inputId = "time_window",
                  label   = "Window length (Minutes/after)",
                  value   = 20,
                  min     = 1
                )
              ),

              conditionalPanel(
                condition = "input.two_periods == 'before'",

                numericInput(
                  inputId = "time_windowb4",
                  label   = "Window length (Minutes/before)",
                  value   = 30,
                  min     = 1
                )
              ),

              conditionalPanel(
                condition = "input.two_periods == 'after'",

                numericInput(
                  inputId = "time_window",
                  label   = "Window length (Minutes/after)",
                  value   = 20,
                  min     = 1
                )
              ),
              
              
            ),
        

            # h4("7. Run Analysis"),
            # hr(),
            # actionButton("run_analysis", "Run Analysis", icon = icon("cogs"), class = "btn-primary"),
            # br(), br(),
            # actionButton("clear_all", "🔄 Clear All & Start Fresh", icon = icon("refresh"), 
            #             class = "btn-warning btn-sm", 
            #             style = "width: 100%;"),
        
        
            actionButton("run_analysis", "Run Analysis", icon = icon("cogs"), class = "btn-primary"),
            br(), br(),
            actionButton("update_view", "Update Plot / Stats", icon = icon("refresh"), class = "btn-info"),
            br(), br(),
            actionButton("clear_all", "🔄 Clear All & Start Fresh", icon = icon("refresh"), 
                         class = "btn-warning btn-sm", 
                         style = "width: 100%;"),
        
        
        
        
            
            hr(),
            
            # h4("8. Download Results"),
            uiOutput("download_buttons")
        ),
        
        mainPanel(
            tabsetPanel(
                tabPanel("📊 Status & Preview", 
                    h4("📊 Data Overview"),
                    verbatimTextOutput("file_info"),
                    hr(),
                    h4("📈 Analysis Status"),
                    verbatimTextOutput("status"),
                    hr(),
                    h4("🎨 Plot Preview"),
                    # plotOutput("plot_preview", height = "600px")
                    # plotOutput("pairwise_diff_plot", height = "600px")
                    
                    fluidRow(
                      column(6, plotOutput("plot_preview")),
                      column(6, plotOutput("ancova_plot"))
                    ),
                    
                    hr(),
                    h4("📊 Statistical Results"),
                    div(
                      style = "border: 1px solid #ddd; border-radius: 4px; padding: 10px; 
                 background-color: #f9f9f9; height: 300px; overflow-y: auto; 
                 font-family: monospace;",
                      verbatimTextOutput("stats_text")
                    ),
                    
                    h4("📊 Statistical Results"),
                    DT::dataTableOutput("stats_table"),
                    
                    h4("🔍 Pairwise Comparisons"),
                    DT::dataTableOutput("pairwise_table")
                    
                    
                ),
                tabPanel("📝 Message Log",
                    h4("System Messages & Debugging Info"),
                    div(style = "margin-bottom: 10px;",
                        actionButton("clear_messages", "Clear Messages", 
                                   class = "btn-sm btn-outline-secondary")
                    ),
                    div(
                        style = "border: 1px solid #ddd; border-radius: 4px; padding: 10px; 
                                 background-color: #f9f9f9; height: 400px; overflow-y: auto; 
                                 font-family: monospace; font-size: 12px; white-space: pre-wrap;",
                        verbatimTextOutput("message_log", placeholder = TRUE)
                    )
                ),
                tabPanel("📋 Data Table",
                    conditionalPanel(
                        condition = "input.data_type == 'dsi'",
                        h4("Processed Data Preview"),
                        DT::dataTableOutput("data_table")
                    ),
                    conditionalPanel(
                        condition = "input.data_type == 'clams'",
                        h4("Data Table View"),
                        p("Data table view is only available for DSI Telemetry data."),
                        p("For CLAMS data, please use the Status & Preview tab.")
                    )
                )
            )
        )
    ),
    
    tags$hr(),
    tags$div(
      style = "text-align:center; color:#777; font-size:12px;",
      HTML("&copy; 2026 RMPC@UPenn | CLAMS Data Processing App v2.3.9")
    )
)

# Define the server logic
server_app1e <- function(input, output, session) {
  
  update_trigger <- reactiveVal(0)
  
  observeEvent(input$update_view, {
    req(analysis_results())
    update_trigger(update_trigger() + 1)
    add_message("Plot/statistics updated", "success")
  })
  
  metadata <-  reactiveVal(NULL)
  
  # r2.2.7
  # Helper: find first existing column among candidates (case-insensitive)
  find_col <- function(df, candidates) {
    ln <- tolower(names(df))
    for (c in tolower(candidates)) {
      hit <- which(ln == c)
      if (length(hit) > 0) return(names(df)[hit[1]])
    }
    return(NULL)
  }
  
  
    
    # Set maximum file upload size to 50GB
    options(shiny.maxRequestSize = 50 * 1024^3)  # 50GB in bytes
    
    # Message logging system
    message_log <- reactiveValues(messages = character(0))
    
    # Analysis state tracking (for clearing outputs)
    analysis_state <- reactiveValues(
        has_run = FALSE,
        results = NULL,
        plot = NULL,
        status_text = NULL,
        clams_folder_selected = FALSE,
        clams_folder_path = NULL,
        metadata_loaded = FALSE,
        excel_loaded = FALSE
    )
    
    # Function to add messages to the log
    add_message <- function(message, type = "info") {
        timestamp <- format(Sys.time(), "%H:%M:%S")
        icon <- switch(type,
                      "error" = "❌",
                      "warning" = "⚠️", 
                      "success" = "✅",
                      "info" = "ℹ️",
                      "📝")
        
        formatted_message <- paste0("[", timestamp, "] ", icon, " ", message)
        message_log$messages <- c(message_log$messages, formatted_message)
        
        # Keep only last 100 messages to prevent memory issues
        if (length(message_log$messages) > 100) {
            message_log$messages <- tail(message_log$messages, 100)
        }
        
        # Also show as notification for immediate feedback
        showNotification(message, type = if(type == "error") "error" else if(type == "warning") "warning" else "message")
    }
    
    # Clear messages handler
    observeEvent(input$clear_messages, {
        message_log$messages <- character(0)
        add_message("Message log cleared", "info")
    })
    
    # Clear all data and reset app handler
    observeEvent(input$clear_all, {
        showModal(modalDialog(
            title = "🔄 Clear All Data?",
            "This will reset the app and clear all uploaded files, analysis results, and settings. Are you sure?",
            footer = tagList(
                modalButton("Cancel"),
                actionButton("confirm_clear", "Yes, Clear Everything", class = "btn-danger")
            )
        ))
    })
    
    # Confirm clear all
    observeEvent(input$confirm_clear, {
        # Clear message log
        message_log$messages <- character(0)
        
        # Clear analysis state
        analysis_state$has_run <- FALSE
        analysis_state$results <- NULL
        analysis_state$plot <- NULL
        analysis_state$status_text <- NULL
        analysis_state$clams_folder_selected <- FALSE
        analysis_state$clams_folder_path <- NULL
        analysis_state$metadata_loaded <- FALSE
        analysis_state$excel_loaded <- FALSE
        
        # Reset file inputs by updating their values (this triggers re-rendering)
        shinyjs::reset("metadata_file")
        shinyjs::reset("excel_file")
        
        # Reset CLAMS directory selection by resetting the shinyDirButton
        # This is done by re-initializing the volumes
        # 10/8/2026
        # shinyDirChoose(input, 'data_dir', roots = getVolumes()(), session = session)
        shinyjs::reset("clams_files")
        
        
        # Reset experiment name
        updateTextInput(session, "experiment_name", value = "MetabolicExperiment")
        
        # Reset plot options
        updateRadioButtons(session, "plot_color_by", selected = "Conditions")
        updateCheckboxInput(session, "use_faceting", value = FALSE)
        updateRadioButtons(session, "dsi_plot_group_by", selected = "mouse")
        updateCheckboxInput(session, "dsi_plot_facet", value = FALSE)
        
        # Reset Y-axis selections
        updateSelectInput(session, "clams_y_var", selected = "si13c")
        updateSelectInput(session, "dsi_y_var", choices = c(), selected = NULL)
        
        # Reset duration inputs to defaults
        updateNumericInput(session, "duration_of_experiment", value = 150)
        updateNumericInput(session, "duration_before_experiment", value = 60)
        updateNumericInput(session, "dsi_duration_of_experiment", value = 150)
        updateNumericInput(session, "dsi_duration_before_experiment", value = 60)
        
        # Close the modal
        removeModal()
        
        # Show success message
        add_message("✅ App reset successfully! All data, folders, and outputs cleared. Ready for new analysis.", "success")
        showNotification("App reset! All files and outputs cleared.", type = "message", duration = 3)
    })
    
    # Message log output
    output$message_log <- renderText({
        if (length(message_log$messages) == 0) {
            return("No messages yet. Select data type and upload files to see status updates here...")
        }
        paste(message_log$messages, collapse = "\n")
    })
    
    # Data type change handler
    observeEvent(input$data_type, {
        data_type_name <- if(input$data_type == "clams") "CLAMS Metabolic Cage" else "DSI Telemetry"
        add_message(paste0("Switched to ", data_type_name, " analysis mode"), "info")
    })
    
    # DSI data processing function adapted from the original DSI app
    process_dsi_data_unified <- function(file_path, metadata = NULL, experiment_duration = 90, baseline_duration = 0, time_anchor = "experiment") {
        # Helper to parse per-mouse start datetime robustly (12h or 24h, ISO or US date)
        parse_mouse_datetime <- function(date_str, time_str) {
            candidates <- c(
                "%Y-%m-%d %I:%M:%S %p",
                "%Y-%m-%d %H:%M:%S",
                "%m/%d/%Y %I:%M:%S %p",
                "%m/%d/%Y %H:%M:%S",
                "%Y-%m-%d %I:%M %p",
                "%Y-%m-%d %H:%M",
                "%m/%d/%Y %I:%M %p",
                "%m/%d/%Y %H:%M",
                "%m/%d/%Y %H:%M:%S",  # Added for 24-hour US date format
                "%Y-%m-%d %H:%M:%S"   # Added for 24-hour ISO date format
            )
            dt <- NA
            for (fmt in candidates) {
                dt <- suppressWarnings(as.POSIXct(paste(date_str, time_str), format = fmt, tz = "UTC"))
                if (!is.na(dt)) {
                    add_message(paste0("Successfully parsed datetime: ", date_str, " ", time_str, " using format: ", fmt), "info")
                    break
                }
            }
            
            if (is.na(dt)) {
                error_msg <- paste0(
                    "DATE/TIME ERROR: Could not parse '", date_str, " ", time_str, "'\n\n",
                    "Expected formats:\n",
                    "• 2024-02-22 11:52:31 AM (ISO date + 12-hour time)\n",
                    "• 2/22/2024 11:52:31 AM (US date + 12-hour time)\n",
                    "• 2024-02-22 14:17:00 (ISO date + 24-hour time)\n",
                    "• 2/22/2024 14:17:00 (US date + 24-hour time)\n\n",
                    "Please check your metadata file format!"
                )
                add_message(error_msg, "error")
            }
            dt
        }
        
        tryCatch({
            add_message("Checking Excel file sheets...", "info")
            # Get available sheets
            sheets <- excel_sheets(file_path)
            add_message(paste0("Available sheets: ", paste(sheets, collapse = ", ")), "info")
            
            # Try to find "Parameters_All Subjects" sheet or use first sheet
            sheet_to_use <- if ("Parameters_All Subjects" %in% sheets) {
                "Parameters_All Subjects"
            } else {
                sheets[1]
            }
            add_message(paste0("Using sheet: ", sheet_to_use), "info")
            
            add_message("Attempting to read Excel file headers...", "info")
            # Read the data structure to understand headers
            raw_data_no_headers <- read_excel(file_path, sheet = sheet_to_use, col_names = FALSE)
            add_message(paste0("Read raw data: ", nrow(raw_data_no_headers), " rows × ", ncol(raw_data_no_headers), " columns"), "info")
            
            # Extract header information (first 3 rows)
            mouse_row <- as.character(raw_data_no_headers[1, ])
            measurement_type_row <- as.character(raw_data_no_headers[2, ])
            measurement_detail_row <- as.character(raw_data_no_headers[3, ])
            
            add_message("Attempting to read Excel data (skip 3 rows)...", "info")
            # Skip first 3 header rows for actual data
            data <- read_excel(file_path, sheet = sheet_to_use, skip = 3, col_names = TRUE)
            add_message(paste0("Read data: ", nrow(data), " rows × ", ncol(data), " columns"), "info")
            
        }, error = function(e) {
            add_message(paste("Error reading Excel structure:", e$message), "error")
            # Fallback: use simple read
            tryCatch({
                add_message("Trying fallback Excel read...", "info")
                data <- read_excel(file_path, sheet = sheet_to_use, skip = 2, col_names = TRUE)
                mouse_row <- rep("Unknown", ncol(data))
                measurement_type_row <- rep("Unknown", ncol(data))
                measurement_detail_row <- names(data)
                add_message(paste0("Fallback read successful: ", nrow(data), " rows × ", ncol(data), " columns"), "info")
            }, error = function(e2) {
                add_message(paste("Fallback Excel read also failed:", e2$message), "error")
                # Try reading from beginning without skipping
                tryCatch({
                    add_message("Trying to read from beginning (no skip)...", "info")
                    data <- read_excel(file_path, sheet = sheet_to_use, col_names = TRUE)
                    mouse_row <- rep("Unknown", ncol(data))
                    measurement_type_row <- rep("Unknown", ncol(data))
                    measurement_detail_row <- names(data)
                    add_message(paste0("No-skip read successful: ", nrow(data), " rows × ", ncol(data), " columns"), "info")
                }, error = function(e3) {
                    add_message(paste("All Excel read methods failed:", e3$message), "error")
                    stop("Cannot read Excel file with any method")
                })
            })
        })
        

                
        # Clean up the data
        clean_data <- data
        
        # Get general columns (first 6) and parse timestamps
        general_cols <- clean_data[, 1:6]
        names(general_cols) <- c("StartDate_Raw", "StartTime_Raw", "EndDate_Raw", "EndTime_Raw", "TimeZone", "ElapsedTime")
        
        # Parse the ISO format timestamps
        general_cols$StartDateTime <- as.POSIXct(general_cols$StartDate_Raw, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
        general_cols$EndDateTime <- as.POSIXct(general_cols$EndDate_Raw, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
        
        # Extract readable date and time components
        general_cols$StartDate <- format(general_cols$StartDateTime, "%Y-%m-%d")
        general_cols$StartTime <- format(general_cols$StartDateTime, "%I:%M:%S %p")
        general_cols$EndDate <- format(general_cols$EndDateTime, "%Y-%m-%d")
        general_cols$EndTime <- format(general_cols$EndDateTime, "%I:%M:%S %p")
        
        # Process mouse-specific columns (starting from column 7)
        if (ncol(clean_data) > 6) {
            mouse_data_cols <- clean_data[, 7:ncol(clean_data)]
            
            # Create proper column names for mouse data
            new_col_names <- c()
            mouse_counter <- 1
            current_mouse <- ""
            
            # Track which mice we've seen for logging
            seen_mice <- character()
            
            for (i in 7:length(mouse_row)) {
                mouse_id <- mouse_row[i]
                measurement_type <- measurement_type_row[i]
                measurement_detail <- measurement_detail_row[i]
                
                # Extract mouse number from header like 'HD-XG(1312228)P12R12' -> Mouse_12
                # Try multiple patterns
                extracted_num <- NA_integer_
                if (!is.na(mouse_id) && nzchar(mouse_id)) {
                    # Pattern 1: Look for R followed by number (R12, R1, etc.)
                    if (grepl("R\\s*[0-9]+", mouse_id)) {
                        num_str <- sub(".*R\\s*([0-9]+).*", "\\1", mouse_id)
                        extracted_num <- suppressWarnings(as.integer(num_str))
                    }
                    # Pattern 2: Look for P followed by number (P12, P1, etc.) if R pattern fails
                    if (is.na(extracted_num) && grepl("P\\s*[0-9]+", mouse_id)) {
                        num_str <- sub(".*P\\s*([0-9]+).*", "\\1", mouse_id)
                        extracted_num <- suppressWarnings(as.integer(num_str))
                    }
                }
                
                if (!is.na(extracted_num)) {
                    current_mouse <- paste0("Mouse_", extracted_num)
                    if (!(current_mouse %in% seen_mice)) {
                        seen_mice <- c(seen_mice, current_mouse)
                    }
                } else {
                    # Fallback: sequential numbering when header changes
                    if (!is.na(mouse_id) && mouse_id != "") {
                        if (i == 7 || mouse_id != mouse_row[i-1]) {
                            current_mouse <- paste0("Mouse_", mouse_counter)
                            if (!(current_mouse %in% seen_mice)) {
                                seen_mice <- c(seen_mice, current_mouse)
                            }
                            mouse_counter <- mouse_counter + 1
                        }
                    } else {
                        current_mouse <- paste0("Mouse_", mouse_counter)
                        if (!(current_mouse %in% seen_mice)) {
                            seen_mice <- c(seen_mice, current_mouse)
                        }
                        mouse_counter <- mouse_counter + 1
                    }
                }
                
                # Clean up measurement type
                if (!is.na(measurement_type) && measurement_type != "") {
                    type_clean <- gsub("[^A-Za-z0-9_]", "_", measurement_type)
                } else {
                    type_clean <- "Unknown"
                }
                
                # Clean up measurement detail
                if (!is.na(measurement_detail) && measurement_detail != "") {
                    detail_clean <- gsub("[^A-Za-z0-9_]", "_", measurement_detail)
                    detail_clean <- gsub("_+", "_", detail_clean)
                    detail_clean <- gsub("^_|_$", "", detail_clean)
                } else {
                    detail_clean <- paste0("Col", i-6)
                }
                
                # Combine into clean column name
                col_name <- paste(current_mouse, type_clean, detail_clean, sep = "_")
                col_name <- gsub("_+", "_", col_name)
                col_name <- gsub("^_|_$", "", col_name)
                
                new_col_names <- c(new_col_names, col_name)
            }
            
            # Apply new column names
            names(mouse_data_cols) <- new_col_names
            
            # Create unified time-based dataset
            unique_mice <- unique(gsub("_.*", "", new_col_names))
            unique_mice <- unique_mice[unique_mice != ""]
            
            # Log detected mice
            add_message(paste0("Detected mice from Excel headers: ", paste(unique_mice, collapse = ", ")), "info")
            
            # Start with helper StartDateTime for filtering
            processed_data <- data.frame(StartDateTime = general_cols$StartDateTime, stringsAsFactors = FALSE)
            
            # Add mouse-specific columns
            for (mouse in unique_mice) {
                mouse_cols <- grep(paste0("^", mouse, "_"), names(mouse_data_cols), value = TRUE)
                
                if (length(mouse_cols) > 0) {
                    # Add mouse-specific time columns
                    processed_data[[paste0(mouse, "_StartDate")]] <- general_cols$StartDate
                    processed_data[[paste0(mouse, "_StartTime")]] <- general_cols$StartTime
                    processed_data[[paste0(mouse, "_EndDate")]] <- general_cols$EndDate
                    processed_data[[paste0(mouse, "_EndTime")]] <- general_cols$EndTime
                    processed_data[[paste0(mouse, "_TimeZone")]] <- general_cols$TimeZone
                    processed_data[[paste0(mouse, "_ElapsedTime")]] <- general_cols$ElapsedTime
                    
                    # Add measurement columns
                    for (col in mouse_cols) {
                        processed_data[[col]] <- mouse_data_cols[[col]]
                    }
                }
            }
        } else {
            processed_data <- general_cols
        }
        
        # Filter data based on metadata if provided
        if (!is.null(metadata)) {
            # Validate metadata
            required_cols <- c("mouse_id", "start_date", "start_time")
            missing_cols <- setdiff(required_cols, names(metadata))
            if (length(missing_cols) > 0) {
                add_message(paste0("Metadata missing required columns: ", paste(missing_cols, collapse = ", ")), "error")
                metadata <- NULL
            } else {
                metadata$mouse_id <- trimws(as.character(metadata$mouse_id))
                metadata$start_date <- trimws(as.character(metadata$start_date))
                metadata$start_time <- trimws(as.character(metadata$start_time))
                
                # Debug: Show available mice in data
                unique_mice <- unique(gsub("_.*", "", names(processed_data)))
                unique_mice <- unique_mice[unique_mice != "StartDateTime" & unique_mice != ""]
                add_message(paste0("Available mice in data: ", paste(unique_mice, collapse = ", ")), "info")
                add_message(paste0("Mice in metadata: ", paste(metadata$mouse_id, collapse = ", ")), "info")
            }
        }
        
        if (!is.null(metadata)) {
            # Process with metadata filtering (complex time windows)
            mouse_dfs <- list()
            failed_mice <- character(0)
            successful_mice <- character(0)
            
            for (i in seq_len(nrow(metadata))) {
                mouse_id <- metadata$mouse_id[i]
                start_date_str <- metadata$start_date[i]
                start_time_str <- metadata$start_time[i]
                
                add_message(paste0("Processing metadata row ", i, ": ", mouse_id, " | ", start_date_str, " | ", start_time_str), "info")
                
                # Check if mouse exists in data
                mouse_cols <- grep(paste0("^", mouse_id, "_"), names(processed_data), value = TRUE)
                if (length(mouse_cols) == 0) {
                    # Try alternative patterns
                    alt_patterns <- c(
                        paste0("^", gsub("_", "", mouse_id), "_"),  # Remove underscores
                        paste0("^", tolower(mouse_id), "_"),        # Try lowercase
                        paste0("^", toupper(mouse_id), "_")         # Try uppercase
                    )
                    for (pattern in alt_patterns) {
                        mouse_cols <- grep(pattern, names(processed_data), value = TRUE, ignore.case = TRUE)
                        if (length(mouse_cols) > 0) {
                            add_message(paste0("Found mouse using alternative pattern: ", pattern), "info")
                            break
                        }
                    }
                }
                
                if (length(mouse_cols) == 0) {
                    failed_mice <- c(failed_mice, paste0(mouse_id, " (not found in data)"))
                    add_message(paste0("ERROR: Mouse '", mouse_id, "' not found in data. Searched for columns starting with '", mouse_id, "_'"), "error")
                    next
                } else {
                    add_message(paste0("Found ", length(mouse_cols), " columns for ", mouse_id), "info")
                }
                
                # Parse mouse-specific date/time
                experiment_start <- parse_mouse_datetime(start_date_str, start_time_str)
                if (is.na(experiment_start)) {
                    failed_mice <- c(failed_mice, paste0(mouse_id, " (invalid date/time)"))
                    add_message(paste0("ERROR: Failed to parse date/time for ", mouse_id, ": ", start_date_str, " ", start_time_str), "error")
                    next
                }
                
                baseline_start <- experiment_start - (baseline_duration * 60)
                experiment_end <- experiment_start + (experiment_duration * 60)
                
                add_message(paste0(mouse_id, " - Experiment window: ", 
                                 format(experiment_start, "%Y-%m-%d %H:%M:%S"), " to ",
                                 format(experiment_end, "%Y-%m-%d %H:%M:%S")), "info")
                add_message(paste0(mouse_id, " - Baseline start: ", 
                                 format(baseline_start, "%Y-%m-%d %H:%M:%S")), "info")
                
                # Debug: Show data time range
                if ("StartDateTime" %in% names(processed_data)) {
                    data_times <- processed_data$StartDateTime[!is.na(processed_data$StartDateTime)]
                    if (length(data_times) > 0) {
                        add_message(paste0("Data time range: ", 
                                         format(min(data_times), "%Y-%m-%d %H:%M:%S"), " to ",
                                         format(max(data_times), "%Y-%m-%d %H:%M:%S")), "info")
                    }
                }
                
                # Filter data for this mouse's time window
                if (identical(time_anchor, "first_sample")) {
                    idx <- which(processed_data$StartDateTime >= baseline_start)
                    if (length(idx) == 0) {
                        failed_mice <- c(failed_mice, paste0(mouse_id, " (no data after baseline start)"))
                        add_message(paste0("ERROR: No data found after baseline start for ", mouse_id), "error")
                        next
                    }
                    tvec <- processed_data$StartDateTime[idx]
                    ord <- order(tvec)
                    idx <- idx[ord]
                    tvec <- tvec[ord]
                    
                    t0_pos <- which(tvec >= experiment_start)[1]
                    if (is.na(t0_pos)) {
                        failed_mice <- c(failed_mice, paste0(mouse_id, " (no data after experiment start)"))
                        next
                    }
                    t0 <- tvec[t0_pos]
                    end_limit <- t0 + (experiment_duration * 60)
                    keep <- tvec <= (end_limit + 1)
                    idx <- idx[keep]
                    tvec <- tvec[keep]
                } else {
                    idx <- which(processed_data$StartDateTime >= baseline_start &
                                processed_data$StartDateTime <= (experiment_end + 1))
                    if (length(idx) == 0) {
                        failed_mice <- c(failed_mice, paste0(mouse_id, " (no data in time window)"))
                        add_message(paste0("ERROR: No data found in time window for ", mouse_id, 
                                         ". Window: ", format(baseline_start, "%Y-%m-%d %H:%M:%S"), 
                                         " to ", format(experiment_end, "%Y-%m-%d %H:%M:%S")), "error")
                        next
                    }
                    tvec <- processed_data$StartDateTime[idx]
                    t0 <- experiment_start
                    add_message(paste0("Found ", length(idx), " data points in time window for ", mouse_id), "info")
                }
                
                successful_mice <- c(successful_mice, mouse_id)
                
                # Calculate relative time
                rel_minutes <- as.numeric(difftime(tvec, t0, units = "mins"))
                
                # Create mouse-specific data frame
                df <- data.frame(RelativeTime_Minutes = rel_minutes, stringsAsFactors = FALSE)
                
                # Add mouse DateTime column for time checking
                df[[paste0(mouse_id, "_DateTime")]] <- tvec
                
                # Add mouse measurement columns
                for (col in mouse_cols) {
                    df[[col]] <- processed_data[[col]][idx]
                }
                
                mouse_dfs[[mouse_id]] <- df
            }
            
            if (length(mouse_dfs) == 0) {
                add_message("❌ No mice could be processed with the given metadata", "error")
                if (length(failed_mice) > 0) {
                    add_message(paste0("Failed mice: ", paste(failed_mice, collapse = ", ")), "error")
                }
                return(data.frame())
            } else {
                add_message(paste0("✅ Successfully processed ", length(successful_mice), " mice: ", 
                                 paste(successful_mice, collapse = ", ")), "success")
                if (length(failed_mice) > 0) {
                    add_message(paste0("⚠️ Failed to process ", length(failed_mice), " mice: ", 
                                     paste(failed_mice, collapse = ", ")), "warning")
                }
            }
            
            # Align all mice by relative time
            min_rows <- min(sapply(mouse_dfs, nrow))
            first_mouse <- names(which.min(sapply(mouse_dfs, function(df) min(df$RelativeTime_Minutes, na.rm = TRUE))))
            final_df <- data.frame(RelativeTime_Minutes = mouse_dfs[[first_mouse]]$RelativeTime_Minutes[1:min_rows], stringsAsFactors = FALSE)
            
            for (m in names(mouse_dfs)) {
                dfm <- mouse_dfs[[m]][1:min_rows, , drop = FALSE]
                dfm$RelativeTime_Minutes <- NULL  # Remove duplicate time column
                final_df <- cbind(final_df, dfm)
            }
            
            processed_data <- final_df
            add_message(paste0("Successfully processed ", length(successful_mice), " mice with metadata filtering"), "success")
        } else {
            # No metadata: simple RelativeTime_Minutes based on first timestamp
            if ("StartDateTime" %in% names(processed_data)) {
                idx0 <- which(!is.na(processed_data$StartDateTime))[1]
                if (!is.na(idx0)) {
                    t0 <- processed_data$StartDateTime[idx0]
                    processed_data$RelativeTime_Minutes <- as.numeric(difftime(processed_data$StartDateTime, t0, units = "mins"))
                    # Put RelativeTime_Minutes first
                    ord <- c("RelativeTime_Minutes", setdiff(names(processed_data), c("RelativeTime_Minutes")))
                    processed_data <- processed_data[, ord, drop = FALSE]
                }
                processed_data$StartDateTime <- NULL
            }
            add_message("DSI data processed with basic time indexing (no metadata filtering)", "info")
        }
        
        return(processed_data)
    }
    
    # ================ CLAMS-SPECIFIC COMPONENTS ================
    # 10/8/2026
    # # shinyFiles server logic for directory selection (CLAMS)
    # volumes <- getVolumes()()
    # shinyDirChoose(input, 'data_dir', roots = volumes, session = session)
    # 
    # # Reactive path for CLAMS
    # data_path <- reactive({
    #     # If folder was cleared, return NULL
    #     if (!analysis_state$clams_folder_selected) {
    #         return(NULL)
    #     }
    #     return(parseDirPath(volumes, input$data_dir))
    # })
    # 
    # # Display the selected path for CLAMS
    # output$dir_path_display <- renderPrint({
    #     # Check if cleared
    #     if (!analysis_state$clams_folder_selected) {
    #         return("No folder selected")
    #     }
    #     req(data_path())
    #     data_path()
    # })
    # 
    # # Track directory selection to add message (CLAMS)
    # observeEvent(input$data_dir, {
    #     path <- parseDirPath(volumes, input$data_dir)
    #     if (!is.null(path) && length(path) > 0) {
    #         analysis_state$clams_folder_selected <- TRUE
    #         analysis_state$clams_folder_path <- path
    #         csv_files <- list.files(path = path, pattern = "\\.csv$", full.names = FALSE)
    #         add_message(paste0("CLAMS data folder selected: ", basename(path), " (", length(csv_files), " CSV files found)"), "success")
    #     }
    # })
    
    output$clams_files_display <- renderPrint({
      if (is.null(input$clams_files)) {
        return("No CLAMS files selected")
      }
      input$clams_files$name
    })
    
    observeEvent(input$clams_files, {
      if (!is.null(input$clams_files)) {
        analysis_state$clams_folder_selected <- TRUE
        analysis_state$clams_folder_path <- NULL
        add_message(
          paste0("CLAMS files uploaded: ", nrow(input$clams_files), " CSV file(s)"),
          "success"
        )
      }
    })
    
    
    
    
    # ================ DSI-SPECIFIC COMPONENTS ================
    
    # DSI Excel file processing
    available_sheets <- reactive({
        if (is.null(input$excel_file)) return(NULL)
        
        tryCatch({
            sheets <- excel_sheets(input$excel_file$datapath)
            return(sheets)
        }, error = function(e) {
            add_message(paste("Error reading Excel sheets:", e$message), "error")
            return(NULL)
        })
    })
    
    # Get available Excel sheets for DSI
    available_sheets <- reactive({
        req(input$excel_file)
        tryCatch({
            readxl::excel_sheets(input$excel_file$datapath)
        }, error = function(e) {
            add_message(paste("Error reading Excel sheets:", e$message), "error")
            return(NULL)
        })
    })
    
    # Sheet selector for DSI
    output$sheet_selector <- renderUI({
        sheets <- available_sheets()
        if (is.null(sheets)) {
            return(p("Upload an Excel file to see available sheets"))
        }
        
        # Default to "Parameters_All Subjects" if available
        default_sheet <- if ("Parameters_All Subjects" %in% sheets) "Parameters_All Subjects" else sheets[1]
        
        selectInput("sheet_name", "Select Excel Sheet:", 
                   choices = sheets, 
                   selected = default_sheet)
    })
    
    # DSI raw data processing
    raw_data <- reactive({
        req(input$excel_file)
        
        sheet_to_use <- if (!is.null(input$sheet_name) && input$sheet_name != "") {
            input$sheet_name
        } else {
            tryCatch({
                sheets <- excel_sheets(input$excel_file$datapath)
                if ("Parameters_All Subjects" %in% sheets) {
                    "Parameters_All Subjects"
                } else {
                    sheets[1]
                }
            }, error = function(e) 1)
        }
        
        tryCatch({
            data <- read_excel(input$excel_file$datapath, 
                             sheet = sheet_to_use,
                             skip = input$start_row - 1,
                             col_names = TRUE)
            
            if (nrow(data) == 0) {
                data <- read_excel(input$excel_file$datapath, 
                                 sheet = sheet_to_use,
                                 col_names = TRUE)
            }
            
            return(data)
        }, error = function(e) {
            add_message(paste("Error reading Excel sheet:", e$message), "error")
            return(data.frame())
        })
    })
    
    # Track DSI file uploads
    observeEvent(input$excel_file, {
        if (!is.null(input$excel_file)) {
            analysis_state$excel_loaded <- TRUE
            add_message(paste0("DSI Excel file uploaded: ", input$excel_file$name), "success")
        }
    })
    
    # ================ SHARED COMPONENTS ================
    
    # Track metadata uploads
    observeEvent(input$metadata_file, {
        if (!is.null(input$metadata_file)) {
            analysis_state$metadata_loaded <- TRUE
            tryCatch({
                req(input$metadata_file, input$start_time_col)
                md <- read.csv(input$metadata_file$datapath, stringsAsFactors = FALSE)
                
              #r2.2.13
              # Original "start_time" → renamed to "start_time0"
              # Selected column → renamed to "start_time"
              req(input$start_time_col)
              md$start_time <- md[[input$start_time_col]]
              
              metadata(md)
              
              add_message(
                paste(
                  "Metadata loaded.",
                  "Using", input$start_time_col,
                  "as start_time."
                ),
                "success"
              )
              
              
                # Validate metadata structure based on data type
                if (input$data_type == "clams") {
                    required_cols <- c("file_index", "start_time", "mouse_id", "cage_num", "sample_name", "sample_dataset", "condition")
                } else {
                    # DSI metadata uses lowercase column names: mouse_id, start_date, start_time
                    required_cols <- c("mouse_id", "start_date", "start_time")
                }
                
                missing_cols <- setdiff(required_cols, colnames(metadata))
                
                if (length(missing_cols) > 0) {
                    add_message(paste0("Metadata loaded but missing columns: ", paste(missing_cols, collapse = ", ")), "warning")
                } else {
                    condition_col <- if (input$data_type == "clams") "condition" else "Condition"
                    num_conditions <- if (condition_col %in% colnames(metadata)) length(unique(metadata[[condition_col]])) else 0
                    add_message(paste0("Metadata file loaded: ", input$metadata_file$name, " (", nrow(metadata), " entries, ", num_conditions, " conditions)"), "success")
                }
                
            }, error = function(e) {
                add_message(paste("Failed to read metadata:", e$message), "error")
            })
        }
    })
    
    # Enhanced file info display function
    get_file_info_static <- function() {
        if (input$data_type == "clams") {
            # CLAMS file info logic
            # Check if folder was cleared
          # 10/8/2026
            # if (!analysis_state$clams_folder_selected || is.null(data_path()) || length(data_path()) == 0) {
            #     return(list(files = c(), file_info = "No CLAMS data folder selected", cage_numbers = c()))
            # }
            
          if (is.null(input$clams_files) || nrow(input$clams_files) == 0) {
            return(list(files = c(), file_info = "No CLAMS CSV files uploaded", cage_numbers = c()))
          }
          
          
            tryCatch({
              # 10/8/2026
                # csv_files <- list.files(path = data_path(), pattern = "\\.csv$", full.names = FALSE)
                csv_files <- input$clams_files$name
                
                if (length(csv_files) == 0) {
                    return(list(files = csv_files, file_info = "No CSV files found in selected directory", cage_numbers = c()))
                }
                
                # Analyze CLAMS cage structure
                all_cage_numbers <- c()
                sample_cage_cols <- c()
                # 10/8/2026
                # for (csv_file in csv_files) {
                #     file_path <- file.path(data_path(), csv_file)
                #     file_data <- tryCatch({
                #         read.csv(file_path, nrows = 1)
                #     }, error = function(e) NULL)
                    
                    for (i in seq_len(nrow(input$clams_files))) {
                      csv_file <- input$clams_files$name[i]
                      file_path <- input$clams_files$datapath[i]
                      file_data <- tryCatch({
                        read.csv(file_path, nrows = 1)
                      }, error = function(e) NULL)    
                    
                    
                    if (!is.null(file_data)) {
                        file_cage_cols <- grep("_\\d+$", colnames(file_data), value = TRUE)
                        if (csv_file == csv_files[1]) {
                            sample_cage_cols <- file_cage_cols
                        }
                        
                        if (length(file_cage_cols) > 0) {
                            file_cage_numbers <- unique(gsub(".*_(\\d+)$", "\\1", file_cage_cols))
                            all_cage_numbers <- c(all_cage_numbers, file_cage_numbers)
                        }
                    }
                }
                
                cage_numbers <- sort(as.numeric(unique(all_cage_numbers)))
                
                file_info <- paste0(
                    "CLAMS Data Analysis\n",
                    "Found ", length(csv_files), " CSV files in directory\n",
                    "Sample file: ", csv_files[1], "\n"
                )
                
                if (length(cage_numbers) > 0) {
                    file_info <- paste0(file_info, "Detected cage numbers: ", paste(cage_numbers, collapse = ", "), "\n")
                    file_info <- paste0(file_info, "Sample columns: ", paste(head(sample_cage_cols, 6), collapse = ", "))
                    if (length(sample_cage_cols) > 6) {
                        file_info <- paste0(file_info, " ... (", length(sample_cage_cols), " total)")
                    }
                } else {
                    file_info <- paste0(file_info, "⚠️ No cage columns detected")
                }
                
                return(list(files = csv_files, file_info = file_info, cage_numbers = cage_numbers))
                
            }, error = function(e) {
                return(list(files = c(), file_info = paste("Error reading CLAMS files:", e$message), cage_numbers = c()))
            })
            
        } else {
            # DSI file info logic
            if (is.null(input$excel_file)) {
                return(list(files = c(), file_info = "No DSI Excel file uploaded", mouse_ids = c()))
            }
            
            tryCatch({
                data <- raw_data()
                if (nrow(data) == 0) {
                    return(list(files = c(), file_info = "Excel file uploaded but no data found", mouse_ids = c()))
                }
                
                # Extract mouse IDs from DSI data
                excel_file_path <- input$excel_file$datapath
                sheet_name <- if (!is.null(input$sheet_name) && input$sheet_name != "") {
                    input$sheet_name
                } else {
                    "Parameters_All Subjects"
                }
                
                raw_data_no_headers <- read_excel(excel_file_path, 
                                                sheet = sheet_name,
                                                col_names = FALSE)
                
                mouse_row <- as.character(raw_data_no_headers[1, ])
                mouse_ids <- c()
                
                for (i in 7:length(mouse_row)) {
                    mouse_id <- mouse_row[i]
                    extracted_num <- NA_integer_
                    if (!is.na(mouse_id) && nzchar(mouse_id)) {
                        if (grepl("R\\s*[0-9]+", mouse_id)) {
                            num_str <- sub(".*R\\s*([0-9]+).*", "\\1", mouse_id)
                            extracted_num <- suppressWarnings(as.integer(num_str))
                        }
                    }
                    if (!is.na(extracted_num)) {
                        mouse_ids <- c(mouse_ids, paste0("Mouse_", extracted_num))
                    }
                }
                mouse_ids <- unique(mouse_ids)
                
                # Get date range
                first_date <- data[[1]][1]
                last_date <- data[[1]][nrow(data)]
                
                file_info <- paste0(
                    "DSI Telemetry Data Analysis\n",
                    "Excel file: ", input$excel_file$name, "\n",
                    "Sheet: ", sheet_name, "\n",
                    "Data rows: ", nrow(data), "\n"
                )
                
                if (length(mouse_ids) > 0) {
                    file_info <- paste0(file_info, "Detected mice: ", paste(mouse_ids, collapse = ", "), "\n")
                } else {
                    file_info <- paste0(file_info, "⚠️ No mouse IDs detected\n")
                }
                
                file_info <- paste0(file_info, "Date range: ", first_date, " to ", last_date)
                
                return(list(files = c(input$excel_file$name), file_info = file_info, mouse_ids = mouse_ids))
                
            }, error = function(e) {
                return(list(files = c(), file_info = paste("Error reading DSI file:", e$message), mouse_ids = c()))
            })
        }
    }
    
    # File info display
    output$file_info <- renderText({
        file_info <- get_file_info_static()
        
        text_lines <- c()
        text_lines <- c(text_lines, file_info$file_info)
        
        # Add metadata comparison if available AND not cleared
        if (!is.null(input$metadata_file) && analysis_state$metadata_loaded) {
            tryCatch({
                req(input$metadata_file, input$start_time_col)
                md <- read.csv(input$metadata_file$datapath, stringsAsFactors = FALSE)
                
                #r2.2.13
                # Selected column → let "start_time" = input$start_time_col
                md$start_time <- md[[input$start_time_col]]  
                metadata(md)
               
                if (input$data_type == "clams") {
                    # CLAMS metadata validation
                    if ("file_index" %in% colnames(metadata)) {
                        unique_files <- unique(metadata$file_index)
                        csv_count <- length(file_info$files)
                        
                        text_lines <- c(text_lines, paste("\nMetadata references", length(unique_files), "file(s), found", csv_count, "CSV file(s)"))
                        
                        if (length(unique_files) != csv_count) {
                            text_lines <- c(text_lines, "⚠️ File count mismatch!")
                        } else {
                            text_lines <- c(text_lines, "✅ File counts match")
                        }
                        
                        if ("cage_num" %in% colnames(metadata) && length(file_info$cage_numbers) > 0) {
                            metadata_cages <- unique(gsub("^_", "", metadata$cage_num))
                            data_cages <- as.character(file_info$cage_numbers)
                            
                            missing_in_data <- setdiff(metadata_cages, data_cages)
                            if (length(missing_in_data) > 0) {
                                text_lines <- c(text_lines, paste("⚠️ Cages in metadata but NOT in data:", paste(missing_in_data, collapse = ", ")))
                            } else {
                                text_lines <- c(text_lines, "✅ All cage numbers match")
                            }
                        }
                    }
                } else {
                    # DSI metadata validation - comprehensive checks
                    text_lines <- c(text_lines, "\n=== DSI METADATA VALIDATION ===")
                    
                    # Check 1: Required columns
                    required_cols <- c("mouse_id", "start_date", "start_time")
                    missing_cols <- setdiff(required_cols, colnames(metadata))
                    
                    if (length(missing_cols) > 0) {
                        text_lines <- c(text_lines, paste("❌ Missing required columns:", paste(missing_cols, collapse = ", ")))
                        text_lines <- c(text_lines, "   Required: mouse_id, start_date, start_time")
                    } else {
                        text_lines <- c(text_lines, "✅ All required columns present")
                    }
                    
                    # Check 2: Mouse ID matching
                    if ("mouse_id" %in% colnames(metadata) && length(file_info$mouse_ids) > 0) {
                        metadata_mice <- unique(trimws(as.character(metadata$mouse_id)))
                        excel_mice <- file_info$mouse_ids
                        
                        text_lines <- c(text_lines, paste("\nMetadata mice:", paste(metadata_mice, collapse = ", ")))
                        text_lines <- c(text_lines, paste("Excel mice:   ", paste(excel_mice, collapse = ", ")))
                        
                        missing_in_excel <- setdiff(metadata_mice, excel_mice)
                        missing_in_metadata <- setdiff(excel_mice, metadata_mice)
                        
                        if (length(missing_in_excel) > 0) {
                            text_lines <- c(text_lines, paste("❌ Mice in metadata but NOT in Excel:", paste(missing_in_excel, collapse = ", ")))
                        }
                        if (length(missing_in_metadata) > 0) {
                            text_lines <- c(text_lines, paste("⚠️ Mice in Excel but NOT in metadata:", paste(missing_in_metadata, collapse = ", ")))
                        }
                        if (length(missing_in_excel) == 0 && length(missing_in_metadata) == 0) {
                            text_lines <- c(text_lines, "✅ All mice match between Excel and metadata")
                        }
                    }
                    
                    # Check 3: Date/time format validation
                    if (all(c("start_date", "start_time") %in% colnames(metadata))) {
                        date_time_valid <- TRUE
                        invalid_entries <- c()
                        
                        for (i in 1:min(nrow(metadata), 10)) {  # Check first 10 entries
                            date_str <- trimws(as.character(metadata$start_date[i]))
                            time_str <- trimws(as.character(metadata$start_time[i]))
                            
                            # Try to parse the datetime
                            parsed <- parse_mouse_datetime_check(date_str, time_str)
                            if (is.na(parsed)) {
                                date_time_valid <- FALSE
                                mouse_name <- if("mouse_id" %in% colnames(metadata)) metadata$mouse_id[i] else paste0("Row ", i)
                                invalid_entries <- c(invalid_entries, paste0(mouse_name, " (", date_str, " ", time_str, ")"))
                            }
                        }
                        
                        if (date_time_valid) {
                            text_lines <- c(text_lines, "✅ Date/time formats valid")
                            # Show sample format
                            sample_date <- trimws(as.character(metadata$start_date[1]))
                            sample_time <- trimws(as.character(metadata$start_time[1]))
                            text_lines <- c(text_lines, paste0("   Sample: ", sample_date, " ", sample_time))
                        } else {
                            text_lines <- c(text_lines, "❌ Invalid date/time formats detected:")
                            for (entry in head(invalid_entries, 3)) {
                                text_lines <- c(text_lines, paste0("   ", entry))
                            }
                            if (length(invalid_entries) > 3) {
                                text_lines <- c(text_lines, paste0("   ... and ", length(invalid_entries) - 3, " more"))
                            }
                        }
                    }
                    
                    # Check 4: Date range comparison with Excel data
                    if ("start_date" %in% colnames(metadata) && !is.null(input$excel_file)) {
                        tryCatch({
                            # Get data time range from Excel
                            data <- raw_data()
                            if (nrow(data) > 0 && ncol(data) > 0) {
                                # Parse first timestamp from Excel
                                first_excel_time <- data[[1]][1]
                                last_excel_time <- data[[1]][nrow(data)]
                                
                                # Try to parse Excel timestamps
                                excel_start <- suppressWarnings(as.POSIXct(first_excel_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
                                excel_end <- suppressWarnings(as.POSIXct(last_excel_time, format = "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
                                
                                if (!is.na(excel_start) && !is.na(excel_end)) {
                                    excel_date_range <- paste(format(excel_start, "%Y-%m-%d"), "to", format(excel_end, "%Y-%m-%d"))
                                    text_lines <- c(text_lines, paste("\nExcel data date range:", excel_date_range))
                                    
                                    # Check if metadata dates fall within Excel date range
                                    metadata_dates <- unique(trimws(as.character(metadata$start_date)))
                                    dates_outside_range <- c()
                                    
                                    for (meta_date in metadata_dates) {
                                        # Try parsing metadata date
                                        meta_parsed <- suppressWarnings(as.Date(meta_date, tryFormats = c("%Y-%m-%d", "%m/%d/%Y", "%m-%d-%Y")))
                                        if (!is.na(meta_parsed)) {
                                            excel_start_date <- as.Date(excel_start)
                                            excel_end_date <- as.Date(excel_end)
                                            
                                            if (meta_parsed < excel_start_date || meta_parsed > excel_end_date) {
                                                dates_outside_range <- c(dates_outside_range, meta_date)
                                            }
                                        }
                                    }
                                    
                                    if (length(dates_outside_range) > 0) {
                                        text_lines <- c(text_lines, paste("⚠️ Metadata dates OUTSIDE Excel data range:", paste(dates_outside_range, collapse = ", ")))
                                        text_lines <- c(text_lines, "   This may cause 'no data in time window' errors!")
                                    } else {
                                        text_lines <- c(text_lines, "✅ Metadata dates within Excel data range")
                                    }
                                }
                            }
                        }, error = function(e) {
                            # Silent error - this is just validation feedback
                        })
                    }
                    
                    # Check 5: Summary
                    text_lines <- c(text_lines, "\n=== VALIDATION SUMMARY ===")
                    total_mice_meta <- if("mouse_id" %in% colnames(metadata)) length(unique(metadata$mouse_id)) else 0
                    total_mice_excel <- length(file_info$mouse_ids)
                    text_lines <- c(text_lines, paste0("Metadata entries: ", nrow(metadata), " (", total_mice_meta, " unique mice)"))
                    text_lines <- c(text_lines, paste0("Excel columns: ", total_mice_excel, " mice detected"))
                }
                
            }, error = function(e) {
                text_lines <- c(text_lines, paste("\n⚠️ Error reading metadata:", e$message))
            })
        } else {
            text_lines <- c(text_lines, "\nNo metadata file uploaded")
        }
        
        paste(text_lines, collapse = "\n")
    })
    
    # Enhanced time parsing function
    parse_datetime <- function(time_str, context = "", data_type = "clams") {
        if (is.na(time_str) || time_str == "") {
            add_message(paste0("Empty or missing time value", if(context != "") paste0(" for ", context) else ""), "error")
            return(NA)
        }
        
        if (data_type == "clams") {
            formats_to_try <- c(
                "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %I:%M:%S %p", "%m/%d/%Y %H:%M", "%m/%d/%Y %I:%M %p",
                "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %I:%M:%S %p", "%Y-%m-%d %H:%M",
                "%m-%d-%Y %H:%M:%S", "%d/%m/%Y %H:%M:%S", "%m/%d/%y %H:%M:%S", "%m/%d/%y %H:%M"
            )
        } else {
            formats_to_try <- c(
                "%Y-%m-%dT%H:%M:%SZ", "%m/%d/%Y %I:%M:%S %p", "%Y-%m-%d %H:%M:%S",
                "%m/%d/%Y %H:%M:%S", "%Y-%m-%d %I:%M:%S %p", "%m-%d-%Y %H:%M:%S"
            )
        }
        
        for (fmt in formats_to_try) {
            result <- tryCatch({
                parsed <- as.POSIXct(time_str, format = fmt, tz = "UTC")
                if (!is.na(parsed)) return(parsed)
                NA
            }, error = function(e) NA)
            
            if (!is.na(result)) {
                return(result)
            }
        }
        
        add_message(paste0("Failed to parse time '", time_str, "'", if(context != "") paste0(" for ", context) else ""), "error")
        return(NA)
    }
    
    # Helper function for validating DSI metadata date/time formats (used in validation checks)
    parse_mouse_datetime_check <- function(date_str, time_str) {
        candidates <- c(
            "%Y-%m-%d %I:%M:%S %p",
            "%Y-%m-%d %H:%M:%S",
            "%m/%d/%Y %I:%M:%S %p",
            "%m/%d/%Y %H:%M:%S",
            "%Y-%m-%d %I:%M %p",
            "%Y-%m-%d %H:%M",
            "%m/%d/%Y %I:%M %p",
            "%m/%d/%Y %H:%M",
            "%m/%d/%Y %H:%M:%S",
            "%Y-%m-%d %H:%M:%S"
        )
        
        for (fmt in candidates) {
            dt <- suppressWarnings(as.POSIXct(paste(date_str, time_str), format = fmt, tz = "UTC"))
            if (!is.na(dt)) {
                return(dt)
            }
        }
        return(NA)
    }
    
    # Enable/disable Run Analysis button based on readiness
    observe({
        # if (input$data_type == "clams") {
        #     folder_ready <- !is.null(data_path()) && length(data_path()) > 0
        #     metadata_ready <- !is.null(input$metadata_file)
        #     ready <- folder_ready && metadata_ready
            # 10/8/2026
            if (input$data_type == "clams") {
              files_ready <- !is.null(input$clams_files) && nrow(input$clams_files) > 0
              metadata_ready <- !is.null(input$metadata_file)
              ready <- files_ready && metadata_ready
            
        } else {
            excel_ready <- !is.null(input$excel_file) && tryCatch({
                dat <- raw_data()
                is.data.frame(dat) && nrow(dat) > 0
            }, error = function(e) FALSE)
            metadata_ready <- if (input$use_metadata) !is.null(input$metadata_file) else TRUE
            ready <- excel_ready && metadata_ready
        }
        
        if (ready) {
            shinyjs::enable("run_analysis")
        } else {
            shinyjs::disable("run_analysis")
        }
    })
    
    # Download handler for metadata template
    output$download_metadata_template <- downloadHandler(
        filename = function() {
            if (input$data_type == "clams") {
                "clams_metadata_template.csv"
            } else {
                "dsi_metadata_template.csv"
            }
        },
        content = function(file) {
            if (input$data_type == "clams") {
                template <- data.frame(
                    file_index = c(1, 1, 1, 1, 2, 2, 2, 2),
                    start_time = c("05/15/2024 13:00:22", "05/15/2024 13:10:22", "05/15/2024 13:20:22", "05/15/2024 13:30:22",
                                  "05/16/2024 13:00:11", "05/16/2024 13:10:11", "05/16/2024 13:20:11", "05/16/2024 13:30:11"),
                    mouse_id = c("1", "11", "6", "16", "2", "7", "12", "17"),
                    cage_num = c("_3", "_4", "_11", "_12", "_1", "_2", "_3", "_4"),
                    sample_name = c("HFD_Ctrl_1", "Chow_Ctrl_11", "HFD_Bima_6", "Chow_Bima_16", 
                                   "HFD_Ctrl_2", "HFD_Bima_7", "Chow_Ctrl_12", "Chow_Bima_17"),
                    sample_dataset = c("HFD_Ctrl_1_D1", "Chow_Ctrl_11_D1", "HFD_Bima_6_D1", "Chow_Bima_16_D1",
                                      "HFD_Ctrl_2_D2", "HFD_Bima_7_D2", "Chow_Ctrl_12_D2", "Chow_Bima_17_D2"),
                    condition = c("HFD_CTRL", "CHOW_CTRL", "HFD_BIMA", "CHOW_BIMA", 
                                 "HFD_CTRL", "HFD_BIMA", "CHOW_CTRL", "CHOW_BIMA"),
                    group_order = c(1, 5, 2, 6, 3, 4, 7, 8),
                    stringsAsFactors = FALSE
                )
                add_message("CLAMS metadata template downloaded", "success")
            } else {
                template <- data.frame(
                    MouseID = c("Mouse_1", "Mouse_2", "Mouse_3", "Mouse_4", "Mouse_5", "Mouse_6"),
                    StartDate = c("2024-10-23", "2024-10-23", "2024-10-23", "2024-10-23", "2024-10-23", "2024-10-23"),
                    StartTime = c("17:17:52", "17:27:52", "17:37:52", "17:47:52", "17:57:52", "18:07:52"),
                    ExperimentDuration = c(150, 150, 150, 150, 150, 150),
                    BaselineDuration = c(60, 60, 60, 60, 60, 60),
                    Treatment = c("Control", "Treatment", "Control", "Treatment", "Control", "Treatment"),
                    Condition = c("Group_A", "Group_B", "Group_A", "Group_B", "Group_A", "Group_B"),
                    stringsAsFactors = FALSE
                )
                add_message("DSI metadata template downloaded", "success")
            }
            write.csv(template, file, row.names = FALSE)
        }
    )
    
    # Main analysis function - implements both CLAMS and DSI processing
    analysis_results <- eventReactive(input$run_analysis, {
        # Mark that analysis is running
        analysis_state$has_run <- TRUE
        
        add_message(paste0("Starting ", toupper(input$data_type), " data analysis..."), "info")
        
        if (input$data_type == "clams") {
            # CLAMS Analysis Logic
          # 10/8/2026
            # req(data_path(), input$metadata_file)
          req(input$clams_files, input$metadata_file)
            
            # Show progress indicator
            withProgress(message = 'Processing CLAMS data...', value = 0, {
                
                # Step 1: Read uploaded data and metadata
                incProgress(0.1, detail = "Reading input files...")
                add_message("Reading CSV files from selected directory...", "info")
                
                # Read data files from selected directory, sorting them for consistency
                # 10/8/2026
                # replace 
                # file_paths <- list.files(path = data_path(), pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
                # file_paths <- sort(file_paths)
                # with
                file_paths <- input$clams_files$datapath
                file_names <- input$clams_files$name
                ord <- order(file_names)
                file_paths <- file_paths[ord]
                file_names <- file_names[ord]
        
                if (length(file_paths) == 0) {
                    error_msg <- "No CSV files found in the selected directory."
                    add_message(error_msg, "error")
                    return(NULL)
                }
                
                add_message(paste0("Found ", length(file_paths), " CSV files to process"), "success")
                
                
                #r2.1
                # The target metric prefixes to guarantee
                target_prefixes <- c("si13c", "vo2", "vco2", "rq", "kcal_hr", "vh2o",
                                     "foodupa", "waterupa", "waterupb", "bodymass",
                                     "xbreak", "ybreak", "zbreak", "pedmeters", "allmeters", "bodytemp")
                
 
                data_list <- tryCatch({
                    result <- lapply(file_paths, function(fp) {
                        tryCatch({
                            # read.csv(fp)
                          
                          #r2.1 3/14/2026
                          df <- read.csv(fp)
                          
                          # --- Detect cages from any columns ending with _<num>, e.g., vo2_3, rq_12
                          cage_cols    <- grep("_\\d+$", names(df), value = TRUE)
                          cage_numbers <- sort(as.integer(unique(sub(".*_(\\d+)$", "\\1", cage_cols))))
                          
                          # If we couldn't detect any cages, keep df as is (a later metadata pass can fix it)
                          if (length(cage_numbers) == 0) return(df)
                          
                          # --- Add any missing metric_<cage> columns with -1
                          df <- ensure_metric_cage_cols(
                            df,
                            prefixes     = target_prefixes,
                            cage_numbers = cage_numbers,
                            fill         = -1
                          )
                          
                          df
                          
                          
                        }, error = function(e) {
                            add_message(paste0("Error reading file ", basename(fp), ": ", e$message), "error")
                            return(NULL)
                        })
                    })
                    # 10/8/2026
                    # names(result) <- basename(file_paths)
                    names(result) <- file_names
                    
                    # Remove any NULL entries from failed file reads
                    result[!sapply(result, is.null)]
                }, error = function(e) {
                    error_msg <- paste("Error reading CSV files:", e$message)
                    add_message(error_msg, "error")
                    return(NULL)
                })
                
                if (length(data_list) == 0) {
                    add_message("No valid CSV files could be read", "error")
                    return(NULL)
                }
                
 
                # Read metadata file
                # mouse_metadata <- tryCatch({
                #      read.csv(input$metadata_file$datapath, stringsAsFactors = FALSE)
                # }, error = function(e) {
                #     error_msg <- paste("Error reading metadata file:", e$message)
                #     add_message(error_msg, "error")
                #     return(NULL)
                # })
                
                
                # Access normalized metadata (do NOT re-read CSV)
                mouse_metadata <- tryCatch({
                  req(metadata())
                  metadata()
                }, error = function(e) {
                  error_msg <- paste("Error accessing metadata:", e$message)
                  add_message(error_msg, "error")
                  return(NULL)
                })
                
                
                
                if (is.null(mouse_metadata)) {
                    return(NULL)
                }
                
                add_message(paste0("Metadata loaded: ", nrow(mouse_metadata), " samples across ", length(unique(mouse_metadata$file_index)), " files"), "success")
                
                # Get parameters from UI
                duration_of_experiment <- input$duration_of_experiment
                duration_before_experiment <- input$duration_before_experiment
                
                add_message(paste0("Experiment parameters: ", duration_of_experiment, " min experiment + ", duration_before_experiment, " min baseline"), "info")
                
                # Group metadata by file index
                metadata_by_file <- split(mouse_metadata, mouse_metadata$file_index)
                
                Times <- list()
                Mouse <- list()
                Cages <- list()
                Samples <- list()
                Sample_Dataset <- list()
                Conditions <- list()
                GroupOrder <- NULL  # Will be initialized if group_order exists in metadata
                
                for (i in seq_along(data_list)) {
                    meta <- metadata_by_file[[as.character(i)]]
                    if (!is.null(meta) && nrow(meta) > 0) {
                        Times[[i]] <- meta$start_time
                        Mouse[[i]] <- meta$mouse_id
                        Cages[[i]] <- meta$cage_num
                        Samples[[i]] <- meta$sample_name
                        Sample_Dataset[[i]] <- meta$sample_dataset
                        Conditions[[i]] <- meta$condition
                        
                        # Add group_order if it exists in metadata
                        if ("group_order" %in% names(meta)) {
                            if (is.null(GroupOrder)) GroupOrder <- list()
                            GroupOrder[[i]] <- meta$group_order
                        }
                    } else {
                        add_message(paste0("No metadata found for file index ", i), "warning")
                        Times[[i]] <- character(0)
                        Mouse[[i]] <- character(0)
                        Cages[[i]] <- character(0)
                        Samples[[i]] <- character(0)
                        Sample_Dataset[[i]] <- character(0)
                        Conditions[[i]] <- character(0)
                    }
                }
                
                # Step 2: Process the data
                incProgress(0.3, detail = "Processing datasets...")
                add_message("Processing individual samples and time windows...", "info")
                
                ExtractedData <- list()
                OrderedGroups <- list()
                FinalData <- list()
                OrderedData <- list()
                
                # Standardize column names
                for (i in c(1:length(data_list))) {
                    if (ncol(data_list[[i]]) >= 7) {
                        names(data_list[[i]])[c(1,2,3,4,5,6,7)] <- c("Date_Time", "DurationMin", "envirolightlux", "envirotemp", "envirorh", "envirooccupancy", "envirosound")
                    } else {
                        add_message(paste0("File ", i, " has fewer than 7 columns, column naming may be incorrect"), "warning")
                    }
                }
                
                total_samples <- sum(sapply(Samples, length))
                progress_counter <- 0
                samples_processed <- 0
                samples_skipped <- 0
                
                for (Dataset in seq_along(data_list)) {
                    if (length(Samples[[Dataset]]) == 0) {
                        add_message(paste0("No samples defined for dataset ", Dataset, ", skipping"), "warning")
                        next
                    }
                    
                    for (Sample in seq_along(Samples[[Dataset]])){
                        
                        progress_counter <- progress_counter + 1
                        incProgress(0.5 / total_samples, detail = paste("Processing Sample", progress_counter, "of", total_samples))
                        
                        sample_name <- Samples[[Dataset]][Sample]
                        time_str <- Times[[Dataset]][Sample]
                        
                        
                        
                        # Find matching time in dataset
                        # match_index <- which(data_list[[Dataset]]$Date_Time == time_str)
                        
                        
                        # Convert Date_Time vector
                        dt <- dmy_hms(data_list[[Dataset]]$Date_Time)
                        
                        # Convert target time (no seconds)
                        target_time <- dmy_hm(time_str)
                        
                        # Match by minute
                        match_index <- which(floor_date(dt, "minute") == target_time)
                        
                        
                        # print(target_time)
                        # print(dt)
                        # print(floor_date(dt, "minute"))
                        
                        
                        if (length(match_index) == 0) {
                            # Try fuzzy matching for time
                            data_times <- data_list[[Dataset]]$Date_Time
                            partial_matches <- grep(gsub("[^0-9]", "", time_str), gsub("[^0-9]", "", data_times), fixed = TRUE)
                            
                            if (length(partial_matches) > 0) {
                                match_index <- partial_matches[1]
                                add_message(paste0("Using fuzzy time match for sample '", sample_name, "' at row ", match_index), "warning")
                            } else {
                                error_msg <- paste0("Start time '", time_str, "' not found for sample '", sample_name, "' in dataset ", Dataset)
                                add_message(error_msg, "error")
                                samples_skipped <- samples_skipped + 1
                                ExtractedData[[Sample]] <- NULL
                                next()
                            }
                        }
                        
                        time_before_start <- match_index[1] - duration_before_experiment
                        end_exp_time <- match_index[1] + duration_of_experiment
                        
                        # Calculate expected number of rows: baseline + experiment + 1 (for inclusive range)
                        expected_rows <- duration_before_experiment + duration_of_experiment + 1
                        
                        # Check if the requested time range is valid
                        if (time_before_start < 1 || end_exp_time > nrow(data_list[[Dataset]])) {
                            error_msg <- paste0("Time window out of range for sample '", sample_name, "' - needed rows ", 
                                              time_before_start, " to ", end_exp_time, " but data only has ", nrow(data_list[[Dataset]]), " rows")
                            add_message(error_msg, "warning")
                            samples_skipped <- samples_skipped + 1
                            ExtractedData[[Sample]] <- NULL
                            next()
                        }
                        
                        ExtractedData[[Sample]] <- data_list[[Dataset]] %>% slice(time_before_start:end_exp_time)
                        
                        if(is.null(ExtractedData[[Sample]])){ 
                            samples_skipped <- samples_skipped + 1
                            next 
                        }

                        Standard_Cols_To_Extract <- colnames(ExtractedData[[Sample]])[1:7]
                        
                        # Enhanced cage column detection for CLAMS data
                        cage_number <- gsub("^_", "", Cages[[Dataset]][Sample])  # Remove leading underscore
                        cage_pattern <- paste0("_", cage_number, "$")  # Look for columns ending with _3, _4, etc.
                        all_cage_cols <- grep(cage_pattern, colnames(ExtractedData[[Sample]]), value = TRUE)
                        
                        # Exclude standard columns from animal columns (avoid duplicates)
                        # Standard columns already have generic names, we only want the metric columns
                        standard_col_names <- c("Date_Time", "DurationMin", "envirolightlux", "envirotemp", "envirorh", "envirooccupancy", "envirosound")
                        Animal_Cols_To_Extract <- all_cage_cols[!sub("_\\d+$", "", all_cage_cols) %in% standard_col_names]
                        
                        if (length(Animal_Cols_To_Extract) == 0) {
                            add_message(paste0("No metric columns found for cage '", Cages[[Dataset]][Sample], "' in sample '", sample_name, "'"), "warning")
                        } else {
                            add_message(paste0("Found ", length(Animal_Cols_To_Extract), " metric columns for cage ", cage_number, " (", length(all_cage_cols), " total cage-specific columns)"), "info")
                        }
                        
                        # Debug: Log the columns being extracted
                        add_message(paste0("Extracting ", length(Standard_Cols_To_Extract), " standard + ", 
                                         length(Animal_Cols_To_Extract), " metric columns = ", 
                                         length(c(Standard_Cols_To_Extract, Animal_Cols_To_Extract)), " total columns"), "info")
                        
                        tryCatch({
                            ExtractedData[[Sample]] <- ExtractedData[[Sample]][,c(Standard_Cols_To_Extract,Animal_Cols_To_Extract)]
                            add_message(paste0("Successfully subsetted to ", ncol(ExtractedData[[Sample]]), " columns, ", 
                                             nrow(ExtractedData[[Sample]]), " rows"), "info")
                        }, error = function(e) {
                            add_message(paste0("ERROR during column subsetting for '", sample_name, "': ", e$message), "error")
                            add_message(paste0("Available columns: ", paste(colnames(ExtractedData[[Sample]]), collapse=", ")), "error")
                            add_message(paste0("Trying to extract: ", paste(c(Standard_Cols_To_Extract, Animal_Cols_To_Extract), collapse=", ")), "error")
                            stop(e)
                        })
                        
                        # Calculate actual number of rows extracted
                        actual_rows <- nrow(ExtractedData[[Sample]])
                        expected_rows <- duration_before_experiment + duration_of_experiment + 1
                        
                        add_message(paste0("Sample '", sample_name, "': Extracted ", actual_rows, " rows (expected ", expected_rows, ")"), "info")
                        
                        if (actual_rows != expected_rows) {
                            add_message(paste0("WARNING: Sample '", sample_name, "' has ", actual_rows, 
                                             " rows but expected ", expected_rows, " rows based on duration settings."), "warning")
                        }
                        
                        # Create Minutes sequence based on actual extracted rows
                        # Start from negative baseline duration, count up to match actual rows
                        minutes_start <- -duration_before_experiment
                        minutes_end <- minutes_start + actual_rows - 1
                        
                        add_message(paste0("Adding metadata columns for '", sample_name, "'..."), "info")
                        
                        # Add metadata columns with detailed error handling
                        tryCatch({
                            ExtractedData[[Sample]]$MouseCage <- rep(paste0("Mouse_Cage", Cages[[Dataset]][Sample]), actual_rows)
                            add_message(paste0("  ✓ MouseCage column added"), "info")
                        }, error = function(e) {
                            add_message(paste0("  ✗ ERROR adding MouseCage: ", e$message), "error")
                            stop(e)
                        })
                        
                        tryCatch({
                            ExtractedData[[Sample]]$Minutes <- seq(minutes_start, minutes_end, by = 1)
                            add_message(paste0("  ✓ Minutes column added (", minutes_start, " to ", minutes_end, ")"), "info")
                        }, error = function(e) {
                            add_message(paste0("  ✗ ERROR adding Minutes: ", e$message), "error")
                            stop(e)
                        })
                        
                        tryCatch({
                            ExtractedData[[Sample]]$Samples <- rep(Samples[[Dataset]][Sample], actual_rows)
                            add_message(paste0("  ✓ Samples column added"), "info")
                        }, error = function(e) {
                            add_message(paste0("  ✗ ERROR adding Samples: ", e$message), "error")
                            stop(e)
                        })
                        
                        tryCatch({
                            ExtractedData[[Sample]]$Conditions <- rep(Conditions[[Dataset]][Sample], actual_rows)
                            add_message(paste0("  ✓ Conditions column added"), "info")
                        }, error = function(e) {
                            add_message(paste0("  ✗ ERROR adding Conditions: ", e$message), "error")
                            stop(e)
                        })
                        
                        tryCatch({
                            ExtractedData[[Sample]]$Sample_Dataset <- rep(Sample_Dataset[[Dataset]][Sample], actual_rows)
                            add_message(paste0("  ✓ Sample_Dataset column added"), "info")
                        }, error = function(e) {
                            add_message(paste0("  ✗ ERROR adding Sample_Dataset: ", e$message), "error")
                            stop(e)
                        })
                        
                        # Add group_order if it exists
                        if (!is.null(GroupOrder) && !is.null(GroupOrder[[Dataset]])) {
                            tryCatch({
                                ExtractedData[[Sample]]$group_order <- rep(GroupOrder[[Dataset]][Sample], actual_rows)
                                add_message(paste0("  ✓ group_order column added"), "info")
                            }, error = function(e) {
                                add_message(paste0("  ✗ ERROR adding group_order: ", e$message), "error")
                                stop(e)
                            })
                        }
                        
                        time_to_minutes <- function(time_str) {
                            parts <- unlist(strsplit(time_str, ":"))
                            return(as.numeric(parts[1]) * 60 + as.numeric(parts[2]) + as.numeric(parts[3]) / 60)
                        }
                        
                        # Time_Numeric no longer needed - we use Minutes column for numeric time values
                        add_message(paste0("  ✓ Sample '", sample_name, "' processing complete"), "success")
                        OrderedGroups[[Sample]] <- ExtractedData[[Sample]]
                        
                        # Rename measurement columns intelligently based on original column names
                        # The animal columns were extracted and should be in positions 8 through (7 + num_animal_cols)
                        num_animal_cols <- length(Animal_Cols_To_Extract)
                        
                        if (num_animal_cols > 0) {
                            # Get the original column names of the extracted animal columns
                            original_animal_col_names <- Animal_Cols_To_Extract
                            
                            # Extract the metric name from each column (e.g., "vo2_1" -> "vo2", "si13c_4" -> "si13c")
                            extracted_metrics <- sapply(original_animal_col_names, function(col_name) {
                                # Remove the cage number suffix (e.g., "_1", "_12")
                                sub("_\\d+$", "", col_name)
                            })
                            
                            # Assign these metric names to columns 8 through (7 + num_animal_cols)
                            names(ExtractedData[[Sample]])[8:(7+num_animal_cols)] <- extracted_metrics
                            
                            add_message(paste0("Sample '", sample_name, "': Renamed ", num_animal_cols, 
                                             " metric columns: ", paste(extracted_metrics, collapse=", ")), "info")
                        } else {
                            add_message(paste0("Sample '", sample_name, "' has no animal-specific measurement columns"), "warning")
                        }
                        
                        samples_processed <- samples_processed + 1
                    }
                    
                    FinalData[[Dataset]] <- ExtractedData
                    OrderedData[[Dataset]] <- OrderedGroups
                }
                
                add_message(paste0("Sample processing complete: ", samples_processed, " processed, ", samples_skipped, " skipped"), 
                           if(samples_skipped > 0) "warning" else "success")
                
                incProgress(0.8, detail = "Combining data...")
                add_message("Combining processed data...", "info")
                
                CombinedDataList <- list()
                for (i in seq_along(data_list)){
                    var_name <- paste0("Exp",i)
                    # Filter out NULL or NA list elements before joining
                    if (!is.null(FinalData[[i]])) {
                        valid_final_data <- FinalData[[i]][!sapply(FinalData[[i]], function(x) is.null(x) || all(is.na(x)))]
                        if(length(valid_final_data) > 0) {
                            CombinedDataList[[var_name]] <- Reduce(function(x,y) full_join(x,y), valid_final_data)
                        }
                    }
                }
                
                if (length(CombinedDataList) == 0) {
                    add_message("No valid data to combine - all samples were skipped", "error")
                    return(NULL)
                }
                
                CombinedData <- do.call(rbind, CombinedDataList)
                
                #r2,2,5   
                temp <- CombinedData
                if(input$baseline_mode == "none"){ CombinedData<-temp}

                processed_data <- reactive({
                  df <- CombinedData
                  
                  if (input$clams_y_var == "si13c" && input$baseline_mode == "adjust") {
                    df <- df %>%
                      group_by(Samples) %>%
                      mutate(
                        baseline_mean = mean(si13c[Minutes >= input$baseline_st & Minutes <= 0], na.rm = TRUE),
                        si13c = si13c - baseline_mean
                      ) %>%
                      ungroup()
                  }
                  
                  df
                })
                
                CombinedData <- processed_data()
                
                
                
                #cutoff
                apply_cutoff <- function(x, method, time){
                  
                  if(method == "none") return(x)
                  if(is.null(time)) return(x)
                  
                  # v4.2              capping the data, remove the outlier
                  for (i in seq_along(x)) {
                    Z <- x[[i]]
                    
                    # Hard cutoff
                    Z[Z < input$lower_cutoff] <- input$lower_cutoff
                    Z[Z > input$upper_cutoff] <- input$upper_cutoff
                    
                    x[[i]] <- Z
                  }
                  
                  return(x)
                }
                
                
                if(input$boundary_cutoff == "cutoff"){
                  
                  metabolic_vars <- c("si13c")
                  
                  #   "vo2","vco2","vh2o","kcal_hr","rq",
                  #  "foodupa","waterupa","waterupb",
                  # "bodymass","xbreak","ybreak","zbreak",
                  #  "pedmeters","allmeters","bodytemp"
                  # )
                  
                  for(var in metabolic_vars){
                    
                    if(var %in% names(CombinedData)){
                      
                      CombinedData[[var]] <-
                        apply_cutoff(
                          CombinedData[[var]],
                          input$boundary_cutoff,
                          time = CombinedData$Minutes
                        )                  
                    }
                  }
                }               
                
                
                
                
                #r3, by DS, 3/10/2026
                #unified filter function
                apply_filter <- function(x, method, time,
                                         t_start_pre = NULL,
                                         t_split = 0,
                                         k = 5,
                                         t0 = 3
                                         
                ){
                  
                  if(method == "none") return(x)
                  if(is.null(time)) return(x)
                  
                  y <- x
                  
                  # --- define regions ---
                  idx_pre  <- !is.null(t_start_pre) & time >= t_start_pre & time <= t_split
                  idx_post <- time >= t_split
                  
                  # --- helper ---
                  filter_segment <- function(x_sub){
                    if(method == "hampel"){
                      if(all(is.na(x_sub))) return(x_sub)
                      
                      x_filled <- zoo::na.approx(x_sub, na.rm = FALSE)
                      x_filled <- zoo::na.locf(x_filled, na.rm = FALSE)
                      x_filled <- zoo::na.locf(x_filled, fromLast = TRUE, na.rm = FALSE)
                      
                      return(pracma::hampel(x_filled, k = k,t0=t0)$y)
                    }
                    
                    return(x_sub)
                  }
                  
                  # --- apply baseline filter ---
                  if(any(idx_pre)){
                    y[idx_pre] <- filter_segment(x[idx_pre])
                  }
                  
                  # --- apply experiment filter ---
                  if(any(idx_post)){
                    y[idx_post] <- filter_segment(x[idx_post])
                  }
                  
                  return(y)
                }
                
                
                if(input$signal_filter != "none"){
                  
                  metabolic_vars <- c("si13c")
                  
                  #   "vo2","vco2","vh2o","kcal_hr","rq",
                  #  "foodupa","waterupa","waterupb",
                  # "bodymass","xbreak","ybreak","zbreak",
                  #  "pedmeters","allmeters","bodytemp"
                  # )
                  
                  for(var in metabolic_vars){
                    
                    if(var %in% names(CombinedData)){
                      
                      CombinedData[[var]] <-
                        apply_filter(
                          CombinedData[[var]],
                          input$signal_filter,
                          time = CombinedData$Minutes,
                          t_start_pre = input$t_start_pre,   # where baseline filtering starts
                          t_split = input$t_split,         # split point
                          k = input$hampel_k,
                          t0 = input$hampel_t0
                        )                  
                    }
                    
                  }
                  
                }
                
                
                  if (input$bodymass_source == "zeromass") {
                    CombinedData$bodymass <- 0
                  }
                
                # r2.2.7
                # === Replace bodymass from metadata (if requested) ===
                if (!is.null(mouse_metadata)) {

                  
                  
                  # Try to detect the mass columns with a few common aliases
                  md <- mouse_metadata
                  lean_col  <- find_col(md, c("BodyMass_Lean","Lean_Mass","Lean","Lean_g","LeanMass"))
                  fat_col   <- find_col(md, c("BodyMass_Fat","Fat_Mass","Fat","Fat_g","FatMass"))
                  total_col <- find_col(md, c("BodyMass_Total","Total_Mass","Total","Body_Mass","BodyMass","Total_g"))
                  
                  # Require at least one of them to proceed
                  if (input$bodymass_source != "original") {
                    # Key column in your metadata that matches CombinedData$Samples is "sample_name"
                    sample_key <- find_col(md, c("sample_name", "Samples", "Sample", "sample"))
                    if (is.null(sample_key)) {
                      add_message("Body mass replacement skipped: could not find 'sample_name' in metadata", "warning")
                    } else {
                      # Build a minimal mapping frame
                      mass_map <- md %>%
                        dplyr::select(
                          !!sample_key,
                          !!lean_col  := dplyr::all_of(if (!is.null(lean_col))  lean_col  else sample_key),
                          !!fat_col   := dplyr::all_of(if (!is.null(fat_col))   fat_col   else sample_key),
                          !!total_col := dplyr::all_of(if (!is.null(total_col)) total_col else sample_key)
                        ) %>%
                        dplyr::distinct()
                      
                      # Coerce to numeric where present
                      for (nm in c(lean_col, fat_col, total_col)) {
                        if (!is.null(nm) && nm %in% names(mass_map)) {
                          mass_map[[nm]] <- suppressWarnings(as.numeric(mass_map[[nm]]))
                        }
                      }
                      
                      # Join per-sample onto CombinedData
                      CombinedData <- CombinedData %>%
                        dplyr::left_join(mass_map, by = c("Samples" = sample_key))
                      
                      # Ensure a bodymass column exists (if missing, create it)
                      if (!("bodymass" %in% names(CombinedData))) {
                        CombinedData$bodymass <- NA_real_
                      }
                      
                      # Choose which column to use
                      chosen_col <- switch(
                        input$bodymass_source,
                        "lean"  = lean_col,
                        "fat"   = fat_col,
                        "total" = total_col,
                        NULL
                      )
                      
                      if (!is.null(chosen_col) && chosen_col %in% names(CombinedData)) {
                        # Replace bodymass
                        n_before <- sum(!is.na(CombinedData$bodymass))
                        CombinedData$bodymass <- CombinedData[[chosen_col]]
                        n_after  <- sum(!is.na(CombinedData$bodymass))
                        
                        add_message(
                          paste0(
                            "Replaced 'bodymass' from metadata (", input$bodymass_source, "). ",
                            "Non-NA entries: ", n_before, " → ", n_after
                          ),
                          if (n_after > 0) "success" else "warning"
                        )
                      } else {
                        add_message(
                          paste0(
                            "Body mass replacement requested (", input$bodymass_source,
                            ") but corresponding column not found in metadata."
                          ),
                          "warning"
                        )
                      }
                    }
                  } # end if not original
                } # end if metadata available
                
                
                
                
                
                
                # Detect available CLAMS metrics in the combined data
                all_metric_names <- c("si13c", "vo2", "vco2", "rq", "kcal_hr", "vh2o", 
                                     "foodupa", "waterupa", "waterupb", "bodymass", 
                                     "xbreak", "ybreak", "zbreak", "pedmeters", "allmeters", "bodytemp")
                available_metrics <- all_metric_names[all_metric_names %in% names(CombinedData)]
                
                add_message(paste0("Available CLAMS metrics in combined data: ", paste(available_metrics, collapse=", ")), "info")
                
                if (length(available_metrics) == 0) {
                    add_message("WARNING: No standard CLAMS metric columns found in data!", "warning")
                    add_message(paste0("Available columns: ", paste(names(CombinedData), collapse=", ")), "info")
                }
                
                # Generate plots
                incProgress(0.9, detail = "Generating plots...")
                add_message("Generating preview plots...", "info")
                
                plot_list <- list()
                
                # Get the selected Y-axis variable with smart defaulting
                y_var <- if (!is.null(input$clams_y_var) && input$clams_y_var %in% names(CombinedData)) {
                    input$clams_y_var
                } else if ("si13c" %in% names(CombinedData)) {
                    "si13c"
                } else if (length(available_metrics) > 0) {
                    # Use first available metric
                    available_metrics[1]
                } else {
                    "si13c"  # Fallback default
                }
                
                add_message(paste0("Using Y-axis variable: ", y_var), "info")
                
                # Define Y-axis labels for each metric
                y_labels <- list(
                    si13c = expression(paste(text=''^13,'C in exhaled CO2 (per mil)')),
                    vo2 = "VO₂ - Oxygen Consumption (ml/kg/hr)",
                    vco2 = "VCO₂ - CO₂ Production (ml/kg/hr)",
                    rq = "Respiratory Quotient (RQ)",
                    kcal_hr = "Energy Expenditure (kcal/hr)",
                    vh2o = "Water Vapor - VH₂O (ml/kg/hr)",
                    foodupa = "Food Intake (grams)",
                    waterupa = "Water Intake A (ml)",
                    waterupb = "Water Intake B (ml)",
                    bodymass = "Body Mass (grams)",
                    xbreak = "X-axis Activity (beam breaks)",
                    ybreak = "Y-axis Activity (beam breaks)",
                    zbreak = "Z-axis Activity (beam breaks)",
                    pedmeters = "Pedometer Counts",
                    allmeters = "Total Activity Counts",
                    bodytemp = "Body Temperature (°C)"
                )
                
                if (nrow(CombinedData) > 0 && y_var %in% colnames(CombinedData)) {
                    
                    # Create base plots for both sample and condition coloring
                    create_base_plot <- function(color_var, title_suffix) {
                        p <- ggplot(data = CombinedData, aes_string(x = "Minutes", y = y_var, group = "Samples", color = color_var)) +
                            geom_line(linetype = "solid", linewidth = 1) + 
                            ylab(label = y_labels[[y_var]]) +
                            xlab(label = "Time (minutes)") +
                            ggpubr::theme_pubr() +
                            geom_vline(xintercept = 0, linetype = 2) 
                    
 
                        
                        # Add faceting if enabled
                        if (input$use_faceting && input$facet_by %in% colnames(CombinedData)) {
                            facet_cols <- if (!is.null(input$facet_cols) && input$facet_cols > 0) input$facet_cols else 2
                            p <- p + facet_wrap(as.formula(paste("~", input$facet_by)), ncol = facet_cols) +
                                 theme(legend.position = "bottom", 
                                       strip.text = element_text(size = 10, face = "bold"))
                            title <- paste("Faceted by", input$facet_by, "-", title_suffix)
                        } else {
                            p <- p + theme(legend.position = "right")
                            title <- title_suffix
                        }
                        
                        p <- p + ggtitle(label = title)
                        
                        # Add appropriate color scale
                        if (color_var == "Samples") {
                            p <- p + scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 8, name = "Accent"),
                                                                 RColorBrewer::brewer.pal(n = 12, name = "Set3")))
                        } else {
                            p <- p + scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                                                                 RColorBrewer::brewer.pal(n = 9, name = "Pastel1")))
                        }
                        
                        return(p)
                    }
                    
                    # Generate plots
                    plot_list$p1 <- create_base_plot("Samples", "Samples Colored by Samples")
                    plot_list$p2 <- create_base_plot("Conditions", "Samples Colored by Conditions")
                    
                    add_message(paste0("Generated ", length(plot_list), " plots for CLAMS data using '", y_var, "' as Y-axis"), "success")
                    
                } else {
                    add_message(paste0("Cannot generate plots - data missing or '", y_var, "' column not found in data"), "warning")
                    add_message(paste0("Available metric columns: ", paste(available_metrics, collapse = ", ")), "info")
                    if (length(available_metrics) > 0) {
                        add_message(paste0("Suggestion: Try using '", available_metrics[1], "' instead"), "info")
                    }
                    add_message(paste0("All columns in data: ", paste(names(CombinedData), collapse = ", ")), "info")
                }
                
                add_message(paste0("CLAMS analysis completed successfully! Processed ", nrow(CombinedData), " data points"), "success")
                
                # Return all the generated data and plots
                return(list(
                    CombinedData = CombinedData,
                    plots = plot_list,
                    individual_samples = if(nrow(CombinedData) > 0) CombinedData %>% group_by(Samples) %>% group_split() else list()
                ))
            })
            
        } else {
            # DSI Analysis Logic - Using sophisticated processing from original DSI app
            req(input$excel_file)
            
            withProgress(message = 'Processing DSI data...', value = 0, {
                
                incProgress(0.1, detail = "Reading DSI Excel file...")
                add_message("Reading DSI Excel file...", "info")
                
                # Process DSI data using sophisticated function from original DSI app
                add_message("Calling DSI processing function...", "info")
                
                # Get parameters with debugging - use DSI-specific input IDs
                exp_duration <- if (!is.null(input$dsi_duration_of_experiment)) as.numeric(input$dsi_duration_of_experiment) else 150
                base_duration <- if (!is.null(input$dsi_duration_before_experiment)) as.numeric(input$dsi_duration_before_experiment) else 60
                add_message(paste0("DSI Parameters: exp_duration=", exp_duration, ", base_duration=", base_duration), "info")
                
                processed_data <- tryCatch({
                    # Check file exists
                    if (!file.exists(input$excel_file$datapath)) {
                        stop("Excel file path does not exist")
                    }
                    add_message(paste0("Processing Excel file: ", input$excel_file$datapath), "info")
                    
                    result <- process_dsi_data_unified(
                        file_path = input$excel_file$datapath,
                        metadata = if (input$use_metadata && !is.null(input$metadata_file)) {
                            tryCatch({
                                meta_data <- read.csv(input$metadata_file$datapath, stringsAsFactors = FALSE)
                                add_message(paste0("Metadata columns: ", paste(names(meta_data), collapse = ", ")), "info")
                                add_message(paste0("Expected columns: mouse_id, start_date, start_time"), "info")
                                meta_data
                            }, error = function(e) {
                                add_message(paste("Error reading metadata:", e$message), "error")
                                NULL
                            })
                        } else {
                            add_message("No metadata provided, using basic processing", "info")
                            NULL
                        },
                        experiment_duration = exp_duration,
                        baseline_duration = base_duration,
                        time_anchor = "experiment"
                    )
                    add_message(paste0("DSI processing function returned data with ", nrow(result), " rows"), "info")
                    result
                }, error = function(e) {
                    add_message(paste("Error processing DSI data:", e$message), "error")
                    add_message(paste("Stack trace:", paste(capture.output(traceback()), collapse = "\n")), "error")
                    return(data.frame())
                })
                
                if (is.null(processed_data) || nrow(processed_data) == 0) {
                    add_message("DSI processing failed or returned no data", "error")
                    return(NULL)
                }
                
                incProgress(0.9, detail = "Finalizing results...")
                
                # Count mouse-specific columns for validation
                per_mouse_cols <- grep("^Mouse_\\d+_", names(processed_data), value = TRUE)
                if (length(per_mouse_cols) > 0) {
                    base_metrics <- sort(unique(sub("^Mouse_\\d+_", "", per_mouse_cols)))
                    unique_mice <- sort(unique(sub("_.*", "", per_mouse_cols)))
                    sample_metrics <- head(base_metrics, 5)
                    add_message(paste0("DSI analysis completed successfully! Found ", length(unique_mice), " mice, ", length(base_metrics), " metrics"), "success")
                    add_message(paste0("Sample metrics: ", paste(sample_metrics, collapse = ", ")), "info")
                    add_message(paste0("Detected mice: ", paste(unique_mice, collapse = ", ")), "info")
                    add_message(paste0("Data dimensions: ", nrow(processed_data), " rows × ", ncol(processed_data), " columns"), "info")
                    if ("RelativeTime_Minutes" %in% names(processed_data)) {
                        add_message(paste0("RelativeTime_Minutes range: ", min(processed_data$RelativeTime_Minutes, na.rm = TRUE), " to ", max(processed_data$RelativeTime_Minutes, na.rm = TRUE)), "info")
                    }
                } else {
                    add_message("Warning: No mouse-specific columns were created. Check Excel file structure and start row setting.", "warning")
                }
                
                # Return DSI-specific data structure
                return(list(
                    processed_data = processed_data,
                    metadata = if (input$use_metadata && !is.null(input$metadata_file)) {
                        tryCatch({
                          # read.csv(input$metadata_file$datapath, stringsAsFactors = FALSE)
                          # req(input$start_time_col)
                          # metadata$start_time <- metadata[[input$start_time_col]]
                          
                          req(metadata())
                          metadata()
                          
                        }, error = function(e) NULL)
                    } else NULL,
                    plots = list(),
                    individual_samples = list()
                ))
            })
        }
    })
    
    # ==== End of analysis_results
    
    
    # ====    

    # Display status message
    output$status <- renderText({
        data_type_name <- if(input$data_type == "clams") "CLAMS" else "DSI"
        
        # If app was reset, show initial prompts
        if (!analysis_state$has_run) {
            if (input$data_type == "clams") {
              # 10/8/2026
                # if (!analysis_state$clams_folder_selected) {
                #     return("📁 Please select a CLAMS data folder containing CSV files")
                  if (is.null(input$clams_files)) {
                    return("📁 Please upload one or more CLAMS CSV files")
                } else if (!analysis_state$metadata_loaded) {
                    return("📋 Please upload a CLAMS metadata CSV file")
                } else {
                    return("⚙️ Ready to run CLAMS analysis - click 'Run Analysis' button")
                }
            } else {
                if (!analysis_state$excel_loaded) {
                    return("📁 Please upload a DSI Excel file")
                } else if (input$use_metadata && !analysis_state$metadata_loaded) {
                    return("📋 Please upload a DSI metadata CSV file")
                } else {
                    return("⚙️ Ready to run DSI analysis - click 'Run Analysis' button")
                }
            }
        }
        
        # If analysis has run, show success message
        if (analysis_state$has_run && !is.null(analysis_results())) {
            return(paste0("✅ ", data_type_name, " Analysis complete!\n\nResults are ready for download."))
        }
        
        # Default fallback
        return("⚙️ Ready to run analysis - click 'Run Analysis' button")
    })
    
    # Update DSI Y-variable choices when DSI data is analyzed
    observe({
        if (input$data_type == "dsi") {
            results <- analysis_results()
            if (!is.null(results) && !is.null(results$processed_data)) {
                data <- results$processed_data
                # Find per-mouse columns (Mouse_X_MetricName pattern)
                per_mouse_cols <- grep("^Mouse_\\d+_", names(data), value = TRUE)
                if (length(per_mouse_cols) > 0) {
                    # Extract unique metric names
                    base_metrics <- sort(unique(sub("^Mouse_\\d+_", "", per_mouse_cols)))
                    
                    # Update the Y-variable choices
                    current_selection <- input$dsi_y_var
                    default_sel <- if (!is.null(current_selection) && current_selection %in% base_metrics) {
                        current_selection
                    } else if ("Glucose_Gavg_mg_dL" %in% base_metrics) {
                        "Glucose_Gavg_mg_dL"
                    } else {
                        base_metrics[1]
                    }
                    
                    updateSelectInput(session, "dsi_y_var", 
                                    choices = base_metrics, 
                                    selected = default_sel)
                    
                    # Log success
                    isolate({
                        add_message(paste0("Found ", length(base_metrics), " DSI metrics for plotting"), "info")
                    })
                } else {
                    updateSelectInput(session, "dsi_y_var", choices = c(), selected = NULL)
                    # Log warning with proper isolation
                    isolate({
                        add_message("No mouse-specific columns found in DSI data", "warning")
                    })
                }
            } else {
                updateSelectInput(session, "dsi_y_var", choices = c(), selected = NULL)
            }
        }
    })
    

    
    # ==== Start of Statistics Analyses ====
    # r2.3.2
  
    gam_results <- reactive({
      req(analysis_results())
      update_trigger()
      
      # req(input$analysis_model == "gam")
      
      results <- analysis_results()
      req(results)
      
      plot_data <- results$CombinedData
      req(plot_data, nrow(plot_data) > 0)
      
      y_var <- input$clams_y_var
      
      run_gam_test <- (
        y_var == "si13c" &&
          input$facet_by == "Conditions"
      )
      
      if (!run_gam_test) {
        return(NULL)
      }
      
      # Ensure variable exists
      req(y_var %in% colnames(plot_data))
      
      # Build model dataframe
      df_model <- plot_data[, c(
        "Minutes",
        "Samples",
        "Conditions",
        "bodymass",
        y_var
      )]
      
      names(df_model)[names(df_model) == y_var] <- "y"
      
      df_model <- df_model[!is.na(df_model$y), ]
      
      df_model$Conditions <- as.factor(df_model$Conditions)
      df_model$Samples <- as.factor(df_model$Samples)
      
      # Optional filter
      df_model <- df_model[df_model$Minutes >= 0, ]
      
      # ===== Main plotting model =====
      gam_model <- mgcv::gam(
        y ~ Conditions +
          bodymass +
          s(Minutes,
            by = Conditions,
            bs = "fs",
            k = input$gam_k) +
          s(Samples, bs = "re"),
        data = df_model,
        method = "REML"
      )
      
      # ===== Full model =====
      gam_full <- mgcv::gam(
        y ~ Conditions +
          bodymass +
          s(Minutes,
            by = Conditions,
            bs = "fs",
            k = input$gam_k),
        data = df_model,
        method = "REML"
      )
      
      # ===== Null model =====
      gam_null <- mgcv::gam(
        y ~ Conditions +
          bodymass +
          s(Minutes,
            bs = "fs",
            k = input$gam_k),
        data = df_model,
        method = "REML"
      )
      
      # ===== Compare models =====
      comp <- anova(gam_null, gam_full, test = "Chisq")
      
      gam_p_value <- comp$`Pr(>Chi)`[2]
      
      # ===== Pairwise comparisons =====
      conds <- unique(df_model$Conditions)
      pairs <- combn(conds, 2, simplify = FALSE)
      
      pairwise_results <- lapply(pairs, function(p) {
        
        sub_df <- df_model[df_model$Conditions %in% p, ]
        
        sub_df$Conditions <- droplevels(sub_df$Conditions)
        
        m_full <- mgcv::gam(
          y ~ Conditions +
            bodymass +
            s(Minutes,
              by = Conditions,
              bs = "fs",
              k = input$gam_k),
          data = sub_df,
          method = "REML"
        )
        
        m_reduced <- mgcv::gam(
          y ~ Conditions +
            bodymass +
            s(Minutes,
              bs = "fs",
              k = input$gam_k),
          data = sub_df,
          method = "REML"
        )
        
        test <- anova(m_reduced, m_full, test = "Chisq")
        
        data.frame(
          Comparison = paste(p, collapse = " vs "),
          p_value = test$`Pr(>Chi)`[2]
        )
      })
      
      pairwise_df <- do.call(rbind, pairwise_results)
      
      pairwise_df$label <- paste0(
        pairwise_df$Comparison,
        ": p = ",
        signif(pairwise_df$p_value, 3)
      )
      
      list(
        df_model = df_model,
        gam_model = gam_model,
        gam_full = gam_full,
        gam_null = gam_null,
        gam_p_value = gam_p_value,
        pairwise_pvals = pairwise_df
      )
    })
    
    # Create pairwise p-value reactive
    gam_pairwise_pvals <- reactive({
      
      req(input$analysis_model == "gam")
      req(input$clams_y_var == "si13c")
      
      df <- gam_results()$df_model
      req(!is.null(df), nrow(df) > 0)
      
      conds <- unique(df$Conditions)
      pairs <- combn(conds, 2, simplify = FALSE)
      
      results <- lapply(pairs, function(p) {
        
        sub_df <- df[df$Conditions %in% p, ]
        
        # Drop unused levels
        sub_df$Conditions <- droplevels(sub_df$Conditions)
        
        # Full model (separate smooths)
        m_full <- mgcv::gam(
          y ~ Conditions + bodymass + s(Minutes, by = Conditions, bs = "fs", k = input$gam_k),# + s(Samples, bs = "re"),
          # y ~ Conditions + s(Minutes, by = Conditions),
          data = sub_df,
          method = "REML"
        )
        
        # Reduced model (shared smooth)
        m_reduced <- mgcv::gam(
          y ~ Conditions + bodymass + s(Minutes, bs = "fs", k = input$gam_k),# + s(Samples, bs = "re"),
          # y ~ Conditions + s(Minutes),
          data = sub_df,
          method = "REML"
        )
        
        # LRT
        test <- anova(m_reduced, m_full, test = "Chisq")
        
        data.frame(
          Comparison = paste(p, collapse = " vs "),
          p_value = test$`Pr(>Chi)`[2]
        )
      })
      
      do.call(rbind, results)
    })
    
    # Format for display
    formatted_pvals <- reactive({
      
      df <- gam_pairwise_pvals()
      req(!is.null(df), nrow(df) > 0)
      
      df$label <- paste0(
        df$Comparison,
        ": p = ",
        signif(df$p_value, 3)
      )
      
      df
    })
    
    # ===== End of GAM MODEL ====
    
    
    # ==== New two-period One data reactive ====
    
    two_period_results <- reactive({
      req(analysis_results())
      update_trigger()
      
      req(analysis_results())
      
      results <- analysis_results()
      
      plot_data <- results$CombinedData
      req(plot_data, nrow(plot_data) > 0)
      
      y_var <- input$clams_y_var
      
      df <- plot_data[, c("Minutes", "Samples", "Conditions", y_var, "bodymass")]
      
      # ----- RECALCULATE BASELINE (only for si13c) r2.3.11-----
      if (y_var == "si13c" && input$start_time_col=="start_time") {
        df <- df %>%
          dplyr::group_by(Samples) %>%
          dplyr::mutate(
            baseline_mean = mean(
              .data[[y_var]][Minutes >= input$baseline_st & Minutes <= 0],
              na.rm = TRUE
            )
          ) %>%
          dplyr::ungroup()
      } else {
        df$baseline_mean <- NA_real_
      } 
    
      
      t0  <- input$cut_time
      w   <- input$time_window
      wb4 <- input$time_windowb4
      
      # ----- FILTER WINDOW -----
      df2 <- dplyr::case_when(
        
        input$two_periods == "before" ~
          df$Minutes >= (t0 - wb4) & df$Minutes <= t0,
        
        input$two_periods == "after" ~
          df$Minutes > t0 & df$Minutes <= (t0 + w),
        
        input$two_periods %in% c("both","bothc") ~
          (df$Minutes >= (t0 - wb4) & df$Minutes <= t0) |
          (df$Minutes > t0 & df$Minutes <= (t0 + w)),
        
        TRUE ~ FALSE
      )
      
      df2 <- df[df2, ]
      
      df2$Period <- ifelse(df2$Minutes <= t0, "Before", "After")
      
      # ----- SUMMARY -----
      summary_df <- df2 %>%
        dplyr::group_by(Samples, Conditions, Period) %>%
        dplyr::summarise(
          mean_y = mean(.data[[y_var]], na.rm = TRUE),
          bodymass = mean(bodymass, na.rm = TRUE),
          baseline_mean = mean(baseline_mean, na.rm = TRUE),
          .groups = "drop"
        )
      
      
      # ----- WIDE FORMAT -----
      wide_df <- summary_df %>%
        tidyr::pivot_wider(names_from = Period, values_from = mean_y)
      
      
      # ----- FILTER VALID CASES -----
      if (input$two_periods == "before") {
        wide_df <- wide_df %>% dplyr::filter(!is.na(Before))
      } else if (input$two_periods == "after") {
        wide_df <- wide_df %>% dplyr::filter(!is.na(After))
      } else {
        wide_df <- wide_df %>% dplyr::filter(!is.na(Before) & !is.na(After))
      }
      
      wide_df
    })
    
    # one model reactive version
    # Two-period + ANCOVA (core upgrade)
    ancova_model <- NULL
    ancova_p_value <- NA
    interaction_p <- NA
    
    ancova_model <- reactive({
      req(analysis_results())
      update_trigger()
      
      df <- two_period_results()
      req(df, nrow(df) > 0)
      
      if (input$analysis_model == "twoperiod") {
        
        if (input$two_periods == "before") {
          lm(Before ~ Conditions + bodymass, data = df)
          
        } else if (input$two_periods == "after") {
          lm(After ~ Conditions + bodymass, data = df)
          
        } else if (input$two_periods == "both") {
          lm(After ~ Before + Conditions + bodymass, data = df)
          
        } else if (input$two_periods == "bothc") {
          lm(After - Before ~ Conditions + bodymass, data = df)
        }
        
      } else if (input$analysis_model == "twoperiodglm") {
        
        if (input$two_periods == "before") {
          lm(Before ~ Conditions + bodymass + Conditions:bodymass, data = df)
          
        } else if (input$two_periods == "after") {
          lm(After ~ Conditions + bodymass + Conditions:bodymass, data = df)
          
        } else if (input$two_periods == "both") {
          lm(After ~ Before + Conditions + bodymass + Conditions:bodymass, data = df)
          
        } else if (input$two_periods == "bothc") {
          lm(After - Before ~ Conditions + bodymass + Conditions:bodymass, data = df)
        }
      }
    })
    
    ancova_p_value <- reactive({
      fit <- ancova_model()
      req(fit)
      
      if (input$analysis_model == "twoperiod") {
        anova(fit)["Conditions", "Pr(>F)"]
      } else {
        car::Anova(fit, type = "II")["Conditions", "Pr(>F)"]
      }
    })
    
    bodymass_p_value <- reactive({
      fit <- ancova_model()
      req(fit)
      
      coefs <- coef(summary(fit))
      
      # Exclude interaction terms
      rows <- grep("bodymass", rownames(coefs))
      rows <- rows[!grepl(":", rownames(coefs)[rows])]
      
      if (length(rows) == 0) return(NA)
      
      coefs[rows, "Pr(>|t|)"][1]
    })

    interaction_p <- reactive({
      fit <- ancova_model()
      coefs <- coef(summary(fit))
      
      rows <- grep("Conditions.*:bodymass", rownames(coefs))
      if (length(rows) == 0) return(NULL)
      
      coefs[rows, "Pr(>|t|)", drop = TRUE]
    })
    
    # Create pairwise comparison reactive
    pairwise_results <- reactive({
      fit <- ancova_model()
      req(fit)
      
      emm <- emmeans(fit, ~ Conditions)
      pw <- pairs(emm, adjust = "tukey")
      
      as.data.frame(pw)
    })
    
    # Clean and format results
    pairwise_table <- reactive({
      req(analysis_results())
      update_trigger()
      
      df <- pairwise_results()
      req(df)
      
      df %>%
        dplyr::mutate(
          Comparison = contrast,
          Estimate = round(estimate, 2),
          p_value = signif(p.value, 3)
        ) %>%
        dplyr::select(Comparison, Estimate, p_value)
    })
    
    # Convert to annotation text
    pairwise_text <- reactive({
      df <- pairwise_table()
      req(df)
      
      paste(
        apply(df, 1, function(row) {
          paste0(row["Comparison"], ": ",
                 row["Estimate"], ", p=", row["p_value"])
        }),
        collapse = "\n"
      )
    })
    
    # --- Shared summary ---
    period_summary <- reactive({
      df <- two_period_results()
      req(df)
      
      # --- Decide which columns to use ---
      period_cols <- switch(
        input$two_periods,
        "before" = "Before",
        "after"  = "After",
        c("Before", "After")   # for both / bothc
      )
      
      df %>%
        tidyr::pivot_longer(
          cols = all_of(period_cols),
          names_to = "Period",
          values_to = "value"
        ) %>%
        dplyr::group_by(Conditions, Period) %>%
        dplyr::summarise(
          mean = mean(value, na.rm = TRUE),
          sd   = sd(value, na.rm = TRUE),
          n    = dplyr::n(),
          se   = sd / sqrt(n),
          ci   = 1.96 * se,
          .groups = "drop"
        )
    })
    
    # --- Window logic ---
    window_df <- reactive({
      t0  <- input$cut_time
      w   <- input$time_window
      wb4 <- input$time_windowb4
      
      period_summary() %>%
        dplyr::mutate(
          xmin = dplyr::case_when(
            Period == "Before" ~ t0 - wb4,
            Period == "After"  ~ t0,
            TRUE ~ t0
          ),
          xmax = dplyr::case_when(
            input$two_periods == "before" ~ t0,
            input$two_periods == "after"  ~ t0 + w,
            Period == "Before" ~ t0,
            Period == "After"  ~ t0 + w,
            TRUE ~ t0
          )
        )
    })
    
    
    # ==== End of Statistics Analyses
    
    
    # This replaces your text-based output with a structured table of results.
    stats_table <- reactive({
      req(analysis_results())
      update_trigger()
      
      model_type <- input$analysis_model
      
      # --- Condition effect ---
      cond_p <- tryCatch(ancova_p_value(), error = function(e) NA)
      
      cond_row <- data.frame(
        Effect = "Condition",
        P_value = cond_p
      )
      
      # --- Bodymass effect ---
      bm_p <- tryCatch(bodymass_p_value(), error = function(e) NA)
      
      bm_row <- data.frame(
        Effect = "Bodymass",
        P_value = bm_p
      )
      
      # --- Interaction effect ---
      int_p <- tryCatch(interaction_p(), error = function(e) NULL)
      
      if (!is.null(int_p)) {
        int_df <- data.frame(
          Effect = paste0("Interaction (", names(int_p), ")"),
          P_value = as.numeric(int_p)
        )
      } else {
        int_df <- NULL
      }
      
      # Combine all
      df <- dplyr::bind_rows(cond_row, bm_row, int_df)
      
      # Format p-values
      df <- df %>%
        dplyr::mutate(
          P_value = ifelse(
            is.na(P_value),
            NA,
            signif(P_value, 3)
          ),
          Significance = dplyr::case_when(
            is.na(P_value) ~ "",
            P_value < 0.001 ~ "***",
            P_value < 0.01  ~ "**",
            P_value < 0.05  ~ "*",
            TRUE ~ ""
          )
        )
      
      df
    })
    
    # Put p-values directly on your plots:
    annot_text <- reactive({
      
      df <- stats_table()
      
      paste(
        apply(df, 1, function(r) {
          paste0(r["Effect"], ": p=", r["P_value"], " ", r["Significance"])
        }),
        collapse = "\n"
      )
    })
    
    # Create a reactive for model fit
    model_fit_info <- reactive({
      m <- ancova_model()
      req(m)
      
      # For lm models (all your cases are lm)
      s <- summary(m)
      
      out <- list(
        R2 = s$r.squared,
        Adj_R2 = s$adj.r.squared,
        AIC = AIC(m),
        BIC = BIC(m)
      )
      
      # Add residual SE if available
      if (!is.null(s$sigma)) {
        out$Residual_SE <- s$sigma
      }
      
      return(out)
    })
    
    # Format the model fit text
    model_fit_text <- reactive({
      fit <- model_fit_info()
      req(fit)
      
      paste(
        paste0("R² = ", signif(fit$R2, 3)),
        paste0("Adj R² = ", signif(fit$Adj_R2, 3)),
        paste0("AIC = ", signif(fit$AIC, 4)),
        paste0("BIC = ", signif(fit$BIC, 4)),
        if (!is.null(fit$Residual_SE)) {
          paste0("Residual SE = ", signif(fit$Residual_SE, 3))
        },
        sep = "\n"
      )
    })
    
    
    
    # disply statistics results 
    output$stats_text <- renderText({
      req(analysis_results())
      update_trigger()
      
      y_var <- input$clams_y_var   #r2.2.4, added
      
      # Check if app was reset
      if (!analysis_state$has_run) {
        return(NULL)
      }
      
      results <- analysis_results()
      if (is.null(results)) return(NULL)
      
      # Always generate plot dynamically to reflect faceting changes
      plot_data <- results$CombinedData
      req(plot_data, nrow(plot_data) > 0)
      
      #r2.2.4
      
      # ===== GAM MODEL COMPARISON (FILTERED si13c ONLY) =====
      gam_model <- NULL
      gam_p_value <- NA
      
      y_var <- input$clams_y_var
      
      run_gam_test <- (
        y_var == "si13c" && input$facet_by=="Conditions" &&
          input$analysis_model == "gam"
      )
      
      
      if (run_gam_test) {
        gam_res <- gam_results()
        if (!is.null(gam_res)) {
          df_model<-gam_res$df_model
          gam_model <- gam_res$gam_model
          gam_p_value <- gam_res$gam_p_value
          pairwise_pvals = gam_res$pairwise_pvals
        }
        
        label_text <- if (is.na(gam_p_value)) {
          "Curve test: NA"
        } else if (gam_p_value < 1e-4) {
          "Curve diff: p < 1e-4"
        } else {
          paste0("Curve diff: p = ", signif(as.numeric(gam_p_value), 3))
        }
        
        pvals <- formatted_pvals()
        # gam_annotation_text <- paste(pvals$label, collapse = "\n")
        
        
        p_val <- gam_p_value
        pairwise_text <- pvals$label
        
        print(pvals)
        print(pvals$label)
        
        
        # paste(
        #   "Statistics Results:",
        #   paste0("Condition Effect: p = ", signif(gam_p_value, 3)),
        #   "Pairwise Comparison:\n",
        #   pvals$label,
        #   
        #   sep = "\n"
        #   
        # )
      }
      
      
      # ✅ Add ANCOVA p-value
      p_val <- ancova_p_value()  # <-- call the reactive

      if (input$two_periods %in% c("both", "bothc", "before", "after") && input$analysis_model == "twoperiodglm") {
        int_p_val <- interaction_p()   #interaction p-val
      } else {int_p_val <- NA}

      label_text <- if (is.null(p_val) || !is.finite(p_val) && input$analysis_model == "twoperiod" ) {
        "ANCOVA p = NA"
      } else {
        paste0("ANCOVA p = ", signif(as.numeric(p_val), 3))
      }


      # Conditionally add interaction p-value if it exists

      if (!is.null(int_p_val) && length(int_p_val) > 0 && all(is.finite(int_p_val)) &&  input$analysis_model == "twoperiodglm" ) {

        if (length(int_p_val) == 1) {
          label_text <- paste(
            paste0("Condition:Bodaymass", ": p = ", signif(int_p_val, 3))
          )
        } else if (length(int_p_val) > 1) {

          label_text <- paste(
            paste0(names(int_p_val), ": p = ", signif(int_p_val, 3)),
            collapse = "\n"
          )

          label_text <- paste("Interaction:\n", label_text)
        }
        }
      
      
      bm_p <- bodymass_p_value()
 
      bm_text <- if (!is.null(bm_p) && is.finite(bm_p)) {
        paste0("Bodymass Effect: p = ", signif(bm_p, 3))
      } else {
        "Bodymass Effect: p = NA"
      }
      
      paste(
        # "Statistics Results:",

        paste0("Condition Effect: p = ", signif(p_val, 3)),
        bm_text,
        label_text, 
        paste0("Pairwise Comparison:\n", pairwise_text()),

        sep = "\n"

       )

      # Plug into your existing renderText
      paste(
        paste0("Condition Effect: p = ", signif(p_val, 3)),
        bm_text,
        label_text,
        
        # ✅ ADD THIS LINE
        "Model Fit:",
        model_fit_text(),
        
        paste0("Pairwise Comparison:\n", pairwise_text()),
        sep = "\n"
      )
      
    })
  
    # Render as a first table (DT)
    output$stats_table <- DT::renderDataTable({
      
      df <- stats_table()
      req(df)
      
      DT::datatable(
        df,
        rownames = FALSE,
        options = list(
          dom = "t",
          paging = FALSE
        )
      ) %>%
        DT::formatStyle(
          "P_value",
          color = DT::styleInterval(
            c(0.001, 0.01, 0.05),
            c("darkred", "red", "orange", "black")
          ),
          fontWeight = DT::styleInterval(
            0.05,
            c("bold", "normal")
          )
        )
    })
    
    # Add a second table instead of messy text:
    output$pairwise_table <- DT::renderDataTable({
      
      df <- pairwise_table()
      req(df)
      
      DT::datatable(
        df,
        rownames = FALSE,
        options = list(pageLength = 6)
      ) %>%
        DT::formatStyle(
          "p_value",
          color = DT::styleInterval(
            c(0.001, 0.01, 0.05),
            c("darkred", "red", "orange", "black")
          ),
          fontWeight = DT::styleInterval(
            0.05,
            c("bold", "normal")
          )
        )
    })
    
    # add second plots for regrestion on body mass
    # Plot preview2 - reactive to faceting changes
    ancova_plot_obj <- reactive({
      req(analysis_results())
      update_trigger()
      
      results <- analysis_results()
      # Always generate plot dynamically to reflect faceting changes
      plot_data <- results$CombinedData
      req(plot_data, nrow(plot_data) > 0)
      
    
      # Fallback plot
      ggplot() + 
        geom_text(aes(x = 0.5, y = 0.5, label = "No data available for plotting"), 
                  size = 6, hjust = 0.5, vjust = 0.5) + 
        theme_void()    
      
      
      # Step 1 — Decide the response variable (critical)
      y_var <- reactive({
        switch(
          input$two_periods,
          "before" = "Before",
          "after"  = "After",
          "both"   = "After",
          "bothc"  = "diff"   # we'll create this
        )
      })
      
      # Step 2 — Prepare data for plotting
      
      df <- two_period_results()
      req(df)
      
      # Handle bothc case
      if (input$two_periods == "bothc") {
        df$diff <- df$After - df$Before
      }
      
      df$Conditions <- factor(df$Conditions)
      
      
      # Step 3 — Create prediction grid
      
      model <- ancova_model()
      req(model)
      
      # sequence of bodymass
      bm_seq <- seq(
        min(df$bodymass, na.rm = TRUE),
        max(df$bodymass, na.rm = TRUE),
        length.out = 100
      )
      
      newdata <- expand.grid(
        bodymass = bm_seq,
        Conditions = levels(df$Conditions)
      )
      
      # Handle "both" model (needs Before)
      if (input$two_periods == "both") {
        newdata$Before <- mean(df$Before, na.rm = TRUE)
      }
      
      # Predictions
      newdata$fit <- predict(model, newdata)
      
      
      # Step 4 — Plot (this is what you want)
      
      y_col <- y_var()
      
      
      title_text <- if (input$bodymass_source == "zeromass") {
        "ANOVA fitted lines by condition"
      } else if (input$analysis_model == "twoperiodglm") {
        "GLM fitted lines by condition"
      } else if (input$analysis_model == "twoperiod") {
        "ANCOVA fitted lines by condition"
      } else {
        "Fitted lines by condition"
      }
      
      
      x_text <- if (input$bodymass_source == "zeromass") {
        "no Body Mass (g)"
      } else if (input$bodymass_source == "lean") {
        "Lean Body Mass (g)"
      } else if (input$bodymass_source == "fat") {
        "Fat Body Mass (g)"
      } else {
        "Total Body Mass (g)"
      }
      
      
      ggplot(df, aes(x = bodymass, y = .data[[y_col]], color = Conditions)) +
        
        # raw data
        geom_point(alpha = 0.5) +
        
        # fitted lines
        geom_line(
          data = newdata,
          aes(x = bodymass, y = fit, color = Conditions),
          linewidth = 1.5
        ) +
        
        labs(
          x = x_text,
          y = y_col,
          title = title_text
        ) +
        
        theme(
          axis.line = element_line(color = "black", linewidth = 1),
          axis.text = element_text(color = "black", size = 12),
          axis.title = element_text(color = "black", size = 12,face=NULL),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
        ) +
        
        annotate(
          "text",
          x = Inf, y = Inf,
          label = annot_text(),
          hjust = 1.1, vjust = 1.5,
          size = 4
        )
        # theme_minimal()
    
  })
    output$ancova_plot <- renderPlot({ ancova_plot_obj() })      
    
    # main Plot preview - reactive to faceting changes
    main_plot_obj <- reactive({
      req(analysis_results())
      update_trigger() 
      
        y_var <- input$clams_y_var   #r2.2.4, added
      
        # Check if app was reset
        if (!analysis_state$has_run) {
            return(NULL)
        }
        
        results <- analysis_results()
        if (is.null(results)) return(NULL)
        
        # Make plot reactive to UI inputs based on data type
        if (input$data_type == "clams") {
            # CLAMS-specific reactive inputs
            input$plot_color_by
            input$use_faceting
            input$facet_by
            input$facet_cols
            input$clams_y_var  # Make reactive to Y-axis selection
        } else {
            # DSI-specific reactive inputs
            input$dsi_y_var
            input$dsi_plot_group_by
            input$dsi_plot_facet
        }
        
        if (input$data_type == "clams") {
            # Define Y-axis labels for each metric
            y_labels <- list(
                si13c = expression(paste(text=''^13,'C in exhaled CO2 (per mil)')),
                vo2 = "VO₂ - Oxygen Consumption (ml/kg/hr)",
                vco2 = "VCO₂ - CO₂ Production (ml/kg/hr)",
                rq = "Respiratory Quotient (RQ)",
                kcal_hr = "Energy Expenditure (kcal/hr)",
                vh2o = "Water Vapor - VH₂O (ml/kg/hr)",
                foodupa = "Food Intake (grams)",
                waterupa = "Water Intake A (ml)",
                waterupb = "Water Intake B (ml)",
                bodymass = "Body Mass (grams)",
                xbreak = "X-axis Activity (beam breaks)",
                ybreak = "Y-axis Activity (beam breaks)",
                zbreak = "Z-axis Activity (beam breaks)",
                pedmeters = "Pedometer Counts",
                allmeters = "Total Activity Counts",
                bodytemp = "Body Temperature (°C)"
            )
            
            # Always generate plot dynamically to reflect faceting changes
            plot_data <- results$CombinedData
            req(plot_data, nrow(plot_data) > 0)
            
            #r2.2.4
           
            # ===== GAM MODEL COMPARISON (FILTERED si13c ONLY) =====
            gam_model <- NULL
            gam_p_value <- NA
            
            y_var <- input$clams_y_var
            
            run_gam_test <- (
              y_var == "si13c" && input$facet_by=="Conditions" &&
                input$analysis_model == "gam"
              #   !is.null(input$signal_filter)  &&
              #   tolower(input$signal_filter) != "none"
            )
            

            if (run_gam_test) {
              gam_res <- gam_results()
              if (!is.null(gam_res)) {
                df_model<-gam_res$df_model
                gam_model <- gam_res$gam_model
                gam_p_value <- gam_res$gam_p_value
                pairwise_pvals = gam_res$pairwise_pvals
              }
            }
            
            
            if (nrow(plot_data) > 0) {
              # Detect available metrics in the data
              metric_names <- c("vo2", "vco2", "vh2o", "kcal_hr", "rq", "foodupa", "waterupa", 
                                "waterupb", "bodymass", "xbreak", "ybreak", "zbreak", "pedmeters", 
                                "allmeters", "si13c", "bodytemp")
              available_metrics <- metric_names[metric_names %in% colnames(plot_data)]
              # Smart Y-variable selection: user choice -> si13c -> first available
              y_var <- if (!is.null(input$clams_y_var)) input$clams_y_var else "si13c"
              if (!y_var %in% available_metrics) {
                    if ("si13c" %in% available_metrics) {
                        y_var <- "si13c"
                        message("Preview: Selected metric not available, defaulting to si13c")
                    } else if (length(available_metrics) > 0) {
                        y_var <- available_metrics[1]
                        message(paste0("Preview: si13c not available, defaulting to ", y_var))
                    } else {
                        return(NULL)
                    }
                }
                
                if (y_var %in% colnames(plot_data)) {
                
                # Determine color variable based on user selection
                color_var <- if (!is.null(input$plot_color_by) && input$plot_color_by == "Samples") "Samples" else "Conditions"
                color_title <- if (color_var == "Samples") "Sample" else "Condition"
                

                cutoff_layers <- NULL
                
                if (input$boundary_cutoff == "cutoff"){
                     
                  cutoff_layers <- list(
                  ggplot2::geom_hline(yintercept = input$lower_cutoff, linetype = 3),
                  ggplot2::geom_hline(yintercept = input$upper_cutoff, linetype = 3),
                  
                  ggplot2::annotate(
                    "text",
                    x = Inf,
                    y = input$lower_cutoff,
                    label = paste("Lower Cutoff:", input$lower_cutoff),
                    hjust = 1.1,
                    vjust = 1.1,
                    size = 4
                  ),
                  
                  ggplot2::annotate(
                    "text",
                    x = Inf,
                    y = input$upper_cutoff,
                    label = paste("Upper Cutoff:", input$upper_cutoff),
                    hjust = 1.1,
                    vjust = 1.1,
                    size = 4
                  )
                  
                  )
                }
                
                
                
                # r2.2.1, DS, 3/23/2026
                
                hampel_layers <- NULL
                
                if (input$signal_filter == "hampel" && input$clams_y_var == "si13c") {
                  
                  hampel_layers <- list(
                    
                    # ggplot2::geom_hline(yintercept = input$lower_cutoff, linetype = 3),
                    # ggplot2::geom_hline(yintercept = input$upper_cutoff, linetype = 3),
                    
                    ggplot2::annotate(
                      "text",
                      x = 0,
                      # y = (input$lower_cutoff + input$upper_cutoff) / 2,
                      y = input$upper_cutoff,
                      label = paste(
                        "Hampel window size(k):", input$hampel_k,
                        "\nHampel threshold(t0):", input$hampel_t0
                      ),
                      hjust = 1.1,
                      vjust = 1.1,
                      size = 4
                    )
                  )
                }
                
                
                
                
                # Create the plot with dynamic coloring and Y-axis
                p <- ggplot(plot_data, aes_string(x = "Minutes", y = y_var, group = "Samples", color = color_var)) +
                    geom_line(linewidth = 1) +
                    ylab(y_labels[[y_var]]) +
                    xlab("Time (minutes)") +
                    ggpubr::theme_pubr() +
                    geom_vline(xintercept = 0, linetype = 2) + 
                    geom_hline(yintercept = 0, linetype = 1) + 
                  
                  hampel_layers +
                  cutoff_layers
                
                #r2.2.4  
                
                if (run_gam_test) {
                  
                  label_text <- if (is.na(gam_p_value)) {
                    "Curve test: NA"
                  } else if (gam_p_value < 1e-4) {
                    "Curve diff: p < 1e-4"
                  } else {
                    paste0("Curve diff: p = ", signif(gam_p_value, 3))
                  }
                  
 
                  p <- p + annotate(
                    "text",
                    x = 0, y = Inf,
                    label = label_text,
                    hjust = 1.1, vjust = 1.5,
                    size = 5,
                    fontface = "bold"
                  )
                   
                }
                

                # --- Two-period mode ---
                if (input$analysis_model == "twoperiod" || input$analysis_model == "twoperiodglm") {

                  # highlight windows
                  
                  if (input$two_periods=="bothc"){
                    p <- p +
                      geom_vline(xintercept = input$cut_time, linetype = "dashed") +
                      
                      annotate("rect",
                               xmin = input$cut_time - input$time_windowb4,
                               xmax = input$cut_time,
                               ymin = -Inf, ymax = Inf,
                               alpha = 0.1) +
                      
                      annotate("rect",
                               xmin = input$cut_time,
                               xmax = input$cut_time + input$time_window,
                               ymin = -Inf, ymax = Inf,
                               alpha = 0.1) +
                      # r2.2.12
                      geom_rect(
                        data = window_df(),
                        aes(
                          xmin = xmin,
                          xmax = xmax,
                          ymin = mean - ci,
                          ymax = mean + ci,
                          fill = Conditions
                        ),
                        alpha = 0.2,
                        inherit.aes = FALSE
                      ) +
                      geom_segment(
                        data = window_df(),
                        aes(
                          x = xmin,
                          xend = xmax,
                          y = mean,
                          yend = mean,
                          color = Conditions
                        ),
                        linewidth = 2,
                        linetype = "solid",
                        inherit.aes = FALSE
                      )
                  }
                  
                  
                  
                  if (input$two_periods=="both"){
                  p <- p +
                    geom_vline(xintercept = input$cut_time, linetype = "dashed") +
                    
                    annotate("rect",
                             xmin = input$cut_time - input$time_windowb4,
                             xmax = input$cut_time,
                             ymin = -Inf, ymax = Inf,
                             alpha = 0.1) +
                    
                    annotate("rect",
                             xmin = input$cut_time,
                             xmax = input$cut_time + input$time_window,
                             ymin = -Inf, ymax = Inf,
                             alpha = 0.1) +
                  # r2.2.12
                    geom_rect(
                      data = window_df(),
                      aes(
                        xmin = xmin,
                        xmax = xmax,
                        ymin = mean - ci,
                        ymax = mean + ci,
                        fill = Conditions
                      ),
                      alpha = 0.2,
                      inherit.aes = FALSE
                    ) +
                    geom_segment(
                      data = window_df(),
                      aes(
                        x = xmin,
                        xend = xmax,
                        y = mean,
                        yend = mean,
                        color = Conditions
                      ),
                      linewidth = 2,
                      linetype = "solid",
                      inherit.aes = FALSE
                    )
                  }
                  
                  if (input$two_periods=="before"){
                    p <- p +
                      geom_vline(xintercept = input$cut_time, linetype = "dashed") +
                      
                      annotate("rect",
                               xmin = input$cut_time - input$time_windowb4,
                               xmax = input$cut_time,
                               ymin = -Inf, ymax = Inf,
                               alpha = 0.1) +
                    # r2.2.12
                    geom_rect(
                      data = window_df(),
                      aes(
                        xmin = xmin,
                        xmax = xmax,
                        ymin = mean - ci,
                        ymax = mean + ci,
                        fill = Conditions
                      ),
                      alpha = 0.2,
                      inherit.aes = FALSE
                    ) +
                      geom_segment(
                        data = window_df(),
                        aes(
                          x = xmin,
                          xend = xmax,
                          y = mean,
                          yend = mean,
                          color = Conditions
                        ),
                        linewidth = 2,
                        linetype = "solid",
                        inherit.aes = FALSE
                      )
                  }
                  
                  if (input$two_periods=="after"){
                    p <- p +
                      geom_vline(xintercept = input$cut_time, linetype = "dashed") +
                      
                      annotate("rect",
                               xmin = input$cut_time,
                               xmax = input$cut_time + input$time_window,
                               ymin = -Inf, ymax = Inf,
                               alpha = 0.1) +
                      # r2.2.12
                      geom_rect(
                        data = window_df(),
                        aes(
                          xmin = xmin,
                          xmax = xmax,
                          ymin = mean - ci,
                          ymax = mean + ci,
                          fill = Conditions
                        ),
                        alpha = 0.2,
                        inherit.aes = FALSE
                      ) +
                      geom_segment(
                        data = window_df(),
                        aes(
                          x = xmin,
                          xend = xmax,
                          y = mean,
                          yend = mean,
                          color = Conditions
                        ),
                        linewidth = 2,
                        linetype = "solid",
                        inherit.aes = FALSE
                      )
                  }
                  
                  
                  
                  # ✅ Add ANCOVA p-value
                  p_val <- ancova_p_value()  # <-- call the reactive
                  
                  if (input$two_periods %in% c("both", "bothc", "before", "after") && input$analysis_model == "twoperiodglm") {
                  int_p_val <- interaction_p()   #interaction p-val
                  } else {int_p_val <- NA}
                  
                  label_text <- if (is.null(p_val) || !is.finite(p_val)) {
                    "ANCOVA p = NA"
                  } else {
                    paste0("ANCOVA p = ", signif(as.numeric(p_val), 3))
                  }

                  p <- p +
                    annotate(
                      "text",
                      x = input$cut_time, y = Inf,
                      label = paste0("Condition Effect: p = ", signif(p_val, 3)),
                      hjust = 1.1,
                      vjust = 1.1,
                      size = 5,
                      fontface = "bold"
                    ) +
                    
                    ggplot2::annotate(
                      "text",
                      x = 10,
                      y = Inf,
                      label = pairwise_text(),
                      hjust = 1.1,
                      vjust = 3.1,
                      size = 4
                    ) 
                    
                    # annotate(
                    #   "text",
                    #   x = input$cut_time, y = Inf,
                    #   label = paste0("Interact p = ", signif(int_p_val, 3)),
                    #   hjust = 1.1,
                    #   vjust = 5.1,
                    #   size = 4,
                    #   fontface = "bold"
                    # )
                  
                  # Conditionally add interaction p-value if it exists
                 
                  if (!is.null(int_p_val) && length(int_p_val) > 0 && all(is.finite(int_p_val))) {
                    
                    if (length(int_p_val) == 1) {
                      label_text <- paste(
                        paste0("Condition:Bodaymass", ": p = ", signif(int_p_val, 3))
                      )
                    } else if (length(int_p_val) > 1) {
                      
                      label_text <- paste(
                        paste0(names(int_p_val), ": p = ", signif(int_p_val, 3)),
                        collapse = "\n"
                      )
                      
                      label_text <- paste("Interaction:\n", label_text)
                      
                    } 
                  
                  p <- p +
                    annotate(
                      "text",
                      # x = input$cut_time,
                      x = 50,
                      y = Inf,
                      
                      label = label_text,
                      hjust = 1.1,
                      vjust = 2.1,
                      size = 4
                    )
                  }
                  
                }
                

                
                # ADD model-based curves:
                if (!is.null(gam_model)) {
                  
                  df_model <- gam_model$model   # safest source
                  
                  newdata <- expand.grid(
                    Minutes = seq(min(df_model$Minutes), max(df_model$Minutes), length.out = 200),
                    Conditions = levels(df_model$Conditions)
                  )
                  
                  # If you HAVE Samples (random effect model)
                  if ("Samples" %in% names(df_model)) {
                    newdata$Samples <- levels(df_model$Samples)[1]
                    newdata$Samples <- factor(newdata$Samples, levels = levels(df_model$Samples))
                  }
                  

                  
                  
                  
                  newdata$Conditions <- factor(newdata$Conditions, levels = levels(df_model$Conditions))
                  newdata$bodymass <- mean(df_model$bodymass)
                  
                  pred <- predict(
                    gam_model,
                    newdata = newdata,
                    se.fit = TRUE,
                    type = "response",
                    exclude = if ("Samples" %in% names(df_model)) "s(Samples)" else NULL
                  )
                  
                  newdata$fit   <- pred$fit
                  newdata$se    <- pred$se.fit
                  newdata$upper <- newdata$fit + 2 * newdata$se
                  newdata$lower <- newdata$fit - 2 * newdata$se
                  
                  
                  # If you HAVE no Samples (no random effect model)
                  if (!("Samples" %in% names(df_model))) {
                    
                  sigma <- sqrt(gam_model$sig2)
                  newdata$upper <- newdata$fit + 2 * sqrt(newdata$se^2 + sigma^2)
                  newdata$lower <- newdata$fit - 2 * sqrt(newdata$se^2 + sigma^2)     
                  }
                  
                  #r2.2.9
                  # Format for display
                  formatted_pvals <- reactive({
                    
                    df <- gam_pairwise_pvals()
                    req(!is.null(df), nrow(df) > 0)
                    
                    df$label <- paste0(
                      df$Comparison,
                      ": p = ",
                      signif(df$p_value, 3)
                    )
                    
                    df
                  })
                  pvals <- formatted_pvals()
                  gam_annotation_text <- paste(pvals$label, collapse = "\n")
                  
                  
                  
                  p <- p +
                    # Confidence bands
                    geom_ribbon(
                      data = newdata,
                      aes(x = Minutes, ymin = lower, ymax = upper, fill = Conditions),
                      alpha = 0.2,
                      inherit.aes = FALSE
                    ) +
                    
                    # Model curves
                    geom_line(
                      data = newdata,
                      aes(x = Minutes, y = fit, color = Conditions),
                      linewidth = 2.5,
                      inherit.aes = FALSE
                    ) +
                    
                    annotate(
                    "text",
                    x = 0,
                    y = input$upper_cutoff,
                    label = gam_annotation_text,
                    hjust = 1.1,
                    vjust = 1.1,
                    size = 4
                  )
                 
                }
                
                
                  
                # Apply appropriate color palette based on what we're coloring by
                if (color_var == "Samples") {
                    p <- p + scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 8, name = "Accent"),
                                                         RColorBrewer::brewer.pal(n = 12, name = "Set3")),
                                               name = color_title)
                } else {
                    p <- p + scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                                                         RColorBrewer::brewer.pal(n = 9, name = "Pastel1")),
                                               name = color_title)
                }
                
                # Apply faceting if enabled
                if (input$use_faceting && !is.null(input$facet_by) && input$facet_by %in% colnames(plot_data)) {
                    facet_cols <- if (!is.null(input$facet_cols) && input$facet_cols > 0) input$facet_cols else 2
                    p <- p + facet_wrap(as.formula(paste("~", input$facet_by)), ncol = facet_cols) +
                         theme(legend.position = "bottom", 
                               strip.text = element_text(size = 10, face = "bold")) +
                         ggtitle(paste("CLAMS Data - Faceted by", input$facet_by, "- Colored by", color_title))
                } else {
                    p <- p + theme(legend.position = "right") +
                         # ggtitle(paste("CLAMS Data Preview - Colored by", color_title))
                    ggtitle(paste(input$experiment_name, " Data Preview - Colored by", color_title))
                    
                }
                
                return(p)
                } else {
                    # No data or Y-variable not found
                    p <- ggplot() + 
                        annotate("text", x = 0.5, y = 0.5, 
                                label = paste0("Data column '", y_var, "' not found.\nAvailable columns:\n", 
                                             paste(head(names(plot_data), 10), collapse = ", ")), 
                                size = 5, hjust = 0.5, vjust = 0.5) +
                        theme_void() +
                        ggtitle("CLAMS Data Plot - Column Not Found")
                    return(p)
                }
            
            
            } 
            
        } else {
            # DSI plot logic - based on original DSI app
            data <- results$processed_data
            metadata <- results$metadata
            
            # Ensure required inputs are available
            if (is.null(input$dsi_y_var) || !nzchar(input$dsi_y_var)) {
                # Return a message plot if no Y-variable is selected
                p <- ggplot() + 
                    annotate("text", x = 0.5, y = 0.5, label = "Please select a Y-axis variable from the dropdown above", 
                             size = 6, hjust = 0.5, vjust = 0.5) +
                    theme_void() +
                    ggtitle("DSI Data Plot - Select Y-Variable")
                return(p)
            }
            
            # Ensure required time column exists
            if (!"RelativeTime_Minutes" %in% names(data)) {
                p <- ggplot() + 
                    annotate("text", x = 0.5, y = 0.5, label = "RelativeTime_Minutes column not found in processed data", 
                             size = 6, hjust = 0.5, vjust = 0.5) +
                    theme_void() +
                    ggtitle("DSI Data Processing Error")
                return(p)
            }
            
            # Discover per-mouse columns and available metrics
            per_mouse_cols <- grep("^Mouse_\\d+_", names(data), value = TRUE)
            if (length(per_mouse_cols) == 0) {
                p <- ggplot() + 
                    annotate("text", x = 0.5, y = 0.5, label = "No mouse-specific columns found in processed data", 
                             size = 6, hjust = 0.5, vjust = 0.5) +
                    theme_void() +
                    ggtitle("DSI Data Processing Error")
                return(p)
            }
            
            base_metrics <- sort(unique(sub("^Mouse_\\d+_", "", per_mouse_cols)))
            y_sel <- input$dsi_y_var
            
            # Columns that match the selected metric
            matching_cols <- per_mouse_cols[sub("^Mouse_\\d+_", "", per_mouse_cols) == y_sel]
            if (length(matching_cols) == 0) {
                p <- ggplot() + 
                    annotate("text", x = 0.5, y = 0.5, 
                             label = paste0("No columns found for metric: ", y_sel, "\nAvailable metrics: ", 
                                          paste(head(base_metrics, 5), collapse = ", ")), 
                             size = 5, hjust = 0.5, vjust = 0.5) +
                    theme_void() +
                    ggtitle("DSI Metric Selection Error")
                return(p)
            }
            
            # Build condition lookup from metadata
            cond_map <- NULL
            if (!is.null(metadata) && "condition" %in% tolower(names(metadata))) {
                names(metadata) <- tolower(trimws(names(metadata)))
                if ("mouse_id" %in% names(metadata)) {
                    metadata$mouse_id <- trimws(as.character(metadata$mouse_id))
                    metadata$condition <- trimws(as.character(metadata$condition))
                    
                    md_unique <- metadata[!duplicated(metadata$mouse_id) & nzchar(metadata$mouse_id), c("mouse_id", "condition")]
                    if (nrow(md_unique) > 0) {
                        cond_map <- setNames(md_unique$condition, md_unique$mouse_id)
                    }
                }
                
                # Fallback: map by order if mouse_id mapping fails
                if (is.null(cond_map)) {
                    mvec <- sub("^(Mouse_\\d+)_.*", "\\1", matching_cols)
                    k <- min(length(mvec), nrow(metadata))
                    if (k > 0) {
                        cond_map <- setNames(metadata$condition[seq_len(k)], mvec[seq_len(k)])
                    }
                }
            }
            
            # Create long-form data for plotting with better error handling
            plot_df_list <- tryCatch({
                lapply(matching_cols, function(colname) {
                    m <- sub("^(Mouse_\\d+)_.*", "\\1", colname)
                    cond_val <- if (!is.null(cond_map)) cond_map[[m]] else NA_character_
                    
                    # Check if column exists and has data
                    if (!colname %in% names(data)) {
                        return(data.frame(RelativeTime_Minutes = numeric(0), value = numeric(0), 
                                        mouse = character(0), condition = character(0), stringsAsFactors = FALSE))
                    }
                    
                    # Extract values safely
                    col_values <- data[[colname]]
                    if (is.null(col_values) || length(col_values) == 0) {
                        return(data.frame(RelativeTime_Minutes = numeric(0), value = numeric(0), 
                                        mouse = character(0), condition = character(0), stringsAsFactors = FALSE))
                    }
                    
                    data.frame(
                        RelativeTime_Minutes = data$RelativeTime_Minutes,
                        value = suppressWarnings(as.numeric(col_values)),
                        mouse = m,
                        condition = ifelse(is.na(cond_val) | !nzchar(cond_val), "Unknown", as.character(cond_val)),
                        stringsAsFactors = FALSE
                    )
                })
            }, error = function(e) {
                p <- ggplot() + 
                    annotate("text", x = 0.5, y = 0.5, 
                             label = paste0("Error creating plot data: ", e$message), 
                             size = 5, hjust = 0.5, vjust = 0.5) +
                    theme_void() +
                    ggtitle("DSI Plotting Error")
                return(p)
            })
            
            plot_df <- do.call(rbind, plot_df_list)
            plot_df <- plot_df[!is.na(plot_df$value), , drop = FALSE]
            
            if (nrow(plot_df) == 0) return(NULL)
            
            # Determine grouping
            group_by <- if (is.null(input$dsi_plot_group_by)) "mouse" else input$dsi_plot_group_by
            
            if (group_by == "condition") {
                grp_levels <- sort(unique(plot_df$condition))
                plot_df$group_label <- factor(plot_df$condition, levels = grp_levels)
                color_lab <- "Condition"
            } else {
                grp_levels <- sort(unique(plot_df$mouse))
                plot_df$group_label <- factor(plot_df$mouse, levels = grp_levels)
                color_lab <- "Mouse"
            }
            
            # Create the plot
            p <- ggplot(plot_df, aes(x = RelativeTime_Minutes, y = value, color = group_label)) +
                geom_line(linewidth = 1.2, alpha = 0.9) +
                labs(title = paste("DSI Data:", y_sel), 
                     x = "Relative Time (minutes)", 
                     y = y_sel, 
                     color = color_lab) +
                ggpubr::theme_pubr() +
                geom_vline(xintercept = 0, linetype = 2)
            
            # Apply faceting if requested
            if (isTRUE(input$dsi_plot_facet)) {
                if (group_by == "condition") {
                    p <- p + facet_wrap(~ condition, scales = "free_y")
                } else {
                    p <- p + facet_wrap(~ mouse, scales = "free_y")
                }
                p <- p + theme(legend.position = "bottom")
            } else {
                p <- p + theme(legend.position = "right")
            }
            
            # Apply color palette
            n <- length(grp_levels)
            if (n <= 9) {
                cols <- RColorBrewer::brewer.pal(max(3, n), "Set1")[seq_len(n)]
            } else {
                base_cols <- RColorBrewer::brewer.pal(8, "Dark2")
                cols <- base_cols[(seq_len(n) - 1) %% length(base_cols) + 1]
            }
            names(cols) <- grp_levels
            p <- p + scale_color_manual(values = cols, drop = FALSE)
            
            return(p)
        }
 
        # Fallback plot
        ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No data available for plotting"), 
                     size = 6, hjust = 0.5, vjust = 0.5) + 
            theme_void()
    })
    output$plot_preview <- renderPlot({ main_plot_obj() })
    
    model_text <- reactive({
      if (input$bodymass_source == "zeromass") {
        "ANOVA"
      } else if (input$analysis_model == "twoperiod") {
        "ANCOVA"
      } else if (input$analysis_model == "twoperiodglm") {
        "GLM"
      } else if (input$analysis_model == "gam") {
        "GAM"
      } else {
        "Unknown"
      }
    })
    
    
    # Data table for DSI
    output$data_table <- DT::renderDataTable({
        if (input$data_type == "dsi" && !is.null(analysis_results())) {
            results <- analysis_results()
            if (!is.null(results$processed_data)) {
                return(DT::datatable(head(results$processed_data, 100), 
                                   options = list(pageLength = 10, scrollX = TRUE)))
            }
        }
        return(DT::datatable(data.frame(Message = "No data to display")))
    })
    
    # Download buttons
    output$download_buttons <- renderUI({
        req(analysis_results())
        results <- analysis_results()
        if (is.null(results)) {
            return(p("Analysis failed - no downloads available"))
        }
        
        tagList(
            downloadButton("download_plots", "Download Plots (PDF)", class = "btn-primary"),
            br(), br(),
            downloadButton("download_two_periods", "Download Selected Data (CSV)", class = "btn-success"),
            br(), br(),
            downloadButton("download_summary", "Download CalR Compatible (CSV)", class = "btn-info"),
            br(), br(),
            downloadButton("download_data", "Download Processed Data (CSV)", class = "btn-success"),
        )
    })
    
    # v2.3.9, v2.3.12
    output$download_two_periods <- downloadHandler(
      filename = function() {
        
        if (input$two_periods == "before") {
          paste0(input$clams_y_var,"_",input$bodymass_source,"_b",input$time_windowb4,"_c",input$cut_time,"_",Sys.Date(), ".csv")
        } else if (input$two_periods == "after") {
          paste0(input$clams_y_var,"_",input$bodymass_source,"_c",input$cut_time,"_a",input$time_window,"_",Sys.Date(), ".csv")
        } else if (input$two_periods == "both") {
          paste0(input$clams_y_var,"_",input$bodymass_source,"_b",input$time_windowb4,"_c",input$cut_time,"_a",input$time_window,"_",Sys.Date(), ".csv")
        }
      },
      content = function(file) {
        df <- two_period_results()   # call reactive
        
        # Optional: check for empty data
        req(df, nrow(df) > 0)
        
        write.csv(df, file, row.names = FALSE)
      }
    )
    
    
    # Download handlers (simplified for demo)
    
    # Download plots/tables
    output$download_plots <- downloadHandler(
        filename = function() {
          if (input$two_periods == "before") {
            paste0(input$experiment_name, "_", input$clams_y_var,"_",input$bodymass_source,"_b",input$time_windowb4,"_c",input$cut_time,"_",model_text(), "_Plots.pdf")
          } else if (input$two_periods == "after") {
            paste0(input$experiment_name, "_", input$clams_y_var,"_",input$bodymass_source,"_c",input$cut_time,"_a",input$time_window,"_",model_text(), "_Plots.pdf")
          } else if (input$two_periods == "both") {
            paste0(input$experiment_name, "_", input$clams_y_var,"_",input$bodymass_source,"_b",input$time_windowb4,"_c",input$cut_time,"_a",input$time_window,"_",model_text(), "_Plots.pdf")
          }
          },
        
        content = function(file) {
            results <- analysis_results()
            if (is.null(results)) {
                add_message("No analysis results available to download", "error")
                return()
            }
        
            # Open PDF device
            pdf(file, width = 10, height = 7)
            
            tryCatch({
                if (input$data_type == "clams") {
                    # CLAMS plot rendering
                    y_var <- input$clams_y_var
                    
                    # Define Y-axis labels
                    y_labels <- list(
                        si13c = expression(paste(text=''^13,'C in exhaled CO2 (per mil)')),
                        vo2 = "VO₂ - Oxygen Consumption (ml/kg/hr)",
                        vco2 = "VCO₂ - CO₂ Production (ml/kg/hr)",
                        rq = "Respiratory Quotient (RQ)",
                        kcal_hr = "Energy Expenditure (kcal/hr)",
                        vh2o = "Water Vapor - VH₂O (ml/kg/hr)",
                        foodupa = "Food Intake (grams)",
                        waterupa = "Water Intake A (ml)",
                        waterupb = "Water Intake B (ml)",
                        bodymass = "Body Mass (grams)",
                        xbreak = "X-axis Activity (beam breaks)",
                        ybreak = "Y-axis Activity (beam breaks)",
                        zbreak = "Z-axis Activity (beam breaks)",
                        pedmeters = "Pedometer Counts",
                        allmeters = "Total Activity Counts",
                        bodytemp = "Body Temperature (°C)"
                    )
                    
                    plot_data <- results$CombinedData
                    color_var <- if (!is.null(input$plot_color_by) && input$plot_color_by == "Samples") "Samples" else "Conditions"
                    color_title <- if (color_var == "Samples") "Sample" else "Condition"
                    facet_cols <- if (!is.null(input$facet_cols)) input$facet_cols else 2

                    #r2.2.4, Replace your current geom_line() section with:
                    
                    # Create plot
                    p <- ggplot(plot_data, aes_string(x = "Minutes", y = y_var, group = "Samples", color = color_var)) +
                        geom_line(size = 0.7) +
                        ylab(y_labels[[y_var]]) +
                        geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
                        ggpubr::theme_pubr() +
                        theme(legend.position = "right")


                    # Apply color scheme (matching preview)
                    if (color_var == "Samples") {
                        p <- p + scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 8, name = "Accent"),
                                                             RColorBrewer::brewer.pal(n = 12, name = "Set3")),
                                                   name = color_title)
                    } else {
                        p <- p + scale_color_manual(values = c(RColorBrewer::brewer.pal(n = 9, name = "Set1"),
                                                             RColorBrewer::brewer.pal(n = 9, name = "Pastel1")),
                                                   name = color_title)
                    }
                    
                    # Apply faceting if enabled
                    if (input$use_faceting && !is.null(input$facet_by) && input$facet_by %in% colnames(plot_data)) {
                        p <- p + facet_wrap(as.formula(paste("~", input$facet_by)), ncol = facet_cols) +
                            theme(legend.position = "bottom",
                                  strip.text = element_text(size = 10, face = "bold"))
                    }
                    
                    p <- p + ggtitle(paste0(input$experiment_name, " - CLAMS Analysis"))
                    # print(p)
                    
                    # Plot 1
                    # print(main_plot_obj())
                    
                    # Plot 2
                    # print(ancova_plot_obj())
                    
                    # grid::grid.newpage()
                    # gridExtra::grid.table(stats_table())
                    # gridExtra::grid.table(pairwise_table())
                    
                    # one page for all
                    gridExtra::grid.arrange(
                      main_plot_obj(),
                      ancova_plot_obj(),
                      gridExtra::tableGrob(stats_table()),
                      gridExtra::tableGrob(pairwise_table()),
                      ncol = 2,
                      heights = c(2, 2)
                    )
                    
                } else {
                    # DSI plot rendering
                    y_var <- input$dsi_y_var
                    if (is.null(y_var) || y_var == "") {
                        plot.new()
                        text(0.5, 0.5, "No Y-variable selected", cex = 1.5)
                    } else {
                        plot_data <- results$processed_data
                        
                        # Get all mouse columns for this variable - match pattern like Mouse_1_Glucose_Gavg_mg_dL
                        # The y_var might be the base name like "Glucose_Gavg_mg_dL"
                        mouse_cols <- grep(paste0("^Mouse_\\d+_", y_var), names(plot_data), value = TRUE)
                        
                        if (length(mouse_cols) > 0) {
                            add_message(paste0("Generating DSI plot for ", y_var, " with ", length(mouse_cols), " mice"), "info")
                            
                            # Reshape data for plotting
                            long_data <- plot_data %>%
                                select(RelativeTime_Minutes, all_of(mouse_cols)) %>%
                                pivot_longer(cols = -RelativeTime_Minutes, names_to = "Mouse_Metric", values_to = "Value") %>%
                                mutate(Mouse = gsub("_.*", "", Mouse_Metric))
                            
                            # Add condition if metadata available
                            if (!is.null(results$metadata)) {
                                metadata <- results$metadata
                                names(metadata) <- tolower(names(metadata))  # Normalize to lowercase
                                
                                if ("mouse_id" %in% names(metadata) && "condition" %in% names(metadata)) {
                                    long_data <- long_data %>%
                                        left_join(metadata %>% select(mouse_id, condition), 
                                                by = c("Mouse" = "mouse_id")) %>%
                                        rename(Condition = condition)  # Rename for consistent use
                                }
                            }
                            
                            # Create plot
                            if (input$dsi_plot_group_by == "mouse" || !"Condition" %in% names(long_data)) {
                                p <- ggplot(long_data, aes(x = RelativeTime_Minutes, y = Value, color = Mouse, group = Mouse)) +
                                    geom_line(size = 0.7) +
                                    labs(x = "Time (minutes)", y = gsub("_", " ", y_var), 
                                         title = paste0(input$experiment_name, " - DSI Analysis: ", y_var))
                            } else {
                                p <- ggplot(long_data, aes(x = RelativeTime_Minutes, y = Value, color = Condition, group = Mouse)) +
                                    geom_line(size = 0.7) +
                                    labs(x = "Time (minutes)", y = gsub("_", " ", y_var),
                                         title = paste0(input$experiment_name, " - DSI Analysis: ", y_var))
                            }
                            
                            p <- p + 
                                geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
                                ggpubr::theme_pubr() +
                                theme(legend.position = "right")
                            
                            if (input$dsi_plot_facet && "Condition" %in% names(long_data)) {
                                p <- p + facet_wrap(~Condition, ncol = 2)
                            }
                            
                            print(p)
                        } else {
                            add_message(paste0("No columns found matching pattern: ^Mouse_\\d+_", y_var), "warning")
                            plot.new()
                            text(0.5, 0.5, paste("No data found for:", y_var), cex = 1.5)
                        }
                    }
                }
                
                add_message("Plots downloaded successfully", "success")
                
            }, error = function(e) {
                add_message(paste("Error generating plot PDF:", e$message), "error")
                plot.new()
                text(0.5, 0.5, paste("Error:", e$message), cex = 1)
            })
            
            dev.off()
        } )
    
    # download ProcessedData.csv
    output$download_data <- downloadHandler(
        filename = function() { paste0(input$experiment_name, "_", toupper(input$data_type), "_ProcessedData.csv") },
        content = function(file) {
            results <- analysis_results()
            if (is.null(results)) {
                add_message("No analysis results to download", "error")
                return()
            }
            
            if (input$data_type == "clams" && !is.null(results$CombinedData)) {
                data_to_export <- results$CombinedData
                add_message(paste0("Exporting CLAMS data: ", nrow(data_to_export), " rows × ", ncol(data_to_export), " columns"), "info")
                write.csv(data_to_export, file, row.names = FALSE)
                add_message("CLAMS processed data downloaded successfully", "success")
                
            } else if (input$data_type == "dsi" && !is.null(results$processed_data)) {
                data_to_export <- results$processed_data
                
                # Count mice in the data
                mouse_cols <- grep("^Mouse_\\d+_", names(data_to_export), value = TRUE)
                unique_mice <- unique(sub("_.*", "", mouse_cols))
                n_mice <- length(unique_mice)
                
                add_message(paste0("Exporting DSI data: ", nrow(data_to_export), " rows × ", ncol(data_to_export), " columns"), "info")
                add_message(paste0("Data includes ", n_mice, " mice: ", paste(unique_mice, collapse = ", ")), "info")
                
                write.csv(data_to_export, file, row.names = FALSE)
                add_message(paste0("DSI processed data downloaded successfully (", n_mice, " mice included)"), "success")
            } else {
                add_message("No data available to download", "error")
            }
        } )
    
    # download CalR compatible csv
    output$download_summary <- downloadHandler(
        filename = function() { paste0(input$experiment_name, "_CalR_Compatible.csv") },
        content = function(file) {
            results <- analysis_results()
            
            if (input$data_type == "clams" && !is.null(results) && !is.null(results$CombinedData)) {
                # CLAMS: Create variable comparison file organized by variable with all samples
                data <- results$CombinedData
                
                # Get all available metric columns
                all_metric_names <- c("vo2", "vco2", "vh2o", "kcal_hr", "rq", "foodupa", "waterupa", 
                                     "waterupb", "bodymass", "xbreak", "ybreak", "zbreak", "pedmeters", 
                                     "allmeters", "si13c", "bodytemp")
                available_metrics <- all_metric_names[all_metric_names %in% names(data)]
                
                if (length(available_metrics) == 0) {
                    add_message("No metric columns found for variable comparison", "warning")
                    # Fallback to simple summary
                    summary_data <- data.frame(
                        Info = "No metric data available",
                        Value = "Check analysis results"
                    )
                    write.csv(summary_data, file, row.names = FALSE)
                    return()
                }
                
                # Get unique samples and conditions for ordering
                all_samples <- unique(data$Samples)
                
                add_message(paste0("Creating CLAMS variable comparison file with ", length(available_metrics), 
                                 " variables across ", length(all_samples), " samples"), "info")
                
                # Create ordering based on group_order (if available) or conditions
                sample_order <- all_samples  # Default order
                
                # Check if group_order exists in the data
                if ("group_order" %in% names(data)) {
                    # Create a lookup of sample to group_order
                    sample_group_order <- data %>%
                        select(Samples, group_order) %>%
                        distinct() %>%
                        arrange(group_order)
                    
                    sample_order <- sample_group_order$Samples
                    
                    add_message(paste0("Organizing columns by custom group_order from metadata"), "info")
                    
                } else if ("Conditions" %in% names(data)) {
                    # Fallback to condition-based ordering
                    sample_condition <- data %>%
                        select(Samples, Conditions) %>%
                        distinct() %>%
                        arrange(Conditions, Samples)
                    
                    sample_order <- sample_condition$Samples
                    
                    add_message(paste0("Organizing columns by condition groups: ", 
                                     paste(unique(data$Conditions), collapse = ", ")), "info")
                }
                
                # Build the organized data frame
                # Start with Minutes column (can be negative for baseline period)
                time_data <- data %>%
                    select(Minutes) %>%
                    distinct() %>%
                    arrange(Minutes)
                
                organized_data <- data.frame(Minutes = time_data$Minutes)
                
                # Add individual Date_Time columns for each sample first (for easy time checking)
                for (sample in sample_order) {
                    sample_data <- data %>%
                        filter(Samples == sample) %>%
                        arrange(Minutes) %>%
                        select(Minutes, Date_Time) %>%
                        distinct(Minutes, .keep_all = TRUE) #r2.1
                    
                    # Get condition for this sample
                    sample_condition <- unique(data$Conditions[data$Samples == sample])
                    if (length(sample_condition) == 0) sample_condition <- "Unknown"
                    
                    # Create Date_Time column name: DateTime__Condition__Sample
                    datetime_col_name <- paste0("DateTime__", sample_condition[1], "__", sample)
                    
                    # Match by Minutes and add Date_Time column
                    if (nrow(sample_data) > 0) {
                      organized_data <- organized_data %>%
                        left_join(
                          sample_data %>% rename(!!datetime_col_name := Date_Time),
                          by = "Minutes"
                        )
                    }
                }
                
                # For each variable, add columns for all samples in the specified order
                for (var_name in available_metrics) {
                    for (sample in sample_order) {
                        # Get data for this sample and variable
                        sample_data <- data %>%
                            filter(Samples == sample) %>%
                            arrange(Minutes) %>%
                            select(Minutes, all_of(var_name)) %>%
                          distinct(Minutes, .keep_all = TRUE) #r2.1
                        
                        # Get condition for this sample
                        sample_condition <- unique(data$Conditions[data$Samples == sample])
                        if (length(sample_condition) == 0) sample_condition <- "Unknown"
                        
                        # Create descriptive column name: Variable__Condition__Sample
                        new_col_name <- paste0(var_name, "__", sample_condition[1], "__", sample)
                        
                        # Match by Minutes and add column
                        if (nrow(sample_data) > 0) {
                            organized_data <- organized_data %>%
                                left_join(
                                    sample_data %>% rename(!!new_col_name := all_of(var_name)),
                                    by = "Minutes"
                                )
                        }
                    }
                }
                
                
                # # Write to CSV
                # write.csv(organized_data, file, row.names = FALSE)
                # 
                # add_message(paste0("CLAMS variable comparison file created: ", nrow(organized_data), 
                #                    " timepoints × ", ncol(organized_data), " columns"), "success")
                # add_message("Columns organized by: Minutes (negative = baseline) → DateTime per sample → Variables per sample", "info")
                # add_message("Columns organized by: Variable → Condition → Sample", "info")
                
                
                
                # Convert to clean data frame for CalR format
                df1 <- organized_data
            
                # First simplify column names
                names(df1) <- sub("^([^_]+).*_([^_]+)$", "\\1_\\2", names(df1))
                
                # fix kcal names
                names(df1) <- sub("^kcal_", "kcal_hr_", names(df1))
                
                # Rename DateTime_1 -> Date_Time_1
                #names(df1)[names(df1) == "DateTime_1"] <- "Date_Time_1"
                
                # 6/18/2026, some time the above code does not work!
                # Rename second col to -> Date_Time_1
                colnames(df1)[2] <- "Date_Time_1"
                
                # Remove all other DateTime columns
                df1 <- df1[, !(grepl("^DateTime_", names(df1)) & names(df1) != "Date_Time_1")]
                
                # Move Date_Time_1 to first column
                if ("Date_Time_1" %in% names(df1)) {
                  df1 <- df1[, c("Date_Time_1", setdiff(names(df1), "Date_Time_1"))]
                }
                
                # Remove Minutes column
                df1 <- df1[, names(df1) != "Minutes"]
                
                # Ensure numeric columns
                df1[-1] <- lapply(df1[-1], as.numeric)
                
                write.csv(df1, file, row.names = FALSE, na = "")
                
            } 
            else if (input$data_type == "dsi" && !is.null(results) && !is.null(results$processed_data)) {
                # DSI: Create variable comparison file organized by variable with all mice
                add_message("=== STARTING DSI VARIABLE COMPARISON GENERATION ===", "info")
                
                data <- results$processed_data
                metadata <- results$metadata
                
                # Comprehensive validation
                if (is.null(data)) {
                    add_message("CRITICAL ERROR: results$processed_data is NULL", "error")
                    write.csv(data.frame(Error = "Processed data is NULL"), file, row.names = FALSE)
                    return()
                }
                
                # Validate data has rows and columns
                if (nrow(data) == 0 || ncol(data) == 0) {
                    add_message(paste0("ERROR: Processed data is empty - ", nrow(data), " rows, ", ncol(data), " columns"), "error")
                    write.csv(data.frame(Error = "No processed data available"), file, row.names = FALSE)
                    return()
                }
                
                add_message(paste0("✓ Processed data validated: ", nrow(data), " rows, ", ncol(data), " columns"), "success")
                add_message(paste0("DEBUG: First 20 column names: ", paste(head(names(data), 20), collapse = ", ")), "info")
                
                # Get all mouse-specific columns
                mouse_cols <- grep("^Mouse_\\d+_", names(data), value = TRUE)
                
                add_message(paste0("DEBUG: Found ", length(mouse_cols), " mouse columns in processed data"), "info")
                add_message(paste0("DEBUG: Sample columns: ", paste(head(mouse_cols, 5), collapse = ", ")), "info")
                
                if (length(mouse_cols) == 0) {
                    add_message("No mouse columns found for variable comparison", "warning")
                    add_message(paste0("Available columns: ", paste(head(names(data), 10), collapse = ", ")), "info")
                    # Fallback to simple summary
                    summary_data <- data.frame(
                        Info = "No mouse-specific data available",
                        Value = "Check analysis results"
                    )
                    
                    write.csv(summary_data, file, row.names = FALSE)
                    return()                    
                }
                
                # Extract unique variables and mice
                # Mouse columns are like: Mouse_1_DateTime, Mouse_2_Glucose_Num, etc.
                # We want mice to be: Mouse_1, Mouse_2, etc.
                # We want variables to be: DateTime, Glucose_Num, etc.
                all_variables <- unique(sub("^Mouse_\\d+_", "", mouse_cols))
                all_mice <- unique(sub("^(Mouse_\\d+)_.*", "\\1", mouse_cols))
                
                # Filter out DateTime from variables since we handle it separately
                all_variables <- all_variables[all_variables != "DateTime"]
                
                add_message(paste0("Creating DSI variable comparison file with ", length(all_variables), 
                                 " variables across ", length(all_mice), " mice"), "info")
                add_message(paste0("DEBUG: all_mice = ", paste(all_mice, collapse=", ")), "info")
                add_message(paste0("DEBUG: Sample variables: ", paste(head(all_variables, 10), collapse=", ")), "info")
                
                if (length(all_variables) == 0) {
                    add_message("ERROR: No variables extracted from column names", "error")
                    write.csv(data.frame(Error = "No variables found", 
                                        AvailableColumns = paste(head(names(data), 20), collapse = "; ")), 
                             file, row.names = FALSE)
                    return()
                }
                
                if (length(all_mice) == 0) {
                    add_message("ERROR: No mice extracted from column names", "error")
                    write.csv(data.frame(Error = "No mice found"), file, row.names = FALSE)
                    return()
                }
                
                # Create ordering based on metadata (group_order preferred, then condition)
                mouse_order <- all_mice  # Default order
                
                if (!is.null(metadata)) {
                    # Normalize column names
                    meta_names <- tolower(names(metadata))
                    names(metadata) <- meta_names
                    
                    # Prefer group_order if available
                    if ("mouse_id" %in% meta_names && "group_order" %in% meta_names) {
                        metadata$mouse_id <- trimws(as.character(metadata$mouse_id))
                        metadata$group_order <- as.numeric(metadata$group_order)
                        
                        # Order by group_order
                        metadata_ordered <- metadata[order(metadata$group_order, metadata$mouse_id), ]
                        
                        # Match mice in data with metadata ordering
                        ordered_mice <- metadata_ordered$mouse_id[metadata_ordered$mouse_id %in% all_mice]
                        remaining_mice <- setdiff(all_mice, ordered_mice)
                        mouse_order <- c(ordered_mice, remaining_mice)
                        
                        add_message(paste0("Organizing columns by custom group_order from metadata"), "info")
                        
                    } else if ("mouse_id" %in% meta_names && "condition" %in% meta_names) {
                        # Fallback to condition-based ordering
                        metadata$mouse_id <- trimws(as.character(metadata$mouse_id))
                        metadata$condition <- trimws(as.character(metadata$condition))
                        
                        # Order by condition, then by mouse_id
                        metadata_ordered <- metadata[order(metadata$condition, metadata$mouse_id), ]
                        
                        # Match mice in data with metadata ordering
                        ordered_mice <- metadata_ordered$mouse_id[metadata_ordered$mouse_id %in% all_mice]
                        remaining_mice <- setdiff(all_mice, ordered_mice)
                        mouse_order <- c(ordered_mice, remaining_mice)
                        
                        add_message(paste0("Organizing columns by condition groups: ", 
                                         paste(unique(metadata$condition), collapse = ", ")), "info")
                    }
                }
                
                # Build the organized data frame
                # Start with RelativeTime_Minutes (can be negative for baseline period)
                organized_data <- data.frame(RelativeTime_Minutes = data$RelativeTime_Minutes)
                add_message(paste0("DEBUG: Starting organized_data with ", nrow(organized_data), " rows, ", 
                                 ncol(organized_data), " columns"), "info")
                
                # Add individual DateTime columns for each mouse first (for easy time checking)
                datetime_cols_added <- 0
                for (mouse in mouse_order) {
                    # Look for DateTime column for this mouse
                    datetime_col <- paste0(mouse, "_DateTime")
                    
                    # Get condition for this mouse
                    mouse_condition <- "Unknown"
                    if (!is.null(metadata) && "condition" %in% names(metadata)) {
                        cond <- metadata$condition[metadata$mouse_id == mouse]
                        if (length(cond) > 0 && !is.na(cond[1])) {
                            mouse_condition <- cond[1]
                        }
                    }
                    
                    # Add DateTime column if available
                    if (datetime_col %in% names(data)) {
                        dt_col_name <- paste0("DateTime__", mouse_condition, "__", mouse)
                        organized_data[[dt_col_name]] <- data[[datetime_col]]
                        datetime_cols_added <- datetime_cols_added + 1
                    }
                }
                add_message(paste0("DEBUG: Added ", datetime_cols_added, " DateTime columns. Now have ", 
                                 ncol(organized_data), " total columns"), "info")
                
                # For each variable, add columns for all mice in the specified order
                variable_cols_added <- 0
                for (var_name in sort(all_variables)) {
                    for (mouse in mouse_order) {
                        col_name <- paste0(mouse, "_", var_name)
                        if (col_name %in% names(data)) {
                            # Create descriptive column name
                            new_col_name <- paste0(var_name, "__", mouse)
                            
                            # Add condition to column name if available
                            if (!is.null(metadata) && "condition" %in% names(metadata)) {
                                mouse_condition <- metadata$condition[metadata$mouse_id == mouse]
                                if (length(mouse_condition) > 0 && !is.na(mouse_condition[1])) {
                                    new_col_name <- paste0(var_name, "__", mouse_condition[1], "__", mouse)
                                }
                            }
                            
                            organized_data[[new_col_name]] <- data[[col_name]]
                            variable_cols_added <- variable_cols_added + 1
                        }
                    }
                }
                
                add_message(paste0("DEBUG: Added ", variable_cols_added, " variable columns. Final: ", 
                                 ncol(organized_data), " total columns"), "info")
                
                # Write to CSV
                write.csv(organized_data, file, row.names = FALSE)
                
                add_message(paste0("DSI variable comparison file created: ", nrow(organized_data), 
                                 " rows × ", ncol(organized_data), " columns"), "success")
                add_message("Columns organized by: Variable → Condition → Mouse", "info")
                
            } else {
                # Fallback
                summary_data <- data.frame(
                    Parameter = c("Analysis Type", "Experiment Name", "Analysis Date"),
                    Value = c(toupper(input$data_type), input$experiment_name, as.character(Sys.Date()))
                )
                write.csv(summary_data, file, row.names = FALSE)
                add_message("Summary report downloaded", "success")
            }
        } )
}

# ====================================================================
# LAUNCH THE SHINY APP
# ====================================================================
# If you're in RStudio, click the "Run App" button at the top right
# Or run this entire script with source() or by pressing Ctrl+Shift+S

cat("\n🚀 Launching Rodent Metabolic Phenotyping Analysis Platform...\n")
cat("📊 Ready to analyze CLAMS and DSI data!\n\n")

# Run the Shiny app
shinyApp(ui = ui_app1e, server = server_app1e)
