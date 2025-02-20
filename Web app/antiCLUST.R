################################################################################
#                               Load Libraries
################################################################################

library(shiny)
library(anticlust)
library(tableone)
library(knitr)
library(shinycssloaders)
library(readr)
library(readxl)
library(DT)

################################################################################
#                     Global Helper & Wrapper Functions
################################################################################

# ------------------------------------------------------------------------------
# categories_to_binary(df)
# ------------------------------------------------------------------------------
# Converts categorical features in a data frame into binary (dummy) variables.
# For each factor/character column, creates dummy columns (0/1) for each level.
# Non-factor columns remain unchanged.
#
# Args:
#   df: A data frame containing one or more columns to be transformed.
# 
# Returns:
#   A data frame (actually a matrix) where categorical variables have been
#   replaced with dummy columns.
# ------------------------------------------------------------------------------
categories_to_binary <- function(df) {
  result_list <- list()
  
  for (col_name in names(df)) {
    # If the column is character or factor, convert to factor first, then dummies
    if (is.character(df[[col_name]]) || is.factor(df[[col_name]])) {
      df[[col_name]] <- as.factor(df[[col_name]])
      dummies <- model.matrix(~ . - 1, data = data.frame(x = df[[col_name]]))
      colnames(dummies) <- paste0(col_name, "_", levels(df[[col_name]]))
      result_list[[col_name]] <- dummies
    } else {
      # If already numeric, just keep as-is
      result_list[[col_name]] <- df[[col_name]]
    }
  }
  
  # Combine all processed columns back into one data structure
  do.call(cbind, result_list)
}

# ------------------------------------------------------------------------------
# anticlustering_shiny(numeric_vars, categorical_vars, must_link_constraints, K)
# ------------------------------------------------------------------------------
# Main wrapper function that decides how to perform anticlustering based on:
#   1) The size of the data (N),
#   2) Whether must-link constraints are provided,
#   3) Whether the data is "large" (> 500 items).
#
# Args:
#   numeric_vars: A data frame of numeric features.
#   categorical_vars: A data frame of categorical features (converted to dummies).
#   must_link_constraints: Factor or vector for must-link grouping.
#   K: Number of batches to form.
#
# Returns:
#   A vector of length N specifying each row's batch assignment.
# ------------------------------------------------------------------------------
anticlustering_shiny <- function(
    numeric_vars = NULL, 
    categorical_vars = NULL, 
    must_link_constraints = NULL, 
    K
) {
  # If must-link constraints are specified, call the dedicated must-link function
  if (!is.null(must_link_constraints)) {
    return(must_link_anticlustering_shiny(numeric_vars, categorical_vars, must_link_constraints, K))
  }
  
  # Compute number of rows in either numeric_vars or categorical_vars
  N <- max(
    ifelse(is.null(numeric_vars), 0, nrow(numeric_vars)),
    ifelse(is.null(categorical_vars), 0, nrow(categorical_vars))
  )
  
  # If the dataset is considered "large," switch to a faster anticlustering method
  if (N > 500) {
    # must-link constraints are not supported in the 'fast_anticlustering_shiny' approach here
    if (!is.null(must_link_constraints)) {
      stop("Deal with must link constraints and large N")
    }
    return(fast_anticlustering_shiny(numeric_vars, categorical_vars, K))
  }
  
  # For smaller data, choose the number of 'repetitions' based on size
  repetitions <- 10
  reps <- "10 repetition"
  if (N <= 200) {
    repetitions <- 100
    reps <- "100 repetitions"
  }
  if (N <= 100) {
    repetitions <- 1000
    reps <- "1000 repetitions"
  }
  
  # Convert categorical variables to binary (dummy) columns if needed
  if (!is.null(categorical_vars)) {
    categorical_vars <- categories_to_binary(categorical_vars)
  }
  
  # Combine numeric and dummy-coded categorical features
  input <- cbind(numeric_vars, categorical_vars)
  
  # Provide a console message (for debugging/logging) indicating chosen method
  message("N = ", N, 
          ". Using anticlustering with diversity criterion, method is 'local-maximum' with ", 
          reps, ".\n")
  
  # Perform the anticlustering with local-maximum method and standardization
  anticlustering(
    input, 
    K = K, 
    method = "local-maximum", 
    repetitions = repetitions, 
    standardize = TRUE
  )
}

# ------------------------------------------------------------------------------
# must_link_anticlustering_shiny(numeric_vars, categorical_vars, must_link_constraints, K)
# ------------------------------------------------------------------------------
# Dedicated function to handle must-link constraints using the "2PML" approach.
#
# Args:
#   numeric_vars: A data frame of numeric features.
#   categorical_vars: A data frame of categorical features (converted to dummies).
#   must_link_constraints: Factor or vector for must-link grouping.
#   K: Number of batches to form.
#
# Returns:
#   A vector of length N specifying each row's batch assignment (with must-link).
# ------------------------------------------------------------------------------
must_link_anticlustering_shiny <- function(numeric_vars, categorical_vars, must_link_constraints, K) {
  # Convert any categorical vars to dummy columns
  if (!is.null(categorical_vars)) {
    categorical_vars <- categories_to_binary(categorical_vars)
  }
  
  # Combine numeric and dummy-coded categorical features
  input <- cbind(numeric_vars, categorical_vars)
  
  N <- nrow(input)
  repetitions <- 2  # default, will be overwritten based on data size
  
  # For extremely large data sets, throw an error (not supported here)
  if (N > 15000) {
    stop("Sorry, we do not offer the must-link for data sets with more than 15000 elements. Try out the R package anticlust directly.")
  }
  
  # Adjust number of repetitions based on data size
  if (N <= 1000) {
    repetitions <- 10
  }
  if (N <= 400) {
    repetitions <- 100
  }
  if (N <= 100) {
    repetitions <- 1000
  }
  reps <- paste0(
    repetitions, 
    " repetitions (", 
    repetitions / 2, "x Phase 1, ", 
    repetitions / 2, "x Phase 2)"
  )
  
  # Log the chosen approach for debugging
  message("N = ", N, 
          ". Using anticlustering with diversity criterion, method is '2PML' (for must-link constraints) with ", 
          reps, ".\n")
  
  # Perform the 2-phase must-link anticlustering
  anticlustering(
    input, 
    K = K, 
    method = "2PML", 
    repetitions = repetitions, 
    standardize = TRUE, 
    must_link = must_link_constraints
  )
}

# ------------------------------------------------------------------------------
# fast_anticlustering_shiny(numeric_vars, categorical_vars, K)
# ------------------------------------------------------------------------------
# Uses a "fast" approach for large data sets, employing k-plus or k-means logic
# for anticlustering. This function is invoked when N > 500.
#
# Args:
#   numeric_vars: A data frame of numeric features.
#   categorical_vars: A data frame of categorical features (converted to dummies).
#   K: Number of batches to form.
#
# Returns:
#   A vector of length N specifying each row's batch assignment.
# ------------------------------------------------------------------------------
fast_anticlustering_shiny <- function(numeric_vars, categorical_vars, K) {
  # Convert categorical variables to dummy columns
  if (!is.null(categorical_vars)) {
    categorical_vars <- categories_to_binary(categorical_vars)
  }
  
  # For numeric variables, use 'kplus_moment_variables' to generate features
  # for the k-plus criterion, otherwise use k-means criterion
  if (!is.null(numeric_vars)) {
    numeric_vars <- kplus_moment_variables(numeric_vars, 2, FALSE)
    criterion <- "k-plus"
  } else {
    criterion <- "k-means"
  }
  
  # Combine numeric and dummy-coded categorical features
  input <- cbind(numeric_vars, categorical_vars)
  
  # Standardize the input for better clustering performance
  input <- scale(input)
  N <- nrow(input)
  
  # For extremely large data sets (> 10000), we reduce overhead with special logic
  if (N > 10000) {
    nn_method <- ifelse(N > 100000, "random", "RANN")
    message(
      "N = ", N, ". Using ", criterion, 
      " anticlustering with fast exchange method (100 exchange partners, selection of exchange partners using method ", 
      nn_method, ").\n"
    )
    return(
      fast_anticlustering(
        input, 
        K = K,
        exchange_partners = generate_exchange_partners(
          100, 
          N = N, 
          features = numeric_vars, 
          method = nn_method
        )
      )
    )
  } else {
    message("N = ", N, ". Using ", criterion, " anticlustering with exchange method.")
    return(fast_anticlustering(input, K = K))
  }
}

################################################################################
#                                   UI
################################################################################
ui <- fluidPage(
  # --------------------------------------------------------------------------
  # Custom CSS & JavaScript
  # --------------------------------------------------------------------------
  tags$head(
    # Inline CSS for styling various elements
    tags$style(HTML("
      /* Tabs Styling */
      .nav-tabs > li > a {
        font-weight: bold;
        color: #333;
        border-radius: 6px 6px 0 0;
        transition: 0.3s ease-in-out;
      }
      
      .nav-tabs > li.active > a, 
      .nav-tabs > li.active > a:focus, 
      .nav-tabs > li.active > a:hover {
        background: linear-gradient(to bottom, #e3f2fd, #bbdefb);
        border-color: #90caf9;
        color: #1565c0;
        font-weight: bold;
      }
      
      .nav-tabs > li > a:hover {
        background: #f1f8ff;
        border-color: #90caf9;
      }
      
      /* Numeric Inputs Styling */
      .form-control {
        border-radius: 6px;
        border: 1px solid #90caf9;
        box-shadow: 1px 1px 4px rgba(0, 0, 0, 0.1);
      }
      
      /* Feature Selection Styling */
      .selectize-input {
        background: #f8f9fa;
        border-radius: 6px;
        border: 1px solid #90caf9;
      }
      
      /* Headers and Summary Tables */
      h4 {
        text-align: center;
        font-weight: bold;
        color: #1f5bbf;
        margin-bottom: 10px;
      }
      
      /* Center align all table headers */
      .table thead th {
        text-align: center !important;
      }

      .summary-container {
        text-align: center;
        padding: 10px;
      }
      
      /* Divider for better structure */
      .divider {
        border-top: 2px solid #90caf9;
        margin-top: 10px;
        margin-bottom: 10px;
      }
      
      /* Footer */
      .footer {
        background-color: #f8f9fa; 
        padding: 10px; 
        border-top: 1px solid #dee2e6; 
        text-align: center;
        font-size: 14px;
        color: #555;
      }
      
      /* Custom Gradient Button for Run Analysis */
      .btn-custom {
        background: linear-gradient(to bottom, #78d1eb, #1f5bbf);
        color: white;
        border-radius: 8px;
        font-size: 18px;
        padding: 14px 28px;
        border: none;
        transition: 0.3s;
        font-weight: bold;
        box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.2);
      }
      
      .btn-custom:hover {
        background: linear-gradient(to bottom, #5bb6d9, #154a9b);
        transform: scale(1.05);
      }

      /* Bold headings for dropdowns, checkboxes, radio buttons */
      .control-label,
      .checkbox label,
      .radio label {
        font-weight: bold !important;
        color: #333 !important;
        font-size: 14px !important;
      }
    ")),
    # JavaScript snippet to activate tooltips
    tags$script(HTML("
      $(function () {
        $('[title]').tooltip({
          placement: 'auto',
          boundary: 'viewport'
        });
      });
    "))
  ),
  
  # --------------------------------------------------------------------------
  # Application Header with Logo
  # --------------------------------------------------------------------------
  div(class = "top-section", 
      img(src = "anticlust logo.png", width = "100%", height = "500")
  ),
  
  # --------------------------------------------------------------------------
  # Main Tabset Panel
  # --------------------------------------------------------------------------
  tabsetPanel(
    
    # ========================================================================
    #                           Synthetic Data Tab
    # ========================================================================
    tabPanel("Synthetic Data",
             br(),
             # -----------------------------------------------------------------
             # Inputs for Synthetic Data
             # -----------------------------------------------------------------
             fluidRow(
               column(3, div(
                 numericInput("batches", "Number of Batches (K)", 
                              value = 20, min = 2, step = 1)
               )),
               column(4, div(
                 selectInput("features", "Select Features", 
                             choices = c("endo", "site", "phase", "stage", "age"), 
                             selected = c("endo", "site", "phase", "stage", "age"), 
                             multiple = TRUE)
               )),
               column(5, div(
                 tags$div("Select Feature Type", 
                          style = "color: #333; font-size: 14px; font-weight: bold; margin-bottom: 10px;"),
                 uiOutput("varTypeSelection")
               ))
             ),
             br(),
             fluidRow(
               column(3, div(
                 title = "Check to activate must link constraints (forces items with the same value to be in the same batch).",
                 checkboxInput("useMustLink", "Use Must Link Feature", value = FALSE)
               )),
               column(4, 
                      # Only show must-link variable selection if 'Use Must Link' is checked
                      conditionalPanel(
                        condition = "input.useMustLink == true",
                        div(
                          title = "Select the variable(s) to be used for must link grouping (e.g., patient_id).",
                          selectInput("must_link_vars", "Select Must Link Variable(s)", 
                                      choices = c("patient_id"), 
                                      selected = "patient_id", multiple = TRUE)
                        )
                      )
               ),
               # Button to run analysis
               column(3, offset = 2, div(
                 actionButton("runBtn", "Run Analysis", icon = icon("play"), 
                              class = "btn-custom", width = "100%")
               ))
             ),
             br(),
             
             # -----------------------------------------------------------------
             # Outputs for Synthetic Data
             # -----------------------------------------------------------------
             conditionalPanel(
               condition = "input.useMustLink == false",
               fluidRow(
                 column(12, h4("Anticlust Summary Table"), 
                        # Loading spinner around the summary table
                        withSpinner(uiOutput("summaryTable")))
               )
             ),
             conditionalPanel(
               condition = "input.useMustLink == true",
               fluidRow(
                 column(12, h4("Must-Link Constrained Anticlust Summary Table"), 
                        withSpinner(uiOutput("mlSummaryTable")))
               )
             ),
             fluidRow(
               column(12, h4("Random Assignment Summary Table"), 
                      withSpinner(uiOutput("randomSummaryTable")))
             ),
             
             # -----------------------------------------------------------------
             # Bold GitHub Repository Link
             # -----------------------------------------------------------------
             fluidRow(
               column(12,
                      br(),
                      tags$p(
                        tags$strong("For further details and to access the complete code and dataset, please visit our GitHub repository: "),
                        tags$a(
                          href = "https://github.com/m-Py/must-link-anticlustering/tree/main/Web%20app",
                          "https://github.com/m-Py/must-link-anticlustering/tree/main/Web%20app",
                          target = "_blank"
                        )
                      )
               )
             )
    ),
    
    # ========================================================================
    #                           Upload Dataset Tab
    # ========================================================================
    tabPanel("Upload Dataset",
             br(),
             # -----------------------------------------------------------------
             # Inputs for Upload Dataset
             # -----------------------------------------------------------------
             fluidRow(
               column(3, div(
                 fileInput("file", "Upload CSV/Excel File", accept = c(".csv", ".xlsx"))
               )),
               column(4, div(
                 uiOutput("feature_select")
               )),
               column(5, div(
                 tags$div("Select Feature Type", 
                          style = "color: #333; font-size: 14px; font-weight: bold; margin-bottom: 10px;"),
                 uiOutput("varTypeSelectionUpload")
               ))
             ),
             br(),
             fluidRow(
               column(3, div(
                 numericInput("batches_upload", "Number of Batches (K)", 
                              value = 20, min = 2, step = 1)
               )),
               column(3, div(
                 title = "Check to use must link constraints on the uploaded dataset.",
                 checkboxInput("useMustLink_upload", "Use Must Link Feature", value = FALSE)
               )),
               column(4, conditionalPanel(
                 condition = "input.useMustLink_upload == true",
                 div(
                   title = "Select the variable(s) to be used for must linking.",
                   uiOutput("must_link_selection_upload")
                 )
               )),
               column(2, div(
                 actionButton("runUploadBtn", "Run Analysis", icon = icon("play"), 
                              class = "btn-custom", width = "100%")
               ))
             ),
             br(),
             
             # -----------------------------------------------------------------
             # Outputs for Upload Dataset
             # -----------------------------------------------------------------
             conditionalPanel(
               condition = "input.useMustLink_upload == false",
               fluidRow(
                 column(12, h4("Assigned Batches Table (No Must-Link)"), 
                        withSpinner(DT::dataTableOutput("batchTable")))
               )
             ),
             conditionalPanel(
               condition = "input.useMustLink_upload == true",
               fluidRow(
                 column(12, h4("Assigned Batches Table (Must-Link)"), 
                        withSpinner(DT::dataTableOutput("batchTableML")))
               )
             ),
             br(),
             div(
               downloadButton("downloadData", "Download Results", class = "btn-custom")
             )
    )
  ),
  
  # --------------------------------------------------------------------------
  # Footer
  # --------------------------------------------------------------------------
  div(class = "footer", p("Anticlust Analysis WebApp by UCSF and University of Düsseldorf"))
)

################################################################################
#                                   SERVER
################################################################################
server <- function(input, output, session) {
  
  ##############################################################################
  #                           Synthetic Data Tab
  ##############################################################################
  
  # --------------------------------------------------------------------------
  # dataset_synthetic() - Reactive that loads the synthetic dataset from a CSV
  # --------------------------------------------------------------------------
  dataset_synthetic <- reactive({
    # Replace with the path to your local CSV, if necessary
    read.csv("~/Desktop/UCSF Lab Work/Anticlust App/synthetic_dataset_20250211.csv", 
             stringsAsFactors = FALSE)
  })
  
  # --------------------------------------------------------------------------
  # varTypeSelection - Dynamically generated radio buttons for each selected feature
  # --------------------------------------------------------------------------
  output$varTypeSelection <- renderUI({
    req(input$features)
    tagList(
      lapply(input$features, function(var) {
        fluidRow(
          style = "margin-bottom: 5px;",
          column(6, tags$label(var, style = "font-weight: normal;")),
          column(6,
                 radioButtons(
                   inputId = paste0("varType_", var),
                   label = NULL,
                   choices = c("Categorical", "Numerical"),
                   inline = TRUE,
                   selected = if (var == "age") "Numerical" else "Categorical"
                 )
          )
        )
      })
    )
  })
  
  # --------------------------------------------------------------------------
  # analysisResult - Event triggered by the 'Run Analysis' button (Synthetic Data)
  #                  Performs anticlustering with or without must-link constraints.
  # --------------------------------------------------------------------------
  analysisResult <- eventReactive(input$runBtn, {
    dataset <- dataset_synthetic()
    req(dataset, input$features)
    
    # Validate user-chosen variable types vs. actual data
    mismatches <- c()
    for (v in input$features) {
      actual_is_numeric <- is.numeric(dataset[[v]])
      user_selection <- input[[paste0("varType_", v)]]
      if (actual_is_numeric && user_selection == "Categorical") {
        mismatches <- c(mismatches, paste0("'", v, "' is numeric in the data, but you selected it as Categorical."))
      } else if (!actual_is_numeric && user_selection == "Numerical") {
        mismatches <- c(mismatches, paste0("'", v, "' is non-numeric in the data, but you selected it as Numerical."))
      }
    }
    # If mismatches exist, show a validation error in the UI
    if (length(mismatches) > 0) {
      validate(
        need(FALSE, paste("Type selection mismatch for the following feature(s):", 
                          paste(mismatches, collapse = " | ")))
      )
    }
    
    # Identify which features are categorical
    catVars <- sapply(input$features, function(v) input[[paste0("varType_", v)]])
    catVars <- names(catVars)[catVars == "Categorical"]
    
    # Prepare data for summary tables by converting categorical to dummy columns
    processed_data <- list()
    for (v in input$features) {
      if (v %in% catVars) {
        processed_data[[v]] <- categories_to_binary(dataset[, v, drop = FALSE])
      } else {
        processed_data[[v]] <- dataset[[v]]
      }
    }
    combined_data <- do.call(cbind, processed_data)
    
    # Compute squared Euclidean distance for local-maximum anticlustering
    input_dist <- dist(combined_data)^2
    
    # Perform normal anticlustering (no must-link)
    dataset$BatchAnticlust <- anticlustering(
      input_dist, K = input$batches, 
      method = "local-maximum", 
      repetitions = 50
    )
    
    # Perform must-link anticlustering if user checked 'Use Must Link'
    if (isTRUE(input$useMustLink)) {
      req(input$must_link_vars)
      must_link_data <- dataset[, input$must_link_vars, drop = FALSE]
      # If multiple must-link columns, combine them via interaction()
      if (ncol(must_link_data) > 1) {
        must_link_factor <- interaction(must_link_data, drop = TRUE)
      } else {
        must_link_factor <- must_link_data[[1]]
      }
      dataset$BatchAnticlustML <- anticlustering(
        input_dist, K = input$batches, 
        method = "2PML", 
        repetitions = 50, 
        must_link = must_link_factor
      )
    }
    
    # Create a random assignment for baseline comparison
    dataset$BatchRandom <- sample(dataset$BatchAnticlust)
    
    # Generate summary tables with tableone
    tabAnticlust <- CreateTableOne(
      vars = input$features,
      strata = "BatchAnticlust",
      data = dataset, 
      factorVars = catVars
    )
    summary_output_anticlust <- knitr::kable(
      print(tabAnticlust, smd = TRUE), 
      format = "html", 
      table.attr = "class='table table-striped'"
    )
    
    summary_output_anticlustML <- NULL
    if (isTRUE(input$useMustLink)) {
      tabAnticlustML <- CreateTableOne(
        vars = input$features,
        strata = "BatchAnticlustML",
        data = dataset, 
        factorVars = catVars
      )
      summary_output_anticlustML <- knitr::kable(
        print(tabAnticlustML, smd = TRUE), 
        format = "html", 
        table.attr = "class='table table-striped'"
      )
    }
    
    tabRandom <- CreateTableOne(
      vars = input$features,
      strata = "BatchRandom",
      data = dataset, 
      factorVars = catVars
    )
    summary_output_random <- knitr::kable(
      print(tabRandom, smd = TRUE), 
      format = "html", 
      table.attr = "class='table table-striped'"
    )
    
    # Return a list with the summary tables for use in the UI
    list(
      summary_anticlust   = summary_output_anticlust,
      summary_anticlustML = summary_output_anticlustML,
      summary_random      = summary_output_random
    )
  })
  
  # --------------------------------------------------------------------------
  # UI Outputs for Synthetic Data
  # --------------------------------------------------------------------------
  
  # Render the normal anticlustering summary table
  output$summaryTable <- renderUI({
    req(analysisResult())
    HTML(analysisResult()$summary_anticlust)
  })
  
  # Render the must-link anticlustering summary table
  output$mlSummaryTable <- renderUI({
    req(analysisResult())
    HTML(analysisResult()$summary_anticlustML)
  })
  
  # Render the random assignment summary table
  output$randomSummaryTable <- renderUI({
    req(analysisResult())
    HTML(analysisResult()$summary_random)
  })
  
  ##############################################################################
  #                           Upload Dataset Tab
  ##############################################################################
  
  # --------------------------------------------------------------------------
  # dataset_reactive() - Loads the user-uploaded CSV or Excel file
  # --------------------------------------------------------------------------
  dataset_reactive <- reactive({
    req(input$file)
    ext <- tools::file_ext(input$file$datapath)
    
    if (ext == "csv") {
      read_csv(input$file$datapath)
    } else if (ext == "xlsx") {
      read_excel(input$file$datapath)
    } else {
      return(NULL)
    }
  })
  
  # --------------------------------------------------------------------------
  # feature_select - Renders a UI selectInput() for choosing features from
  #                  the uploaded dataset
  # --------------------------------------------------------------------------
  output$feature_select <- renderUI({
    req(dataset_reactive())
    selectInput("features_upload", "Select Features", 
                choices = names(dataset_reactive()), multiple = TRUE)
  })
  
  # --------------------------------------------------------------------------
  # varTypeSelectionUpload - Renders radio buttons for each selected feature
  #                          in the uploaded data, to specify Cat vs Num
  # --------------------------------------------------------------------------
  output$varTypeSelectionUpload <- renderUI({
    req(input$features_upload)
    tagList(
      lapply(input$features_upload, function(var) {
        fluidRow(
          style = "margin-bottom: 5px;",
          column(6, tags$label(var, style = "font-weight: normal;")),
          column(6,
                 radioButtons(
                   inputId = paste0("varType_upload_", var),
                   label = NULL,
                   choices = c("Categorical", "Numerical"),
                   inline = TRUE,
                   selected = "Numerical"
                 )
          )
        )
      })
    )
  })
  
  # --------------------------------------------------------------------------
  # must_link_selection_upload - UI for must-link variable selection in 
  #                              uploaded data
  # --------------------------------------------------------------------------
  output$must_link_selection_upload <- renderUI({
    req(dataset_reactive())
    selectInput("must_link_vars_upload", "Select Must Link Variable(s)", 
                choices = names(dataset_reactive()), selected = NULL, multiple = TRUE)
  })
  
  # --------------------------------------------------------------------------
  # analysisResultUpload - Event triggered by 'Run Analysis' button (Upload Tab)
  #                        Uses anticlustering_shiny() with or without must-link.
  # --------------------------------------------------------------------------
  analysisResultUpload <- eventReactive(input$runUploadBtn, {
    dataset <- dataset_reactive()
    req(dataset, input$features_upload)
    
    # Validate user-chosen variable types vs. actual data
    mismatches <- c()
    for (v in input$features_upload) {
      actual_is_numeric <- is.numeric(dataset[[v]])
      user_selection <- input[[paste0("varType_upload_", v)]]
      if (actual_is_numeric && user_selection == "Categorical") {
        mismatches <- c(mismatches, paste0("'", v, "' is numeric in the data, but you selected it as Categorical."))
      } else if (!actual_is_numeric && user_selection == "Numerical") {
        mismatches <- c(mismatches, paste0("'", v, "' is non-numeric in the data, but you selected it as Numerical."))
      }
    }
    if (length(mismatches) > 0) {
      validate(
        need(FALSE, paste("Type selection mismatch for the following feature(s):", 
                          paste(mismatches, collapse = " | ")))
      )
    }
    
    # Separate numeric and categorical variables based on user selection
    feature_types <- sapply(input$features_upload, function(v) input[[paste0("varType_upload_", v)]])
    numeric_vars <- dataset[, names(feature_types)[feature_types == "Numerical"], drop = FALSE]
    categorical_vars <- dataset[, names(feature_types)[feature_types == "Categorical"], drop = FALSE]
    
    # Perform anticlustering using the new wrapper function
    if (isTRUE(input$useMustLink_upload)) {
      # If must-link is checked, gather the must-link variables
      req(input$must_link_vars_upload)
      must_link_data_upload <- dataset[, input$must_link_vars_upload, drop = FALSE]
      if (ncol(must_link_data_upload) > 1) {
        must_link_factor_upload <- interaction(must_link_data_upload, drop = TRUE)
      } else {
        must_link_factor_upload <- must_link_data_upload[[1]]
      }
      # Must-link assignment
      batch_assignment <- anticlustering_shiny(
        numeric_vars = numeric_vars,
        categorical_vars = categorical_vars,
        must_link_constraints = must_link_factor_upload,
        K = input$batches_upload
      )
      dataset$BatchAnticlustML <- batch_assignment
    } else {
      # No must-link assignment
      batch_assignment <- anticlustering_shiny(
        numeric_vars = numeric_vars,
        categorical_vars = categorical_vars,
        must_link_constraints = NULL,
        K = input$batches_upload
      )
      dataset$BatchAnticlust <- batch_assignment
    }
    
    # Return the dataset with new batch assignment columns
    dataset
  })
  
  # --------------------------------------------------------------------------
  # UI Outputs for the Uploaded Data Tab
  # --------------------------------------------------------------------------
  
  # Show assignment table if 'Use Must Link' is NOT checked
  output$batchTable <- DT::renderDataTable({
    req(analysisResultUpload())
    req(!isTRUE(input$useMustLink_upload))
    data_noML <- analysisResultUpload()
    # Remove the must-link column if it exists
    if ("BatchAnticlustML" %in% names(data_noML)) {
      data_noML$BatchAnticlustML <- NULL
    }
    data_noML
  })
  
  # Show assignment table if 'Use Must Link' IS checked
  output$batchTableML <- DT::renderDataTable({
    req(analysisResultUpload())
    req(isTRUE(input$useMustLink_upload))
    data_ML <- analysisResultUpload()
    # Remove the non-must-link column if it exists
    if ("BatchAnticlust" %in% names(data_ML)) {
      data_ML$BatchAnticlust <- NULL
    }
    data_ML
  })
  
  # --------------------------------------------------------------------------
  # downloadData - Allows user to download the final data with batch assignments
  # --------------------------------------------------------------------------
  output$downloadData <- downloadHandler(
    filename = function() { "assigned_batches.csv" },
    content = function(file) {
      final_data <- analysisResultUpload()
      
      # If must-link was used, remove the normal batch column, and vice versa
      if (isTRUE(input$useMustLink_upload)) {
        if ("BatchAnticlust" %in% names(final_data)) {
          final_data$BatchAnticlust <- NULL
        }
      } else {
        if ("BatchAnticlustML" %in% names(final_data)) {
          final_data$BatchAnticlustML <- NULL
        }
      }
      
      # Write to CSV
      write.csv(final_data, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)

