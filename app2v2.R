# https://upenn.shinyapps.io/rmpc_app2_model_v10/
# 7/2/2026, by Daohang Sha, ver2026.07.02.1.0
# 7/3/2026, by Daohang Sha, ver2026.07.02.1.0.1 added AIC/BIC
# 7/13/2026, to publish this app
# 7/30/2026, updated adjusted means plots - replaced box plots with color regression lines w/raw data points for ANVOVA/GLM, added visit count ver2.0

# install.packages('rsconnect')
# rsconnect::setAccountInfo(name='upenn', token='AAB7CC53CDC0E1F588FF9F04414C5911', secret='cekxUWs1sUpSlwqGpPvqMTGLwBFztZI9fHBbxLAc')
# 
# library(rsconnect)
# rsconnect::deployApp('C:/COREY/Apps/RMPC_app2_model_v1.0.1')


library(shiny)
library(bslib)
library(tidyverse)
library(emmeans)
library(broom)
library(DT)



# ── UI ────────────────────────────────────────────────────────────────────────

ui_app2 <- page_navbar(
  title = "RMPC Data Analysis App2",
  theme = bs_theme(bootswatch = "cerulean"),   #all other available themes (e.g., "flatly", "cerulean", "cosmo", "darkly", "cyborg", etc.)
  fillable = FALSE,
  footer = div(
    style = "text-align: center; padding: 10px; font-size: 0.8rem; color: #6c757d;",
    HTML("&copy; 2026 Corey. All rights reserved."),
    br(),
    textOutput("visit_text")
  ),
  
  header = tags$head(
    HTML("
      <!-- Global site tag (gtag.js) - Google Analytics -->
      <script async src='https://www.googletagmanager.com/gtag/js?id=G-RH1P6ZMNBZ'></script>
      <script>
        window.dataLayer = window.dataLayer || [];
        function gtag(){dataLayer.push(arguments);}
        gtag('js', new Date());
        gtag('config', 'G-RH1P6ZMNBZ');
      </script>
    ")
  ),
  
  

  nav_panel(
    "Analysis v2.0",
    layout_sidebar(
      fillable = FALSE,

      # ── Sidebar ─────────────────────────────────────────────────────────────
      sidebar = sidebar(
        width = 310,

        card(
          card_header("1. Load Data"),
          fileInput("file", NULL, accept = ".csv",
                    buttonLabel = "Browse CSV…",
                    placeholder = "No file selected")
        ),

        card(
          card_header("2. Column Mapping"),
          uiOutput("ui_outcome"),
          uiOutput("ui_treatment"),
          uiOutput("ui_baseline_col"),
          uiOutput("ui_gender_col")
        ),

        card(
          card_header("3. Model Type"),
          radioButtons("model_type", NULL,
            choices = c(
              "ANOVA  (no body mass)"             = "anova",
              "ANCOVA  (body mass as covariate)"  = "ancova",
              "GLM  (body mass × treatment)"      = "glm_int"
            )
          ),
          conditionalPanel(
            "input.model_type !== 'anova'",
            uiOutput("ui_bodymass_col")
          ),
          hr(style = "margin: 6px 0;"),
          strong("Additional Covariates"),
          uiOutput("ui_cb_baseline"),
          uiOutput("ui_cb_gender")
        ),

        actionButton("run", "Run Analysis",
                     class = "btn-primary w-100",
                     icon = icon("play-circle"))
      ),

      # ── Main panel ──────────────────────────────────────────────────────────
      uiOutput("formula_banner"),

      navset_card_tab(
        id = "main_tabs",

        nav_panel(
          "Data Preview",
          DTOutput("out_data_preview")
        ),

        nav_panel(
          "Model Summary",
          verbatimTextOutput("out_summary")
        ),

        nav_panel(
          "ANOVA Results",
          DTOutput("out_anova")
        ),
        
        nav_panel(
          "Model Fit Statistics",
          DTOutput("out_fitstats")
        ),
        
        nav_panel(
          "Marginal Means & Contrasts",
          h5("Estimated Marginal Means (adjusted for covariates)"),
          DTOutput("out_emmeans"),
          br(),
          h5("Pairwise Contrasts (Tukey adjustment)"),
          DTOutput("out_contrasts")
        ),

        nav_panel(
          "Diagnostic Plots",
          layout_columns(
            col_widths = c(6, 6),
            card(card_header("Raw Data by Group"),
                 plotOutput("plot_raw",   height = "300px")),
            
            # card(card_header("Adjusted Means ± 95% CI"),
            #      plotOutput("plot_emm",   height = "300px")),
            
            # 7/30/2026, replace above chunk
            card(
              card_header(textOutput("plot_emm_title")),
              plotOutput("plot_emm", height = "300px")
            ),
            
            card(card_header("Residuals vs Fitted"),
                 plotOutput("plot_resid", height = "300px")),
            card(card_header("Normal Q-Q"),
                 plotOutput("plot_qq",    height = "300px"))
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server_app2 <- function(input, output, session) {
  # 7/30/2026 
  # Global-like counter in this R session
  if (!exists("visit_counter", envir = .GlobalEnv)) {
    assign("visit_counter", 0L, envir = .GlobalEnv)
  }
  visit_counter <<- visit_counter + 1L
  this_visit <- visit_counter

  output$visit_text <- renderText({
    paste("Visit count (since last restart):", this_visit)
  })
  
  
  # ── Reactive data ────────────────────────────────────────────────────────
  raw_data <- reactive({
    req(input$file)
    read.csv(input$file$datapath, stringsAsFactors = FALSE)
  })

  all_cols     <- reactive({ req(raw_data()); names(raw_data()) })
  numeric_cols <- reactive({
    req(raw_data())
    names(raw_data())[sapply(raw_data(), is.numeric)]
  })

  # ── Column-mapping UI ────────────────────────────────────────────────────

  output$ui_outcome <- renderUI({
    nc <- numeric_cols()
    detected <- grep("^before$|^outcome$|^response$|^value$|^si13c$", nc,
                     ignore.case = TRUE, value = TRUE)
    selected <- if (length(detected) > 0) detected[1] else nc[length(nc)]
    selectInput("outcome", "Outcome Variable", choices = nc, selected = selected)
  })

  output$ui_treatment <- renderUI({
    ac <- all_cols()
    detected <- grep("condition|treatment|group", ac,
                     ignore.case = TRUE, value = TRUE)
    selected <- if (length(detected) > 0) detected[1] else ac[1]
    selectInput("treatment", "Treatment / Group", choices = ac, selected = selected)
  })

  output$ui_baseline_col <- renderUI({
    nc <- numeric_cols()
    choices  <- c("(None)" = "", nc)
    # Match "baseline" but not "before" to avoid capturing the outcome
    detected <- grep("baseline|^pre$|^t0$", nc, ignore.case = TRUE, value = TRUE)
    selected <- if (length(detected) > 0) detected[1] else ""
    selectInput("baseline_col", "Baseline Column",
                choices = choices, selected = selected)
  })

  output$ui_gender_col <- renderUI({
    ac <- all_cols()
    choices  <- c("(None)" = "", ac)
    detected <- grep("^gender$|^sex$", ac, ignore.case = TRUE, value = TRUE)
    selected <- if (length(detected) > 0) detected[1] else ""
    selectInput("gender_col", "Gender Column",
                choices = choices, selected = selected)
  })

  output$ui_bodymass_col <- renderUI({
    nc <- numeric_cols()
    lean  <- grep("\\blean\\b",  nc, ignore.case = TRUE, value = TRUE)
    fat  <- grep("\\bfat\\b",  nc, ignore.case = TRUE, value = TRUE)
    total <- grep("\\btotal\\b", nc, ignore.case = TRUE, value = TRUE)
    default <- grep("\\bbodymass\\b", nc, ignore.case = TRUE, value = TRUE)
    
    # Restrict to only lean / total columns
    bm_choices <- unique(c(lean, fat, total, default))
    selected   <- if (length(lean) > 0) lean[1] else if (length(bm_choices) > 0) bm_choices[1] else character(0)

    radioButtons("bodymass_col", "Body Mass Variable",
                 choices = if (length(bm_choices) > 0) bm_choices else "(none available)",
                 selected = selected)
  })

  # Checkbox for baseline (only when a baseline column is selected)
  output$ui_cb_baseline <- renderUI({
    if (!is.null(input$baseline_col) && nzchar(input$baseline_col)) {
      checkboxInput("use_baseline",
                    paste("Include:", input$baseline_col), value = FALSE)
    } else {
      p(class = "text-muted small mt-1",
        "Select a baseline column above to enable.")
    }
  })

  # Checkbox for gender (only when a gender column is selected)
  output$ui_cb_gender <- renderUI({
    if (!is.null(input$gender_col) && nzchar(input$gender_col)) {
      checkboxInput("use_gender",
                    paste("Include:", input$gender_col), value = FALSE)
    } else {
      p(class = "text-muted small mt-1",
        "Select a gender column above to enable.")
    }
  })
  
  
  # 7/30/2026, add
  
  output$plot_emm_title <- renderText({
    req(input$model_type)
    if (input$model_type == "anova") {
      "Adjusted Means ± 95% CI"
    } else if (input$model_type == "ancova") {
      "ANCOVA: Outcome vs Body Mass"
    } else if (input$model_type == "glm_int") {
      "GLM: Outcome vs Body Mass (Interaction)"
    } else {
      "Model Plot"
    }
  })

  # ── Formula construction ─────────────────────────────────────────────────
  formula_str <- reactive({
    req(input$outcome, input$treatment, input$model_type)

    outcome   <- input$outcome
    treatment <- input$treatment

    # Collect optional covariates (baseline, gender)
    covs <- c()
    if (isTRUE(input$use_baseline) &&
        !is.null(input$baseline_col) && nzchar(input$baseline_col)) {
      covs <- c(covs, input$baseline_col)
    }
    if (isTRUE(input$use_gender) &&
        !is.null(input$gender_col) && nzchar(input$gender_col)) {
      covs <- c(covs, input$gender_col)
    }

    rhs <- switch(input$model_type,

      # ANOVA: outcome ~ [covs +] treatment
      "anova" = paste(c(covs, treatment), collapse = " + "),

      # ANCOVA: outcome ~ bodymass + [covs +] treatment
      "ancova" = {
        req(input$bodymass_col)
        paste(c(input$bodymass_col, covs, treatment), collapse = " + ")
      },

      # GLM interaction: outcome ~ [covs +] bodymass * treatment
      "glm_int" = {
        req(input$bodymass_col)
        bm_x_trt <- paste0(input$bodymass_col, " * ", treatment)
        paste(c(covs, bm_x_trt), collapse = " + ")
      }
    )

    paste(outcome, "~", rhs)
  })

  # ── Formula banner ───────────────────────────────────────────────────────
  output$formula_banner <- renderUI({
    tryCatch({
      fs <- formula_str()
      card(
        style = "margin-bottom: 12px; background: #eef2f7; border-left: 4px solid #2c8c99;",
        card_body(
          padding = "8px 14px",
          tags$span(tags$b("Model formula: "), tags$code(fs))
        )
      )
    }, error = function(e) NULL)
  })

  # ── Prepared data (factors applied) ──────────────────────────────────────
  prepped_data <- reactive({
    req(raw_data(), input$treatment)
    df <- raw_data()
    df[[input$treatment]] <- as.factor(df[[input$treatment]])
    if (!is.null(input$gender_col) && nzchar(input$gender_col) &&
        isTRUE(input$use_gender)) {
      df[[input$gender_col]] <- as.factor(df[[input$gender_col]])
    }
    df
  })

  # ── Model fit ────────────────────────────────────────────────────────────
  fitted_model <- eventReactive(input$run, {
    req(prepped_data(), formula_str())
    tryCatch(
      lm(as.formula(formula_str()), data = prepped_data()),
      error = function(e) {
        showNotification(paste("Model error:", e$message), type = "error", duration = 8)
        NULL
      }
    )
  })

  # Estimated marginal means
  emm_obj <- reactive({
    req(fitted_model(), input$treatment)
    tryCatch(
      emmeans(fitted_model(), specs = input$treatment),
      error = function(e) {
        showNotification(paste("emmeans error:", e$message), type = "warning")
        NULL
      }
    )
  })

  # ── Output: data preview ─────────────────────────────────────────────────
  output$out_data_preview <- renderDT({
    req(raw_data())
    datatable(raw_data(),
              options  = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE,
              class    = "compact stripe")
  })

  # ── Output: model summary ────────────────────────────────────────────────
  output$out_summary <- renderPrint({
    req(fitted_model())
    summary(fitted_model())
  })

  # ── Output: ANOVA table ──────────────────────────────────────────────────
  output$out_anova <- renderDT({
    req(fitted_model())

    tbl <- tidy(anova(fitted_model())) |>
      rename(
        Term      = term,
        df        = df,
        SS        = sumsq,
        MS        = meansq,
        F         = statistic,
        `p-value` = p.value
      ) |>
      mutate(
        across(c(SS, MS, F), ~round(., 3)),
        `p-value` = case_when(
          is.na(`p-value`) ~ "—",
          `p-value` < 0.001 ~ "< 0.001",
          TRUE ~ as.character(round(`p-value`, 4))
        )
      )

    datatable(tbl,
              options  = list(dom = "t", pageLength = 20),
              rownames = FALSE,
              class    = "compact stripe") |>
      formatStyle("p-value",
                  color = styleInterval(c(0), c("black", "black")))
  })

  
  # v1.0.1, 7/3/2026, added Model Fit Statistics Table
  
  output$out_fitstats <- renderDT({
    req(fitted_model())
    
    fm <- fitted_model()
    
    datatable(
      tibble(
        Statistic = c("AIC", "BIC"),
        Value = round(c(AIC(fm), BIC(fm)), 3)
      ),
      options = list(dom = "t"),
      rownames = FALSE
    )
  })
  
  
  # ── Output: marginal means ───────────────────────────────────────────────
  output$out_emmeans <- renderDT({
    req(emm_obj())
    as.data.frame(emm_obj()) |>
      mutate(across(where(is.numeric), ~round(., 3))) |>
      datatable(options  = list(dom = "t"),
                rownames = FALSE,
                class    = "compact stripe")
  })

  # ── Output: pairwise contrasts ───────────────────────────────────────────
  output$out_contrasts <- renderDT({
    req(emm_obj())
    pairs(emm_obj(), adjust = "tukey") |>
      as.data.frame() |>
      mutate(
        across(where(is.numeric), ~round(., 4)),
        p.value = case_when(
          p.value < 0.001 ~ "< 0.001",
          TRUE ~ as.character(round(p.value, 4))
        )
      ) |>
      datatable(options  = list(dom = "t", pageLength = 30),
                rownames = FALSE,
                class    = "compact stripe")
  })

  # ── Output: raw data plot ────────────────────────────────────────────────
  # output$plot_raw <- renderPlot({
  #   req(fitted_model(), input$outcome, input$treatment)
  #   df <- prepped_data()
  #   ggplot(df, aes(x = .data[[input$treatment]],
  #                  y = .data[[input$outcome]])) +
  #     geom_boxplot(alpha = 0.35, outlier.shape = NA) +
  #     geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
  #     labs(x = input$treatment, y = input$outcome) +
  #     theme_bw(base_size = 13)
  # })
  
  # 7/30/2026, added color
  output$plot_raw <- renderPlot({
    req(fitted_model(), input$outcome, input$treatment)
    df <- prepped_data()
    
    ggplot(
      df,
      aes(
        x = .data[[input$treatment]],
        y = .data[[input$outcome]],
        color = .data[[input$treatment]]
      )
    ) +
      geom_boxplot(alpha = 0.35, outlier.shape = NA) +
      geom_jitter(width = 0.15, size = 2.5, alpha = 0.85) +
      labs(
        x = input$treatment,
        y = input$outcome,
        color = input$treatment
      ) +
      theme_bw(base_size = 13) +
      theme(legend.position = "none")  # remove if you want the legend
  })
  
  

  # ── Output: adjusted means plot ──────────────────────────────────────────
  # output$plot_emm <- renderPlot({
  #   req(emm_obj())
  #   emm_df  <- as.data.frame(emm_obj())
  #   trt_col <- names(emm_df)[1]
  # 
  #   ggplot(emm_df, aes(x = .data[[trt_col]], y = emmean)) +
  #     geom_point(size = 4) +
  #     geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
  #                   width = 0.2, linewidth = 0.9) +
  #     labs(x = trt_col, y = "Estimated Marginal Mean") +
  #     theme_bw(base_size = 13)
  # })
  
  # 7/30/2026, replace above chunk
  # ── Output: adjusted means plot / body mass fitted lines ─────────────────
  output$plot_emm <- renderPlot({
    req(fitted_model())
    
    # If model is ANOVA, keep the original adjusted means plot
    if (input$model_type == "anova") {
      req(emm_obj())
      emm_df  <- as.data.frame(emm_obj())
      trt_col <- names(emm_df)[1]
      
      ggplot(emm_df, aes(
        x = .data[[trt_col]],
        y = emmean,
        color = .data[[trt_col]])) +
        
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                      width = 0.2, linewidth = 0.9) +
        labs(x = trt_col, y = "Estimated Marginal Mean") +
        theme_bw(base_size = 13)
      
    } else {
      # ANCOVA or GLM: plot raw data + fitted lines vs body mass
      req(input$bodymass_col, input$outcome, input$treatment)
      
      df <- prepped_data()
      
      # Create a grid of body mass values for each treatment to draw fitted lines
      # Use range of observed body mass within each treatment
      bm_var  <- input$bodymass_col
      trt_var <- input$treatment
      y_var   <- input$outcome
      
      # Basic check that body mass exists and is numeric
      validate(
        need(bm_var %in% names(df), "Body mass column not found in data."),
        need(is.numeric(df[[bm_var]]),
             "Selected body mass column must be numeric.")
      )
      
      # New data frame for predictions: sequence over body mass per treatment
      pred_grid <- df |>
        group_by(.data[[trt_var]]) |>
        summarize(
          bm_min = min(.data[[bm_var]], na.rm = TRUE),
          bm_max = max(.data[[bm_var]], na.rm = TRUE),
          .groups = "drop"
        ) |>
        # build sequences for each treatment
        rowwise() |>
        mutate(
          bm_seq = list(seq(bm_min, bm_max, length.out = 50))
        ) |>
        unnest(bm_seq) |>
        rename(!!bm_var := bm_seq) |>
        select(-bm_min, -bm_max)
      
      # Add other covariates (baseline, gender) if used.
      # For plotting fitted lines, we hold them at some representative value:
      #   numeric covariates at their overall mean,
      #   factor covariates at reference (first) level.
      # This keeps the predictions consistent across plotting.
      covars <- c()
      if (isTRUE(input$use_baseline) &&
          !is.null(input$baseline_col) && nzchar(input$baseline_col)) {
        covars <- c(covars, input$baseline_col)
      }
      if (isTRUE(input$use_gender) &&
          !is.null(input$gender_col) && nzchar(input$gender_col)) {
        covars <- c(covars, input$gender_col)
      }
      
      if (length(covars) > 0) {
        for (cv in covars) {
          if (is.numeric(df[[cv]])) {
            # numeric covariate: set to overall mean
            pred_grid[[cv]] <- mean(df[[cv]], na.rm = TRUE)
          } else {
            # factor/character: set to most frequent level
            ref_level <- names(sort(table(df[[cv]]), decreasing = TRUE))[1]
            pred_grid[[cv]] <- ref_level
          }
        }
      }
      
      # Ensure treatment is factor in prediction grid as in model
      pred_grid[[trt_var]] <- factor(pred_grid[[trt_var]],
                                     levels = levels(df[[trt_var]]))
      
      # Get fitted values from the model
      pred_grid$fit <- predict(fitted_model(), newdata = pred_grid)
      
      # Plot raw data + fitted lines
      ggplot(df, aes(x = .data[[bm_var]], y = .data[[y_var]],
                     color = .data[[trt_var]])) +
        # raw data points
        geom_point(alpha = 0.7, size = 2) +
        # fitted lines
        geom_line(data = pred_grid,
                  aes(x = .data[[bm_var]], y = fit,
                      color = .data[[trt_var]]),
                  linewidth = 1) +
        labs(
          x = bm_var,
          y = y_var,
          color = trt_var
        ) +
        theme_bw(base_size = 13)
    }
  })
  
  
  # ── Output: residuals vs fitted ──────────────────────────────────────────
  output$plot_resid <- renderPlot({
    req(fitted_model())
    data.frame(
      Fitted    = fitted(fitted_model()),
      Residuals = residuals(fitted_model())
    ) |>
      ggplot(aes(x = Fitted, y = Residuals)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", se = FALSE,
                  color = "steelblue", linewidth = 0.8) +
      labs(x = "Fitted Values", y = "Residuals") +
      theme_bw(base_size = 13)
  })

  # ── Output: Q-Q plot ─────────────────────────────────────────────────────
  output$plot_qq <- renderPlot({
    req(fitted_model())
    data.frame(resid = residuals(fitted_model())) |>
      ggplot(aes(sample = resid)) +
      stat_qq(size = 2.5, alpha = 0.8) +
      stat_qq_line(color = "red", linewidth = 0.9) +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      theme_bw(base_size = 13)
  })

}

shinyApp(ui=ui_app2, server=server_app2)
