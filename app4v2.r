# ============================================================
# Shiny App: Linear Mixed-Effects Model — Metabolic Analysis
# Updated: Dynamic condition column + group selection
# ============================================================

library(shiny)
library(shinythemes)
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(DT)
library(colourpicker)

# ── Global helpers ────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

KNOWN_LABELS <- c(
  Fcirc_C_perMouse_umol = "Fcirc C per Mouse (\u00b5mol)",
  Vox_perMouse_umol     = "Vox per Mouse (\u00b5mol)",
  si13c                 = "si\u00b9\u00b3C",
  vco2                  = "VCO\u2082",
  vo2                   = "VO\u2082",
  rq                    = "Respiratory Quotient (RQ)"
)

# All recognised condition/treatment column names
CONDITION_CANDIDATES <- c("Condition", "condition", "Treatment", "treatment")

# Columns excluded from outcome variable detection
META_COLS <- c( "file_index","Sample..",
  "mouse_id", "MouseID", "Subject", "ID",
  "Time.Pt", "time_pt", "TimePoint", "timepoint",
  CONDITION_CANDIDATES,
  "Lean", "Fat", "Total",
  "Lean_c", "Fat_c", "Total_c"
)

FALLBACK_PALETTE <- c(
  "#4DBBD5", "#E64B35", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000"
)

# ============================================================
# UI
# ============================================================
ui_app4 <- fluidPage(
  theme = shinytheme("flatly"),

  div(
    style = "padding:15px 0 5px 0;",
    h3(icon("flask"),
       " Linear Mixed-Effects Model \u2014 Metabolic Analysis",
       style = "margin:0;"),
    p("Upload CSV \u2192 configure options \u2192 click Run Analysis",
      style = "color:#666; font-size:13px; margin-top:4px;")
  ),
  hr(style = "margin:5px 0 15px 0;"),

  sidebarLayout(

    # ── Sidebar ───────────────────────────────────────────────
    sidebarPanel(
      width = 3,

      # Data upload
      h5(strong("\U0001F4C1  Data")),
      fileInput("file", NULL, accept = ".csv",
                buttonLabel = "Browse CSV\u2026",
                placeholder = "No file selected"),
      hr(),

      # Outcome variable
      h5(strong("\U0001F4CA  Outcome Variable")),
      uiOutput("outcome_ui"),
      hr(),

      # Time points
      h5(strong("\U0001F551  Time Points")),
      radioButtons(
        "phase", NULL,
        choices = c(
          "After Injection  (CL-0 / CL-10 / CL-20)" = "after",
          "Before Injection (90 min / 120 min)"     = "before",
          "Custom selection"                        = "custom"
        ),
        selected = "after"
      ),
      conditionalPanel(
        condition = "input.phase == 'custom'",
        uiOutput("custom_tp_ui")
      ),
      hr(),

      # ── NEW: Condition column + group selector ─────────────
      h5(strong("\U0001F9EC  Condition Groups")),

      # Step 1: which column is the grouping variable?
      uiOutput("cond_col_ui"),

      # Step 2: which groups to include?
      uiOutput("condition_ui"),
      hr(),

      # Body-mass covariate
      h5(strong("\u2696  Body Mass Covariate")),
      radioButtons(
        "covariate", NULL,
        choices  = c("None"  = "none",
                     "Lean"  = "Lean",
                     "Fat"   = "Fat",
                     "Total" = "Total"),
        selected = "none"
      ),
      helpText(style = "font-size:11px; color:#888;",
               "Fat / Total columns must exist in the CSV."),
      hr(),

      # Per-group color pickers
      h5(strong("\U0001F3A8  Group Colors")),
      uiOutput("color_ui"),
      hr(),

      # Run
      actionButton(
        "run",
        label = tagList(icon("play"), " Run Analysis"),
        class = "btn-success btn-block",
        style = "font-weight:bold; font-size:14px;"
      )
    ),

    # ── Main panel ────────────────────────────────────────────
    mainPanel(
      width = 9,

      uiOutput("status_banner"),

      tabsetPanel(
        id = "tabs",

        tabPanel(
          title = tagList(icon("table"), " Summary Stats"),
          br(),
          downloadButton("dl_sum", "Download CSV",
                         class = "btn-sm btn-default"),
          br(), br(),
          DTOutput("sum_table")
        ),

        tabPanel(
          title = tagList(icon("chart-bar"), " Bar Plot"),
          br(),
          fluidRow(
            column(3, sliderInput("bar_ptsize", "Point size",
                                  0.5, 6, 2.4, 0.1)),
            column(3, sliderInput("bar_alpha",  "Transparency",
                                  0.1, 1, 0.85, 0.05)),
            column(3, sliderInput("bar_jitter", "Jitter width",
                                  0, 0.3, 0.08, 0.01)),
            column(3, sliderInput("dodge_w",    "Dodge width",
                                  0.3, 1.0, 0.75, 0.05))
          ),
          downloadButton("dl_bar", "Download PNG",
                         class = "btn-sm btn-default"),
          br(), br(),
          plotOutput("bar_plot", height = "520px")
        ),

        tabPanel(
          title = tagList(icon("chart-line"), " Line Plot"),
          br(),
          downloadButton("dl_line", "Download PNG",
                         class = "btn-sm btn-default"),
          br(), br(),
          plotOutput("line_plot", height = "520px")
        ),

        tabPanel(
          title = tagList(icon("project-diagram"), " EMM Plot"),
          br(),
          plotOutput("emm_plot", height = "460px")
        ),

        tabPanel(
          title = tagList(icon("cogs"), " Model Summary"),
          br(),
          verbatimTextOutput("model_summary")
        ),

        tabPanel(
          title = tagList(icon("calculator"), " ANOVA"),
          br(),
          verbatimTextOutput("anova_table")
        ),

        tabPanel(
          title = "EMM: Cond | Time",
          br(),
          verbatimTextOutput("emm_cond_time")
        ),

        tabPanel(
          title = "EMM: Time (overall)",
          br(),
          verbatimTextOutput("emm_time_overall")
        ),

        tabPanel(
          title = "EMM: Time | Cond",
          br(),
          verbatimTextOutput("emm_time_cond")
        )
      )
    )
  )
)

# ============================================================
# Server
# ============================================================
server_app4 <- function(input, output, session) {

  # ── 1. Load CSV ──────────────────────────────────────────
  raw_data <- reactive({
    req(input$file)
    tryCatch(
      read.csv(input$file$datapath, stringsAsFactors = FALSE),
      error = function(e) {
        showNotification(paste("Error reading file:", e$message),
                         type = "error", duration = 10)
        NULL
      }
    )
  })

  # ── 2. Outcome variable (auto-detected numeric columns) ──
  output$outcome_ui <- renderUI({
    df <- raw_data()
    if (is.null(df)) {
      choices <- setNames(names(KNOWN_LABELS), KNOWN_LABELS)
    } else {
      num_cols <- names(df)[sapply(df, is.numeric)]
      num_cols <- setdiff(num_cols, META_COLS)
      labels   <- ifelse(num_cols %in% names(KNOWN_LABELS),
                         KNOWN_LABELS[num_cols], num_cols)
      choices  <- setNames(num_cols, labels)
    }
    selectInput("outcome", NULL, choices = choices)
  })

  # ── 3. Condition column selector (NEW) ───────────────────
  # Finds whichever of the four candidate names exist in the CSV.
  # Falls back to all character/factor columns if none match.
  output$cond_col_ui <- renderUI({
    req(raw_data())
    df <- raw_data()

    # Check for candidate column names
    found <- intersect(CONDITION_CANDIDATES, names(df))

    if (length(found) == 0) {
      # Fallback: any character or factor column that isn't an ID/time col
      found <- names(df)[sapply(df, function(x) is.character(x) | is.factor(x))]
      found <- setdiff(found,
                       c("mouse_id", "MouseID", "Subject", "ID",
                         "Time.Pt", "time_pt", "TimePoint", "timepoint"))
    }

    validate(
      need(length(found) > 0,
           "No suitable condition/treatment column detected in this CSV.")
    )

    tagList(
      selectInput(
        "cond_col",
        label   = "Condition column:",
        choices = found,
        selected = found[1]
      )
    )
  })

  # ── 4. Group checkboxes (dynamic, from selected column) ──
  output$condition_ui <- renderUI({
    req(raw_data(), input$cond_col)
    df     <- raw_data()
    groups <- sort(unique(as.character(df[[input$cond_col]])))

    checkboxGroupInput(
      "conditions",
      label    = "Include groups:",
      choices  = groups,
      selected = groups          # all groups selected by default
    )
  })

  # ── 5. Custom time-point checkboxes ─────────────────────
  output$custom_tp_ui <- renderUI({
    req(raw_data())
    tps <- sort(unique(raw_data()$Time.Pt))
    checkboxGroupInput("custom_tp", "Select time points:",
                       choices = tps, selected = tps)
  })

  # ── 6. Per-group color pickers (no hardcoded names) ──────
  output$color_ui <- renderUI({
    req(input$conditions)
    tagList(
      lapply(seq_along(input$conditions), function(i) {
        grp <- input$conditions[i]
        # Cycle through fallback palette; no hardcoded group names
        def <- FALLBACK_PALETTE[((i - 1) %% length(FALLBACK_PALETTE)) + 1]
        colourInput(
          inputId    = paste0("col_", make.names(grp)),
          label      = grp,
          value      = def,
          showColour = "background"
        )
      })
    )
  })

  # ── 7. Process data (fires on Run) ──────────────────────
  proc_data <- eventReactive(input$run, {
    req(raw_data(), input$conditions, input$outcome,
        input$cond_col)

    df       <- raw_data()
    cond_col <- input$cond_col   # e.g. "Treatment", "Condition", etc.

    # Validate required columns
    req_cols <- c("mouse_id", "Time.Pt", cond_col, input$outcome)
    miss     <- setdiff(req_cols, names(df))
    validate(need(
      length(miss) == 0,
      paste("Missing columns:", paste(miss, collapse = ", "))
    ))

    # ── Filter time points by phase ──────────────────────
    if (input$phase == "after") {
      df        <- df %>% filter(Time.Pt %in% c("120/CL-0", "CL-10", "CL-20"))
      df$Time.Pt[df$Time.Pt == "120/CL-0"] <- "CL-0"
      tp_levels <- c("CL-0", "CL-10", "CL-20")

    } else if (input$phase == "before") {
      df        <- df %>% filter(Time.Pt %in% c("90", "120/CL-0"))
      df$Time.Pt[df$Time.Pt == "120/CL-0"] <- "120"
      tp_levels <- c("90", "120")

    } else {
      req(input$custom_tp)
      df        <- df %>% filter(Time.Pt %in% input$custom_tp)
      tp_levels <- input$custom_tp
    }

    # ── Filter selected groups from the chosen column ────
    df <- df %>%
      filter(as.character(.data[[cond_col]]) %in% input$conditions)

    validate(need(
      nrow(df) > 0,
      "No rows remain after filtering. Adjust time points or group selection."
    ))

    # ── Rename condition column → internal "condition" ───
    # This keeps all downstream model/plot code unchanged
    # regardless of the original column name.
    if (cond_col != "condition") {
      df <- df %>% rename(condition = all_of(cond_col))
    }

    # ── Center body-mass covariates (if present) ─────────
    for (cov_name in c("Lean", "Fat", "Total")) {
      if (cov_name %in% names(df)) {
        df[[paste0(cov_name, "_c")]] <-
          as.numeric(scale(df[[cov_name]], center = TRUE, scale = FALSE))
      }
    }

    # ── Set factor levels ─────────────────────────────────
    df$mouse_id  <- factor(df$mouse_id)
    df$condition <- factor(df$condition,
                           levels = intersect(input$conditions,
                                              unique(df$condition)))
    df$Time.Pt   <- factor(df$Time.Pt, levels = tp_levels)
    df
  })

  # ── 8. Color map ─────────────────────────────────────────
  color_map <- eventReactive(input$run, {
    req(input$conditions)
    cm <- vapply(input$conditions, function(grp) {
      val <- input[[paste0("col_", make.names(grp))]]
      if (is.null(val) || !nzchar(val)) "#888888" else val
    }, character(1))
    setNames(cm, input$conditions)
  })

  # ── 9. Y-axis label ──────────────────────────────────────
  y_label <- reactive({
    req(input$outcome)
    lbl <- KNOWN_LABELS[input$outcome]
    if (is.na(lbl)) input$outcome else as.character(lbl)
  })

  # ── 10. Summary statistics ───────────────────────────────
  sum_stats <- reactive({
    req(proc_data(), input$outcome)
    out <- input$outcome
    proc_data() %>%
      group_by(Condition = condition, `Time Point` = Time.Pt) %>%
      summarise(
        N    = n(),
        Mean = round(mean(.data[[out]], na.rm = TRUE), 4),
        SD   = round(sd(.data[[out]],   na.rm = TRUE), 4),
        SE   = round(SD / sqrt(N),                     4),
        .groups = "drop"
      )
  })

  # ── 11. Model formula ────────────────────────────────────
  model_formula <- reactive({
    req(input$outcome, input$covariate)
    out <- input$outcome
    cov <- input$covariate
    if (cov == "none") {
      as.formula(
        sprintf("`%s` ~ condition * Time.Pt + (1 | mouse_id)", out)
      )
    } else {
      as.formula(
        sprintf("`%s` ~ condition * Time.Pt + %s_c + (1 | mouse_id)",
                out, cov)
      )
    }
  })

  # ── 12. Fit LME model ────────────────────────────────────
  fit_model <- reactive({
    req(proc_data(), model_formula())
    df  <- proc_data()
    cov <- input$covariate

    if (cov != "none") {
      cov_c <- paste0(cov, "_c")
      validate(need(
        cov_c %in% names(df),
        sprintf("Column '%s' not found. Cannot adjust for %s.", cov_c, cov)
      ))
    }

    withCallingHandlers(
      tryCatch(
        lmer(model_formula(), data = df, REML = TRUE),
        error = function(e) {
          validate(need(FALSE, paste("Model failed:", conditionMessage(e))))
        }
      ),
      warning = function(w) {
        showNotification(paste("Model warning:", conditionMessage(w)),
                         type = "warning", duration = 8)
        invokeRestart("muffleWarning")
      }
    )
  })

  # ── 13. Status banner ────────────────────────────────────
  output$status_banner <- renderUI({
    if (input$run == 0) {
      div(class = "alert alert-info", style = "margin-bottom:15px;",
          icon("info-circle"),
          " Upload a CSV file and click ",
          strong("Run Analysis"), " to begin.")
    }
  })

  # ══ Rendered outputs ══════════════════════════════════════

  # Summary table
  output$sum_table <- renderDT({
    req(sum_stats())
    datatable(sum_stats(),
              options  = list(pageLength = 20, dom = "tip"),
              rownames = FALSE,
              class    = "stripe hover compact")
  })

  output$dl_sum <- downloadHandler(
    filename = function() paste0("summary_stats_", Sys.Date(), ".csv"),
    content  = function(f) write.csv(sum_stats(), f, row.names = FALSE)
  )

  # ── Bar plot ─────────────────────────────────────────────
  bar_plot_obj <- reactive({
    req(sum_stats(), proc_data(), color_map(), input$outcome,
        input$dodge_w, input$bar_ptsize, input$bar_alpha, input$bar_jitter)

    ss  <- sum_stats() %>%
             rename(condition = Condition, Time.Pt = `Time Point`)
    df  <- proc_data()
    out <- input$outcome
    cm  <- color_map()
    dw  <- input$dodge_w

    ggplot(ss, aes(x = Time.Pt, y = Mean, fill = condition)) +
      geom_col(position  = position_dodge(dw),
               width     = dw * 0.87,
               color     = "black",
               linewidth = 0.4) +
      geom_errorbar(
        aes(ymin = Mean - SE, ymax = Mean + SE),
        position  = position_dodge(dw),
        width     = 0.15,
        linewidth = 0.5
      ) +
      geom_point(
        data        = df,
        aes(x = Time.Pt, y = .data[[out]], fill = condition),
        position    = position_jitterdodge(jitter.width = input$bar_jitter,
                                            dodge.width  = dw),
        shape       = 21,
        size        = input$bar_ptsize,
        stroke      = 0.3,
        color       = "black",
        alpha       = input$bar_alpha,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = cm) +
      labs(x = "Time Point", y = y_label(), fill = "Condition") +
      theme_classic(base_size = 13) +
      theme(
        axis.line       = element_line(color = "black", linewidth = 0.5),
        axis.ticks      = element_line(color = "black", linewidth = 0.5),
        axis.text       = element_text(color = "black"),
        axis.title      = element_text(color = "black"),
        legend.position = "top",
        legend.title    = element_text(face = "bold")
      )
  })

  output$bar_plot <- renderPlot({ req(bar_plot_obj()); bar_plot_obj() })

  output$dl_bar <- downloadHandler(
    filename = function() paste0("bar_plot_", Sys.Date(), ".png"),
    content  = function(f) ggsave(f, bar_plot_obj(),
                                  width = 8, height = 5.5,
                                  dpi = 300, bg = "white")
  )

  # ── Line plot ────────────────────────────────────────────
  line_plot_obj <- reactive({
    req(proc_data(), color_map(), input$outcome)
    df  <- proc_data()
    out <- input$outcome
    cm  <- color_map()

    ggplot(df, aes(x = Time.Pt, y = .data[[out]],
                   color = condition, group = condition)) +
      stat_summary(fun      = mean, geom = "line",
                   aes(group = condition),
                   linewidth = 0.9,
                   position  = position_dodge(0.15)) +
      stat_summary(fun = mean, geom = "point",
                   aes(fill = condition),
                   shape    = 21, size = 3.2,
                   stroke   = 0.4, color = "black",
                   position = position_dodge(0.15)) +
      stat_summary(fun.data  = mean_se,
                   geom      = "errorbar",
                   width     = 0.08, linewidth = 0.5,
                   position  = position_dodge(0.15)) +
      geom_point(aes(fill = condition),
                 shape    = 21, size  = 2.2,
                 stroke   = 0.25, color = "black",
                 alpha    = 0.6,
                 position = position_jitterdodge(jitter.width = 0.06,
                                                  dodge.width  = 0.15)) +
      scale_color_manual(values = cm) +
      scale_fill_manual(values  = cm) +
      labs(x = "Time Point", y = y_label(),
           color = "Condition", fill = "Condition") +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "top",
        legend.title    = element_text(face = "bold"),
        axis.text       = element_text(color = "black"),
        axis.title      = element_text(color = "black")
      )
  })

  output$line_plot <- renderPlot({ req(line_plot_obj()); line_plot_obj() })

  output$dl_line <- downloadHandler(
    filename = function() paste0("line_plot_", Sys.Date(), ".png"),
    content  = function(f) ggsave(f, line_plot_obj(),
                                  width = 8, height = 5.5,
                                  dpi = 300, bg = "white")
  )

  # ── EMM interaction plot ─────────────────────────────────
  output$emm_plot <- renderPlot({
    req(fit_model(), color_map())
    emmip(fit_model(), condition ~ Time.Pt, CIs = TRUE) +
      scale_color_manual(values = color_map()) +
      labs(
        x     = "Time Point",
        y     = paste0("EMM: ", y_label()),
        title = "Estimated Marginal Means \u00b1 95% CI",
        color = "Condition"
      ) +
      theme_classic(base_size = 13) +
      theme(
        legend.position = "top",
        legend.title    = element_text(face = "bold"),
        axis.text       = element_text(color = "black"),
        axis.title      = element_text(color = "black")
      )
  })

  # ── Model summary ────────────────────────────────────────
  output$model_summary <- renderPrint({
    req(fit_model())
    cat("Model formula:\n  ", deparse(formula(fit_model())), "\n")
    cat("Condition column used: '", input$cond_col, "'\n\n", sep = "")
    print(summary(fit_model()))
  })

  # ── ANOVA ─────────────────────────────────────────────────
  output$anova_table <- renderPrint({
    req(fit_model())
    cat("Type III ANOVA (Satterthwaite df):\n\n")
    print(anova(fit_model()))
  })

  # ── EMM: condition | Time.Pt ─────────────────────────────
  output$emm_cond_time <- renderPrint({
    req(fit_model())
    emm <- emmeans(fit_model(), ~ condition | Time.Pt)
    cat("=== EMM: Condition within each Time Point ===\n\n")
    print(summary(emm))
    cat("\n=== Pairwise comparisons (Tukey) ===\n\n")
    print(pairs(emm, adjust = "tukey"))
  })

  # ── EMM: Time.Pt marginal ────────────────────────────────
  output$emm_time_overall <- renderPrint({
    req(fit_model())
    emm_t <- emmeans(fit_model(), ~ Time.Pt)
    cat("=== EMM: Time Points (marginal over conditions) ===\n\n")
    print(summary(emm_t))
    cat("\n=== Pairwise comparisons (Tukey) ===\n\n")
    print(pairs(emm_t, adjust = "tukey"))
  })

  # ── EMM: Time.Pt | condition ─────────────────────────────
  output$emm_time_cond <- renderPrint({
    req(fit_model())
    emm_tc <- emmeans(fit_model(), ~ Time.Pt | condition)
    cat("=== EMM: Time Points within each Condition ===\n\n")
    print(summary(emm_tc))
    cat("\n=== Pairwise comparisons (Tukey) ===\n\n")
    print(pairs(emm_tc, adjust = "tukey"))
  })
}

# ============================================================
shinyApp(ui=ui_app4, server=server_app4)
