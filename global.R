
cat("🔍 Checking for required packages...\n")

required_packages <- c(
  "shiny",
  "tidyverse",
  "ggpubr",
  "RColorBrewer",
  "shinyjs",
  "readxl",
  "DT"
)

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

for (pkg in required_packages) {
  install_if_missing(pkg)
}

cat("\n🎉 All required packages are installed and loaded!\n")
cat("========================================\n\n")

library(shiny)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(shinyjs)
library(readxl)
library(DT)

# for app2
# library(shiny)
library(bslib)
# library(tidyverse)
library(emmeans)
library(broom)
# library(DT)



# library(pracma)
# library(mgcv)
# library(gratia)
# library(emmeans)
# library(lubridate)

ensure_metric_cage_cols <- function(df, prefixes, cage_numbers, fill = -1) {
  prefixes <- unique(trimws(as.character(prefixes)))
  prefixes <- prefixes[nzchar(prefixes)]
  cage_numbers <- sort(unique(suppressWarnings(as.integer(cage_numbers))))
  cage_numbers <- cage_numbers[!is.na(cage_numbers)]
  if (length(prefixes) == 0 || length(cage_numbers) == 0) return(df)
  
  expected <- as.vector(outer(prefixes, cage_numbers, function(p, k) paste0(p, "_", k)))
  missing  <- setdiff(expected, names(df))
  
  if (length(missing) > 0) {
    for (col in missing) df[[col]] <- fill
  }
  df
}

options(shiny.maxRequestSize = 50 * 1024^3)