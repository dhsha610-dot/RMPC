source("global.R")

# Source all UIs/Servers
# source("app1v2.R")
source("app2v2.R")
source("app3v5.R")
source("app1e_Rodent_Metabolic_Phenotyping_App_v2.3.13.R")
source("app4v2.R")

ui <- navbarPage("RMPC Data Anlysis Platform",
                 # tabPanel("🔵 Cage data alignment", ui_app1), 
                 tabPanel("🟠 Cage data alignment & summary", ui_app1e),
                 tabPanel("🔵 nmolar Calculator", ui_app3), 
                 tabPanel("🟢 Statistics analysis", ui_app2),
                 tabPanel("🟢 Linear Mixed-Effects Model", ui_app4)
)

server <- function(input, output, session) {
  # server_app1(input, output, session)
  server_app1e(input, output, session)
  server_app3(input, output, session)
  server_app2(input, output, session)
  server_app4(input, output, session)
}

shinyApp(ui, server)