####################
##### Shiny UI #####
####################
shinyUI(pageWithSidebar(
  
  # Application title.
  headerPanel(""),
  
  sidebarPanel(
    uiOutput("choose_patient"),

    div(style="display:inline-block", actionButton("renew_list", "Renew List")),    
    div(style="display:inline-block", actionButton("plotSurv", "Plot Survivl Curves")),

    br(),

    uiOutput("time_slider"),

    div(style="display:inline-block", actionButton("varSelect", "Find Best Moedel")),    
    div(style="display:inline-block", actionButton("plotHazard", "Plot Fitted Cumulative Hazards"))
  ), 
  
  mainPanel(
    tabsetPanel(
      tabPanel("Hazard Predicted", plotOutput("haz_predicted")),
      tabPanel("Survival", plotOutput("plot_predicted")),
      tabPanel("Model Summary", verbatimTextOutput("summary")),
      tabPanel("Data", tableOutput("raw_data")),
      tabPanel("Best Model", tableOutput("best_model")),
      tabPanel("Hazard Fitted",htmlOutput("haz_fitted")),
      tabPanel("Test",textOutput("test")),
      id = "tabs"
    ) 
  )
))
