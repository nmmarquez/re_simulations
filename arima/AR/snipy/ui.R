rm(list=ls())

library(shiny)
library(shinydashboard)
library(leaflet)
library(data.table)


header <- dashboardHeader(
    title = 'AR Estimation Comparsion'
)

body <- dashboardBody(
    fluidRow(
        column(width=12,
               tabBox(id='tabvals', width=NULL,
                      tabPanel('Estimation', plotOutput('estplot'), value=1),
                      tabPanel('Data', plotOutput('dataplot'), value=2)
               )
        ) 
    ),
    tags$head(tags$style(HTML('
                              section.content {
                              height: 2500px;
                              }
                              ')))
)



sidebar <- dashboardSidebar(
    sliderInput('N', 'Observations', 10, 10000, 100),
    sliderInput('mu', 'Constant', 0, 5, 0, .2),
    selectInput("d", "Differences", 0:4, 0),
    selectInput("M", "Lag Terms", 1:5),
    uiOutput("sliders")
)

dashboardPage(
    header,
    sidebar,
    body
)
