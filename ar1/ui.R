rm(list=ls())
library(shiny)
library(shinydashboard)

source("./sim_ar1.R")

header <- dashboardHeader(
    title = 'AR1 simulation'
)

body <- dashboardBody(
    fluidRow(
        column(width=12,
               tabBox(
                   id='tabvals',
                   width=NULL,
                   tabPanel('plots',
                            plotOutput('d1'), value=1)
               )
        ) 
    ),
    tags$head(tags$style(HTML('
                              section.content {
                              height: 2500px;
                              }
                              ')))
)

start <- "20160822_all-cause-tests_LDI-only_Global-secular"

sidebar <- dashboardSidebar(
    sliderInput('N_obs', 'In Sample Observations', min=0, max=100, step=1, value=20),
    sliderInput('N_forecast', 'Forecasted Years', min=0, max=100, step=1, value=20),
    sliderInput('rho', 'Rho', min=0, max=1, step=.01, value=.99),
    textInput('sigma', 'Sigma', value=1),
    textInput('trend', 'Trend', value=0)
)

dashboardPage(
    header,
    sidebar,
    body
)