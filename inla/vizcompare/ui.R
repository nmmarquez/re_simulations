rm(list=ls())

library(shiny)
library(shinydashboard)
library(data.table)


header <- dashboardHeader(
    title = 'Spatiotemporal Field Estimation'
)

body <- dashboardBody(
    fluidRow(
        column(width=12,
               tabBox(id='tabvals', width=NULL,
                      tabPanel('Latent', plotOutput('latent'), value=1),
                      tabPanel('Parameters', plotOutput('param'), value=2)#,
                      #tabPanel('Metadata', plotOutput('meta'), value=3)
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
    selectInput('sigma', 'Sigma', seq(.01, 2, .25)),
    selectInput('range', 'Range', seq(.01, 1, .15)),
    selectInput('rho', 'Rho', seq(.01, 1, .15)),
    selectInput('size', 'Size', c(50, 100, 500, 1000, 10000)),
    selectInput('type', 'Type', c('Mean', 'Variance'))
)

dashboardPage(
    header,
    sidebar,
    body
)
