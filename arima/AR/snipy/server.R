library(shiny)
library(shinydashboard)
library(ggplot2)
source("./utilities.R")

shinyServer( function(input, output, session) {
    
    output$sliders <- renderUI({
        numIndividuals <- as.integer(input$M)
        lapply(1:numIndividuals, function(i) {
            sliderInput(paste0("p", i), paste0("Lag Term ", i), -2, 2, 0, .05)
        })
    })
    
    data <- reactive({sim_ar_data(input$N, 
                                  sapply(1:input$M, function(x) input[[paste0("p",x)]]), 
                                  input$mu)})
    
    output$dataplot <- renderPlot({
        ggplot(data(), aes(x=time, y=obs)) + geom_line()
    })
    
    output$estplot <- renderPlot({
        print(head(data()$obs))
        print(as.integer(input$M))
        plot_models(data()$obs, as.integer(input$M))
    })
})