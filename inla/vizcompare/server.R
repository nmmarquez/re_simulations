library(shiny)
source("./utilities.R")

shinyServer(function(input,output){
    output$latent <- renderPlot({
        if(input$type == "mean"){
            res <-plot_latent(input$sigma, input$range, input$rho, input$size)
        }
        else{
            res <- plot_var(input$sigma, input$range, input$rho, input$size)
        }
        res
    })
    output$param <- renderPlot({
        plot_params(input$sigma, input$range, input$rho, input$size)
    })
})