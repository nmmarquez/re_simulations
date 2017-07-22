library(shiny)
source("./utilities.R")

shinyServer(function(input,output){
    output$latent <- renderPlot({
        if(input$type == "mean"){
            plot_latent(input$sigma, input$range, input$rho, input$N)
        }
        else{
            plot_var(input$sigma, input$range, input$rho, input$N)
        }
    })
    output$param <- renderPlot({
        plot_params(input$sigma, input$range, input$rho, input$N)
    })
})