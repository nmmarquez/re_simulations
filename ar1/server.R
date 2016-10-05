rm(list=ls())
library(shiny)
library(ggplot2)
source("./sim_ar1.R")


shinyServer(function(input,output){
    
    output$d1 <- renderPlot({
        df <- ar1_forecast(input$N_obs, input$N_forecast, input$rho, 
                           as.numeric(input$sigma), as.numeric(input$trend))
        ggplot(df, aes(x=t, y=y)) + geom_line() + 
            geom_ribbon(aes(ymin=lower_bound, ymax=upper_bound), alpha=.4)
        
    })
})