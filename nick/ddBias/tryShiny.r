rm(list=ls())

#
library(shiny)
library(shinydashboard)
library(ggplot2)

#
data <-  data.frame(x=c(1,2,3,4),y=c(10,11,12,13))

#
ui <- dashboardPage(
  dashboardHeader(),
  dashboardSidebar(sliderInput("sliderA","A", min=1, max=3, step=0.5, value=1),
               sliderInput("sliderK","K", min=1, max=10, step=1, value=1)),
  dashboardBody(
    fluidRow(column(6,plotOutput('waveplot')))
  ))

#
server <- function(input, output, session) { 
  output$waveplot <- renderPlot({
    x <- seq(0,10,0.1)
    yfxn <- function(x) { input$sliderA*sin(input$sliderK*x) }
    y <- yfxn(x)
    df <- data.frame(x,y)
    ggplot(df,aes_string(x=x,y=y))+geom_point(size=2)+geom_line()+ 
         scale_x_continuous()
  })
}

#
shinyApp(ui, server)
