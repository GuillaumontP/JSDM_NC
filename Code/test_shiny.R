
server <- function(input, output, session) {    
  data <- reactive({
    subset(
      iris,
      Species == input$species
    )
  })    
  output$plot <- renderPlot({
    ggplot(data = data()) + 
      geom_point(aes(x = Sepal.Length, y = Sepal.Width))
  })    
}
ui <- fluidPage(
  titlePanel("iris"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "species",
        label = "Species",
        choices = unique(iris$Species)
      )
    ),
    mainPanel(  
      plotOutput("plot")
    )
  )
)
shinyApp(ui, server)


