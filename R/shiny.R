shiny <- function(jSDM_bin_pro){
  #' Create shiny output
  #'
  #' @description
  #' Create a shiny output where you can see variability of each variable for each species
  #'
  #' @import shiny
  #' @import stars
  #' @import here
  #' @export

  # load(jSDM_bin_pro_path)

  ui <- fluidPage(
    titlePanel(title = h1("Visualization effect of explanatories variables on each species", align = "center")),
    fluidRow(
      column(3,
             selectInput(inputId = "spec",
                         label = "Please select a species : ",
                         choices = sort(colnames(read.csv2(here("data_raw", "NCpippn", "Presence_Absence.csv"), sep = ","))[2:878]))
      ),
      column(3,
             selectInput("variable",
                         "Please select an explanatory variable : ",
                         choices = colnames(jSDM_bin_pro$mcmc.sp[[1]])[3:7])
      )
    ),
    mainPanel(
      plotOutput("distPlot")
    )
  )

  server <- shinyServer(function(input, output) {
    data <- reactive({

      library(here)
      library(stars)
      library(matrixStats)
      library(stringr)
      load(here("jSDM_bin_pro.RData"))
      beta <- matrix(0, nrow = 12, ncol = dim(jSDM_bin_pro$model_spec$presence_data)[2])
      colnames(beta) <- colnames(jSDM_bin_pro$model_spec$presence_data)
      for (i in 1:dim(jSDM_bin_pro$model_spec$presence_data)[2]) {
        beta[, i] <- colMeans(jSDM_bin_pro$mcmc.sp[[i]])[1:12]
      }
      mins <- rowMins(beta) # min of each var
      maxs <- rowMaxs(beta) # max of each var
      names(mins) <- names(maxs) <- colnames(jSDM_bin_pro$mcmc.sp[[i]])[1:12]
      rownames(beta) <- names(maxs)
      input$variable
      input$spec
      values <- seq(mins[input$variable], maxs[input$variable], length = 100)
      x <- values
      y <- values * beta[input$variable, input$spec] + values^2 * beta[paste0("beta_`", str_split(str_extract(input$variable, "_.*"), "_")[[1]][2], "^2`"), input$spec]
      df <- data.frame(x = x, y = y)
      return(df)
    })

    output$distPlot <- renderPlot({
      library(ggplot2)
      ggplot(data = data()) +
        geom_point(aes(x = x, y = y), color = "salmon") +
        theme_bw() +
        labs(x = paste(input$variable, " values"),
             y = "",
             title = paste0(input$spec, " with only \n", input$variable, " moving")) +
        theme(plot.title = element_text(hjust = 0.5))
    })
  })

  shinyApp(ui = ui, server = server)
}
