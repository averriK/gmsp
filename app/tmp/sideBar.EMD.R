sideBar.EMD <- function(ID="EMD"){
  
  # boundary <- req(input$boundary)
  # max.imf <- req(input$max.imf) |> as.integer()
  # stop.rule <- req(input$stop.rule)
  
  
  fluidRow(
    # Number of IMFs ----
    fluidRow(
      helpText("Max IMFs:"),
      column(width = 10,
             sliderTextInput(
               inputId = NS(id=ID,"max.imf"),
               label = NULL, 
               choices = seq(1,15),
               selected = 12,
               grid = TRUE
             )
      )
    ),
    
    # Stop rule ----
    fluidRow(
      helpText("Stop rule:"),
      column(width = 10,
             prettyRadioButtons(
               inputId = NS(id=ID,"stop.rule"),
               label = NULL,
               choices=c("type1","type2","type3","type4","type5"),
               selected="type5",
               inline = FALSE, 
               status = "danger",
               fill = TRUE
             )
      )
    ),
    
    
    
    # Boundary ----
    fluidRow(
      helpText("Boundary type:"),
      column(width = 10,
             
             prettyRadioButtons(
               inputId = NS(id=ID,"boundary"),
               label = NULL,
               choices=c("none","wave","symmetric","periodic","evenodd"),
               selected="wave",
               inline = FALSE, 
               status = "danger",
               fill = TRUE
             )
      )
    )
  )
}