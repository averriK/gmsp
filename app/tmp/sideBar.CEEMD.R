
sideBar.CEEMD <- function(ID="CEEMD"){
  
  # boundary <- req(input$boundary)
  # max.imf <- req(input$max.imf) |> as.integer()
  # stop.rule <- req(input$stop.rule)
  # noise.type <- req(input$noise.type)
  # noise.amp <- req(input$noise.amp) |> as.double()
  # trials <- req(input$trials) |> as.integer()
  
  fluidRow(
    # Number of IMFs ----
    fluidRow(
      helpText("Max IMFs:"),
      column(width = 10,
             sliderTextInput(
               inputId = NS(id=ID,"max.imf"),
               label = NULL, 
               choices = seq(1,15),
               selected = 5,
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
               inline = TRUE, 
               status = "danger",
               fill = TRUE
             )
      )
    ),
    
    # Number of Iterations ----
    fluidRow(
      helpText("Trials"),
      column(width = 10,
             sliderTextInput(
               inputId = NS(id=ID,"trials"),
               label = NULL, 
               choices = c(1,3,5,10,25,50,100,200),
               selected = 5,
               grid = TRUE
             )
      )
    ),
    
    # Noise ----
    fluidRow(
      helpText("Noise type:"),
      column(width = 10,
             prettyRadioButtons(
               inputId = NS(id=ID,"noise.type"),
               label = NULL,
               choices=c("uniform","gaussian"),
               selected="gaussian",
               inline = TRUE, 
               status = "danger",
               fill = TRUE
             )
      ),
      helpText("Noise amplitude:"),
      column(width = 10,
             numericInput(
               inputId = NS(id=ID,"noise.amp"),
               label = NULL, 
               value = 0.5e-7,
               min =0.5e-12,               
               max =1
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