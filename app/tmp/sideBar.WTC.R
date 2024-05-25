sideBar.WTC <- function(ID="WTC"){
  
  
  fluidRow(
    # Noise ----
    fluidRow(
      helpText("Wavelet:"),
      column(width = 10,
             prettyRadioButtons(
               inputId = NS(id=ID,"mother"),
               label = NULL,
               choices=c("morlet","dog","paul"),
               selected="morlet",
               inline = TRUE, 
               status = "danger",
               fill = TRUE
             )
      )
    ),
    
    
    
    
    # Scale ----
    fluidRow(
      helpText("Separación entre escalas j=1/dj:"),
      column(width = 12,
             sliderTextInput(
               inputId = NS(id=ID,"j"),
               label = NULL, 
               choices = c(2,4,8,12,16,24,32,48,64,72,128,256),
               selected = 32,
               grid = TRUE
             )
      )
    ),
    # Phase ----
    fluidRow(
      helpText("Fases:"),
      column(width = 10,
             # prettySwitch(
             #   inputId = NS(id=ID,"phase"),
             #   label = "Incluir fase", 
             #   value = default,
             #   status = "success",
             #   fill = TRUE)
             prettyRadioButtons(
               inputId = NS(id=ID,"phase"),
               label = NULL,
               choiceNames = c("Incluir","Ocultar"),
               choiceValues = c(TRUE,FALSE),
               # selected=TRUE,
               inline = TRUE, 
               status = "danger",
               fill = TRUE
             )
             
      )
    ),
  )
}
