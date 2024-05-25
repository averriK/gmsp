sideBar.STFT <- function(ID="STFT"){
  fluidRow(
    
    # Windows length ----
    fluidRow(
      helpText("Windows Length k=N/NW:"),
      column(width = 10,
             sliderTextInput(
               inputId = NS(id=ID,"k"),
               label = NULL, 
               choices = c(4,8,16,32,64,128),
               selected = 32,
               grid = TRUE)
      )
    ),
    
    # Windows Overlap ----
    fluidRow(
      helpText("Windows Overlap [%]:"),
      column(width = 10,
             sliderTextInput(
               inputId = NS(id=ID,"ovlp"),
               label = NULL, 
               choices = c(10,25,50,75,90),
               selected = 75,
               grid = TRUE
             )
      )
    ),
    
    # Nyqist ----
    fluidRow(
      helpText("Frecuencia máxima fmax/fNyquist:"),
      column(width = 10,
             sliderTextInput(
               inputId = NS(id=ID,"kf"),
               label = NULL, 
               choices = c(seq(0.01,0.14,by=0.01),seq(0.15,0.95,by=0.05)),
               selected = 0.25,
               grid = FALSE)
      ),
    )
  )
}
