source("setup.R",local = TRUE)
source("global.R",local = TRUE)
source("R/helpers.R",local = TRUE)

source("R/GMDB.R",local = TRUE)
source("R/TS.R",local = TRUE)
source("R/SDOF.R",local = TRUE)
source("R/FFT.R",local = TRUE)
source("R/STFT.R",local = TRUE)
source("R/WTC.R",local = TRUE)
source("R/EMD.R",local = TRUE)

server <- function(input, output, session) {

  TS.server(id="TS",.data=TSL,series="AT",yAxis.legend="A(t)")
  TS.server(id="TS",.data=TSL,series="VT",yAxis.legend="V(t)")
  TS.server(id="TS",.data=TSL,series="DT",yAxis.legend="D(t)")
  
  FFT.server(id="FFT",.data=TSL,series="AT")
  FFT.server(id="FFT",.data=TSL,series="VT")
  FFT.server(id="FFT",.data=TSL,series="DT")
  # 
  
  SDOF.server(id="SDOF",.data=TSL)
  
  STFT.server(id="STFT",.data=TSL,series="AT")
  STFT.server(id="STFT",.data=TSL,series="VT")
  STFT.server(id="STFT",.data=TSL,series="DT")
  # 
  
  WTC.server(id="WTC",.data=TSL,series="AT")
  WTC.server(id="WTC",.data=TSL,series="VT")
  WTC.server(id="WTC",.data=TSL,series="DT")
  # 
  
  
  EMD.server(id="EMD",.data=TSL,series="AT",yAxis.legend="A(t)")
  EMD.server(id="EMD",.data=TSL,series="VT",yAxis.legend="V(t)")
  EMD.server(id="EMD",.data=TSL,series="DT",yAxis.legend="D(t)")
  # 
}
