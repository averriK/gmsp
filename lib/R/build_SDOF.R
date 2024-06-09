build_SDOF<- function(TSL,xi=0.05,Units){
  
  PSW <- TSL[ID=="AT",gmsp::build_Spectra(.x=.SD,Units = Units,xi=xi),by=.(RecordSN,OCID,DIR)]
  
}
