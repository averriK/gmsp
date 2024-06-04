message("Reading AT2 data")
DATASET <- readRDS(file.path("data/AT2.Rds"))

UNITS <- "mm"
RSN <- 4000 
RAW <- DATASET[[RSN]]
FMAX <- 25
RESAMPLE <- FALSE
LOWPASS <- TRUE
IMF1 <- 0
IMFN <-0
message("Building TimeSeries")

AUX <- gmsp::buildTS(
  x=RAW$AT,
  dt=RAW$dt,
  Units=RAW$SourceUnits,
  OrderTS =2,
  OrderEMD = NULL,
  Fmax=FMAX,
  Resample = RESAMPLE,
  LowPass = LOWPASS,
  TargetUnits=UNITS,
  Astop.AT=1e-4,
  Apass.AT=5e-4,
  removeIMF1 = IMF1,
  removeIMFn = IMFN)

TSL <- AUX$TSL[ OCID  %in%  c("H1","H2","UP")]
rm(AUX)
