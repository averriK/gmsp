devtools::load_all()
root_folder <- file.path("inst/test")
DIR <- "HHE" # HHE/HHN
ID <- "167_SRK03" # 167_SRK03, 235_SRK03
record_folder <- file.path(root_folder, ID, DIR)
record_file <- "R_167_SRK03_HHE_D0-99.txt"
record_path <- file.path(record_folder, record_file)
# Read header, dt, NPTS
header_line <- readLines(record_path, n = 2)[2]
# read file, ignore header
VT <- fread(record_path, skip = 2)
# extract "dt" from "dt=0.005, npts=25418" in header_line
dt <- sub(".*dt=([0-9.]+).*", "\\1", header_line) |> as.numeric()
# build time vector
npts <- sub(".*npts=([0-9]+).*", "\\1", header_line) |> as.numeric()
ts <- seq(0, (npts - 1) * dt, by = dt)
# build TimeHistory with build_TS()


build_TS_VT(
    x = VT,
    dt = dt,
    ts = ts,
    Fmax = 25,
    Resample = TRUE,
    Units = "m",
    TargetUnits = "mm",
)

build_TS(
    x = VT,
    dt = dt,
    ts = ts,
    Fmax = 25,
    Resample = TRUE,
    Units = "m",
    TargetUnits = "mm",
)
