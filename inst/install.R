file.remove("NAMESPACE")
# usethis::use_data_raw()
usethis::use_proprietary_license(copyright_holder = "Alejandro Verri Kozlowski")
devtools::document()
devtools::check()
remove.packages("gmsp")
devtools::install()
## Commit Push
## remotes::install_github("averrik/gmsp",auth_token = Sys.getenv("PAT"))


# generar un m*.R por cada serie patologica
#
# generar plots IMF
# generar plots FFT.IMF
#
# generar plot de las tres series para una direccion y registro dados
#
# identificar frecuencias en desplazamientos y eliminarlas
# identificar frecuencias en aceleraciones y eliminarlas
# definir si se hace un paper corto o un dashboard

