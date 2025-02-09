file.remove("NAMESPACE") |> suppressWarnings()
# usethis::use_data_raw()
usethis::use_proprietary_license(copyright_holder = "Alejandro Verri Kozlowski")
devtools::check(document = TRUE)
remove.packages("gmsp") |> suppressWarnings()
# devtools::build()
## Commit Push
devtools::install()
# remotes::install_github("averriK/gmsp")


