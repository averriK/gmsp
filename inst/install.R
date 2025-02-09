file.remove("NAMESPACE") |> suppressWarnings()
# usethis::use_data_raw()
usethis::use_proprietary_license(copyright_holder = "Alejandro Verri Kozlowski")
devtools::check(document = TRUE)
remove.packages("gmsp") |> suppressWarnings()
# devtools::build()
# devtools::install()
## Commit Push
# remotes::install_github("averrik/gmsp",auth_token = Sys.getenv("PAT"))
devtools::install()
remotes::install_github("averriK/gmsp")


