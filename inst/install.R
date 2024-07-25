file.remove("NAMESPACE")
# usethis::use_data_raw()
usethis::use_proprietary_license(copyright_holder = "Alejandro Verri Kozlowski")
devtools::document()
devtools::check()
remove.packages("gmsp")
# devtools::build()
devtools::install()
## Commit Push
# remotes::install_github("averrik/gmsp",auth_token = Sys.getenv("PAT"))
# remotes::install_github("averriK/gmsp")


