file.remove("NAMESPACE")
usethis::use_proprietary_license(copyright_holder = "Alejandro Verri Kozlowski")
devtools::document()
devtools::check()
# remove.packages("deeplr")
# devtools::install()
## Commit Push
## remotes::install_github("averrik/deeplr",auth_token = Sys.getenv("PAT"))
