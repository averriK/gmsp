file.remove("NAMESPACE")
usethis::use_proprietary_license(copyright_holder = "Alejandro Verri Kozlowski")
devtools::document()
devtools::check()
remove.packages("gmsp")
devtools::install()
## Commit Push
## remotes::install_github("averrik/gmsp",auth_token = Sys.getenv("PAT"))
