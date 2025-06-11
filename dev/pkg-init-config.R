####################################################################################################
# pkg-init-config.R
# Initial R package setup
#
####################################################################################################
library(devtools)
library(usethis)
library(testthat)
# library(desc)

# #--------------------------------------------------
# #----- Initial config, file initialization
# #--------------------------------------------------
# #--- Create package project skeleton
# usethis::create_package("../criteval")

# # --- Add entries to .Rbuildignore
# # Add files to .Rbuildignore
# usethis::use_build_ignore(c("LICENSE.md", ".Rhistory", ".Rproj.user"))
# # To add an entire directory and all its contents (recursively) to .Rbuildignore
# # include the line ^path/to/dir/ (or possibly ^path/to/dir/.*$ might be safer?).
# usethis::edit_r_buildignore() # add lines ^.*\.Rproj$ ^ignore/ and ^dev/
# usethis::use_build_ignore("man-roxygen")

# #----- Create skeleton {packageName}-package.R file
# usethis::use_package_doc()

#----- Create DESCRIPTION file
per1 <- utils::person(
  given = "David",
  family = "Fleischer",
  email = "david.p.fleischer@gmail.com",
  role = c("cre", "aut")
)

usethis::use_description(fields = list(
  Title = "Simulations and Comparisons of Splitting Criteria for Generalized Random Forests",
  Description = paste(
    "An evaluation framework for the simulation and comparison of splitting criteria used",
    "by GRF in its recursive partitioning mechanism, particularly for causal forests and",
    "multi-dimensional varying-coefficient models."
  ),
  `Authors@R` = c(per1)
))

usethis::use_mit_license()

#----- Setup testthat file structure
# usethis::use_testthat()

#----- Set some package dependencies
usethis::use_package("dplyr", type = "Imports")
usethis::use_package("tidyr", type = "Imports")
usethis::use_package("bench", type = "Imports")
usethis::use_package("ggplot2", type = "Imports")
usethis::use_import_from("mvtnorm", "rmvnorm")
usethis::use_import_from("ggbeeswarm", c("geom_beeswarm", "geom_quasirandom"))
