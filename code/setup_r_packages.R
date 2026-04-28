options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c("dplyr", "magrittr", "data.table"))
install.packages(
  "cmdstanr",
  repos = c("https://stan-dev.r-universe.dev", getOption("repos"))
)

library(cmdstanr)
install_cmdstan(cores = 20)

install.packages("brms")

library(brms)
options(brms.backend = "cmdstanr", brms.debug = TRUE)

cat("R version:", R.version.string, "\n")
cat("cmdstan path:", cmdstan_path(), "\n")
cat("cmdstan version:", as.character(cmdstan_version()), "\n")
cat("brms version:", as.character(packageVersion("brms")), "\n")