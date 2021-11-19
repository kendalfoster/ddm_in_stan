Function to load the chains from file (argument `nums` is a vector of the numbers/indices of the chains)
```
read_stan_chain_data_csv <- function(filename, nums) {
  if (tolower(substr(filename, nchar(filename)-3, nchar(filename))) == ".csv") {
    filename <- substr(filename, 1, nchar(filename)-4)
  }
  if (substr(filename, nchar(filename), nchar(filename)) == .Platform$file.sep) {
    filename <- paste0(filename, "chain_data")
  }
  filename <- paste0(filename, "_", nums, ".csv")

  stanfit <- rstan::read_stan_csv(filename)
  return(stanfit)
}
```

Example extraction from file
```
chains <- read_stan_chain_data_csv(file.path(getwd(), "fits", "completed_non-hierarchical_chains", "chain_data"), c(1, 2, 4))
```
