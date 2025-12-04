#' @title Simulation Study for LDART (Scenario A)
#'
#' @description
#' Runs a simulation study for longitudinal dynamic treatment regimes using the GcompBART package.
#' This script simulates data under Scenario A and evaluates bias, RMSE, and coverage.
#'
#' @details
#' - Scenario A uses \eqn{\lambda = (0.25, 0.5, 0.75, 1)} and true mean \eqn{\text{mtm} = 36.42767}.
#' - Predictors are generated using \code{sim_lfried()} with 4 time points and 5 variables per time point.
#' - Correlation structure is defined by \code{lags = c(0.4, 0.2, 0.1)}.
#' - Outcome is generated with \eqn{\sigma = 10}.
#'
#' The script fits:
#' 1. Block models using \code{BMfits()}.
#' 2. G-computation using \code{gcompbart()}.
#'
#' Metrics computed:
#' - Absolute relative bias
#' - RMSE
#' - 95% coverage
#'
#' @param n_rep Integer. Number of simulation repetitions (default: 500).
#' @param sample_sizes Integer vector. Sample sizes to evaluate (default: c(250, 1000, 2000)).
#' @param n_cores Integer. Number of cores for parallel processing (default: 16).
#'
#' @return Writes results to a text file with columns:
#' \code{sample_size | abs_bias | RMSE | coverage}.
#'
#' @examples
#' # Run Scenario A simulation
#' # Source this script and execute:
#' # Rscript sim_ldart_A.R
#'
#' @seealso \code{\link{sim_lfried}}, \code{\link{BMfits}}, \code{\link{gcompbart}}

############################################################
# Simulation study: Longitudinal DART using GcompBART
# Scenario A: lambda = c(0.25, 0.5, 0.75, 1), mtm = 36.42767
############################################################

library(GcompBART)
library(foreach)
library(doParallel)

source("sim_data.R") # contains helper functions

########################################################
# --- Simulation settings ---
n_tree   <- 50
n_burn   <- 300
n_thin   <- 4
n_save   <- 500
n_rep    <- 500          # number of repetitions
n_cores  <- 16           # parallel cores

# Variable types and time grouping
var_type <- c(rep("X0", 5), rep("X", 15), "Y")
t_group  <- c(rep(1:4, each = 5), 4)

# Hyperparameters for BART
hypers_sim <- BaseHypers(phi = rep(1, max(t_group) - 1),
                         alpha_vec = rep(1, max(t_group) - 1),
                         num_tree = n_tree)


opts_sim <- Opts(num_burn = n_burn,
                 num_thin = n_thin,
                 num_save = n_save,
                 update_s = FALSE,
                 update_alpha = FALSE,
                 update_tvp = TRUE,
                 update_alpha_vec = TRUE,
                 update_eta = TRUE,
                 update_phi = TRUE)

# --- Define scenarios ---
scenarios <- list(
  A = list(lambda = c(0.25, 0.5, 0.75, 1), mtm = 36.42767),
  B = list(lambda = c(0, 0, 0, 1), mtm = 14.57107),
  C = list(lambda = c(0.25, 0.25, 0.25, 0.25), mtm = 14.57107)
)

# --- Scenario A ---
scenario <- scenarios$A
lambda   <- scenario$lambda
mtm      <- scenario$mtm

# --- sample_sizes ---
sample_sizes <- c(250, 1000, 2000)
n <- sample_sizes[1]

# Output file
fna <- "ldart_results_A.txt"

########################################################
# Parallel setup
myCluster <- makeCluster(n_cores)
registerDoParallel(myCluster)

foreach(i = 1:n_rep, .packages = c("GcompBART")) %dopar% {

    # Simulate longitudinal data
    training_data <- sim_lfried(N = n, n_time = 4, n_var = 5,
                                lags = c(0.4, 0.2, 0.1),
                                sigma = 10,
                                lambda = lambda)

    # Fit BM model
    BM <- BMfits(training_data[, 1:21],
                 var.type = var_type,
                 tgroup = t_group,
                 opts = opts_sim,
                 base_hypers = hypers_sim)

    # Fit gcompbart
    out_bart <- gcompbart(training_data[, 1:21],
                          var.type = var_type,
                          J = 10000,
                          tgroup = t_group,
                          BModels = BM)

    # Compute metrics
    out <- c(abs(out_bart$summary_out[1] - mtm) / mtm, # abs relative bias
             rmse(out_bart$y_hat, mtm),                # RMSE
             coverage(mtm, out_bart$y_hat))            # 95% coverage

    # Write results
    write(c(n, out), fna, ncolumns = 4, append = TRUE)
  }

stopCluster(myCluster)

