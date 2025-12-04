library(GcompBART)
library(foreach)
library(doParallel)

betula_data <- read.table("betula_data")

vartype_nc <- c(rep("X0",6),rep("X",3),"Y","S",
                rep("X",5),"X","Y","S",
                rep("X",4),"X","Y","S",
                rep("X",4),"X","Y")

vartype_inc <- c(rep("X0",6),rep("X",3),"Y","S",
                 rep("X",5),"Ra","Y","S",
                 rep("X",4),"Ra","Y","S",
                 rep("X",4),"Ra","Y")

tgroup <- c(rep(0,3),rep(1,8),rep(2,8),rep(3,7),rep(4,6))

# --- Threshold ---
threshold_value <- 140  # change to 130, 140, or 150
sb_threshold <- rep(threshold_value, 3)

# --- Shift ---
delta <- cbind(-2*sd(betula_data$BP_systolic.1),
               0,
               -2*sd(betula_data$BP_systolic.1) + 0.1)

n_burn <- 300
n_thin <- 100
n_save <- 150
n_tree <- 50

opts <- Opts(num_burn = n_burn, num_thin = n_thin, num_save = n_save,
             update_s = FALSE, update_alpha = FALSE,
             update_tvp = TRUE, update_alpha_vec = TRUE,
             update_eta = TRUE, update_phi = TRUE)


myCluster <- makeCluster(16)
registerDoParallel(myCluster)

foreach(i=1:16, .packages = c("GcompBART") %dopar% {
  
  BM <- BMfits(betula_data,
               var.type = vartype_nc,
               drop_param = drop_par,
               opts = opts,
               tgroup = tgroup,
               base_hypers = BaseHypers(num_tree = n_tree))
  
  ##################################################
  
  ages <- seq(35, 70, by = 5)
  
  for (age in ages) {
    sub_dat <- subset(betula_data, Age_cohort_T1 == age)
    
    frs <- gcompbart(sub_dat,
                     var.type = vartype_nc,
                     J = 10000,
                     opts = opts,
                     weighted = FALSE,
                     tgroup = tgroup,
                     BModels = BM)
    
    frs_inc1 <- gcompbart(sub_dat,
                          var.type = vartype_inc,
                          random.regime = rep("triangular", 3),
                          param = list(delta, delta, delta),
                          nat_value = TRUE,
                          above = TRUE,
                          threshold = sb_threshold,
                          incremental = TRUE,
                          J = 10000,
                          opts = opts,
                          tgroup = tgroup,
                          BModels = BM)
    
    
    out <- t(rbind(frs$y_hat,
                   frs_inc1$y_hat,
                   frs$s_hat,
                   frs_inc1$s_hat))
    
    fna <- paste("MAR_lbart_age", age, sep = "_")
    
    write.table(out, fna, append=TRUE, col.names = FALSE, row.names = FALSE )
    
  }
}
stopCluster(myCluster)

