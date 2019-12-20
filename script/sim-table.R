# load package
require(foreach)
require(doParallel)
require(peer)

# set rereplication number
repetition <- 2
# set cluster number
cl.num <- 2

# set work dir
dir.create("~/table-result/", showWarnings = FALSE, recursive = TRUE)
setwd("~/table-result/")

set.seed(10)

for (snr in c(0.25)) {
  cat("starting with snr: ", snr, "\n", sep = "")
  time <- Sys.time()
  filepath <- paste("./snr=", snr, sep = "")
  dir.create(filepath, showWarnings = FALSE, recursive = TRUE)
  for (p in c(100)) {
    for (r in c(3)) {
      cat("p = ", p, ", r = ", r, '\n', sep = "")
      # define simulation setup
      setting <- sim.control(n = 100, p = p, q = 100, nrank = r, s = 4, snr = snr, erank = r + 1)
      # begin simulation

      # parallel computing
      cl <- makeCluster(cl.num)
      registerDoParallel(cl)
      table <- foreach(i = 1:repetition, .combine = 'rbind', .packages = 'peer') %dopar% {
        sim.table(control = setting)
      }
      # end parallel computing
      stopImplicitCluster()
      stopCluster(cl)

      filename <- paste("p=", p, "--r=", r, ".Rdata", sep = "")
      path <- paste(filepath, "/", filename, sep = "")
      save(table, file = path)
    }
  }
  time <- Sys.time() - time
  print(time)
  cat("\n")
}
