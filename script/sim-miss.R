require(foreach)
require(doParallel)
require(peer)

dir.create("~/miss-result", showWarnings = FALSE)
setwd("~/miss-result")

set.seed(10)

p <- 200
q <- n <- 100
s <- 4
nrank <- 3
erank <- nrank + 1

repetition <- 2
cl.num <- 2
cl <- makeCluster(cl.num)
registerDoParallel(cl)

for (snr in c(0.25)) {
  filepath <- paste("./snr-", snr, sep = "")
  dir.create(path = filepath, showWarnings = FALSE)
  table <- data.frame()
  for (miss in c(0.1)) {
    # define simulation setup
    setting <- sim.control(n = n, p = p, q = q, nrank = nrank, s = s, snr = snr, erank = erank)

    temp1 <- foreach(i = 1:repetition, .combine = 'rbind', .packages = c("peer")) %dopar% {
      sim.miss(control = setting)
    }

    miss_col <- data.frame(miss = rep(miss, nrow(temp1)))
    temp2 <- cbind(temp1, miss_col)
    table <- rbind(table, temp2)
  }
  save(file = paste(filepath, "/miss.RData", sep = ""), list = "table")
}
stopImplicitCluster()
stopCluster(cl)
