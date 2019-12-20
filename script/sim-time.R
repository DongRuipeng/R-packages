require(peer)
require(foreach)
require(doParallel)

dir.create("~/time-result", showWarnings = FALSE, recursive = TRUE)
setwd("~/time-result")

fix.n <- fix.p <- fix.q <- 100
nrank <- 3
s <- 4
snr <- 0.5
erank <- nrank + 1

start <- 100
end <- 100
step <- 100

times <- 2
# begin parallel computing
cl.num <- 2
cl <- makeCluster(cl.num)
registerDoParallel(cl)

# increase p
vec.p <- seq(from = start, to = end, by = step)
cat("begin to increase p...\n")
block.time <- Sys.time()
T.p <- data.frame()
for (p in vec.p) {
  set.seed(10)
  sim.data <- sim.setup(n = fix.n, p = p, q = fix.q, nrank = nrank, s = s, snr = snr)
  X <- sim.data$X
  Y <- sim.data$Y
  ave.time <- 0
  ave.time <- foreach(k = 1:times, .combine = '+', .packages = c("peer")) %dopar% {
    sim.time(X, Y, erank = erank)
  }
  ave.time <- ave.time / times
  method <- row.names(ave.time)
  temp <- data.frame(method = method, time = ave.time, p = rep(p, length(method)))
  rownames(temp) <- NULL
  T.p <- rbind(T.p, temp)
  cat("p =", p, "is over \n", sep = " ")
}
cat("finish increasing p \n")
block.time <- Sys.time() - block.time
print(block.time)
cat("\n")


# increase q
vec.q <- seq(from = start, to = end, by = step)
cat("begin to increase q...\n")
block.time <- Sys.time()
T.q <- data.frame()
for (q in vec.q) {
  set.seed(10)
  sim.data <- sim.setup(n = fix.n, p = fix.p, q = q, nrank = nrank, s = s, snr = snr)
  X <- sim.data$X
  Y <- sim.data$Y
  ave.time <- 0
  ave.time <- foreach(k = 1:times, .combine = '+', .packages = c("peer")) %dopar% {
    sim.time(X, Y, erank = erank)
  }
  ave.time <- ave.time / times
  method <- row.names(ave.time)
  temp <- data.frame(method = method, time = ave.time, q = rep(q, length(method)))
  rownames(temp) <- NULL
  T.q <- rbind(T.q, temp)
  cat("q =", q, "is over \n", sep = " ")
}
cat("finish increasing q \n")
block.time <- Sys.time() - block.time
print(block.time)
cat("\n")


# increase n
vec.n <- seq(from = start, to = end, by = step)
cat("begin to increase n...\n")
block.time <- Sys.time()
T.n <- data.frame()
for (n in vec.n) {
  set.seed(10)
  sim.data <- sim.setup(n = n, p = fix.p, q = fix.q, nrank = nrank, s = s, snr = snr)
  X <- sim.data$X
  Y <- sim.data$Y
  ave.time <- 0
  ave.time <- foreach(k = 1:times, .combine = '+', .packages = c("peer")) %dopar% {
    sim.time(X, Y, erank = erank)
  }
  ave.time <- ave.time / times
  method <- row.names(ave.time)
  temp <- data.frame(method = method, time = ave.time, n = rep(n, length(method)))
  rownames(temp) <- NULL
  T.n <- rbind(T.n, temp)
  cat("n =", n, "is over \n", sep = " ")
}
cat("finish increasing n \n")
block.time <- Sys.time() - block.time
print(block.time)
cat("\n")


# end parallel computing
stopImplicitCluster()
stopCluster(cl)

# save T.p, T.q and T.n 
save(file = "./time.RData", list = c("T.p", "T.q", "T.n"))