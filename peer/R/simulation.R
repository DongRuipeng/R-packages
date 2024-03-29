#' estimation error
#' @param C true coefficient 
#' @param Chat estimated coefficient 
#' @return estimation error
#' @export
est.err <- function(C, Chat) {
  p <- nrow(C)
  q <- ncol(C)
  err <- (norm(C - Chat, "F") ^ 2 / (p * q)) * 1000
  return(err)
}

#' prediction error 
#' @param X predictors matrix 
#' @param C true coeffiecient 
#' @param Chat estimated coefficient 
#' @return prediction error
#' @export
pred.err <- function(X, C, Chat) {
  p <- nrow(C)
  n <- nrow(X)
  q <- ncol(C)
  err <- (norm(X %*% C - X %*% Chat, "F") ^ 2 / (n * q)) * 1000
  return(err)
}

#' FDR and FNR 
#' @param U left singular vectors 
#' @param Uhat estimated left singular vectors
#' @param dhat estimated singualr values 
#' @return 
#'  \item{rate}{consisting FDR and FNR}
#' @export
frate <- function(U, Uhat, dhat) {
  p <- nrow(U)
  r <- ncol(U)
  #if did't find out all eigenvalues, use 0 to append
  rhat <- length(dhat)
  if (rhat < r) {
    append <- r - rhat
    for (k in 1:append) {
      dhat[rhat + k] <- 0
      Uhat <- cbind(Uhat, rep(0, p))
    }
  }
  #sort d with decreasing order to calculate numerator
  ind <- sort(dhat, decreasing = TRUE, index.return = TRUE)$ix
  fpv <- fnv <- rep(0, r)
  for (i in 1:r) {
    fpv[i] <- sum((U[, i] == 0) & (Uhat[, ind[i]] != 0))
    fnv[i] <- sum((U[, i] != 0) & (Uhat[, ind[i]] == 0))
  }
  fp <- sum(fpv)
  fn <- sum(fnv)
  neg <- sum(U == 0)
  pos <- sum(U != 0)
  fprate <- fp / neg * 100
  fnrate <- fn / pos * 100
  return(list(fprate = fprate, fnrate = fnrate))
}


#' generate a list consisting simulation setup 
#' @param n number of observations
#' @param p dimension of predictor vector
#' @param q dimension of reponse vector 
#' @param nrank rank of true coefficient 
#' @param s sparsity of left singular vectors
#' @param snr singnal to noise rate 
#' @param erank rank for estimator 
#' @return 
#'  \item{control}{a list consisting simulation setup}
#' @export
sim.control <- function(n, p, q, nrank, s, snr, erank, miss = 0) {
  return(list(n = n, p = p, q = q, nrank = nrank, s = s,
  snr = snr, erank = erank, miss = miss))
}

#' generate result table with RRR, SRRR, SEED, PEER(l1) and PEER(l0), which
#' includes ErC, ErXC, FPR, FNR and time. 
#' @param control a list consists simlation setting, which can be generated by sim.control() function
#' @return \item{table}{result table with RRR, SRRR, SEED, PEER(l1) and PEER(l0), which includes ErC, ErXC, FPR, FNR and time.}
#' @export
#' @import rrpack
sim.table <- function(control = control) {
  n <- control$n
  p <- control$p
  q <- control$q
  nrank <- control$nrank
  erank <- control$erank
  s <- control$s
  snr <- control$snr

  sim.data <- sim.setup(n, p, q, nrank, s, snr)
  X <- sim.data$X
  Y <- sim.data$Y
  U <- sim.data$U
  V <- sim.data$V
  d <- sim.data$d
  C <- sim.data$C

  rrr.fit <- srrr.fit <- seed.fit <- peer.fit.l1 <-
  peer.fit.l0 <- data.frame(ErC = 0, ErXC = 0, FPR = 0, FNR = 0, time = 0)

  # rrr (D is diagonal matrix decrease)
  start_time <- Sys.time()
  fit <- rrr(Y, X, penaltySVD = c("rank", "ann")[1], ic.type = "GIC")
  end_time <- Sys.time()
  Chat <- fit$coef
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- diag(fit$D)
  rrr.fit$ErC <- est.err(C, Chat)
  rrr.fit$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  rrr.fit$FPR <- rate$fprate #why not 100 in rrr
  rrr.fit$FNR <- rate$fnrate
  rrr.fit$time <- end_time - start_time

  # srrr (D diagonal matrix, may not increase)
  start_time <- Sys.time()
  fit <- srrr(Y, X, nrank = erank, ic.type = "GIC")
  end_time <- Sys.time()
  Chat <- fit$U %*% fit$D %*% t(fit$V)
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- diag(fit$D)
  srrr.fit$ErC <- est.err(C, Chat)
  srrr.fit$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  srrr.fit$FPR <- rate$fprate
  srrr.fit$FNR <- rate$fnrate
  srrr.fit$time <- end_time - start_time

  # seed (d is a vector)
  start_time <- Sys.time()
  fit <- seed(X, Y, nrank = erank)
  end_time <- Sys.time()
  Chat <- fit$C
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- fit$d
  seed.fit$ErC <- est.err(C, Chat)
  seed.fit$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  seed.fit$FPR <- rate$fprate
  seed.fit$FNR <- rate$fnrate
  seed.fit$time <- end_time - start_time

  # peer(l1) (d is a vector)
  fit <- peer(X, Y, nrank = erank, penalty = "l1")
  Chat <- fit$C
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- fit$d
  peer.fit.l1$ErC <- est.err(C, Chat)
  peer.fit.l1$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  peer.fit.l1$FPR <- rate$fprate
  peer.fit.l1$FNR <- rate$fnrate
  peer.fit.l1$time <- fit$time

  # peer(l0) (d is a vector)
  fit <- peer(X, Y, nrank = erank, penalty = "l0")
  Chat <- fit$C
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- fit$d
  peer.fit.l0$ErC <- est.err(C, Chat)
  peer.fit.l0$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  peer.fit.l0$FPR <- rate$fprate
  peer.fit.l0$FNR <- rate$fnrate
  peer.fit.l0$time <- fit$time

  table <- rbind(rrr.fit, srrr.fit, seed.fit, peer.fit.l1, peer.fit.l0)
  method <- c("RRR", "SRRR", "SEED", "PEER(l1)", "PEER(l0)") 
  table <- cbind(table, method = method)
  rownames(table) <- NULL 
  # print(table)
  return(table)
}

#' generate table for missing data simulation. Then tbale includes ErC, ErXC, FPR, FNR and time with 
#' l1 and l0 penalty. 
#' @param control a list consists simlation setting, which can be generated by sim.control() function 
#' @return \item{table}{including ErC, ErXC, FPR, FNR and time with l1 and l0 penalty.}
#' @export 
sim.miss <- function(control = sim.control()) {
  n <- control$n
  p <- control$p
  q <- control$q
  nrank <- control$nrank
  erank <- control$erank
  s <- control$s
  snr <- control$snr
  miss <- control$miss

  sim.data <- sim.setup(n, p, q, nrank, s, snr, miss)
  X <- sim.data$X
  Y <- sim.data$Y
  U <- sim.data$U
  V <- sim.data$V
  d <- sim.data$d
  C <- sim.data$C

  peer.l1 <- peer.l0 <- data.frame(ErC = 0, ErXC = 0, FPR = 0, FNR = 0, time = 0)

  # peer with l1 penalty 
  fit <- peer(X, Y, nrank = erank, penalty = "l1")
  Chat <- fit$C
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- fit$d
  peer.l1$ErC <- est.err(C, Chat)
  peer.l1$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  peer.l1$FPR <- rate$fprate
  peer.l1$FNR <- rate$fnrate
  peer.l1$time <- fit$time

  # peer with l0 penalty 
  fit <- peer(X, Y, nrank = erank, penalty = "l0")
  Chat <- fit$C
  Uhat <- fit$U
  Vhat <- fit$V
  dhat <- fit$d
  peer.l0$ErC <- est.err(C, Chat)
  peer.l0$ErXC <- pred.err(X, C, Chat)
  rate <- frate(U, Uhat, dhat)
  peer.l0$FPR <- rate$fprate
  peer.l0$FNR <- rate$fnrate
  peer.l0$time <- fit$time

  table <- rbind(peer.l1, peer.l0)
  method <- c("peer(l1)", "peer(l0)")
  table <- cbind(table, data.frame(method = method))
  rownames(table) <- NULL

  return(table)
}

#' conuting time for each methos (SRRR, SEED, PEER(l1), PEER(l0)) 
#' @param X predictors matrix 
#' @param Y responses matrix 
#' @param erank rank of estimator 
#' @return {Time}{times of SRRR, SEED, PEER(l1), PEER(l0), respectively}
#' @import rrpack
#' @export 
sim.time <- function(X, Y, erank) {
  Time <- matrix(0, 2, 1)
  rownames(Time) <- c("PEER(l1)", "PEER(l0)")

  # # SRRR
  # start_time <- proc.time()
  # fit <- srrr(Y, X, nrank = erank, ic.type = "GIC")
  # end_time <- proc.time()
  # Time[1,] <- as.double(end_time[1] - start_time[1])
  # # SEED
  # start_time <- proc.time()
  # fit <- seed(X, Y, nrank = erank)
  # end_time <- proc.time()
  # Time[2,] <- as.double(end_time[1] - start_time[1])
  # PEER(l1)
  fit <- peer(X, Y, nrank = erank, penalty = "l1")
  Time[1,] <- fit$time
  # PEER(l0)
  fit <- peer(X, Y, nrank = erank, penalty = "l0")
  Time[2,] <- fit$time

  return(Time)
}