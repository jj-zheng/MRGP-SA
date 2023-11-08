source("para_mttdl.R")
source("commonMatrix.R")

mttdf.list <- numeric(0)
mttdl.list <- numeric(0)
ava.list <- numeric(0)

mttdf.array <- c(1:100)*10000
for (mttdf in mttdf.array) {
  
  lambda <- 1/mttdf
  source("raid6_mat.R")
  source("raid6_vec.R")

  # Q_{ij}
  Q00 <- -rowSums(G0G1E)
  Q01 <- G0G1E
  Q11 <- G1G1E + diag(-rowSums(G1G1E)) + diag(-rowSums(G1G2E))
  Q12 <- G1G2E
  Q22 <- 0

  tQ00 <- solve(-Q00)

  # P_{i,j}
  P10 <- G1G0P0
  P11 <- G1G1P0
  P20 <- G2G0P1

  # E_i
  E1 <- expm(Q11*t1)
  E2 <- exp(Q22*t2)
  #bE11 <- mexpAm.unif(A=Q11, t=t1, m=Q11, cumulative=TRUE)$cx
  #bE12 <- mexpAm.unif(A=Q11, t=t1, m=Q12, cumulative=TRUE)$cx
  bE1 <- mexpAm.unif(A=Q11, t=t1, cumulative=TRUE)$cx
  bE2 <- t2

  # eP_{i,j}
  eP11 <- E1 %*% P11 + E1 %*% P10 %*% tQ00 %*% Q01 + bE1 %*% Q11
  eP12 <- bE1 %*% Q12
  eP21 <- E2 %*% P20 * tQ00 %*% Q01

  Pemc <- rBind(
    cBind(eP11,eP12),
    cBind(eP21,0))

  # Steady-satte probability
  stres <- ctmc.st(Pemc - eyeM(dim(Pemc)[1]))
  stopifnot(stres$convergence)
  piemc <- stres$x

  ind1 <- 1:dim(eP11)[1]
  ind2 <- max(ind1) + 1:dim(eP12)[2]

  piemc1 <- piemc[ind1]
  piemc2 <- piemc[ind2]

  # sojourn time
  st0 <- as.vector(piemc1 %*% E1 %*% P10 %*% tQ00 + piemc2 %*% E2 %*% P20 %*% tQ00)
  st1 <- as.vector(piemc1 %*% bE1)
  st2 <- as.vector(piemc2 %*% bE2)

  # steady-state probabilities
  pi0 <- as.vector(st0/(sum(st0) + sum(st1) + sum(st2)))
  pi1 <- as.vector(st1/(sum(st0) + sum(st1) + sum(st2)))
  pi2 <- as.vector(st2/(sum(st0) + sum(st1) + sum(st2)))

  # MTTDL
  ava <- pi0+sum(pi1)
  mttdl <- (ava * t2/(1-ava)) / (365 * 24) 
  cat("# ", mttdf, ava, mttdl)
  
  mttdf.list <- c(mttdf.list, mttdf)
  mttdl.list <- c(mttdl.list, mttdl)
  ava.list <- c(ava.list, ava)
}

res <- cbind(mttdf.list, ava.list, mttdl.list)
write.table(res, "~/n10/txt/10res24t1.txt", row.names = F, col.names = F, append = F)

