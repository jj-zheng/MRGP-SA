source("para_sens.R")
source("commonMatrix.R")
source("mexp.sens.unif.R")

comp.mttdl.raid6 <- function(t1, t2, param) {

    tQ00 <- solve(-param$Q00)

    # E_i
    E1 <- expm(param$Q11*t1)
    E2 <- exp(param$Q22*t2)
    #bE11 <- mexpAm.unif(A=Q11, t=t1, m=Q11, cumulative=TRUE)$cx
    #bE12 <- mexpAm.unif(A=Q11, t=t1, m=Q12, cumulative=TRUE)$cx
    bE1 <- mexpAm.unif(A=param$Q11, t=t1, cumulative=TRUE)$cx
    bE2 <- t2

    # eP_{i,j}
    eP11 <- E1 %*% param$P11 + E1 %*% param$P10 %*% tQ00 %*% param$Q01 + bE1 %*% param$Q11
    eP12 <- bE1 %*% param$Q12
    eP21 <- E2 %*% param$P20 * tQ00 %*% param$Q01

    Pemc <- rbind(
      cbind(eP11,eP12),
      cbind(eP21,0))

    # Steady-satte probability
    stres <- ctmc.st(Pemc - eyeM(dim(Pemc)[1]))
    stopifnot(stres$convergence)
    piemc <- stres$x

    ind1 <- 1:dim(eP11)[1]
    ind2 <- max(ind1) + 1:dim(eP12)[2]

    piemc1 <- piemc[ind1]
    piemc2 <- piemc[ind2]

    # sojourn time
    st0 <- as.vector(piemc1 %*% E1 %*% param$P10 %*% tQ00 + piemc2 %*% E2 %*% param$P20 %*% tQ00)
    st1 <- as.vector(piemc1 %*% bE1)
    st2 <- as.vector(piemc2 %*% bE2)

    # steady-state probabilities
    pi0 <- as.vector(st0/(sum(st0) + sum(st1) + sum(st2)))
    pi1 <- as.vector(st1/(sum(st0) + sum(st1) + sum(st2)))
    pi2 <- as.vector(st2/(sum(st0) + sum(st1) + sum(st2)))

    # MTTDL
    ava <- pi0+sum(pi1)
    mttdl <- (ava * t2/(1-ava)) / (365 * 24) 
    #cat("# ", ava, mttdl)

    list(E1=E1,E2=E2,bE1=bE1,bE2=bE2,piemc1=piemc1,piemc2=piemc2,Pemc=Pemc,st0=st0,st1=st1,st2=st2,pi0=pi0,pi1=pi1,pi2=pi2,ava=ava,mttdl=mttdl)

}

#- 20230412 revised on Apr 25th, 2023
comp.mttf.raid6 <- function(t1, param) {
  
  tQ00 <- solve(-param$Q00)
  
  # E_1 and bE1
  E1 <- expm(param$Q11*t1)
  bE1 <- mexpAm.unif(A=param$Q11, t=t1, cumulative=TRUE)$cx
  
  # m0 and m1
  X <- eyeM(dim(E1 %*% param$P11)[1]) - E1 %*% param$P11 - E1 %*% param$P10 %*% tQ00 %*% param$Q01
  Y <- bE1 %*% vone(dim(bE1)[1]) + E1 %*% param$P10 %*% tQ00
  m1 <- solve(X) %*% Y
  m0 <- tQ00 + tQ00 %*% param$Q01 %*% m1
  
  # mttf
  #t2 <- 24
  #res <- comp.mttdl.raid6(t1, t2, param.q)
  #pi0 <- res$pi0
  #mttf <- pi0 %*% m0 / (365 * 24) 
  mttf <- m0 / (365 * 24)
  mttf
  #cat("# ", mttf)
  
}
#-

comp.mttff.raid6 <- function(t1, t2, param) {
  
  tQ00 <- solve(-param$Q00)
  
  # E_i
  E1 <- expm(param$Q11*t1)
  E2 <- exp(param$Q22*t2)
  bE1 <- mexpAm.unif(A=param$Q11, t=t1, cumulative=TRUE)$cx
  bE2 <- t2
  
  # eP_{i,j}
  eP11 <- E1 %*% param$P11 + E1 %*% param$P10 %*% tQ00 %*% param$Q01 + bE1 %*% param$Q11
  eP12 <- bE1 %*% param$Q12
  eP21 <- E2 %*% param$P20 * tQ00 %*% param$Q01
  
  Pemc <- rbind(
    cbind(eP11,eP12),
    cbind(eP21,0))
  
  # Steady-satte probability
  stres <- ctmc.st(Pemc - eyeM(dim(Pemc)[1]))
  stopifnot(stres$convergence)
  piemc <- stres$x
  
  ind1 <- 1:dim(eP11)[1]
  piemc1 <- piemc[ind1]
  
  # sojourn time
  st1 <- as.vector(piemc1 %*% bE1)
  
  # initial state probability distribution
  alp <- c(1,0)
  
  # mttff
  tQ <- solve(eyeM(dim(eP11)[1]) - eP11) 
  #temp <- as.vector(alp %*% tQ %*% st1)
  temp <- as.vector(sum(alp %*% tQ %*% bE1))
  mttff <- (temp / (365 * 24)) 
  cat("# ", mttff)
  
#  list(E1=E1,E2=E2,bE1=bE1,bE2=bE2,piemc1=piemc1,piemc2=piemc2,Pemc=Pemc,st0=st0,st1=st1,st2=st2,ava=ava,mttdl=mttdl)
  
}

comp.sens.raid6 <- function(t1, t2, param.q, param.dq) {

    res <- comp.mttdl.raid6(t1, t2, param.q)

    #1 n
    #dE1n <- mexp.sens.unif(param.q$Q11, t1, param.dq$dq11n)$S
    #dbE1n <- cmexp.sens.unif.t(param.q$Q11, t1, param.dq$dq11n)$X
    #dE2n <- 0
    #dbE2n <- 0

    #deP11n <- dE1n %*% (param.q$P11 + 1/(n*lambda) * param.q$P10 %*% param.q$Q01) + res$E1 %*% (1/(n*lambda) * param.q$P10 %*% param.dq$dq01n - (1/{(n*lambda)}^2) * lambda * param.q$P10 %*% param.q$Q01) + dbE1n %*% param.q$Q11 + res$bE1 %*% param.dq$dq11n
    #deP12n <- dbE1n %*% param.q$Q12 + res$bE1 %*% param.dq$dq12n
    #deP21n <- -(1/{(n*lambda)}^2) * lambda * param.q$P20 %*% param.q$Q01 + (1/(n*lambda)) * param.q$P20 %*% param.dq$dq01n
    #deP22n <- 0

    #dPemcn <- rbind(
    #      cbind(deP11n,deP12n),
    #      cbind(deP21n,0))

    # Steady-satte probability
    #stres <- ctmc.st(eyeM(dim(res$Pemc)[1]) - res$Pemc)
    #stopifnot(stres$convergence)
    #dpiemcn <- stres$x

    #ind1 <- 1:dim(deP11n)[1]
    #ind2 <- max(ind1) + 1:dim(deP12n)[2]

    #dpiemc1n <- dpiemcn[ind1]
    #dpiemc2n <- dpiemcn[ind2]

    # sojourn time
    #dst0n <- as.vector(-(1/{(n*lambda)}^2) * lambda * (res$piemc1 %*% res$E1 %*% param.q$P10 + res$piemc2) + (1/(n*lambda)) * (dpiemc1n %*% res$E1 %*% param.q$P10 + res$piemc1 %*% dE1n %*% param.q$P10 + dpiemc2n))
    #dst1n <- as.vector(dpiemc1n %*% res$bE1 + res$piemc1 %*% dbE1n)
    #dst2n <- as.vector(dpiemc2n * t2)

    # steady-state probabilities
    #dpi0n <- as.vector((dst0n * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st0 * (sum(dst0n) + sum(dst1n) + sum(dst2n)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)
    #dpi1n <- as.vector((dst1n * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st1 * (sum(dst0n) + sum(dst1n) + sum(dst2n)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)
    #dpi2n <- as.vector((dst2n * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st2 * (sum(dst0n) + sum(dst1n) + sum(dst2n)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)

    # MTTDL
    #davan <- dpi0n+sum(dpi1n)
    #dmttdln <- (davan * t2)/(1-res$ava)^2
    #cat("# ", davan, dmttdln)

    #2 lambda
    dE1r <- mexp.sens.unif(param.q$Q11, t1, param.dq$dq11r)$S
    dbE1r <- cmexp.sens.unif.t(param.q$Q11, t1, param.dq$dq11r)$X
    dE2r <- 0
    dbE2r <- 0

    deP11r <- dE1r %*% (param.q$P11 + 1/(n*lambda) * param.q$P10 %*% param.q$Q01) + res$E1 %*% (1/(n*lambda) * param.q$P10 %*% param.dq$dq01r - (1/{(n*lambda)}^2) * n * param.q$P10 %*% param.q$Q01) + dbE1r %*% param.q$Q11 + res$bE1 %*% param.dq$dq11r
    deP12r <- dbE1r %*% param.q$Q12 + res$bE1 %*% param.dq$dq12r
    deP21r <- -(1/{(n*lambda)}^2) * n * param.q$P20 %*% param.q$Q01 + (1/(n*lambda)) * param.q$P20 %*% param.dq$dq01r
    deP22r <- 0

    dPemcr <- rbind(
          cbind(deP11r,deP12r),
          cbind(deP21r,0))

    # Steady-satte probability
    stres <- ctmc.st(eyeM(dim(res$Pemc)[1]) - res$Pemc)
    stopifnot(stres$convergence)
    dpiemcr <- stres$x

    ind1 <- 1:dim(deP11r)[1]
    ind2 <- max(ind1) + 1:dim(deP12r)[2]

    dpiemc1r <- dpiemcr[ind1]
    dpiemc2r <- dpiemcr[ind2]

    # sojourn time
    dst0r <- as.vector(-(1/{(n*lambda)}^2) * n * (res$piemc1 %*% res$E1 %*% param.q$P10 + res$piemc2) + (1/(n*lambda)) * (dpiemc1r %*% res$E1 %*% param.q$P10 + res$piemc1 %*% dE1r %*% param.q$P10 + dpiemc2r))
    dst1r <- as.vector(dpiemc1r %*% res$bE1 + res$piemc1 %*% dbE1r)
    dst2r <- as.vector(dpiemc2r * t2)

    # steady-state probabilities
    dpi0r <- as.vector((dst0r * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st0 * (sum(dst0r) + sum(dst1r) + sum(dst2r)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)
    dpi1r <- as.vector((dst1r * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st1 * (sum(dst0r) + sum(dst1r) + sum(dst2r)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)
    dpi2r <- as.vector((dst2r * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st2 * (sum(dst0r) + sum(dst1r) + sum(dst2r)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)

    # MTTDL
    davar <- dpi0r+sum(dpi1r)
    dmttdlr <- (davar * t2)/(1-res$ava)^2
    #cat("# r", davar, dmttdlr)

    #3 t1
    dE1t1 <- param.q$Q11 %*% expm(param.q$Q11*t1)
    dbE1t1 <- expm(param.q$Q11*t1)
    dE2t1 <- 0
    dbE2t1 <- 0

    deP11t1 <- dE1t1 %*% (param.q$P11 + 1/(n*lambda) * param.q$P10 %*% param.q$Q01) + dbE1t1 %*% param.q$Q11
    deP12t1 <- dbE1t1 %*% param.q$Q12
    deP21t1 <- zeroM(1,2)
    deP22t1<- 0

    dPemct1 <- rbind(
          cbind(deP11t1,deP12t1),
          cbind(deP21t1,0))

    # Steady-satte probability
    stres <- ctmc.st(eyeM(dim(res$Pemc)[1]) - res$Pemc)
    stopifnot(stres$convergence)
    dpiemct1 <- stres$x

    ind1 <- 1:dim(deP11t1)[1]
    ind2 <- max(ind1) + 1:dim(deP12t1)[2]

    dpiemc1t1 <- dpiemct1[ind1]
    dpiemc2t1 <- dpiemct1[ind2]

    # sojourn time
    dst0t1 <- as.vector((1/(n*lambda)) * (dpiemc1t1 %*% res$E1 %*% param.q$P10 + res$piemc1 %*% dE1t1 %*% param.q$P10 + dpiemc2t1))
    dst1t1 <- as.vector(dpiemc1t1 %*% res$bE1 + res$piemc1 %*% dbE1t1)
    dst2t1 <- as.vector(dpiemc2t1 * t2)

    # steady-state probabilities
    dpi0t1 <- as.vector((dst0t1 * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st0 * (sum(dst0t1) + sum(dst1t1) + sum(dst2t1)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)
    dpi1t1 <- as.vector((dst1t1 * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st1 * (sum(dst0t1) + sum(dst1t1) + sum(dst2t1)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)
    dpi2t1 <- as.vector((dst2t1 * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st2 * (sum(dst0t1) + sum(dst1t1) + sum(dst2t1)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)

    # MTTDL
    davat1 <- dpi0t1+sum(dpi1t1)
    dmttdlt1 <- (davat1 * t2)/(1-res$ava)^2
    #cat("# t1", davat1, dmttdlt1)

    #4 t2
    dE1t2 <- zeroM(2,2)
    dbE1t2 <- zeroM(2,2)
    dE2t2 <- 0
    dbE2t2 <- 1

    deP11t2 <- zeroM(2,2)
    deP12t2 <- zeroM(2,1) 
    deP21t2 <- zeroM(1,2)
    deP22t2<- 0

    dPemct2 <- rbind(
          cbind(deP11t2,deP12t2),
          cbind(deP21t2,0))

    # Steady-satte probability
    stres <- ctmc.st(eyeM(dim(res$Pemc)[1]) - res$Pemc)
    stopifnot(stres$convergence)
    dpiemct2 <- stres$x

    ind1 <- 1:dim(deP11t2)[1]
    ind2 <- max(ind1) + 1:dim(deP12t2)[2]

    dpiemc1t2 <- dpiemct2[ind1]
    dpiemc2t2 <- dpiemct2[ind2]

    # sojourn time
    dst0t2 <- as.vector((1/(n*lambda)) * (dpiemc1t2 %*% res$E1 %*% param.q$P10 + dpiemc2t2))
    dst1t2 <- as.vector(dpiemc1t2 %*% res$bE1)
    dst2t2 <- as.vector(dpiemc2t2 * t2 + res$piemc2)

    # steady-state probabilities
    dpi0t2 <- as.vector((dst0t2 * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st0 * (sum(dst0t2) + sum(dst1t2) + sum(dst2t2)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)
    dpi1t2 <- as.vector((dst1t2 * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st1 * (sum(dst0t2) + sum(dst1t2) + sum(dst2t2)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)
    dpi2t2 <- as.vector((dst2t2 * (sum(res$st0) + sum(res$st1) + sum(res$st2)) - res$st2 * (sum(dst0t2) + sum(dst1t2) + sum(dst2t2)))/((sum(res$st0) + sum(res$st1) + sum(res$st2)))^2)

    # MTTDL
    davat2 <- dpi0t2+sum(dpi1t2)
    dmttdlt2 <- (davat2 * t2)/(1-res$ava)^2
    #cat("# t2", davat2, dmttdlt2)
    
    list(davar=davar,dmttdlr=dmttdlr,davat1=davat1,dmttdlt1=dmttdlt1,davat2=davat2,dmttdlt2=dmttdlt2,dE1r=dE1r,dbE1r=dbE1r,dE1t1=dE1t1,dbE1t1=dbE1t1,dpi0r=dpi0r,dpi0t1=dpi0t1)
  
}

# added on Mon Apr 17, 2023
comp.sens.mttf.raid6 <- function(t1, param.q, param.dq) {
  
  tQ00 <- solve(-param.q$Q00)
  
  # E_1 and bE1
  t2 <- 24
  res <- comp.mttdl.raid6(t1, t2, param.q)
  
  # dE1r, dbE1r, dE1t1, dbE1t1
  sens <- comp.sens.raid6(t1, t2, param.q, param.dq)
  
  # dXr, dYr, dM0r, dXt1, dYt1, dM0t1
  X <- eyeM(dim(res$E1 %*% param.q$P11)[1]) - res$E1 %*% param.q$P11 - res$E1 %*% param.q$P10 %*% tQ00 %*% param.q$Q01
  Y <- res$bE1 %*% vone(dim(res$bE1)[1]) + res$E1 %*% param.q$P10 %*% tQ00
  dXr <- -sens$dE1r %*% param.q$P11 - sens$dE1r %*% param.q$P10 %*% tQ00 %*% param.q$Q01 - res$E1 %*% param.q$P10 %*% (tQ00 %*% param.dq$dq00r %*% tQ00 %*% param.q$Q01 + tQ00 %*% param.dq$dq01r)
  dYr <- sens$dbE1r %*% vone(dim(sens$dbE1r)[1]) + sens$dE1r %*% param.q$P10 %*% tQ00 + res$E1 %*% param.q$P10 %*% tQ00 %*% param.dq$dq00r %*% tQ00
  dM0r <- tQ00 %*% param.dq$dq00r %*% tQ00 + tQ00 %*% param.dq$dq00r %*% tQ00 %*% param.q$Q01 %*% solve(X) %*% Y + tQ00 %*% (param.dq$dq01r %*% solve(X) %*% Y + param.q$Q01 %*% (solve(X) %*% dYr - solve(X) %*% dXr %*% solve(X) %*% Y))
  
  dXt1 <- -sens$dE1t1 %*% param.q$P11 - sens$dE1t1 %*% param.q$P10 %*% tQ00 %*% param.q$Q01
  dYt1 <- sens$dbE1t1 %*% vone(dim(sens$dbE1t1)[1]) + sens$dE1t1 %*% param.q$P10 %*% tQ00
  dM0t1 <- tQ00 %*% param.q$Q01 %*% (solve(X) %*% dYt1 - solve(X) %*% dXt1 %*% solve(X) %*% Y)
  
  # dmttfr, dmmtft1
  M0 <- tQ00 + tQ00 %*% param.q$Q01 %*% solve(X) %*% Y
  #dmttfr <- sens$dpi0r %*% M0 + res$pi0 %*% dM0r
  #dmttft1 <- sens$dpi0t1 %*% M0 + res$pi0 %*% dM0t1
  dmttfr <- dM0r
  dmttft1 <- dM0t1
  
  # cMTTF
  cMTTF <- M0 / (365 * 24)
  # dAmttf
  dAmttf <- (sens$davar * dmttfr + sens$davat1 * dmttft1) / (dmttfr^2 + dmttft1^2)

  list(cMTTF=cMTTF,dAmttf=dAmttf,davar=sens$davar,dmttdlr=sens$dmttdlr,davat1=sens$davat1,dmttdlt1=sens$dmttdlt1,davat2=sens$davat2,dmttdlt2=sens$dmttdlt2)
  
}

# investigate the effect of t1
comp.mttdl.raid6.t1 <- function(t2, param) {
  
  t1.list <- numeric(0)
  mttdl.list <- numeric(0)
  ava.list <- numeric(0)
  
  t1.array <- c(2:24)
  for (t1 in t1.array) {
    
    tQ00 <- solve(-param$Q00)
  
   # E_i
    E1 <- expm(param$Q11*t1)
    E2 <- exp(param$Q22*t2)
    #bE11 <- mexpAm.unif(A=Q11, t=t1, m=Q11, cumulative=TRUE)$cx
    #bE12 <- mexpAm.unif(A=Q11, t=t1, m=Q12, cumulative=TRUE)$cx
    bE1 <- mexpAm.unif(A=param$Q11, t=t1, cumulative=TRUE)$cx
    bE2 <- t2
  
    # eP_{i,j}
    eP11 <- E1 %*% param$P11 + E1 %*% param$P10 %*% tQ00 %*% param$Q01 + bE1 %*% param$Q11
    eP12 <- bE1 %*% param$Q12
    eP21 <- E2 %*% param$P20 * tQ00 %*% param$Q01
  
    Pemc <- rbind(
      cbind(eP11,eP12),
      cbind(eP21,0))
  
    # Steady-satte probability
    stres <- ctmc.st(Pemc - eyeM(dim(Pemc)[1]))
    stopifnot(stres$convergence)
    piemc <- stres$x
  
    ind1 <- 1:dim(eP11)[1]
    ind2 <- max(ind1) + 1:dim(eP12)[2]
  
    piemc1 <- piemc[ind1]
    piemc2 <- piemc[ind2]
  
    # sojourn time
    st0 <- as.vector(piemc1 %*% E1 %*% param$P10 %*% tQ00 + piemc2 %*% E2 %*% param$P20 %*% tQ00)
    st1 <- as.vector(piemc1 %*% bE1)
    st2 <- as.vector(piemc2 %*% bE2)
  
    # steady-state probabilities
    pi0 <- as.vector(st0/(sum(st0) + sum(st1) + sum(st2)))
    pi1 <- as.vector(st1/(sum(st0) + sum(st1) + sum(st2)))
    pi2 <- as.vector(st2/(sum(st0) + sum(st1) + sum(st2)))
  
    # MTTDL
    ava <- pi0+sum(pi1)
    mttdl <- (ava * t2/(1-ava)) / (365 * 24) 
    cat("# ", ava, mttdl)
  
    t1.list <- c(t1.list, t1)
    mttdl.list <- c(mttdl.list, mttdl)
    ava.list <- c(ava.list, ava)
  }
  
  res <- cbind(t1.list, ava.list, mttdl.list)
  write.table(res, "~/ieeedata/PE23-11/raid6/n100/txt/rest1.txt", row.names = F, col.names = F, append = F)
  
}

