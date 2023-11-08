# transition probability matrix
mexp.sens.unif <- function(Q, t, dQ) {
	P <- unif(Q)$P
	qv <- unif(Q)$qv

	invI <- eyeM(dim(P)[1])
	F <- invI
	bF <- list(F)

	U <- 100
	for (i in 1:U) {
		F <- F %*% P 
		bF <- c(bF, list(F))
	}

	A = dQ
	X <- zeroM(dim(P))
	bX <- list(X)
	for (i in 1:U) {
		X <- bF[[i+1]] %*% A + X %*% P
		bX <- c(bX, list(X))
	}

	S <- zeroM(dim(P))
	for (i in 0:U) {
		S <- S + 1/qv * dpois(i+1,qv * t) * bX[[i+1]]
	}
	
	list(Q=Q, dQ=dQ, S=S)
}


cmexp.sens.unif.pdf.inf <- function(Q, dQ) {
	T <- 100
	X <- zeroM(dim(Q))
	for (t in 1:T) {
		X <- X + mexp.sens.unif(Q, t, dQ)$S * dgamma(t, a, b)
	}
	
	list(T=T, X=X)
}

cmexp.sens.unif.t <- function(Q, t, dQ) {
	T <- t
	X <- zeroM(dim(Q))
	for (t in 1:T) {
		X <- X + mexp.sens.unif(Q, t, dQ)$S
	}
	
	list(T=T, X=X)
}

cmexp.sens.unif.t.inf <- function(Q, dQ) {
	TT <- 10
	Y <- zeroM(dim(Q))
	for (tt in 1:TT) {
		X <- zeroM(dim(Q))
		for (t in 1:tt) {
			X <- X + mexp.sens.unif(Q, t, dQ)$S
		}
		X <- X * dgamma(t, a, b)
		Y <- Y + X
	}

	list(T=TT, Y=Y)	
}

dgammaa <- function(t){
		t^{a-1}*exp(-t)*log(t)
}

pgamma.a <- function(t) {
	((b^a*log(b)*t^{a-1}+b^a*t^{a-1}*log(t))*exp(-b*t)*gamma(a)-b^a*t^{a-1}*exp(-b*t)*integrate(dgammaa, 0, Inf)$value)/(gamma(a))^2
}

pgamma.b <- function(t) {
	((a*b^{a-1}-t*b^a)*t^{a-1}*exp(-b*t))/gamma(a)
}

cmexp.sens.ppdfa.inf <- function(Q) {
	T <- 19
	X <- zeroM(dim(Q))
	for (t in 1:T) {
		X <- X + expm(Q*t) * pgamma.a(t)
	}
	
	list(T=T, X=X)
}

cmexp.sens.ppdfb.inf <- function(Q) {
	T <- 19
	X <- zeroM(dim(Q))
	for (t in 1:T) {
		X <- X + expm(Q*t) * pgamma.b(t)
	}
	
	list(T=T, X=X)
}

cmexp.sens.ppdfa.t.inf <- function(Q) {
	TT <- 19
	Y <- zeroM(dim(Q))
	for (tt in 1:TT) {
		X <- zeroM(dim(Q))
		for (t in 1:tt) {
			X <- X + expm(Q*t)
		}
		X <- X * pgamma.a(t)
		Y <- Y + X
	}

	list(T=TT, Y=Y)	
}

cmexp.sens.ppdfb.t.inf <- function(Q) {
	TT <- 19
	Y <- zeroM(dim(Q))
	for (tt in 1:TT) {
		X <- zeroM(dim(Q))
		for (t in 1:tt) {
			X <- X + expm(Q*t)
		}
		X <- X * pgamma.b(t)
		Y <- Y + X
	}

	list(T=TT, Y=Y)	
}

