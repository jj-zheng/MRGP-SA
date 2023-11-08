library(RMarkov)
library(Matrix)

# failure rate
lambda <- 1/1000000

# disk nubmer
n <- 4

# Gen
Trebuild.dist <- 1
Trecon.dist <- 1

source("raid6_mat.R")
source("raid6_vec.R")

# mean time to disk rebuild
#t1 <- 2

# mean time to storage reconstruction
#t2 <- 24

makeQ <- function(G0G1E,G1G1E,G1G2E,G1G0P0,G1G1P0,G2G0P1) {

	# Q_{ij}
    Q00 <- -rowSums(G0G1E)
    Q01 <- G0G1E
    Q11 <- G1G1E + diag(-rowSums(G1G1E)) + diag(-rowSums(G1G2E))
    Q12 <- G1G2E
    Q22 <- 0

    # P_{i,j}
    P10 <- G1G0P0
    P11 <- G1G1P0
    P20 <- G2G0P1

    list(Q00=Q00,Q01=Q01,Q11=Q11,Q12=Q12,Q22=Q22,P10=P10,P11=P11,P20=P20)
}

param.q <- makeQ(G0G1E=G0G1E,G1G1E=G1G1E,G1G2E=G1G2E,G1G0P0=G1G0P0,G1G1P0=G1G1P0,G2G0P1=G2G0P1)

makedQ <- function(lambda,n) {

	 # first derivatives: n, r, t1, t2
	dq00n <- -lambda
	dq01n <- c(0,lambda)
	dq11n <- rbind(
				c(-lambda,0),
				c(lambda,-lambda))
	dq12n <- rbind(
				c(lambda),
				c(0))

	dq00r <- -n
	dq01r <- c(0,n)
	dq11r <- rbind(
				c(-(n-2),0),
				c((n-1),-(n-1)))
	dq12r <- rbind(
				c(n-2),
				c(0))

	dq00t1 <- 0
	dq01t1 <- c(0,0)
	dq11t1 <- zeroM(2)
	dq12t1 <- rbind(
				0,
				0)

	dq00t2 <- 0
	dq01t2 <- c(0,0)
	dq11t2 <- zeroM(2)
	dq12t2 <- rbind(
				0,
				0)

	list(dq00n=dq00n,dq01n=dq01n,dq11n=dq11n,dq12n=dq12n,dq00r=dq00r,dq01r=dq01r,dq11r=dq11r,dq12r=dq12r,dq00t1=dq00t1,dq01t1=dq01t1,dq11t1=dq11t1,dq12t1=dq12t1,dq00t2=dq00t2,dq01t2=dq01t2,dq11t2=dq11t2,dq12t2=dq12t2)
}

param.dq <- makedQ(lambda=lambda,n=n) 
