source("SMatrix.R")

setClass("kSparseMatrix", representation(x="list", OP="character"), contains="sSparseMatrix")

kronMatrix <- function(..., op = "*") {
  x <- list(...)
  mdims <- c(1,1)
  for (m in x) {
    stopifnot(op == "*" || (op == "+" && dim(m)[1] == dim(m)[2]))
    mdims <- mdims * dim(m)
  }
  km <- new("kSparseMatrix")
  km@x <- x
  km@MDim <- as.integer(mdims)
  km@OP <- op
  km
}

convKtoM <- function(km) {
  mat <- c(1)
  for (m in km@x) {
    mat <- kronecker(X=mat, Y=as(m, "TsparseMatrix"))
  }
  mat
}

convKPtoM <- function(km) {
  stopifnot(km@MDim[1] == km@MDim[2])
  mat <- Matrix(0, km@MDim[1], km@MDim[2])
  left <- 1
  right <- km@MDim[1]
  for (m in km@x) {
    right <- right / dim(m)[1]
    leftI <- sparseMatrix(i=1:left, j=1:left, x=rep.int(x=1, times=left), giveCsparse=TRUE)
    rightI <- sparseMatrix(i=1:right, j=1:right, x=rep.int(x=1, times=right))
    mat <- mat + kronecker(X=kronecker(X=leftI, Y=as(m, "TsparseMatrix")), Y=rightI, giveCsparse=TRUE)
    left <- left * dim(m)[1]
  }
  mat
}

setAs("kSparseMatrix", "TsparseMatrix", function(from, to) {
  if (from@OP == "*") {
    as(convKtoM(from), "TsparseMatrix")
  } else if (from@OP == "+") {
    as(convKPtoM(from), "TsparseMatrix")
  } else {
    stop("Slot OP is neither '*' nor '+'.")
  }
})

setAs("kSparseMatrix", "CsparseMatrix", function(from, to) {
  if (from@OP == "*") {
    as(convKtoM(from), "CsparseMatrix")
  } else if (from@OP == "+") {
    as(convKPtoM(from), "CsparseMatrix")
  } else {
    stop("Slot OP is neither '*' nor '+'.")
  }
})

setAs("kSparseMatrix", "RsparseMatrix", function(from, to) {
  if (from@OP == "*") {
    as(convKtoM(from), "RsparseMatrix")
  } else if (from@OP == "+") {
    as(convKPtoM(from), "RsparseMatrix")
  } else {
    stop("Slot OP is neither '*' nor '+'.")
  }
})
