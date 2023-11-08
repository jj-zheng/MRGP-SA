
source("SMatrix.R")

setClass("fSparseMatrix", representation(x="list", OP="character"), contains="sSparseMatrix")

plusMatrix <- function(...) {
  x <- list(...)
  mdims <- dim(x[[1]])
  for (m in x) {
    stopifnot(mdims[1] == dim(m)[1] && mdims[2] == dim(m)[2])
  }
  fm <- new("fSparseMatrix")
  fm@x <- x
  fm@MDim <- as.integer(mdims)
  fm@OP <- "+"
  fm
}

convFPtoM <- function(fm) {
  mat <- Matrix(0, fm@MDim[1], fm@MDim[2])
  for (m in fm@x) {
    mat <- mat + as(m, "TsparseMatrix")
  }
  mat
}

setAs("fSparseMatrix", "TsparseMatrix", function(from, to) {
  switch(from@OP,
    "+" = as(convFPtoM(from), "TsparseMatrix"),
    stop("Slot OP is neither '*' nor '+'."))
})

setAs("fSparseMatrix", "CsparseMatrix", function(from, to) {
  switch(from@OP,
    "+" = as(convFPtoM(from), "CsparseMatrix"),
    stop("Slot OP is neither '*' nor '+'."))
})

setAs("fSparseMatrix", "RsparseMatrix", function(from, to) {
  switch(from@OP,
    "+" = as(convFPtoM(from), "RsparseMatrix"),
    stop("Slot OP is neither '*' nor '+'."))
})
