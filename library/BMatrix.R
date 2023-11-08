
source("SMatrix.R")

setClass("bSparseMatrix", representation(i="integer", j="integer", Dim="integer", x="list", si="integer", sj="integer", ii="integer", jj="integer"), contains="sSparseMatrix")

blockMatrix <- function(i, j, x, dims, si, sj) {

  stopifnot(length(i) == length(j) && length(i) == length(x))

  if (missing(dims)) {
    if (missing(si)) {
      mi <- max(i)
    } else {
      mi <- length(si)
    }
    if (missing(sj)) {
      mj <- max(j)
    } else {
      mj <- length(sj)
    }
    dims <- c(mi, mj)
  }
  if (missing(si)) {
    si <- integer(dims[1])
  }
  if (missing(sj)) {
    sj <- integer(dims[2])
  }
  stopifnot(length(si) == dims[1], length(sj) == dims[2])

  for (u in 1:length(x)) {
    v <- i[u]
    si[v] <- dim(x[[u]])[1]
  }

  for (u in 1:length(x)) {
    v <- j[u]
    sj[v] <- dim(x[[u]])[2]
  }

  # check si sj
  stopifnot(all(si != 0) && all(sj != 0))

  mdims <- c(sum(si), sum(sj))
  ii <- cumsum(c(1,si))[1:dims[1]]
  jj <- cumsum(c(1,sj))[1:dims[2]]

  bm <- new("bSparseMatrix")
  bm@Dim <- as.integer(dims)
  bm@i <- as.integer(i)
  bm@j <- as.integer(j)
  bm@x <- x
  bm@si <- as.integer(si)
  bm@sj <- as.integer(sj)
  bm@ii <- as.integer(ii)
  bm@jj <- as.integer(jj)
  bm@MDim <- as.integer(mdims)
  bm
}

convBtoM <- function(bm) {
  mat <- Matrix(0, nrow = bm@MDim[1], ncol = bm@MDim[2])
  for (u in 1:length(bm@x)) {
    i <- bm@ii[bm@i[u]]
    j <- bm@jj[bm@j[u]]
    i2 <- i + dim(bm@x[[u]])[1] - 1
    j2 <- j + dim(bm@x[[u]])[2] - 1
    mat[i:i2,j:j2] <- as(bm@x[[u]], "TsparseMatrix")
  }
  mat
}

setAs("bSparseMatrix", "TsparseMatrix", function(from, to) {
    as(convBtoM(from), "TsparseMatrix")
})

setAs("bSparseMatrix", "CsparseMatrix", function(from, to) {
  as(convBtoM(from), "CsparseMatrix")
})

setAs("bSparseMatrix", "RsparseMatrix", function(from, to) {
  as(convBtoM(from), "RsparseMatrix")
})

setMethod("print", signature(x = "bSparseMatrix"), function(x, ...) {
  mat <- matrix(".", x@Dim[1], x@Dim[2])
  for (u in 1:length(x@x)) {
    mat[x@i[u], x@j[u]] <- gettextf(fmt="(%d,%d)", dim(x@x[[u]])[1], dim(x@x[[u]])[2])
  }
  mat <- data.frame(mat)
  rownames(mat) <- sapply(x@si, gettext)
  colnames(mat) <- sapply(x@sj, gettext)
  print(mat)
})
