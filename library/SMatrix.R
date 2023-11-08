setClass("sSparseMatrix", representation(MDim="integer"))

setMethod("dim", signature(x = "sSparseMatrix"), function(x) x@MDim)
setMethod("as.matrix", signature(x = "sSparseMatrix"), function(x) as(x, "TsparseMatrix"))
