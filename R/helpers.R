
#' Internal function to convert (row,col) to vector index
indexing = function(row,col,nrows) {
    (col-1)*nrows+row
}

#' Internal function to generate diagonal matrices with possibly an offset with possibly mirrored diagonal
#' @param x vector
#' @param index integer offset index for the diagonal (can be negative)
diag_ext = function(x, index, mirror=FALSE) {
    stopifnot(is.integer(index))
    if (index==0L) return(diag(x)) # ignores mirror
    if (mirror) return(diag_ext(x,index) + diag_ext(x,-index))
    if (index<0L) return(t(diag_ext(x,-index)))
    n = length(x)
    dim = n+index
    diag = matrix(0, dim, dim)
    diag[indexing(1:n,(1:n)+index,dim)] = x
    ## for (i in 1:n)
    ##     diag[i,i+index-1L] = x[i]
    diag
}

#' To be completed:)
makeH = function(dim) {
    d2 = rep(1, dim-2L)
    d1 = c(-2, rep(-4, dim-3), -2)
    d0 = c(1,5,rep(6,dim-4),5,1)
    diag_ext(d1,1L,TRUE) + diag_ext(d2,2L,TRUE) + diag(d0)
}

#' To be completed:)
simpson = function(dim) {
    stopifnot(dim>2, dim %% 2 == 1)
    tiles = rep(c(4,2),dim %/% 2)
    c(1,tiles[1:(dim-2)],1)/3
}

#' To be completed:)
smooth = function(v) {
    dim = length(v)
    smat = diag(c(1,rep(0.5, dim-1L))) + diag_ext(rep(0.5, dim-1),-1L)
    smat %*% v
}
