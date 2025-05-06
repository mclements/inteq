#' Solve a Volterra equation of the first kind
#' @param k kernel function of two time scales
#' @param f left hand side (free) function with f(a)=0
#' @param a lower bound of the integral
#' @param b upper bound of the integral
#' @param method string for the method
#' @export
volterra_solve =
    function(k, f = function(x) x, a = 0, b = 1, num = 1000L,
             method=c("midpoint","trapezoid")) {
        method = match.arg(method)
        sgrid = seq(a,b,length.out=num+1)[-1]
        ktril = outer(sgrid,sgrid,"k")
        ktril[upper.tri(ktril)] = 0
        if (method=="trapezoid") {
            diag(ktril) = diag(ktril)/2
            ktril[,1] = ktril[,1] + k(sgrid,0)/2 # should the zero be 'a'?
        }
        ggrid = solve(ktril, f(sgrid))
        data.frame(sgrid, ggrid = ggrid * num / (b-a))
}

## a = 0; b = 1; num=10L
## k = function(x,y) x+y
## sgrid = seq(a,b,length.out=num+1)[-1]
## ktril = outer(sgrid, sgrid, k)
## volterra_solve(k,num=10L)

#' Solve a Volterra equation of the second kind
#' @param k kernel function of two time scales
#' @param f left hand side (free) function with f(a)=0
#' @param a lower bound of the integral
#' @param b upper bound of the integral
#' @param method string for the method
#' @export
volterra_solve2 =
    function(k, f = function(x) x,
             a = 0, b = 1, num = 1001L,
             method=c("trapezoid","midpoint")) {
        method = match.arg(method)
        sgrid = seq(a,b,length.out=num)
        h = (b-a)/(num-1L)
        ktril = outer(sgrid,sgrid,k)
        ktril[upper.tri(ktril)] = 0
        if (method=="trapezoid") {
            ktril[,1] = ktril[,1]/2
            diag(ktril) = diag(ktril)/2
        }
        ktril[1,1] = 0
        ggrid = solve(diag(num) - h*ktril, f(sgrid))
        data.frame(sgrid, ggrid)
}
