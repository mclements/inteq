#' Solve a Volterra equation of the first kind
#' @param k kernel function of two time scales
#' @param f left hand side (free) function with f(a)=0
#' @param a lower bound of the integral
#' @param b upper bound of the integral
#' @param num integer for the number of evaluation points
#' @param method string for the method
#' @return data-frame with evaluation points 'sgrid' and calculations 'ggrid'
#' @examples
#' k <- function(s,t) {
#'         cos(t-s)
#' }
#' trueg <- function(s) {
#'     (2+s**2)/2
#' }
#' 
#' res <- volterra_solve(k,a=0,b=1,num=1000)
#' 
#' plot(
#'     res$sgrid, res$ggrid,
#'     type = "l",
#'     col = "blue",
#'     xlim = c(0, 1),
#'     #ylim = c(-1, 1),
#'     xlab = "s",
#'     ylab = "g(s)",
#'     main = "Volterra Equation Solution first kind"
#' )
#' # add the true solution
#' lines(res$sgrid, trueg(res$sgrid), col = "red", lty = 2)
#' legend( 
#'     "topright",
#'     legend = c("Estimated Value", "True Value"),
#'     col = c("blue", "red"),
#'     lty = c(1, 2)
#' )
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
#' @param num integer for the number of evaluation points
#' @param method string for the method
#' @return data-frame with evaluation points 'sgrid' and calculated values 'ggrid' 
#' @examples
#' k <- function(s,t) {
#'         0.5 * (t-s)** 2 * exp(t-s)
#' }
#' free <- function(t) {
#'     0.5 * t**2 * exp(-t)
#' }
#' true <- function(t) {
#'     1/3 * (1 - exp(-3*t/2) * (cos(sqrt(3)/2*t) + sqrt(3) * sin(sqrt(3)/2*t)))
#' }
#' 
#' res <- volterra_solve2(k,free,a=0,b=6,num=100)
#' 
#' plot(
#'     res$sgrid, res$ggrid,
#'     type = "l",
#'     col = "blue",
#'     xlim = c(0, 6),
#'     #ylim = c(-1, 1),
#'     xlab = "s",
#'     ylab = "g(s)",
#'     main = "Volterra Equation Solution second kind"
#' )
#' # add the true solution
#' lines(res$sgrid, true(res$sgrid), col = "red", lty = 2)
#' legend( 
#'     "topright",
#'     legend = c("Estimated Value", "True Value"),
#'     col = c("blue", "red"),
#'     lty = c(1, 2)
#' )
#' 
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
