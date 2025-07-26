#' Solve a Fredholm equation of the first and second kind
#' @param k kernel function of two time scales
#' @param f left hand side function with f(a)=0
#' @param a lower bound of the grid for the integral approximation 
#' @param b upper bound of the grid for the integral approximation 
#' @param num number of points for the grid of the integral approximation 
#' @param smin lower bound of enforcement values for equation
#' @param smax upper bound of enforcement values for equation
#' @param snum number of points for the grid for the equation
#' @param gamma regularization parameter
#' @return data-frame with evaluation points 'ygrid' and calculations 'ggrid'
#' @examples
#' # Define the kernel function
#' k <- function(s, t) {
#'     ifelse(abs(s - t) <= 3, 1 + cos(pi * (t - s) / 3), 0)
#' }
#' 
#' # Define the right-hand side function
#' f <- function(s) {
#'     sp <- abs(s)
#'     sp3 <- sp * pi / 3
#'     ((6 - sp) * (2 + cos(sp3)) + (9 / pi) * sin(sp3)) / 2
#' }
#' 
#' # Define the true solution for comparison
#' trueg <- function(s) {
#'     k(0, s)
#' }
#' 
#' # Solve the Fredholm equation
#' res <- fredholm_solve(
#'     k, f, -3, 3, 1001L,
#'     smin = -6, smax = 6, snum = 2001L,
#'     gamma = 0.01
#' )
#' 
#' # Plot the results on the same graph using base graphics
#' plot(
#'     res$ygrid, res$ggrid,
#'     type = "l",
#'     col = "blue",
#'     xlim = c(-3, 3),
#'     #ylim = c(-1, 1),
#'     xlab = "s",
#'     ylab = "g(s)",
#'     main = "Fredholm Equation Solution"
#' )
#' # add the true solution
#' lines(res$ygrid, trueg(res$ygrid), col = "red", lty = 2)
#' legend( 
#'     "topright",
#'     legend = c("Estimated Value", "True Value"),
#'     col = c("blue", "red"),
#'     lty = c(1, 2)
#' )
#' @export
fredholm_solve =
    function(k, f = function(x) x, a = -0, b = 1, num = 41L, smin = 0, smax = 1, snum = 41L, gamma = 0.001) {

        #num needs to be odd to apply simpson's rule
        if (num %% 2 == 0) {
            warning("num should be odd for Simpson's rule, proceding with num +1")
            num = num + 1L
        }

        # grid for enforcement points
        sgrid = seq(smin,smax,length.out=snum+1)[-1]
        # grid for integration points
        ygrid = seq(a,b,length.out=num+1)[-1]

        weights = simpson(num)

        # quadrature matrix
        ksqur = outer(sgrid, ygrid, "k") %*% diag(weights)

        if(gamma != 0) {
            Hmat = makeH(num)
            AAgH = t(ksqur) %*% ksqur + gamma * Hmat
        } else {
            AAgH = t(ksqur) %*% ksqur
        }

        ggrid = solve(AAgH, t(ksqur) %*% f(sgrid))
        data.frame(ygrid, ggrid = ggrid * num / (b-a))
}

