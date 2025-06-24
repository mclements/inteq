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

