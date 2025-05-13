library(inteq)

expect_eps <- function(expr, value, eps=1e-7)
    expect_lt(max(abs(expr-value)),eps)

context("helpers")
##
test_that("indexing", {
    expect_true(all(indexing(1:3, 1:3, 4) == c(1,6,11)))
})
test_that("diag_ext", {
    expect_true(all(diag_ext(1:3,1L,TRUE) ==
                    c(0,1,0,0,
                      1,0,2,0,
                      0,2,0,3,
                      0,0,3,0)))
})
test_that("makeH", {
    expect_true(all(makeH(5) ==
                    c(1,-2,1,0,0,
                      -2,5,-4,1,0,
                      1,-4,6,-4,1,
                      0,1,-4,5,-2,
                      0,0,1,-2,1)))
})
test_that("simpson", {
    expect_true(all(simpson(7L) == c(1,4,2,4,2,4,1)/3))
})
test_that("smooth", {
    expect_true(all(smooth(1:5) == c(1,1.5,2.5,3.5,4.5)))
})


context("fredholm")
##
test_that("fredholm_solve", {
    k <- function(s, t) {
    ifelse(abs(s-t)<=3,1 + cos(pi*(t-s)/3), 0)
    }
    f <- function(s) {
    sp <- abs(s)
    sp3 <- sp * pi / 3
    ((6 - sp) * (2 + cos(sp3)) + (9 / pi) * sin(sp3)) / 2
    }
    trueg <- function(s) {
        k(0,s)
    }
    res = fredholm_solve(k,f,-3,3,1001L,smin=-6,smax=6,snum=2001L,gamma=0.01)
    expect_true(sum((res$ggrid - sapply(res$ygrid,trueg))**2) < 1e-6)
})

context("volterra")
##
test_that("volterra_solve", {
    k <- function(s,t) {
        cos(t-s)
    }
    trueg <- function(s) {
        (2+s**2)/2
    }
        res = fredholm_solve(k,f,-3,3,1001L,smin=-6,smax=6,snum=2001L,gamma=0.01)
        expect_true(sum((volterra_solve(k,a=0,b=1,num=1000)$ggrid - trueg(seq(0,1,length.out=1000)))**2) < 1e-3
    )
})

test_that("volterra_solve2", {
        k <- function(s,t) {
        0.5 * (t-s)** 2 * exp(t-s)
    }
    free <- function(t) {
        0.5 * t**2 * exp(-t)
    }
    true <- function(t) {
        #like this python cde : 1/3*(1-np.exp(-3*t/2)*(np.cos(np.sqrt(3)/2 * t) + np.sqrt(3) * np.sin(np.sqrt(3)/2 * t)))
        1/3 * (1 - exp(-3*t/2) * (cos(sqrt(3)/2*t) + sqrt(3) * sin(sqrt(3)/2*t)))
    }
    espect_true(sum((volterra_solve2(k,free,a=0,b=6,num=100)$ggrid - true(seq(0,6,length.out=100)))**2) < 1e-11)
})


