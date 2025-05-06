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

