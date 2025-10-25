
test_that("validate_cigars()", {
    expect_null(validate_cigars(c("11M", "11=", "11X")))

    res <- validate_cigars(c("11M", "11= ", "11X"))
    expect_true(isSingleString(res))
    expect_true(grepl("element 2 is invalid", res))
})

test_that("explode_cigar_ops() and explode_cigar_oplens()", {
    cigars <- c("1S66M2H", "3M22N1M", "5M2D3M1X3M", "10M2I8M")

    res <- explode_cigar_ops(cigars)
    expected <- list(c("S", "M", "H"),
                     c("M", "N", "M"),
                     c("M", "D", "M", "X", "M"),
                     c("M", "I", "M"))
    expect_identical(res, expected)
    res <- explode_cigar_oplens(cigars)
    expected <- list(c(1L, 66L, 2L),
                     c(3L, 22L, 1L),
                     c(5L, 2L, 3L, 1L, 3L),
                     c(10L, 2L, 8L))

    res <- explode_cigar_ops(cigars, ops=c("M", "X"))
    expected <- list(c("M"),
                     c("M", "M"),
                     c("M", "M", "X", "M"),
                     c("M", "M"))
    expect_identical(res, expected)
    res <- explode_cigar_oplens(cigars, ops=c("M", "X"))
    expected <- list(c(66L),
                     c(3L, 1L),
                     c(5L, 3L, 1L, 3L),
                     c(10L, 8L))
})

