
test_that("narrow_cigars_along_ref()", {
    cigar <- c("25M4D10M", "6S17M6I3M3S")

    current <- narrow_cigars_along_ref(cigar)
    expected <- c("25M4D10M", "17M6I3M")
    attr(expected, "rshift") <- c(0L, 0L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_ref(cigar, start=3, end=-3)
    expected <- c("23M4D8M", "15M6I1M")
    attr(expected, "rshift") <- c(2L, 2L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_ref(cigar, start=7, end=-4)
    expected <- c("19M4D7M", "11M")
    attr(expected, "rshift") <- c(6L, 6L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_ref(cigar, start=8, end=-5)
    expected <- c("18M4D6M", "9M")
    attr(expected, "rshift") <- c(7L, 7L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_ref(cigar[1], start=26, end=-10)
    expected <- "1M"
    attr(expected, "rshift") <- 29L
    expect_identical(current, expected)
})

test_that("narrow_cigars_along_query()", {
    cigars <- c("25M4D10M", "6S17M6I3M3S")

    current <- narrow_cigars_along_query(cigars)
    expected <- cigars
    attr(expected, "rshift") <- c(0L, 0L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_query(cigars, start=3, end=-3)
    expected <- c("23M4D8M", "4S17M6I3M1S")
    attr(expected, "rshift") <- c(2L, 0L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_query(cigars, start=7, end=-4)
    expected <- c("19M4D7M", "17M6I3M")
    attr(expected, "rshift") <- c(6L, 0L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_query(cigars, start=8, end=-5)
    expected <- c("18M4D6M", "16M6I2M")
    attr(expected, "rshift") <- c(7L, 1L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_query(cigars, start=25)
    expected <- c("1M4D10M", "5I3M3S")
    attr(expected, "rshift") <- c(24L, 17L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_query(cigars, start=26)
    expected <- c("10M", "4I3M3S")
    attr(expected, "rshift") <- c(29L, 17L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_query(cigars, start=26, end=-8)
    expected <- c("3M", "3I")
    attr(expected, "rshift") <- c(29L, 17L)
    expect_identical(current, expected)

    current <- narrow_cigars_along_query(cigars, start=26, end=-10)
    expected <- c("1M", "1I")
    attr(expected, "rshift") <- c(29L, 17L)
    expect_identical(current, expected)
})

