### Subset ranges by CIGAR op.
.subset_by_op <- function(ir, ops) {
    stopifnot(is(ir, "IRanges"), !is.null(names(ir)), is.character(ops))
    ir[names(ir) %in% ops]
}

### Drop empty ranges.
.drop_empty_ranges <- function(ir) {
    stopifnot(is(ir, "IRanges"))
    ir[width(ir) != 0L]
}

### A customized reduce() that collapses the names and the "oplen"
### metadata column of the supplied IRanges object.
.reduce2 <- function(ir) {
    stopifnot(is(ir, "IRanges"))
    ans <- reduce(ir, with.revmap=TRUE)
    revmap <- mcols(ans)[ , "revmap"]
    names(ans) <- unstrsplit(extractList(names(ir), revmap))
    ans_oplen <- extractList(mcols(ir)[ , "oplen"], revmap)
    mcols(ans) <- DataFrame(oplen=ans_oplen)
    ans
}

.do_test_cigars_as_ranges_along_ref <- function(cigars, lmmpos, ops, expected)
{
    expected <- lapply(expected, .subset_by_op, ops)

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, ops=ops,
                                          with.ops=TRUE, with.oplens=TRUE)
    expect_true(is(current, "CompressedIRangesList"))
    expect_equal(length(current), length(expected))
    expect_null(names(current))
    for (i in seq_along(current))
        expect_identical(current[[i]], expected[[i]])

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, ops=ops,
                                          reduce.ranges=TRUE,
                                          with.ops=TRUE, with.oplens=TRUE)
    expect_true(is(current, "CompressedIRangesList"))
    expect_equal(length(current), length(expected))
    expect_null(names(current))
    for (i in seq_along(current))
        expect_identical(current[[i]], .reduce2(expected[[i]]))

    expected <- lapply(expected, .drop_empty_ranges)

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, ops=ops,
                                          drop.empty.ranges=TRUE,
                                          with.ops=TRUE, with.oplens=TRUE)
    expect_true(is(current, "CompressedIRangesList"))
    expect_equal(length(current), length(expected))
    expect_null(names(current))
    for (i in seq_along(current))
        expect_identical(current[[i]], expected[[i]])

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, ops=ops,
                                          drop.empty.ranges=TRUE,
                                          reduce.ranges=TRUE,
                                          with.ops=TRUE, with.oplens=TRUE)
    expect_true(is(current, "CompressedIRangesList"))
    expect_equal(length(current), length(expected))
    expect_null(names(current))
    for (i in seq_along(current))
        expect_identical(current[[i]], .reduce2(expected[[i]]))
}

test_that("cigars_as_ranges_along_ref()", {
    cigars <- c("30M5000N10M", "50M4S", "90=10X5I50M10D40M", "18M10I22M", "99I")
    ir1 <- IRanges(c(M="1-30", N="31-5030", M="5031-5040"),
                   oplen=c(30L, 5000L, 10L))
    ir2 <- IRanges(c(M="1-50", S="51-50"), oplen=c(50L, 4L))
    ir3 <- IRanges(c(`=`="1-90", X="91-100", I="101-100",
                     M="101-150", D="151-160", M="161-200"),
                   oplen=c(90L, 10L, 5L, 50L, 10L, 40L))
    ir4 <- IRanges(c(M="1-18", I="19-18", M="19-40"), oplen=c(18L, 10L, 22L))
    ir5 <- IRanges(c(I="1-0"), oplen=99L)

    offsets <- c(100L, 200L, 1000L, 300L, 1200L)
    ir1 <- shift(ir1, offsets[[1]])
    ir2 <- shift(ir2, offsets[[2]])
    ir3 <- shift(ir3, offsets[[3]])
    ir4 <- shift(ir4, offsets[[4]])
    ir5 <- shift(ir5, offsets[[5]])

    lmmpos <- offsets + 1

    ## Without specifying 'f' (in which case an **unnamed**
    ## CompressedIRangesList object is returned).

    list_of_ops <- list(
        CIGAR_OPS,
        c("M", "=", "X", "I", "D"),
        c("M", "=", "X", "I"),
        c("M", "=", "X", "D"),
        c("M", "N", "I"),
        c("M", "N", "S"),
        c("M", "I"),
        c("N", "I"),
        c("I", "S"),
        "M",
        "=",
        "I",
        character(0)
    )
    expected <- list(ir1, ir2, ir3, ir4, ir5)
    for (ops in list_of_ops)
        .do_test_cigars_as_ranges_along_ref(cigars, lmmpos, ops, expected)

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos,
                                          with.oplens=TRUE)
    expect_true(is(current, "CompressedIRangesList"))
    expect_equal(length(current), length(cigars))
    expect_identical(current[[1]], unname(ir1))
    expect_identical(current[[2]], unname(ir2))
    expect_identical(current[[3]], unname(ir3))
    expect_identical(current[[4]], unname(ir4))
    expect_identical(current[[5]], unname(ir5))

    ## Specifying 'f' (in which case a **named** SimpleIRangesList object
    ## is returned).

    rnames <- factor(c("chr6", "chr6", "chr2", "chr6", "chr2"),
                     levels=c("chr2", "chr6"))

    names(ir1) <- names(ir2) <- names(ir3) <- names(ir4) <- names(ir5) <- NULL
    mcols(ir1) <- mcols(ir2) <- mcols(ir3) <- mcols(ir4) <- mcols(ir5) <- NULL

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, f=rnames)
    expect_true(is(current, "SimpleIRangesList"))
    expect_identical(names(current), levels(rnames))
    expect_identical(current[[1]], c(ir3, ir5))
    expect_identical(current[[2]], c(ir1, ir2, ir4))

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, f=rnames,
                                          reduce.ranges=TRUE)
    expect_true(is(current, "SimpleIRangesList"))
    expect_identical(names(current), levels(rnames))
    expect_identical(current[[1]], c(reduce(ir3), reduce(ir5)))
    expect_identical(current[[2]], c(reduce(ir1), reduce(ir2), reduce(ir4)))

    ops <- c("M", "=", "X", "I", "D")
    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, f=rnames,
                                          ops=ops)
    expect_true(is(current, "SimpleIRangesList"))
    expect_identical(names(current), levels(rnames))
    expect_identical(current[[1]], c(ir3, ir5))
    expect_identical(current[[2]], c(ir1[-2], ir2[1], ir4))

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, f=rnames,
                                          ops=ops,
                                          reduce.ranges=TRUE)
    expect_true(is(current, "SimpleIRangesList"))
    expect_identical(names(current), levels(rnames))
    expect_identical(current[[1]], c(reduce(ir3), reduce(ir5)))
    expect_identical(current[[2]],
                     c(reduce(ir1[-2]), reduce(ir2[1]), reduce(ir4)))

    ops <- c("M", "=", "X", "I")
    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, f=rnames,
                                          ops=ops)
    expect_true(is(current, "SimpleIRangesList"))
    expect_identical(names(current), levels(rnames))
    expect_identical(current[[1]], c(ir3[-5], ir5))
    expect_identical(current[[2]], c(ir1[-2], ir2[1], ir4))

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, f=rnames,
                                          ops=ops,
                                          reduce.ranges=TRUE)
    expect_true(is(current, "SimpleIRangesList"))
    expect_identical(names(current), levels(rnames))
    expect_identical(current[[1]], c(reduce(ir3[-5]), reduce(ir5)))
    expect_identical(current[[2]],
                     c(reduce(ir1[-2]), reduce(ir2[1]), reduce(ir4)))

    current <- cigars_as_ranges_along_ref(cigars, lmmpos=lmmpos, f=rnames,
                                          ops=character(0))
    expect_true(is(current, "SimpleIRangesList"))
    expect_identical(lengths(current),
                     setNames(integer(length(levels(rnames))), levels(rnames)))
})

