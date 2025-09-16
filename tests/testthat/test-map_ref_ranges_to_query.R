
test_that("map_ref_ranges_to_query()", {
    start <- end <- 1:34
    cigars <- c("4M", "5M3I4M", "4M3D5M", "3M", "10M", "25M")
    lmmpos <- c(2, 11, 21, 31, 41, 5)
    df <- map_ref_ranges_to_query(start, end, cigars, lmmpos)
    expect_true(is.data.frame(df))
    expect_identical(colnames(df), c("start", "end", "from_hit", "to_hit"))
    expect_false(is.unsorted(df[ , "from_hit"]))

    revhits <- Hits(df[ , "to_hit"], df[ , "from_hit"],
                    nLnode=length(cigars), nRnode=length(start),
                    start=df[ , "start"], end=df[ , "end"],
                    sort.by.query=TRUE)

    ## Group input ranges per cigar, that is, for each cigar extract the
    ## input range indices that hit the cigar:
    rgidx_per_cigar <- as(revhits, "IntegerList")
    expect_equal(length(rgidx_per_cigar), length(cigars))
    starts_per_cigar <- relist(mcols(revhits)$start, rgidx_per_cigar)
    ends_per_cigar <- relist(mcols(revhits)$end, rgidx_per_cigar)

    ## Closely examine each group:

    ## cigar 4M / lmmpos = 2
    expect_identical(rgidx_per_cigar[[1L]], 2:5)
    expect_identical(starts_per_cigar[[1L]], 1:4)
    expect_identical(ends_per_cigar[[1L]], 1:4)

    ## cigar 5M3I4M / lmmpos = 11
    expect_identical(rgidx_per_cigar[[2L]], 11:19)
    expect_identical(starts_per_cigar[[2L]], c(1:5, 9:12))
    expect_identical(ends_per_cigar[[2L]], c(1:5, 9:12))

    ## cigar 4M3D5M / lmmpos = 21
    expect_identical(rgidx_per_cigar[[3L]], 21:32)
    expect_identical(starts_per_cigar[[3L]], c(1:5, 5L, 5L, 5:9))
    expect_identical(ends_per_cigar[[3L]], c(1:4, 4L, 4L, 4:9))

    ## cigar 3M / lmmpos = 31
    expect_identical(rgidx_per_cigar[[4L]], 31:33)
    expect_identical(starts_per_cigar[[4L]], 1:3)
    expect_identical(ends_per_cigar[[4L]], 1:3)

    ## cigar 10M / lmmpos = 41
    expect_identical(rgidx_per_cigar[[5L]], integer(0))
    expect_identical(starts_per_cigar[[5L]], integer(0))
    expect_identical(ends_per_cigar[[5L]], integer(0))

    ## cigar 25M / lmmpos = 5
    expect_identical(rgidx_per_cigar[[6L]], 5:29)
    expect_identical(starts_per_cigar[[6L]], 1:25)
    expect_identical(ends_per_cigar[[6L]], 1:25)

    ## The data.frame returned by fast_map_ref_ranges_to_query() is
    ## expected to always be the same as the data.frame returned by
    ## map_ref_ranges_to_query() **modulo** the order of the rows.
    ## However, if we use 'strictly.sort.hits=TRUE' the two data.frames
    ## are expected to be identical:
    df2 <- fast_map_ref_ranges_to_query(start, end, cigars, lmmpos,
                                        strictly.sort.hits=TRUE)
    expect_identical(df2, df)
})

