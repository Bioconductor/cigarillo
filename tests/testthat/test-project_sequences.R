
test_that("project_sequences()", {
    PROJECTION_SPACES <- cigarillo:::PROJECTION_SPACES

    cigar_extent_funs <- list(
        function(cigar) cigar_extent_along_ref(cigar),
        function(cigar) cigar_extent_along_ref(cigar,
                                               N.regions.removed=TRUE),
        function(cigar) cigar_extent_along_query(cigar),
        function(cigar) cigar_extent_along_query(cigar,
                                                 before.hard.clipping=TRUE),
        function(cigar) cigar_extent_along_query(cigar,
                                                 after.soft.clipping=TRUE),
        function(cigar) cigar_extent_along_pwa(cigar),
        function(cigar) cigar_extent_along_pwa(cigar,
                                               N.regions.removed=TRUE),
        function(cigar) cigar_extent_along_pwa(cigar, dense=TRUE)
    )

    cigars <- c("3H2S4M1D2M2I1M5N3M6H", "5M1I3M2D4M2S")

    projected_sequences <- list(
        BStringSet(c(A="AAAA-BBC.....DDD", B="AAAAABBB--CCCC")),
        BStringSet(c(A="AAAA-BBCDDD", B="AAAAABBB--CCCC")),
        BStringSet(c(A="++AAAABBiiCDDD", B="AAAAAiBBBCCCC++")),
        BStringSet(c(A="+++++AAAABBiiCDDD++++++", B="AAAAAiBBBCCCC++")),
        BStringSet(c(A="AAAABBiiCDDD", B="AAAAAiBBBCCCC")),
        BStringSet(c(A="AAAA-BBiiC.....DDD", B="AAAAAiBBB--CCCC")),
        BStringSet(c(A="AAAA-BBiiCDDD", B="AAAAAiBBB--CCCC")),
        BStringSet(c(A="AAAABBCDDD", B="AAAAABBBCCCC"))
    )

    for (i in seq_along(PROJECTION_SPACES))
        expect_identical(width(projected_sequences[[i]]),
                         cigar_extent_funs[[i]](cigars))

    project_sequences2 <- function(x, cigars, from, to)
        project_sequences(x, cigars, from=from, to=to, I.letter="i")

    identical_XStringSet_objects <- function(x, y)
    {
        ok1 <- identical(class(x), class(y))
        ok2 <- identical(width(x), width(y))
        ok3 <- all(x == y)
        ok4 <- identical(names(x), names(y))
        ok1 && ok2 && ok3 && ok4
    }

    for (i in seq_along(PROJECTION_SPACES))
        for (j in seq_along(PROJECTION_SPACES)) {
            expected <- projected_sequences[[j]]
            current <- project_sequences2(projected_sequences[[i]], cigars,
                                          from=PROJECTION_SPACES[[i]],
                                          to=PROJECTION_SPACES[[j]])
            expect_true(identical_XStringSet_objects(current, expected))
        }
})

