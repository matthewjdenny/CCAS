test_that("Log space multinomial sampler works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)
    # have to assign variables to the global environment so they will run in
    # automatic tests

    k <- 100
    dist <- runif(k)
    u <- runif(1)
    bins <- cut(dist, k)
    uppers <- as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", levels(bins)))
    j <- findInterval(u, uppers, all.inside = TRUE)

    seed <- 123456
    # get the result
    result <- test_internal_functions(
        Test_Log_Space_Multinomial_Sampler = TRUE,
        distribution = log(dist), seed = seed, u = u)

    expect_that(j, equals(result))
})

