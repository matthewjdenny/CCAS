test_that("Log space multinomial sampler works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)
    # have to assign variables to the global environment so they will run in
    # automatic tests
    dist <- round(runif(100, min = -500, max = -400))
    seed <- 123456
    # get the result
    result <- test_internal_functions(
        Test_Log_Space_Multinomial_Sampler = TRUE,
        distribution = dist,
        seed = seed)

})

