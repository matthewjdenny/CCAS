test_that("That LDA contribution works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    num_topics = 100
    alpha = 5
    alpha_m = rep(alpha/num_topics,num_topics)

    result <- test_internal_functions(
        Test_RDirichlet = TRUE,
        alpha_m = alpha_m)

    plot(result)

    # no errors, will write an analytical test

})
