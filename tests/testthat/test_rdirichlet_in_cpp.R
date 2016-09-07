test_that("That LDA contribution works", {
    skip_on_cran()

    # create an example distribution
    num_topics = 100
    alpha = 5
    alpha_m = rep(alpha/num_topics,num_topics)

    set.seed(12345)
    result <- test_internal_functions(
        Test_RDirichlet = TRUE,
        alpha_m = alpha_m)

    set.seed(12345)
    y <- rgamma(length(alpha_m), alpha_m)
    x <- y / sum(y)

    expect_equal(x, as.numeric(result))
})
