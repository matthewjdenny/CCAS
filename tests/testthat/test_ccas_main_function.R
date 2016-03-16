test_that("main function works", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)
    data(ComNet_data)

    net_covariate <- matrix(runif(100),10,10)

    # specify a formula that we will use for testing.
    formula <- ComNet_data ~ euclidean(d = 2) +
                             sender("age") +
                             nodemix("gender", base = "male") +
                             netcov("net_covariate")

})
