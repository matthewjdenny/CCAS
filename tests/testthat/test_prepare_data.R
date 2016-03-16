test_that("prepare data function works, and that it returns a ComNet object", {
    skip_on_cran()

    # create an example distribution
    set.seed(12345)

    num_documents = 20
    num_terms = 100
    num_actors = 10

    document_edge_matrix <- matrix(round(runif(num_actors * num_documents)),
                                   nrow = num_documents,
                                   ncol = num_actors)
    document_term_matrix <- matrix(round(rpois(num_documents*num_terms, lambda = 0.3)),
                                   nrow = num_documents,
                                   ncol = num_terms)

    author_indexes <- floor(runif(n = num_documents, min = 0, max = num_actors - 0.00001))

    covariate_data <- data.frame(id = LETTERS[1:10],
                                 gender = c(rep("male",5),rep("female",5)),
                                 age = rpois(10,10) + 10)

    # create a comnet object
    ComNet_data <- prepare_data(document_authors = author_indexes,
                             document_edge_matrix = document_edge_matrix,
                             document_term_matrix = document_term_matrix,
                             covariate_data = covariate_data,
                             vocabulary = paste0("word_",1:num_terms))

    # devtools::use_data(ComNet_data, overwrite = TRUE)

    expect_equal(ComNet_data@num_tokens, sum(document_term_matrix))
    expect_equal(ComNet_data@num_edges, sum(document_edge_matrix))
})
