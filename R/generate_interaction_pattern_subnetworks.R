generate_interaction_pattern_subnetworks <- function(CCAS_Object) {

    # assign temporary variables we will operate on.
    num_actors <- CCAS_Object@ComNet_Object@num_actors
    num_IPs <- CCAS_Object@interaction_patterns
    token_topic_assignments <- CCAS_Object@topic_model_results$token_topic_assignments
    topic_IPs <- CCAS_Object@topic_model_results$topic_interaction_patterns
    # get the last iteration
    topic_IPs <- as.numeric(topic_IPs[nrow(topic_IPs),])
    doc_edge_matrix <- CCAS_Object@ComNet_Object@document_edge_matrix
    senders <-  CCAS_Object@ComNet_Object@document_authors

    # generate list of adjacency matrices we will fill in:
    interaction_pattern_networks <- vector(mode = "list",
                                           length = num_IPs)
    # count the total number of emails where some edge weight was contributed
    # to the interaction pattern sociomatrix
    emails_represented <- rep(0,num_IPs)
    # populate the list object with blank networks
    for (i in 1:num_IPs) {
        interaction_pattern_networks[[i]] <- matrix(0,
                                                    nrow = num_actors,
                                                    ncol = num_actors)
    }
    # make sure it has interpretable names
    names(interaction_pattern_networks) <- paste("interaction_pattern_",
                                                 1:num_IPs, sep = "")

    # fill in the adjacency matrices
    for (i in 1:CCAS_Object@ComNet_Object@num_documents) {
        # get the current token topic assignments
        current_tta <- as.numeric(token_topic_assignments[[i]])
        # create a blank vector of proportions of edge weight to be assigned to
        # each cluster
        cluster_proportions <- rep(0,num_IPs)

        # make sure there is atleast one token in the document
        if (!is.null(current_tta)) {
            # loop over tokens in the current document
            for (j in 1:length(current_tta)) {
                # figure out the current token topic cluster assignment
                IP <- topic_IPs[current_tta[j]]
                cluster_proportions[IP] <- cluster_proportions[IP] + 1
            }
            cluster_proportions <- cluster_proportions / length(current_tta)
        }

        # increment the interaction pattern specific networks
        for (k in 1:num_IPs) {
            interaction_pattern_networks[[k]][senders,] <-
                interaction_pattern_networks[[k]][senders,] +
                cluster_proportions[k] * doc_edge_matrix[i,]
            # increment the count of emails represented
            if (cluster_proportions[k] > 0) {
                emails_represented[k] <- emails_represented[k] + 1
            }
        }
    }

    # return the list of networks
    return(list(interaction_pattern_networks = interaction_pattern_networks,
                emails_represented = emails_represented))
}
