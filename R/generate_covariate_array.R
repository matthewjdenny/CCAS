generate_covariate_array <- function(formula,
                                     possible_structural_terms,
                                     possible_covariate_terms,
                                     possible_network_terms,
                                     ComNet_Object){

    node_covariates_list <- parse_specification(
        formula = formula,
        possible_structural_terms = possible_structural_terms,
        possible_covariate_terms = possible_covariate_terms,
        possible_network_terms = possible_network_terms,
        terms_to_parse = "covariate")
    network_covariates_list <- parse_specification(
        formula = formula,
        possible_structural_terms = possible_structural_terms,
        possible_covariate_terms = possible_covariate_terms,
        possible_network_terms = possible_network_terms,
        terms_to_parse = "network")

    # determine whether covariates were provided and the total number of
    # covariates
    if(length(node_covariates_list) < 1){
        node_covariates_provided <- FALSE
    }else{
        node_covariates_provided <- TRUE
    }
    if(length(network_covariates_list) < 1){
        network_covariates_provided <- FALSE
    }else{
        network_covariates_provided <- TRUE
    }

    num_covariates <- 0

    if (!node_covariates_provided) {
        cat("No node level covariates provided.\n")
    } else {
        covariate_data <- ComNet_Object@covariate_data
        node_names <- rownames(covariate_data)
        for (i in 1:length(node_covariates_list)) {
            type <- node_covariates_list[[i]]$term
            if (type == "sender" | type == "receiver" | type == "absdiff" |
                type == "nodecov" | type == "intercept" | type == "nodematch") {
                num_covariates <- num_covariates + 1
            } else if (type == "nodemix") {
                # need to get the number of levels
                covar <- node_covariates_list[[i]]$covariate
                index <- which(colnames(covariate_data) == covar)
                node_covariates_list[[i]]$levels <- unique(covariate_data[,index])
                num_levels <- length(node_covariates_list[[i]]$levels)
                node_covariates_list[[i]]$num_levels <- num_levels
                if(node_covariates_list[[i]]$base == "NULL"){
                    # if the user specified base = NULL then we use all levels for matching
                    num_covariates <- num_covariates + num_levels*num_levels
                    node_covariates_list[[i]]$base_index <- 0
                }else{
                    num_covariates <- num_covariates + num_levels*num_levels - 1
                    base_index <- which(node_covariates_list[[i]]$levels == node_covariates_list[[i]]$base)
                    node_covariates_list[[i]]$base_index <- base_index
                }
            }else{
                stop(paste("You specified a node level covariate term:",node_covariates_list[[i]]$term, "Node level covariate effects must be one of: 'sender', 'receiver', or 'nodemix' , please respecify."))
            }
        }
        cat("You have specified", num_covariates,"node level covariate effects.\n")
    }

    #determine the type and number of user provided network covariates (will not be altered)
    if(!network_covariates_provided){
        cat("No network covariates provided.\n")
        num_additional_covars <- 0
    }else{
        num_covariates <- num_covariates + length(network_covariates_list)
        num_additional_covars <- length(network_covariates_list)
        cat("You have provided",num_additional_covars,"network covariates.\n")
    }

    generate_covariate_effect_matrix <- function(num_nodes,
                                                 node_names,
                                                 covariates,
                                                 covariate_column,
                                                 effect_type,
                                                 level = NA,
                                                 level2 = NA){
        return_matrix <- matrix(0,num_nodes,num_nodes)
        if(effect_type == "intercept"){
            return_matrix <- matrix(1,num_nodes,num_nodes)
            diag(return_matrix) <- 0
        }
        if(effect_type == "sender"){
            for(j in 1:num_nodes){
                for(k in 1:num_nodes){
                    if(j != k){
                        row <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
                        return_matrix[j,k] <- covariates[row,covariate_column]
                    }
                }
            }
        }
        if(effect_type == "receiver"){
            for(j in 1:num_nodes){
                for(k in 1:num_nodes){
                    if(j != k){
                        row <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
                        return_matrix[j,k] <- covariates[row,covariate_column]
                    }
                }
            }
        }
        if(effect_type == "absdiff"){
            for(j in 1:num_nodes){
                for(k in 1:num_nodes){
                    if(j != k){
                        row1 <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
                        row2 <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
                        return_matrix[j,k] <- abs(covariates[row1,covariate_column] - covariates[row2,covariate_column])
                    }
                }
            }
        }
        if(effect_type == "nodecov"){
            for(j in 1:num_nodes){
                for(k in 1:num_nodes){
                    if(j != k){
                        row1 <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
                        row2 <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
                        return_matrix[j,k] <- covariates[row1,covariate_column] + covariates[row2,covariate_column]
                    }
                }
            }
        }
        if(effect_type == "nodematch"){
            for(j in 1:num_nodes){
                for(k in 1:num_nodes){
                    if(j != k){
                        row1 <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
                        row2 <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
                        # handle both numeric and categorical values
                        if(is.numeric(covariates[row1,covariate_column]) & is.numeric(covariates[row2,covariate_column])){
                            check <- abs(covariates[row1,covariate_column] - covariates[row2,covariate_column])
                        }else{
                            check <- 1
                            if(covariates[row1,covariate_column] == covariates[row2,covariate_column]){
                                check <- 0
                            }
                        }
                        if(check == 0){
                            return_matrix[j,k] <- 1
                        }
                    }
                }
            }

        }
        if(effect_type == "nodemix"){
            for(j in 1:num_nodes){
                for(k in 1:num_nodes){
                    if(j != k){
                        col1 <- which(toupper(rownames(covariates)) == toupper(node_names)[k])
                        row1 <- which(toupper(rownames(covariates)) == toupper(node_names)[j])
                        colval <- covariates[col1,covariate_column]
                        rowval <- covariates[row1,covariate_column]
                        if(level == colval & level2 == rowval){
                            return_matrix[j,k] <- 1
                        }
                    }
                }
            }
        }
        return(return_matrix)
    }

    #########################################################
    # Construct transformed_covariates: An array of covariates which parameterize
    # the latent space.

    # Generate array which covaries will be transformed into
    num_nodes <- ComNet_Object@num_actors
    transformed_covariates <- array(0, dim = c(num_nodes, num_nodes,
                                               num_covariates))

    #omit this since we are adding in an intercept term
    #transformed_covariates[,,1] <- 1
    # Set a slice counter to keep track of where we should add the covariates
    # in the resulting covariate array.
    slice_counter <- 1

    # generate vector to store covariate names
    if (!node_covariates_provided & !network_covariates_provided) {
        # do nothing
    } else {
        slice_names <- rep("",length = num_covariates)

        #remove since we are specifying an intercept
        #slice_names[1] <- "network"
    }



    #1. generate sender and reciever effects
    if (node_covariates_provided) {

        # Determine if type_of_effect and covariates_to_use have the same length or
        # if no covariates_to_use is provided, that the number of columns in
        # covariate_data is the same length.
        node_covariates_list[[i]]$num_levels
        # Loop through covariates
        for (i in 1:length(node_covariates_list)) {
            if (node_covariates_list[[i]]$term == "intercept") {
                col_index <- 0
            } else {
                col_index <- which(tolower(colnames(covariate_data)) ==
                                       tolower(node_covariates_list[[i]]$covariate))
            }
            if (length(col_index) == 0) {
                stop(paste("There is no matching column name in covariate_data for:",
                           tolower(node_covariates_list[[i]]$covariate)))
            }

            if (node_covariates_list[[i]]$term == "nodemix") {
                levels_to_include <- node_covariates_list[[i]]$num_levels
                for (j in 1:levels_to_include) {
                    for (k in 1:levels_to_include) {
                        if (node_covariates_list[[i]]$levels[j] ==
                            node_covariates_list[[i]]$base &
                            node_covariates_list[[i]]$levels[k] ==
                            node_covariates_list[[i]]$base) {
                            #do nothing since we do not include a term for the base
                        } else {
                            add <- generate_covariate_effect_matrix(
                                num_nodes = num_nodes,
                                node_names = node_names,
                                covariates = covariate_data,
                                covariate_column = col_index,
                                effect_type = "nodemix",
                                level = node_covariates_list[[i]]$levels[j],
                                level2 = node_covariates_list[[i]]$levels[k])
                            #print(add)
                            transformed_covariates[,,slice_counter] <- add
                            slice_names[slice_counter] <- paste(
                                node_covariates_list[[i]]$covariate,
                                node_covariates_list[[i]]$term,
                                node_covariates_list[[i]]$levels[j],
                                node_covariates_list[[i]]$levels[k],
                                sep = "_")
                            slice_counter <- slice_counter + 1
                        }
                    }
                }
            }else{
                add <- generate_covariate_effect_matrix(
                    num_nodes = num_nodes,
                    node_names = node_names,
                    covariates = covariate_data,
                    covariate_column = col_index,
                    effect_type = node_covariates_list[[i]]$term)
                transformed_covariates[,,slice_counter] <- add
                if (node_covariates_list[[i]]$term == "intercept") {
                    slice_names[slice_counter] <- "intercept"
                } else {
                    slice_names[slice_counter] <- paste(
                        node_covariates_list[[i]]$covariate,
                        node_covariates_list[[i]]$term,
                        sep = "_")
                }
                slice_counter <- slice_counter + 1
            }
        } # End of loop over node level covariates
    } # End of condition that we have node level covariates



    #2. tack on any user supplied effects
    if (network_covariates_provided) {
        for (i in 1:num_additional_covars) {
            transformed_covariates[,,slice_counter] <- network_covariates_list[[i]]$network_matrix_object
            slice_names[slice_counter] <- paste(
                network_covariates_list[[i]]$network,
                network_covariates_list[[i]]$term,
                sep = "_")
            slice_counter <- slice_counter + 1
        }
    }

    dimnames(transformed_covariates) <- list(node_names,
                                             node_names,
                                             slice_names)

    return(list(covariate_array = transformed_covariates,
                number_of_covariates = length(slice_names)))
} # End of function definition.
