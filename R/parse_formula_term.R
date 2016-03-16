# Parse a formula entry
parse_formula_term <- function(term,
                               possible_structural_terms,
                               possible_covariate_terms,
                               possible_network_terms){
    # split up the formula term
    parsed <- stringr::str_split(term,"[\\(\\)]+")[[1]]
    # get rid of trailing spaces
    rm_spaces <- which(parsed == "")
    if(length(rm_spaces) > 0){
        parsed <-  parsed[-rm_spaces]
    }

    # generate return list object
    return_list <- list(term = parsed[1],
                        d = 1,
                        covariate = NA,
                        base = NA,
                        network = NA,
                        threshold = 0,
                        levels = NA,
                        same = NA,
                        parens_no_arg = NA,
                        network_matrix_object = NA,
                        num_levels = NA,
                        base_index = NA)
    possible_fields <- c("term","d","covariate", "base", "network",
                         "threshold", "levels", "same", "")
    # if there is an argument to the term -- this will be a lazy implementation
    # where if you do not get the name right, it will simply not be set and a
    # warning will be thrown.
    if(length(parsed) > 1){

        # take the arguments and further deparse them and turn them into a list
        # object, removing any quotes and further splitting each argument on an
        # equals sign so we can get the key value pair.
        args <- stringr::str_split(parsed[2],",")[[1]]
        args <- as.list(args)
        for(i in 1:length(args)){
            args[[i]] <- stringr::str_split(args[[i]],"=")[[1]]
            for(j in 1:length(args[[i]])){
                args[[i]][j] <- stringr::str_replace_all(args[[i]][j],"\"","")[[1]]
                args[[i]][j] <- stringr::str_replace_all(args[[i]][j],"\'","")[[1]]
            }
        }

        # if we are only supplied one term, without a name, we will default to
        # setting the exponential weight if it is a structural term and the
        # covariate name if it is a covariate term.

        ###### for structural terms ######
        if(return_list$term %in% possible_structural_terms){
            if(length(args[[1]]) == 1){
                # if an argument is supplied without an expicit argument name
                # assignment, then treat it as a weight and check if it is numeric
                # make sure the argument is numeric
                args[[1]] <- as.numeric(args[[1]])
                if(is.numeric(args[[1]])){
                    return_list$d <- as.numeric(args[[1]])
                }else{
                    stop(paste("You must supply a numeric d for latent space specification term:", return_list$term))
                }
            }else{
                which_arg <- which(possible_fields == args[[1]][1])
                if(length(which_arg) == 0){
                    stop(paste("You supplied an argument:",args[[1]][1],"which is not recognized."))
                }else{
                    return_list[[which_arg]] <- as.numeric(args[[1]][2])
                }
            }
            ###### for covariates ######
        }else if(return_list$term %in% possible_covariate_terms){
            if(length(args[[1]]) == 1){
                # if an argument is supplied without an expicit argument name
                # assignment, then treat it as a weight and check if it is numeric
                if(is.character(args[[1]])){
                    return_list$covariate <- args[[1]]
                }else{
                    stop(paste("You must supply a valid name for each node level covariate. You specified:", return_list$term))
                }
            }else{
                which_arg <- which(possible_fields == args[[1]][1])
                if(length(which_arg) == 0){
                    stop(paste("You supplied an argument:",args[[1]][2],"which is not recognized."))
                }else{
                    return_list[[which_arg]] <- args[[1]][2]
                }
            }
            ###### for network covariates ######
        }else if(return_list$term %in% possible_network_terms){
            if(length(args[[1]]) == 1){
                # if an argument is supplied without an expicit argument name
                # assignment, then treat it as a weight and check if it is numeric
                if(is.character(args[[1]])){
                    return_list$network <- args[[1]]
                }else{
                    stop(paste("You must supply a valid name for each network covariate. You specified:", return_list$term))
                }
            }else{
                which_arg <- which(possible_fields == args[[1]][1])
                if(length(which_arg) == 0){
                    stop(paste("You supplied an argument:",args[[1]][2],"which is not recognized."))
                }else{
                    return_list[[which_arg]] <- args[[1]][2]
                }
            }
            # check to make sure a valid network matrix is specified
            if(!is.na(return_list$network) & !is.null(return_list$network)){
                temp_net <- dynGet(as.character(as.character(return_list$network)),
                                   ifnotfound = get(as.character(as.character(return_list$network))))
                return_list$network_matrix_object <- temp_net
                if(class(return_list$network_matrix_object) != "matrix"){
                    stop(paste("You must supply network covariates as matrix objects."))
                }
            }else{
                stop(paste("The network covariate matrix:",return_list$network,"does not appear to exist, please check that it is loaded in your current R session."))
            }
        }else{
            # throw an error because the term was not valid
            stop(paste("You specified term:",return_list$term,"which is not recognized. Please respecify."))
        }

        # if we are given two arguments.
        if(length(args) > 1){
            # loop through the rest of the arguments
            for(i in 2:length(args)){
                ###### for structural terms ######
                if(return_list$term %in% possible_structural_terms){
                    which_arg <- which(possible_fields == args[[i]][1])
                    if(length(which_arg) == 0){
                        stop(paste("You supplied an argument:",args[[i]][2],"which is not recognized."))
                    }else{
                        return_list[[which_arg]] <- args[[i]][2]
                    }
                    ###### for covariates ######
                }else{
                    which_arg <- which(possible_fields == args[[i]][1])
                    if(length(which_arg) == 0){
                        stop(paste("You supplied an argument:",args[[i]][2],"which is not recognized."))
                    }else{
                        return_list[[which_arg]] <- args[[i]][2]
                    }
                }
            }
        }
    }

    # return everything
    return(return_list)
}
