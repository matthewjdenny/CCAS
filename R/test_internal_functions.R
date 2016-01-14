#' @title Test C++ functions
#' @description Allows the user to interface with the individual C++ functions.
#'
#' @param Test_Log_Space_Multinomial_Sampler Defualts to FALSE. If TRUE, then
#' optional arguments 'distribution' and 'seed' must be provided.
#' @param envir Should not be changed by the user, captures the current environment
#' to facilitate testing.
#' @param ... optional arguments necessary to run each of the internal functions.
#' @return Whatever is returned by the internal function being tested
#' @export
test_internal_functions <- function(Test_Log_Space_Multinomial_Sampler = FALSE,
                                    envir = environment(),
                                    ...){

    # set all optional variables potentially used in the function to NULL inside
    # a conditional that is never satisfied so that we avoid R CMD check error.
    dont_preassign <- TRUE
    if(!dont_preassign){
        return_object <- NULL
        seed <- NULL
        distribution <- NULL
    }

    object <- as.list(substitute(list(...)))[-1L]
    if (length(object) > 0) {
        for (i in 1:length(object)){
            # try both direct assignment and get()
            if(typeof(object[[i]]) == "symbol"){
                # have to do this double get trick to make it work in all contexts.
                temp <- dynGet(as.character(object[[i]]),
                                       ifnotfound = get(as.character(object[[i]]),
                                                        envir = envir))
                assign(names(object)[i],temp)
            }else{
                assign(names(object)[i],object[[i]])
            }
        }
    }

    # test the lsms function
    if(Test_Log_Space_Multinomial_Sampler){
        if(is.null(distribution)){
            stop("you must provide an optional variable 'distribution', which is an unnormalized vetor of log probabilities.")
        }
        if(is.null(seed)){
            stop("you must provide an optional variable 'seed', which is an integer.")
        }
        return_object <- lsms(distribution,
                              seed)
    }

    # return whatever needs to be returned
    return(return_object)
}
