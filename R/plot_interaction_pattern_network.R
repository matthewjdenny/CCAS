#' @title Plot interaction-pattern sepcific sub-networks
#' @description Plots a sub network associated with a particular interaction
#' pattern with nodes placed according to their posterior mean latent positions.
#'
#' @param CCAS_Object The object returned by the ccas() main estimation
#' function.
#' @param interaction_pattern_index The index of the cluster specific
#' communication network we wish to plot.
#' @param plot_color_category The name of a categorical covariate the user
#' wishes to color nodes by. For example, nodes might be colored by gender. This
#' function only currently supports plotting by a categroical variable
#' containing five or less unique values. Defaults to NULL, in which case all
#' nodes are plotted in the same color.
#' @param generate_plot Logical indicating whether plot should be printed.
#' Defaults to TRUE, but can be set to FALSE if the user only wishes to access
#' the mean latent positions of actors.
#' @return An data.frame containing the mean latent coordinates of each actor.
#' @export
plot_interaction_pattern_network <- function(CCAS_Object,
                                             interaction_pattern_index,
                                             plot_color_category = NULL,
                                             generate_plot = TRUE){
    # shorthand name to be used in the function
    cluster <- interaction_pattern_index

    # define light colors
    UMASS_BLUE <- rgb(51,51,153,20,maxColorValue = 255)
    UMASS_RED <- rgb(153,0,51,20,maxColorValue = 255)
    UMASS_GREEN <- rgb(0,102,102,20,maxColorValue = 255)
    UMASS_YELLOW <- rgb(255,255,102,20,maxColorValue = 255)
    UMASS_ORANGE <- rgb(255,204,51,20,maxColorValue = 255)
    light_colors <- c(UMASS_BLUE, UMASS_RED, UMASS_GREEN,
                     UMASS_YELLOW, UMASS_ORANGE)
    # define dark colors
    UMASS_BLUE <- rgb(51,51,153,255,maxColorValue = 255)
    UMASS_RED <- rgb(153,0,51,255,maxColorValue = 255)
    UMASS_GREEN <- rgb(0,102,102,255,maxColorValue = 255)
    UMASS_YELLOW <- rgb(255,255,102,255,maxColorValue = 255)
    UMASS_ORANGE <- rgb(255,204,51,255,maxColorValue = 255)
    dark_colors <- c(UMASS_BLUE, UMASS_RED, UMASS_GREEN,
                     UMASS_YELLOW, UMASS_ORANGE)

    # the possible node shapes
    possible_shapes <- c(16,17,15,18,8)

    # get the number of samples after burnin and thinning
    samples <- nrow(CCAS_Object@MCMC_output$intercepts)

    # get the number of latent dimensions, actors, etc.
    lat_dim <- CCAS_Object@latent_space_dimensions
    num_actors <- CCAS_Object@ComNet_Object@num_actors

    # extract the appropriate network from the CCAS Object
    network <- CCAS_Object@model_output$interaction_pattern_subnetworks[[cluster]]

    # extract a count of the number of emails represented by this interaction
    # pattern
    emails_represented <- CCAS_Object@model_output$emails_represented[cluster]

    # extract latent positions
    latent_positions <- CCAS_Object@MCMC_output$latent_positions[cluster,,]

    #generate mean coordinates
    base <- latent_positions[,1:lat_dim]
    mean_coordinates <- matrix(0, nrow = num_actors, ncol = lat_dim)
    for (i in 1:samples) {
        start <- lat_dim * (i - 1) + 1
        end <- lat_dim * i
        temp <- latent_positions[,start:end]
        rotated <- vegan::procrustes(base, temp, scale = F)$Yrot
        mean_coordinates <- mean_coordinates + rotated
    }

    #take the average
    mean_coordinates <- mean_coordinates / samples
    mean_coordinates <- as.data.frame(mean_coordinates)

    # either color all of the nodes the same color, or color them according to
    # the node level covariate categories
    node_colors <- rep(dark_colors[1], num_actors)
    point_colors <- rep(light_colors[1], num_actors)
    node_shapes <- rep(possible_shapes[1], num_actors)
    # set the color indicator to FALSE (Default)
    colored_by_category <- FALSE

    # now determine whether a covariate was provided and how many and what
    # levels it has.
    node_covariates <- CCAS_Object@ComNet_Object@covariate_data
    num_categories <- cats <- NULL

    # if a plot color category was provided
    if (!is.null(plot_color_category)) {
        # if covariate data was provided
        if (!is.null(node_covariates)) {
            # get the column index
            col <- which(colnames(node_covariates) == plot_color_category)
            # if there was a unique entry
            if (length(col) == 1) {
                # get the unique number of different categories
                col_vec <- as.vector(node_covariates[,col])
                cats <- unique(col_vec)

                # only color by categories if there are at least two categories
                # but less than six as we only have 5 colors to work with
                # currently
                if (length(cats) <= length(dark_colors) & length(cats) > 1) {
                    # then for the second to the last category, replace the
                    # colors with the 2 through nth colors in our color vectors
                    for (l in 2:length(cats)) {
                        cur_cat <- which(col_vec == cats[l])
                        node_colors[cur_cat] <- dark_colors[l]
                        point_colors[cur_cat] <- light_colors[l]
                        node_shapes[cur_cat] <- possible_shapes[l]
                    }
                    # set our indicator to TRUE
                    colored_by_category <- TRUE
                    num_categories <- length(cats)
                } else if (length(cats) > length(dark_colors)) {
                    cat("Only five categories are currently supported you provided:",
                         length(cats),"\n")
                } else if (length(cats) < 2) {
                    cat("You provided less than two categories so all nodes will be the same color.\n")
                }
            } else {
                cat("There was not a unique column in the node covariate data with name",
                    plot_color_category,
                    ". Either it does not exist, or there is more than one.\n")
            }
        }
    } # end of conditional to generate node colors if a category was provided


    #now generate edge colors which get darker as there is more weight
    linecolor <- colorRampPalette(c("grey90", "black"))
    cols <- linecolor(21)
    min <- min(network)
    max <- max(network)
    difference <- max - min
    vec <- seq(min,max,difference/20)

    # generate a color for each edge
    edge_colors <- rep("white",length(network))
    for (j in 1:length(network)) {
        already <- F
        for (k in 1:length(vec)) {
            if (!already) {
                if (network[j] == 0) {
                    edge_colors[j] <- cols[k]
                    already <- T
                } else {
                    if (vec[k] >= network[j]){
                        if ((k + 1) > length(vec)) {
                            index <- length(vec)
                        } else {
                            index <- k + 1
                        }
                        edge_colors[j] <- cols[index]
                        already <- T
                    }
                }
            }
        }
    }

    # if we are actaully generating a plot, then do so
    if (generate_plot) {
        # make sure the background for the plot is white
        par(bg = "white")

        # start the plot!
        plot(mean_coordinates,
             main = paste("Interaction Pattern:",cluster,"Fraction of Edges:",
                          round(sum(network)/CCAS_Object@ComNet_Object@num_edges,3),
                          "--- Number of Emails Represented:",emails_represented,
                          " of ",CCAS_Object@ComNet_Object@num_documents, "\n",
                          "Darker edges indicate more communication."),
             pch = 20,
             col = "white",
             axes = F,
             xlab = "",
             ylab = "")
        box(which = "plot")

        # add in points
        base <- latent_positions[,1:lat_dim]
        for (i in 1:samples) {
            start <- lat_dim * (i - 1) + 1
            end <- lat_dim * i
            temp <- latent_positions[,start:end]
            rotated <- vegan::procrustes(base, temp, scale = F)$Yrot
            # fill in points
            for (j in 1:num_actors) {
                points(x = rotated[j,1],
                       y = rotated[j,2],
                       col = point_colors[j],
                       pch = node_shapes,
                       cex = .4)
            }
        }

        #add in lines between actors in order from dimmest to brightest
        total <- length(which(network > 0))

        #iff there are any edges assigned to this cluster
        if(total > 0){
            ordering <- order(network, decreasing = F)
            correct_ordering <- rep(0,length(ordering))
            for (i in 1:length(ordering)) {
                correct_ordering[ordering[i]] <- i
            }
            add_matrix <- matrix(correct_ordering,
                                 ncol = num_actors,
                                 nrow = num_actors)
            colormat <- matrix(edge_colors,
                               ncol = num_actors,
                               nrow = num_actors)

            # fill in the edges in reverse order of darkness (weight) so that the
            # darkest edges are on top
            for (l in 1:length(ordering)) {
                for (i in 1:num_actors) {
                    for (j in 1:num_actors) {
                        if (add_matrix[i,j] == l & network[i,j] > 0) {
                            lines(c(mean_coordinates[j,1],
                                    mean_coordinates[i,1]) ,
                                  c(mean_coordinates[j,2],
                                    mean_coordinates[i,2]),
                                  col = colormat[i,j],
                                  lwd = 1.5)
                        }
                    }
                }
            }
        }

        # now add in points for the nodes themselves
        points(mean_coordinates,
               col = node_colors,
               pch = node_shapes,
               cex = 1.5)

        # if we used
        if (colored_by_category) {
            legend("topright",
                   inset = 0.05,
                   cats,
                   fill = dark_colors[1:num_categories],
                   horiz = FALSE,
                   bg = "white")
        }
    }

    return(mean_coordinates)
}
