# CCAS [![Travis-CI Build Status](https://travis-ci.org/matthewjdenny/CCAS.svg?branch=master)](https://travis-ci.org/matthewjdenny/CCAS)
A statistical model for communication network data (R package)

## Overview
This package provides a number of functions to read in and model communication network (email) data. See `?prepare_data`, and `?ccas` to get started, as well as the example at the bottom of this page.

## Installation

### Requirements for using C++ code with R

Note that if you are using a Mac, you will need to start by making sure you have Xcode + developer tools installed or you will not be able to compile the C++ code that is used in the samplers for this package. You will need to go here: <https://developer.apple.com/xcode/downloads/> and then select the link to the "additional tools" page which will prompt you to enter your apple ID. This will let you download the developer tools. This requirement is not unique to this package, but is necessary to compile all packages that use Rcpp.  
  
If you are using a Windows machine, you will need to make sure you have the latest release of R (3.2.0+) and will also need to install `Rtools` (v33 or higher, available here <http://cran.r-project.org/bin/windows/Rtools/>)  before you can use any packages with C++ code in them. It is also highly advised that you use [RStudio](http://www.rstudio.com/) to download and install the package as it seems to play nicer with Rcpp under Windows. You may also want to visit [this blog post](https://cdrv.wordpress.com/2013/01/12/getting-compilers-to-work-with-rcpp-rcpparmadillo/) which has more information on making C++ work with R under Windows. 
  
If you are using a Linux distro, make sure you have a C++ complier installed, but in general, you should not run into as many issues. 

More generally, I suggest you check out this [tutorial on using C++ with R](http://www.mjdenny.com/Rcpp_Intro.html). It goes over some of the code used in this package and also covers a number of potential problems you might run into when trying to compile C++ code on your computer, so it is a good reference. 

### Installing The Package
  
To install this package from Github, you will need to Hadley Wickham's devtools package installed.

    install.packages("devtools")
    
Now we can install from Github using the following line:

    devtools::install_github("matthewjdenny/CCAS")

I have  had success installing with R 3.2.3+ installed but please email me if you hit any issues.

### Using The Package

The package provides two main functions `prepare_data()` and `ccas()`, which allow the user to read in raw communication network (email) data, and run the statistical model on that data (also generating output). For illustration, we will use some example county government department manager email data, consisting of 121 emails sent among 20 department managers in a North Carolina county during a three month period in 2013.

We begin by loading in the data (which are included with the package), and creating a `ComNet` object which will catalogue and store the data in the correct format so that we can run our model on the data. This object is also a much cleaner way to store data than as a list, or in multiple different data objects. You can learn more about each raw data object by using the following queries: `?author_attributes`, `?document_edge_matrix`, `?document_word_matrix`, and `?vocabulary`. 

	 # load in example county government email data.
	 data(author_attributes)
	 data(document_edge_matrix)
	 data(document_word_matrix)
	 data(vocabulary)

	 # the first column of the doc-edge matrix is the author index. Take it out and
	 # then remove it from the doc-edge matrix.
	 document_authors <- document_edge_matrix[,1]
	 document_edge_matrix <- document_edge_matrix[,-1]
	 
	 # create the ComNet object
	 ComNet <- prepare_data(document_authors = document_authors,
	                        document_edge_matrix = document_edge_matrix,
	                        document_term_matrix = document_word_matrix,
	                        covariate_data = author_attributes,
	                        vocabulary = vocabulary)

Once the data have been read in, we can proceed to run the CCAS model on our ComNet object. This model is specified in a very similar fashion to a model specification using the `latentnet` package. The key component of model specification using the `ccas()` function is a formula object of the form `ComNet ~ euclidean(d = 2)` where `d` is the number of dimensions in the latent space that the user would like to include. The formula may also include optional terms `sender("covariate_name")`, `receiver("covariate_name")`, `nodemix("covariate_name", base = value)` and `netcov("network_covariate")`, which are defined analogously to the arguments in the `latentnet` package, and allow for essentially limitless flexibility in covariate effects (thanks to the `netcov` argument). In this example, we are going to use a two dimensional latent space and gender mixing-matrices, holding out Male-Male interaction as a base case. This is a relatively standard specification, although users will want to increase the `iterations` argument to 2-5,000 depending on their dataset to ensure model convergence (although this number will vary). Fortunately, the `ccas()` function provides automatic convergence diagnostics. The user will also want to increase the `final_metropolis_hastings_iterations` and `final_metropolis_hastings_burnin` to somewhere on the order of 1,000,000 to ensure proper mixing (while making sure to decrease the `thin` argument to somewhere around 1/200 - 1/500). 

	# specify a formula that we will use for testing.
	formula <- ComNet_data ~ euclidean(d = 2) +
	           nodemix("Gender", base = "Male")
	
	# run the model 
	Model_Output <- ccas(formula,
	                    interaction_patterns = 4,
	                    topics = 40,
	                    alpha = 1,
	                    beta = 0.01,
	                    iterations = 50,
	                    metropolis_hastings_iterations = 500,
	                    final_metropolis_hastings_iterations = 10000,
	                    final_metropolis_hastings_burnin = 5000,
	                    thin = 1/10,
	                    target_accept_rate = 0.25,
	                    tolerance = 0.05,
	                    adaptive_metropolis_update_size = 0.05,
	                    LSM_proposal_variance = .5,
	                    LSM_prior_variance = 1,
	                    LSM_prior_mean = 0,
	                    slice_sample_alpha_m = TRUE,
	                    slice_sample_step_size = 1,
	                    generate_plots = TRUE,
	                    output_directory = NULL,
	                    output_name_stem = NULL)
						
Once the model has completed running, it will automatically generate diagnostic plots and output for interpreting estimation results if the `generate_plots` argument is set to `TRUE`. Alternatively, the user may generate output plots themselves by accessing the following functions: `plot_topic_model_log_likelihood()`, `plot_interaction_pattern_log_likelihood()`, `plot_interaction_pattern_network()`, `plot_parameter_estimates()`, and `top_words()` functions, which are each separately documented. Addition summary output can be found in the `@model_output` of the object returned by the `ccas()` function, which is of class `CCAS`.
