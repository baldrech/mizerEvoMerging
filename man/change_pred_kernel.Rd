% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper_functions.R
\name{change_pred_kernel}
\alias{change_pred_kernel}
\title{Changes the predation kernel to allow size-dependent PPMR}
\usage{
change_pred_kernel(params, pred_kernel)
}
\arguments{
\item{params}{A \linkS4class{MizerParams} object}

\item{pred_kernel}{An array (species x predator size x prey size) 
that holds the predation coefficient of each predator at size on 
each prey size. The dimensions are thus no_sp, no_w, no_w_full.}
}
\value{
A \linkS4class{MizerParams} object
}
\description{
One way to set up a model with size-dependent predator/prey mass ratio
(PPMR) is to first use one of the standard set up functions to create
a MizerParams object and then to call this function with that object 
and an array holding the predation kernel.
}
\examples{
\dontrun{
## Set up a MizerParams object
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)

## Create a predation kernel
beta <- params@species_params$beta
sigma <- params@species_params$sigma
w <- params@w
w_full <- params@w_full
pk = array(beta, dim = c(length(beta),length(w), length(w_full)))
pk <- exp(-0.5 * sweep(log(sweep(sweep(pk, 3, w_full, "*") ^ -1, 2, w, "*")),
                       1, sigma, "/") ^ 2)
pk <- sweep(pk, c(2, 3), combn(w_full, 1, function(x, w) x < w, w = w), "*")

## Create a new MizerParamsVariablePPMR object to use the new predation kernel
params_old <- change_pred_kernel(params, pred_kernel = pk)
}
}
