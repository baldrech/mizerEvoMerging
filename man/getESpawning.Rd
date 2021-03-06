% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_methods.R
\name{getESpawning}
\alias{getESpawning}
\title{Alias for getERepro}
\usage{
getESpawning(object, n, n_pp, e = getEReproAndGrowth(object, n = n, n_pp
  = n_pp))
}
\arguments{
\item{object}{A \code{MizerParams} object.}

\item{n}{A matrix of species abundance (species x size).}

\item{n_pp}{A vector of the plankton abundance by size.}

\item{e}{The energy available for reproduction and growth (optional). A
matrix of size no. species x no. size bins. If not supplied, is calculated
internally using the \code{getEReproAndGrowth()} method.}
}
\value{
A two dimensional array (prey species x prey size)
}
\description{
An alias provided for backward compatibility with mizer version <= 1.0
}
\examples{
\dontrun{
data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
# Project with constant fishing effort for all gears for 20 time steps
sim <- project(params, t_max = 20, effort = 0.5)
# Get the energy at a particular time step
getERepro(params,sim@n[21,,],sim@n_pp[21,])
}
}
\seealso{
\code{\link{project}} and \code{\link{getEReproAndGrowth}}.
}
