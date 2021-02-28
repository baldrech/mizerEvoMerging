# Plotting methods for the MizerSim class

# Copyright 2012 Finlay Scott and Julia Blanchard.
# Copyright 2018 Gustav Delius and Richard Southwell.
# Development has received funding from the European Commission's Horizon 2020 
# Research and Innovation Programme under Grant Agreement No. 634495 
# for the project MINOUW (http://minouw-project.eu/).
# Distributed under the GPL 3 or later 
# Maintainer: Gustav Delius, University of York, <gustav.delius@york.ac.uk>

# Hackiness to get past the 'no visible binding ... ' warning when running check
utils::globalVariables(c("time", "value", "Species", "w", "gear", "Age",
                         "x", "y", "Year", "Yield"))

#' Helper function to produce nice breaks on logarithmic axes
#'
#' This is needed when the logarithmic y-axis spans less than one order of
#' magnitude, in which case the ggplot2 default produces no ticks.
#' Thanks to Heather Turner at
#' https://stackoverflow.com/questions/14255533/pretty-ticks-for-log-normal-scale-using-ggplot2-dynamic-not-manual
#'
#' @param n Approximate number of ticks
#'
#' @return A function that can be used as the break argument in calls to
#'   scale_y_continuous() or scale_x_continuous()
log_breaks <- function(n = 6){
  n <- max(1, n)  # Because n=0 could lead to R crash
  function(x) {
    grDevices::axisTicks(log10(range(x, na.rm = TRUE)),
                         log = TRUE, nint = n)
  }
}


#' Display frames
#' 
#' @param f1 Data frame for left plot
#' @param f2 Data frame for right plot
#' @param params A MizerParams object
#' @param y_ticks The approximate number of ticks desired on the y axis
#' 
#' @return ggplot2 object
#' @export
display_frames <- function(f1, f2, params, y_ticks = 6) {
  var_names <- names(f1)
  if (!(length(var_names) == 3)) {
    stop("A frame needs to have three variables.")
  }
  if (!all(names(f2) == var_names)) {
    stop("Both frames need to have the same variable names.")
  }
  f <- rbind(cbind(f1, Simulation = 1), cbind(f2, Simulation = 2))
  p <- ggplot(f, aes_string(x = names(f)[1], y = names(f)[3],
                            colour = names(f)[2], linetype = names(f)[2])) +
    scale_y_log10(breaks = log_breaks(n = y_ticks), labels = prettyNum) +
    geom_line() +
    facet_wrap(~ Simulation) +
    scale_colour_manual(values = params@linecolour) +
    scale_linetype_manual(values = params@linetype)
  return(p)
}


#' Get data frame of spawning stock biomass of species through time, 
#' ready for ggplot2
#'
#' After running a projection, the spawning stock biomass of each species can be
#' plotted against time.
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all foreground species are plotted.
#' @param start_time The first time to be plotted. Default is the beginning
#'   of the time series.
#' @param end_time The last time to be plotted. Default is the end of the
#'   time series.
#' @param ylim A numeric vector of length two providing limits of for the
#'   y axis. Use NA to refer to the existing minimum or maximum. Any values
#'   below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total SSB from
#'   all species is plotted as well. Default is FALSE
#'   
#' @return A data frame that can be used in \code{\link{display_frames}}
#' @export
#' @seealso \code{\link{getSSB}}
getSSBFrame <- function(sim,
                        species = dimnames(sim@n)$sp[!is.na(sim@params@A)],
                        start_time = as.numeric(dimnames(sim@n)[[1]][1]),
                        end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
                        ylim = c(NA, NA), total = FALSE){
  b <- getSSB(sim)
  if (start_time >= end_time) {
    stop("start_time must be less than end_time")
  }
  # Select time range
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) &
           (as.numeric(dimnames(b)[[1]]) <= end_time), , drop = FALSE]
  b_total <- rowSums(b)
  # Include total
  if (total) {
    b <- cbind(b, Total = b_total)
    species <- c("Total", species)
  }
  bm <- reshape2::melt(b)
  # Implement ylim and a minimal cutoff
  min_value <- 1e-20
  bm <- bm[bm$value >= min_value &
             (is.na(ylim[1]) | bm$value >= ylim[1]) &
             (is.na(ylim[2]) | bm$value <= ylim[1]), ]
  names(bm) <- c("Year", "Species", "SSB")
  # Force Species column to be a factor (otherwise if numeric labels are
  # used they may be interpreted as integer and hence continuous)
  bm$Species <- as.factor(bm$Species)
  # Select species
  bm <- bm[bm$Species %in% species, ]
  return(bm)
}


#' Get data frame of biomass of species through time, ready for ggplot2
#'
#' After running a projection, the biomass of each species can be plotted
#' against time. The biomass is calculated within user defined size limits 
#' (min_w, max_w, min_l, max_l, see \code{\link{getBiomass}}). 
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param start_time The first time to be plotted. Default is the beginning
#'   of the time series.
#' @param end_time The last time to be plotted. Default is the end of the
#'   time series.
#' @param ylim A numeric vector of length two providing limits of for the
#'   y axis. Use NA to refer to the existing minimum or maximum. Any values
#'   below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total biomass from
#'   all species is plotted as well. Default is FALSE
#' @param ... Other arguments to pass to \code{getBiomass} method, for example
#'   \code{min_w} and \code{max_w}
#'   
#' @return A data frame that can be used in \code{\link{display_frames}}
#' @export
#' @seealso \code{\link{getBiomass}}
getBiomassFrame <- function(sim,
                            species = dimnames(sim@n)$sp[!is.na(sim@params@A)],
                            start_time = as.numeric(dimnames(sim@n)[[1]][1]),
                            end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
                            ylim = c(NA, NA), total = FALSE, ...){
  b <- getBiomass(sim, ...)
  if (start_time >= end_time) {
    stop("start_time must be less than end_time")
  }
  # Select time range
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) &
           (as.numeric(dimnames(b)[[1]]) <= end_time), , drop = FALSE]
  b_total <- rowSums(b)
  # Include total
  if (total) {
    b <- cbind(b, Total = b_total)
    species <- c("Total", species)
  }
  bm <- reshape2::melt(b)
  # Implement ylim and a minimal cutoff
  min_value <- 1e-20
  bm <- bm[bm$value >= min_value &
             (is.na(ylim[1]) | bm$value >= ylim[1]) &
             (is.na(ylim[2]) | bm$value <= ylim[1]), ]
  names(bm) <- c("Year", "Species", "Biomass")
  # Force Species column to be a factor (otherwise if numeric labels are
  # used they may be interpreted as integer and hence continuous)
  bm$Species <- as.factor(bm$Species)
  # Select species
  bm <- bm[bm$Species %in% species, ]
  
  return(bm)
}


#' Plot the biomass of species through time
#'
#' After running a projection, the biomass of each species can be plotted
#' against time. The biomass is calculated within user defined size limits 
#' (min_w, max_w, min_l, max_l, see \code{\link{getBiomass}}). 
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param start_time The first time to be plotted. Default is the beginning
#'   of the time series.
#' @param end_time The last time to be plotted. Default is the end of the
#'   time series.
#' @param y_ticks The approximate number of ticks desired on the y axis
#' @param ylim A numeric vector of length two providing limits of for the
#'   y axis. Use NA to refer to the existing minimum or maximum. Any values
#'   below 1e-20 are always cut off.
#' @param print_it Display the plot, or just return the ggplot2 object. Default
#'   value is TRUE
#' @param total A boolean value that determines whether the total biomass from
#'   all species is plotted as well. Default is FALSE
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param ... Other arguments to pass to \code{getBiomass} method, for example
#'   \code{min_w} and \code{max_w}
#'   
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getBiomass}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2)
#' plotBiomass(sim)
#' plotBiomass(sim, species = c("Cod", "Herring"), total = TRUE)
#' plotBiomass(sim, min_w = 10, max_w = 1000)
#' plotBiomass(sim, start_time = 10, end_time = 15)
#' plotBiomass(sim, y_ticks = 3)
#' }
plotBiomass <- function(sim,
                        species = dimnames(sim@n)$sp[!is.na(sim@params@A)],
                        start_time = as.numeric(dimnames(sim@n)[[1]][1]),
                        end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
                        y_ticks = 6, print_it = TRUE,
                        ylim = c(NA, NA),
                        total = FALSE, background = TRUE, returnData = FALSE, ...){
  b <- getBiomass(sim, ...)
  if (start_time >= end_time) {
    stop("start_time must be less than end_time")
  }
  # Select time range
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) &
           (as.numeric(dimnames(b)[[1]]) <= end_time), , drop = FALSE]
  b_total <- rowSums(b)
  # Include total
  if (total) {
    b <- cbind(b, Total = b_total)
    species <- c("Total", species)
  }
  names(dimnames(b)) <- c("time", "Species")
  bm <- reshape2::melt(b)
  # Force Species column to be a factor (otherwise if numeric labels are
  # used they may be interpreted as integer and hence continuous)
  bm$Species <- as.factor(bm$Species)
  # Implement ylim and a minimal cutoff
  min_value <- 1e-20
  bm <- bm[bm$value >= min_value &
             (is.na(ylim[1]) | bm$value >= ylim[1]) &
             (is.na(ylim[2]) | bm$value <= ylim[1]), ]
  # Select species
  spec_bm <- bm[bm$Species %in% species, ]
  x_label <- "Year"
  y_label <- "Biomass [g]"
  p <- ggplot(spec_bm, aes(x = time, y = value)) +
    scale_y_continuous(trans = "log10", breaks = log_breaks(n = y_ticks),
                       labels = prettyNum, name = y_label) +
    scale_x_continuous(name = x_label) +
    scale_colour_manual(values = sim@params@linecolour) +
    scale_linetype_manual(values = sim@params@linetype)
  
  if (background) {
    # Add background species in light grey
    back_sp <- dimnames(sim@n)$sp[is.na(sim@params@A)]
    back_bm <- bm[bm$Species %in% back_sp, ]
    p <- p + geom_line(aes(group = Species), data = back_bm,
                       colour = "lightgrey")
  }
  
  if ( (length(species) + total) > 12) {
    p <- p + geom_line(aes(group = Species))
  } else {
    p <- p +
      geom_line(aes(colour = Species, linetype = Species))
  }
  
  if(returnData) return(p) else if (print_it) return(p)
  
}


#' Plot the total yield of species through time
#'
#' After running a projection, the total yield of each species across all 
#' fishing gears can be plotted against time. The yield is obtained with
#' \code{\link{getYield}}.
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param sim2 An optional second object of class \linkS4class{MizerSim}. If
#'   this is provided its yields will be shown on the same plot in bolder lines.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species contained in \code{sim} are plotted.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE.
#' @param total A boolean value that determines whether the total yield from
#'   all species in the system is plotted as well. Default is FALSE.
#' @param log Boolean whether yield should be plotted on a logarithmic axis. 
#'   Defaults to true.
#' @param ... Other arguments to pass to \code{\link{getYield}} method.
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getYield}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2)
#' plotYield(sim)
#' plotYield(sim, species = c("Cod", "Herring"), total = TRUE)
#' 
#' # Comparing with yield from twice the effort
#' sim2 <- project(params, effort=2, t_max=20, t_save = 0.2)
#' plotYield(sim, sim2, species = c("Cod", "Herring"), log = FALSE)
#' }
plotYield <- function(sim, sim2,
                      species = dimnames(sim@n)$sp,
                      print_it = TRUE, total = FALSE, log = TRUE, ...){
  if (missing(sim2)) {
    y <- getYield(sim, ...)
    y_total <- rowSums(y)
    y <- y[, (as.character(dimnames(y)[[2]]) %in% species),
           drop = FALSE]
    if (total) {
      # Include total
      y <- cbind(y, "Total" = y_total)
    }
    ym <- reshape2::melt(y, varnames = c("Year", "Species"),
                         value.name = "Yield")
    ym$Species <- as.factor(ym$Species)
    ym <- subset(ym, ym$Yield > 0)
    if (nlevels(ym$Species) > 12) {
      p <- ggplot(ym) +
        geom_line(aes(x = Year, y = Yield, group = Species))
    } else {
      p <- ggplot(ym) +
        geom_line(aes(x = Year, y = Yield,
                      colour = Species, linetype = Species))
    }
    if (log) {
      p <- p + scale_y_continuous(trans = "log10", name = "Yield [g/year]",
                                  breaks = log_breaks(),
                                  labels = prettyNum)
    } else {
      p <- p + scale_y_continuous(name = "Yield [g/year]")
    }
    p <- p +
      scale_colour_manual(values = sim@params@linecolour) +
      scale_linetype_manual(values = sim@params@linetype)
    if (print_it) {
      print(p)
    }
    return(p)
  } else {
    if (!all(dimnames(sim@n)$time == dimnames(sim2@n)$time)) {
      stop("The two simulations do not have the same times")
    }
    y <- getYield(sim, ...)
    y2 <- getYield(sim2, ...)
    y_total <- rowSums(y)
    y <- y[, (as.character(dimnames(y)[[2]]) %in% species) & colSums(y) > 0,
           drop = FALSE]
    y2_total <- rowSums(y2)
    y2 <- y2[, (as.character(dimnames(y2)[[2]]) %in% species),
             drop = FALSE]
    if (total) {
      # Include total
      y <- cbind(y, Total = y_total)
      y2 <- cbind(y2, Total = y2_total)
    }
    ym <- reshape2::melt(y, varnames = c("Year", "Species"),
                         value.name = "Yield")
    ym2 <- reshape2::melt(y2, varnames = c("Year", "Species"),
                          value.name = "Yield")
    ym$Simulation <- 1
    ym2$Simulation <- 2
    ym <- rbind(ym, ym2)
    ym$Species <- as.factor(ym$Species)
    ym$Simulation <- as.factor(ym$Simulation)
    ym <- subset(ym, ym$Yield > 0)
    if (nlevels(ym$Species) > 12) {
      p <- ggplot(ym) +
        geom_line(aes(x = Year, y = Yield, group = Species))
    } else {
      p <- ggplot(ym) +
        geom_line(aes(x = Year, y = Yield, colour = Species,
                      linetype = Species))
    }
    if (log) {
      p <- p + scale_y_continuous(trans = "log10", name = "Yield [g/year]")
    } else {
      p <- p + scale_y_continuous(name = "Yield [g/year]")
    }
    p <- p + facet_wrap(~ Simulation)
    if (print_it) {
      print(p)
    }
    return(p)
  }
}


#' Plot the total yield of each species by gear through time
#'
#' After running a projection, the total yield of each species by fishing gear
#' can be plotted against time. 
#' 
#' This plot is pretty easy to do by hand. It just
#' gets the biomass using the \code{\link{getYieldGear}} method and plots using
#' the ggplot2 package. You can then fiddle about with colours and linetypes
#' etc. Just look at the source code for details.
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param print_it Display the plot, or just return the ggplot2 object. 
#'   Defaults to TRUE
#' @param total A boolean value that determines whether the total yield
#'   per gear over all species in the system is plotted as well. Default is FALSE
#' @param ... Other arguments to pass to \code{\link{getYieldGear}} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getYieldGear}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 0.2)
#' plotYieldGear(sim)
#' plotYieldGear(sim, species = c("Cod", "Herring"), total = TRUE)
#' }
plotYieldGear <- function(sim,
                          species = dimnames(sim@n)$sp,
                          print_it = TRUE, total = FALSE, ...){
  y <- getYieldGear(sim, ...)
  y_total <- rowSums(y, dims = 2)
  y <- y[, , dimnames(y)$sp %in% species, drop = FALSE]
  names(dimnames(y))[names(dimnames(y)) == "sp"] <- "Species"
  ym <- reshape2::melt(y)
  if (total) {
    yt <- reshape2::melt(y_total)
    yt$Species <- "Total"
    ym <- rbind(ym, yt)
  }
  ym <- subset(ym, ym$value > 0)
  if (length(species) > 12) {
    p <- ggplot(ym) + geom_line(aes(x = time, y = value, group = Species))
  } else {
    p <- ggplot(ym) +
      geom_line(aes(x = time, y = value, colour = Species, linetype = gear))
  }
  p <- p + scale_y_continuous(trans = "log10", name = "Yield [g]") +
    scale_x_continuous(name = "Year") +
    scale_colour_manual(values = sim@params@linecolour)
  if (print_it) {
    print(p)
  }
  return(p)
}


#' Plot the abundance spectra
#' 
#' Plots the number density multiplied by a power of the weight, with the power
#' specified by the \code{power} argument.
#'
#' When called with a \linkS4class{MizerSim} object, the abundance is averaged
#' over the specified time range (a single value for the time range can be used
#' to plot a single time step). When called with a \linkS4class{MizerParams}
#' object the initial abundance is plotted.
#' 
#' @param object An object of class \linkS4class{MizerSim} or \linkS4class{MizerParams}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step. Ignored when called with a \linkS4class{MizerParams}
#'   object.
#' @param min_w Minimum weight to be plotted (useful for truncating the
#'   plankton spectrum). Default value is a hundredth of the minimum size
#'   value of the community.
#' @param ylim A numeric vector of length two providing limits of for the
#'   y axis. Use NA to refer to the existing minimum or maximum. Any values
#'   below 1e-20 are always cut off.
#' @param power The abundance is plotted as the number density times the weight
#' raised to \code{power}. The default \code{power = 1} gives the biomass 
#' density, whereas \code{power = 2} gives the biomass density with respect
#' to logarithmic size bins.
#' @param biomass Obsolete. Only used if \code{power} argument is missing. Then
#'   \code{biomass = TRUE} is equivalent to \code{power=1} and 
#'   \code{biomass = FALSE} is equivalent to \code{power=0}
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param total A boolean value that determines whether the total over all
#'   species in the system is plotted as well. Default is FALSE
#' @param plankton A boolean value that determines whether plankton is included.
#'   Default is TRUE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param ... Other arguments (currently unused)
#'   
#' @return A ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotSpectra(sim)
#' plotSpectra(sim, min_w = 1e-6)
#' plotSpectra(sim, time_range = 10:20)
#' plotSpectra(sim, time_range = 10:20, power = 0)
#' plotSpectra(sim, species = c("Cod", "Herring"), power = 1)
#' }
plotSpectra <- function(object, species = NULL,
                        time_range,
                        min_w, ylim = c(NA, NA),
                        power = 1, biomass = TRUE, print_it = TRUE,
                        total = FALSE, plankton = TRUE,
                        background = TRUE, ...) {
  if (is(object, "MizerSim")) {
    if (missing(time_range)){
      time_range  <- max(as.numeric(dimnames(object@n)$time))
    }
    if (missing(min_w)){
      min_w <- min(object@params@w) / 100
    }
    # to deal with old-type biomass argument
    if (missing(power)) {
      power <- as.numeric(biomass)
    }
    time_elements <- get_time_elements(object,time_range)
    n <- apply(object@n[time_elements, , ,drop = FALSE], c(2, 3), mean)
    n_pp <- apply(object@n_pp[time_elements,,drop = FALSE], 2, mean)
    ##AAsp##
    #n_bb <- apply(object@n_bb[time_elements,,drop = FALSE], 2, mean)
    
    ps <- plot_spectra(object@params, n = n, n_pp = n_pp, 
                       species = species, min_w = min_w, ylim = ylim,
                       power = power, print_it = print_it,
                       total = total, plankton = plankton, 
                       background = background)
    return(ps)
  } else {
    if (missing(power)) {
      power <- as.numeric(biomass)
    }
    if (missing(min_w)) {
      min_w <- min(object@w) / 100
    }
    ps <- plot_spectra(object, n = object@initial_n,
                       n_pp = object@initial_n_pp,
                       species = species, min_w = min_w, ylim = ylim,
                       power = power, print_it = print_it,
                       total = total, plankton = plankton, 
                       background = background)
    return(ps)
  }
}


plot_spectra <- function(params, n, n_pp, 
                         species, min_w, ylim, power, print_it,
                         total, plankton, background) {
  if (total) {
    # Calculate total community abundance
    fish_idx <- (length(params@w_full) - length(params@w) + 1):
      length(params@w_full)
    total_n <- n_pp
    total_n[fish_idx] <- total_n[fish_idx] + colSums(n)
    total_n <- total_n * params@w_full^power
  }
  # Set species if missing to list of all non-background species
  if (is.null(species)) {
    species <- dimnames(params@initial_n)$sp[!is.na(params@A)]
  }
  # Deal with power argument
  if (power %in% c(0, 1, 2)) {
    y_label <- c("Number density [1/g]", "Biomass density",
                 "Biomass density [g]")[power + 1]
  } else {
    y_label <- paste0("Number density * w^", power)
  }
  n <- sweep(n, 2, params@w^power, "*")
  # Select only the desired species and background species
  spec_n <- n[as.character(dimnames(n)[[1]]) %in% species, , drop = FALSE]
  # Make data.frame for plot
  plot_dat <- data.frame(value = c(spec_n),
                         Species = as.factor(dimnames(spec_n)[[1]]),
                         w = rep(params@w,
                                 each = dim(spec_n)[[1]]))
  if (plankton) {
    # Decide where to cut off plankton
    max_w <- min(params@species_params$w_mat)
    if (is.na(max_w)) {
      max_w <- Inf
    }
    plankton_sel <- params@w_full >= min_w &
      params@w_full < max_w
    w_plankton <- params@w_full[plankton_sel]
    plank_n <- n_pp[plankton_sel] * w_plankton^power
    plot_dat <- rbind(plot_dat,
                      data.frame(value = c(plank_n),
                                 Species = "Plankton",
                                 w = w_plankton))
  }
  if (total) {
    plot_dat <- rbind(plot_dat,
                      data.frame(value = c(total_n),
                                 Species = "Total",
                                 w = params@w_full))
  }
  # lop off 0s and apply min_w
  plot_dat <- plot_dat[(plot_dat$value > 0) & (plot_dat$w >= min_w), ]
  # Impose ylim
  if (!is.na(ylim[1])) {
    plot_dat <- plot_dat[plot_dat$value < ylim[1], ]
  }
  if (is.na(ylim[2])) {
    ylim[2] <- 1e-20
  }
  plot_dat <- plot_dat[plot_dat$value > ylim[2], ]
  # Create plot
  p <- ggplot(plot_dat, aes(x = w, y = value)) +
    scale_x_continuous(name = "Size [g]", trans = "log10",
                       breaks = log_breaks()) +
    scale_y_continuous(name = y_label, trans = "log10",
                       breaks = log_breaks()) +
    scale_colour_manual(values = params@linecolour) +
    scale_linetype_manual(values = params@linetype)
  if (background) {
    back_n <- n[is.na(params@A), , drop = FALSE]
    plot_back <- data.frame(value = c(back_n),
                            Species = as.factor(dimnames(back_n)[[1]]),
                            w = rep(params@w,
                                    each = dim(back_n)[[1]]))
    # lop off 0s and apply min_w
    plot_back <- plot_back[(plot_back$value > 0) & (plot_back$w >= min_w), ]
    # Impose ylim
    if (!is.na(ylim[1])) {
      plot_back <- plot_back[plot_back$value < ylim[1], ]
    }
    plot_back <- plot_back[plot_back$value > ylim[2], ]
    # Add background species in grey
    p <- p +
      geom_line(aes(group = Species), colour = "grey",
                data = plot_back)
  }
  if ( (length(species) + plankton + total) > 13) {
    p <- p + geom_line(aes(group = Species))
  } else {
    p <- p + geom_line(aes(colour = Species, linetype = Species))
  }
  if (print_it)
    print(p)
  return(p)
}


#' Plot the feeding level of species by size
#' 
#' After running a projection, plot the feeding level of each species by size. 
#' The feeding level is averaged over the specified time range (a single value
#' for the time range can be used).
#' 
#' @param sim An object of class \linkS4class{MizerSim}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param ... Other arguments to pass to \code{getFeedingLevel} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getFeedingLevel}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotFeedingLevel(sim)
#' plotFeedingLevel(sim, time_range = 10:20)
#' }
plotFeedingLevel <- function(sim,
                             species = dimnames(sim@n)$sp,
                             time_range = max(as.numeric(dimnames(sim@n)$time)),
                             print_it = TRUE, returnData = F, ...) {
  feed_time <- getFeedingLevel(sim, time_range = time_range, ##AA
                               drop = FALSE, ...)
  feed <- apply(feed_time, c(2, 3), mean)
  feed <- feed[as.character(dimnames(feed)[[1]]) %in% species, ,
               drop = FALSE]
  plot_dat <- data.frame(value = c(feed),
                         Species = dimnames(feed)[[1]],
                         w = rep(sim@params@w, each = length(species)))
  if (length(species) > 12) {
    p <- ggplot(plot_dat) +
      geom_line(aes(x = w, y = value, group = Species))
  } else {
    p <- ggplot(plot_dat) +
      geom_line(aes(x = w, y = value, colour = Species, linetype = Species))
  }
  p <- p +
    scale_x_continuous(name = "Size [g]", trans = "log10") +
    scale_y_continuous(name = "Feeding Level", limits = c(0, 1)) +
    scale_colour_manual(values = sim@params@linecolour) +
    scale_linetype_manual(values = sim@params@linetype)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}


#' Plot predation mortality rate of each species against size
#' 
#' After running a projection, plot the predation mortality rate of each species
#' by size. The mortality rate is averaged over the specified time range (a
#' single value for the time range can be used to plot a single time step).
#' 
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param ... Other arguments to pass to \code{getM2} method.
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getM2}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotM2(sim)
#' plotM2(sim, time_range = 10:20)
#' }
plotM2 <- function(sim, species = dimnames(sim@n)$sp,
                   time_range = max(as.numeric(dimnames(sim@n)$time)),
                   print_it = TRUE, ...) {
  m2_time <- getM2(sim, time_range = time_range, intakeScalar = sim@intTempScalar[,,time_range], drop = FALSE, ...)
  m2 <- apply(m2_time, c(2, 3), mean)
  m2 <- m2[as.character(dimnames(m2)[[1]]) %in% species, , 
           drop = FALSE]
  plot_dat <- data.frame(value = c(m2),
                         Species = dimnames(m2)[[1]],
                         w = rep(sim@params@w, each = length(species)))
  if (length(species) > 12) {
    p <- ggplot(plot_dat) +
      geom_line(aes(x = w, y = value, group = Species))
  } else {
    p <- ggplot(plot_dat) +
      geom_line(aes(x = w, y = value, colour = Species, linetype = Species))
  }
  p <- p +
    scale_x_continuous(name = "Size [g]", trans = "log10") +
    scale_y_continuous(name = "Predation mortality [1/year]",
                       limits = c(0, max(plot_dat$value))) +
    scale_colour_manual(values = sim@params@linecolour) +
    scale_linetype_manual(values = sim@params@linetype)
  if (print_it) {
    print(p)
  }
  return(p)
}


#' Plot total fishing mortality of each species by size
#' 
#' After running a projection, plot the total fishing mortality of each species
#' by size. The total fishing mortality is averaged over the specified time
#' range (a single value for the time range can be used to plot a single time
#' step).
#' 
#' @param sim An object of class \linkS4class{MizerSim}.
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the final time step.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param ... Other arguments to pass to \code{getFMort} method
#'
#' @return A ggplot2 object
#' @export
#' @seealso \code{\link{getFMort}}
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotFMort(sim)
#' plotFMort(sim, time_range = 10:20)
#' }
plotFMort <- function(sim, species = dimnames(sim@n)$sp,
                      time_range = max(as.numeric(dimnames(sim@n)$time)),
                      print_it = TRUE, ...){
  f_time <- getFMort(sim, time_range = time_range, drop = FALSE, ...)
  f <- apply(f_time, c(2, 3), mean)
  f <- f[as.character(dimnames(f)[[1]]) %in% species, , drop = FALSE]
  plot_dat <- data.frame(value = c(f),
                         Species = dimnames(f)[[1]],
                         w = rep(sim@params@w, each = length(species)))
  if (length(species) > 12) {
    p <- ggplot(plot_dat) + geom_line(aes(x = w, y = value, group = Species))
  } else {
    p <- ggplot(plot_dat) +
      geom_line(aes(x = w, y = value, colour = Species, linetype = Species))
  }
  p <- p +
    scale_x_continuous(name = "Size [g]", trans = "log10") +
    scale_y_continuous(name = "Fishing mortality [1/Year]",
                       limits = c(0, max(plot_dat$value))) +
    scale_colour_manual(values = sim@params@linecolour) +
    scale_linetype_manual(values = sim@params@linetype)
  if (print_it) {
    print(p)
  }
  return(p)
}


#' Plot growth curves giving weight as a function of age
#' 
#' If given a \linkS4class{MizerSim} object, uses the growth rates at the final
#' time of a simulation to calculate the size at age. If given a
#' \linkS4class{MizerParams} object, uses the initial growth rates instead.
#' 
#' When the growth curve for only a single species is plotted, horizontal
#' lines are included that indicate the maturity size and the maximum size for 
#' that species. If furthermore the species parameters contain the variables
#' a and b for length to weight conversion and the von Bertalanffy parameter
#' k_vb, then the von Bertalanffy growth curve is superimposed in black.
#' 
#' @param object MizerSim or MizerParams object
#' @param species Name or vector of names of the species to be plotted. By
#'   default all species are plotted.
#' @param max_age The age up to which the weight is to be plotted. Default is 20
#' @param percentage Boolean value. If TRUE, the size is shown as a percentage
#'   of the maximal size.
#' @param print_it Display the plot, or just return the ggplot2 object.
#'   Defaults to TRUE
#' @param ... Other arguments (unused)
#' 
#' @return A ggplot2 object
#' @export
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plotGrowthCurves(sim, percentage = TRUE)
#' plotGrowthCurves(sim, species = "Cod", max_age = 24)
#' }
plotGrowthCurves <- function(object, species,
                             max_age = 20, percentage = TRUE, print_it = TRUE) {
  if (is(object, "MizerSim")) {
    sim <- object
    if (missing(species)) {
      species <- dimnames(sim@n)$sp
    }
    # reorder list of species to coincide with order in sim
    idx <- which(dimnames(sim@n)$sp %in% species)
    species <- dimnames(sim@n)$sp[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list("Species" = species, "Age" = age))
    g <- getEGrowth(sim@params, sim@n[dim(sim@n)[1], , ], 
                    sim@n_pp[dim(sim@n)[1], ], sim@n_bb[dim(sim@n)[1], ], sim@n_aa[dim(sim@n)[1], ], sim@intTempScalar[,,1], sim@metTempScalar[,,1]) #AA
    for (j in 1:length(species)) {
      i <- idx[j]
      g_fn <- stats::approxfun(sim@params@w, g[i, ])
      myodefun <- function(t, state, parameters){
        return(list(g_fn(state)))
      }
      ws[j, ] <- deSolve::ode(y = sim@params@species_params$w_min[i],
                              times = age, func = myodefun)[, 2]
      if (percentage) {
        ws[j, ] <- ws[j, ] / sim@params@species_params$w_inf[i] * 100
      }
    }
    plot_dat <- reshape2::melt(ws)
    plot_dat$Species <- as.character(plot_dat$Species)
    if (length(species) > 12) {
      p <- ggplot(plot_dat) +
        geom_line(aes(x = Age, y = value, group = Species))
    } else {
      p <- ggplot(plot_dat) +
        geom_line(aes(x = Age, y = value,
                      colour = Species, linetype = Species))
    }
    y_label <- if (percentage) "Percent of maximum size" else "Size [g]"
    p <- p +
      scale_x_continuous(name = "Age [Years]") +
      scale_y_continuous(name = y_label) +
      scale_colour_manual(values = sim@params@linecolour) +
      scale_linetype_manual(values = sim@params@linetype)
    
    # Extra stuff for single-species case
    if (length(species) == 1 && !percentage) {
      w_inf <- sim@params@species_params$w_inf[idx[1]]
      p <- p + geom_hline(yintercept = w_inf) +
        annotate("text", 0, w_inf, vjust = -1, label = "Maximum")
      w_mat <- sim@params@species_params$w_mat[idx[1]]
      p <- p + geom_hline(yintercept = w_mat) +
        annotate("text", 0, w_mat, vjust = -1, label = "Maturity")
      if (all(c("a", "b", "k_vb") %in% names(sim@params@species_params))) {
        a <- sim@params@species_params$a[idx[1]]
        b <- sim@params@species_params$b[idx[1]]
        k_vb <- sim@params@species_params$k_vb[idx[1]]
        L_inf <- (w_inf/a)^(1/b)
        vb <- a * (L_inf * (1 - exp(-k_vb * age)))^b
        dat <- data.frame("x" = age, "y" = vb)
        p <- p + geom_line(data = dat, aes(x = x, y = y))
      }
    }
    
    if(percentage) p <- p + scale_y_continuous(breaks = seq(0,100,by = 20))
    
    if (print_it) return(p)
  } else {
    # Plot growth curves using a MizerParams object.
    sim <- project(object, t_max = 1) # construct the temperature scalars
    params <- object
    if (missing(species)) {
      species <- dimnames(params@initial_n)$sp
    }
    # reorder list of species to coincide with order in params
    idx <- which(dimnames(params@initial_n)$sp %in% species)
    species <- dimnames(params@initial_n)$sp[idx]
    age <- seq(0, max_age, length.out = 50)
    ws <- array(dim = c(length(species), length(age)),
                dimnames = list(Species = species, Age = age))
    g <- getEGrowth(params, params@initial_n, params@initial_n_pp, params@initial_n_bb, params@initial_n_aa, 
                    intakeScalar = sim@intTempScalar[,,1], metScalar = sim@metTempScalar[,,1]) ##AA
    for (j in 1:length(species)) {
      i <- idx[j]
      g_fn <- stats::approxfun(params@w, g[i, ])
      myodefun <- function(t, state, parameters){
        return(list(g_fn(state)))
      }
      ws[j, ] <- deSolve::ode(y = params@species_params$w_min[i], 
                              times = age, func = myodefun)[, 2]
      if (percentage) {
        ws[j, ] <- ws[j, ] / params@species_params$w_inf[i] * 100
      }
    }
    plot_dat <- reshape2::melt(ws)
    plot_dat$Species <- as.character(plot_dat$Species)
    if (length(species) > 12) {
      p <- ggplot(plot_dat) +
        geom_line(aes(x = Age, y = value, group = Species))
    } else {
      p <- ggplot(plot_dat) +
        geom_line(aes(x = Age, y = value,
                      colour = Species, linetype = Species))
    }
    y_label <- if (percentage) "Percent of maximum size" else "Size [g]"
    p <- p +
      scale_x_continuous(name = "Age [Years]") +
      scale_y_continuous(name = y_label) +
      scale_colour_manual(values = params@linecolour) +
      scale_linetype_manual(values = params@linetype)
    
    # Extra stuff for single-species case
    if (length(species) == 1 && !percentage) {
      w_inf <- params@species_params$w_inf[idx[1]]
      p <- p + geom_hline(yintercept = w_inf) +
        annotate("text", 0, w_inf, vjust = -1, label = "Maximum")
      w_mat <- params@species_params$w_mat[idx[1]]
      p <- p + geom_hline(yintercept = w_mat) +
        annotate("text", 0, w_mat, vjust = -1, label = "Maturity")
      if (all(c("a", "b", "k_vb") %in% names(params@species_params))) {
        a <- params@species_params$a[idx[1]]
        b <- params@species_params$b[idx[1]]
        k_vb <- params@species_params$k_vb[idx[1]]
        L_inf <- (w_inf/a)^(1/b)
        vb <- a * (L_inf * (1 - exp(-k_vb * age)))^b
        dat <- data.frame(x = age, y = vb)
        p <- p + geom_line(data = dat, aes(x = x, y = y))
      }
    }
    
    if (print_it) return(p)
  }
}


#' Plot diets composition of by predator / prey and their sizes 

plotDietComp<-function(object, prey=dimnames(object@diet_comp)$prey, min_w=.001,
                       predator=dimnames(object@diet_comp)$predator, timeaverage=FALSE, print_it = T){
  
  prey_nam<-prey
  pred_nam<-predator
  
  out<-object@diet_comp 
  
  prey<-apply(out, c(1,2,3), FUN=sum) #Sum across size classess with in prey 
  tot<-apply(prey, c(1,2), FUN=sum) #Sum across prey species 
  
  prey_prop<-sweep(prey, c(1,2), tot, "/") # Get proportion of diet for each species
  
  no_pred<- length(dimnames(prey_prop)[[1]])
  no_pred_w<- length(dimnames(prey_prop)[[2]])
  no_prey<- length(dimnames(prey_prop)[[3]])
  
  #Stacked  bar chart 
  plot_dat<-expand.grid(dimnames(prey_prop)[[1]], dimnames(prey_prop)[[2]], dimnames(prey_prop)[[3]])
  colnames(plot_dat)<-c("predator","predsize","prey")
  plot_dat$predsize<-as.numeric(as.character(log10(as.numeric(as.character(plot_dat$predsize)))))
  
  plot_dat$value<- as.vector(prey_prop)
  
  
  species<-object@params@species_params$species
  wmin<-object@params@w[object@params@species_params$w_min_idx]
  wmax<-object@params@w[object@params@species_params$w_max_idx]
  
  for ( i in 1:length(species)){
    plot_dat$value[plot_dat$predator==species[i] & plot_dat$predsize < log10(wmin[i])]<- 0
    plot_dat$value[plot_dat$predator==species[i] & plot_dat$predsize > log10(wmax[i])]<- 0
  }
  
  
  plot_dat<-subset(plot_dat,predsize> log10(min_w))
  
  dsub<-plot_dat[ plot_dat$prey %in% prey_nam, ]
  dsub<-dsub[ dsub$predator %in% pred_nam, ]
  dsub[is.na(dsub)]<-0
  
  p<-  ggplot(data = dsub, aes(x = predsize, y = value, fill = prey)) + geom_area( position = 'stack')  + facet_wrap(~predator, ncol=5) + scale_color_brewer(palette="Set1") +
    scale_x_continuous(name = "log10 predator mass (g)") + scale_y_continuous(name = "Proportion of diet by mass (g)")
  
  if(print_it) return(p)
  
}

#' Plot PPRM values for the selected time period based on diet compositions 

plotPPMR<-function(object=object, grid=T, observed=FALSE, prey=dimnames(object@diet_comp)$prey, 
                   predator=dimnames(object@diet_comp)$predator, timeaverage=FALSE ){
  
  prey_nam<-prey
  pred_nam<-predator
  
  out<-object@diet_comp 
  
  prey<-apply(out, c(1,2,4), FUN=sum) #Sum across size classess with in prey 
  tot<-apply(prey, c(1,2), FUN=sum) #Sum across prey weight size classes 
  
  prop_prey<-sweep(prey, 1:2, tot, "/" ) #proportion of prey weight in each each prey size class
  prop_prey[is.na(prop_prey)]<-0
  
  #make matrix of realized PPMR
  ppmr_mat<-outer(object@params@w, object@params@w_full, FUN="/")
  ppmr_frac<-sweep(prop_prey, 2:3, ppmr_mat, "*")
  ppmr_tot<-apply(ppmr_frac, c(1,2), FUN=sum) #Sum across prey weight size classes 
  
  plot_dat<-expand.grid(dimnames(ppmr_tot)[[1]], dimnames(ppmr_tot)[[2]])
  colnames(plot_dat)<-c("predator","predsize")
  plot_dat$predsize<-as.numeric(as.character(log10(as.numeric(as.character(plot_dat$predsize)))))
  
  plot_dat$value<- as.vector(ppmr_tot)
  
  
  species<-object@params@species_params$species
  wmin<-object@params@w[object@params@species_params$w_min_idx]
  wmax<-object@params@w[object@params@species_params$w_max_idx]
  
  for ( i in 1:length(species)){
    plot_dat$value[plot_dat$predator==species[i] & plot_dat$predsize < log10(wmin[i])]<- 0
    plot_dat$value[plot_dat$predator==species[i] & plot_dat$predsize > log10(wmax[i])]<- 0
  }
  
  plot_dat$value[plot_dat$value==0]<-NA
  
  plot_dat<-plot_dat[!plot_dat$predator=="Benthos",]
  plot_dat<-plot_dat[!plot_dat$predator=="Detritus",]
  
  #Plot realized vs prefferred PPMR; for individuals over 10 g 
  
  plot_dat$preferred
  
  ma<-match(plot_dat$predator, object@params@species_params$species)
  
  plot_dat$preferred<-object@params@species_params$beta[ma]
  plot_dat$obs_PPMR<- object@params@species_params$obs_beta[ma] #originally I assumed that preferred ppmr was 1.7 times higher than realized (following Harvig simulations)
  
  
  #PPMR vs body size
  
  if(observed==FALSE & grid==TRUE){
    p <- ggplot(plot_dat) + geom_line(aes(x=predsize, y = value,  colour = predator, size = 1)) + ggtitle("PPMRbio for predators > 10 g") +
      scale_size(range = c(2))+ scale_x_continuous(name = "log10 predator mass (g)", trans="log10") + 
      scale_y_continuous(name = "log10 realized PPMRbio", trans="log10") + geom_abline(intercept = 0, slope=1) + facet_wrap(~predator, ncol=5)
  } 
  
  if(observed==FALSE & grid==FALSE){
    p <- ggplot(plot_dat) + geom_line(aes(x=predsize, y = value,  colour = predator, size = 1)) + ggtitle("PPMRbio for predators > 10 g") +
      scale_size(range = c(2))+ scale_x_continuous(name = "log10 predator mass (g)", trans="log10") + 
      scale_y_continuous(name = "log10 realized PPMRbio", trans="log10") + geom_abline(intercept = 0, slope=1) 
  }
  
  if(observed==TRUE & grid==TRUE){
    p <- ggplot(plot_dat) + geom_line(aes(x=obs_PPMR, y = value,  colour = predator, size = 1)) + ggtitle("PPMRbio for predators > 10 g") +
      scale_size(range = c(2))+ scale_x_continuous(name = "log10 observed realPPMRbio", trans="log10") + 
      scale_y_continuous(name = "log10 simulated realPPMRbio", trans="log10") + geom_abline(intercept = 0, slope=1) + facet_wrap(~predator, ncol=5)
  }
  
  if(observed==TRUE & grid==FALSE){
    p <- ggplot(plot_dat) + geom_line(aes(x=obs_PPMR, y = value,  colour = predator, size = 1)) + ggtitle("PPMRbio for predators > 10 g") +
      scale_size(range = c(2))+ scale_x_continuous(name = "log10 observed realPPMRbio", trans="log10") + 
      scale_y_continuous(name = "log10 simulated realPPMRbio", trans="log10") + geom_abline(intercept = 0, slope=1) 
  }
  
  
  print(p)
  return(p)
}







#### plot ####
#' Summary plot for \code{MizerSim} objects
#' 
#' After running a projection, produces 5 plots in the same window: feeding
#' level, abundance spectra, predation mortality and fishing mortality of each
#' species by size; and biomass of each species through time. This method just
#' uses the other plotting methods and puts them all in one window.
#' 
#' @param x An object of class \linkS4class{MizerSim}
#' @param y Not used
#' @param ...  For additional arguments see the documentation for
#'   \code{\link{plotBiomass}},
#'   \code{\link{plotFeedingLevel}},\code{\link{plotSpectra}},\code{\link{plotM2}}
#'   and \code{\link{plotFMort}}.
#' @return A viewport object
#' @export
#' @rdname plotMizerSim
#' @examples
#' \dontrun{
#' data(NS_species_params_gears)
#' data(inter)
#' params <- MizerParams(NS_species_params_gears, inter)
#' sim <- project(params, effort=1, t_max=20, t_save = 2)
#' plot(sim)
#' plot(sim, time_range = 10:20) # change time period for size-based plots
#' plot(sim, min_w = 10, max_w = 1000) # change size range for biomass plot
#' }
setMethod("plot", signature(x = "MizerSim", y = "missing"),
          function(x, ...) {
            p1 <- plotFeedingLevel(x, print_it = FALSE, ...)
            p2 <- plotSpectra(x, print_it = FALSE, ...)
            p3 <- plotBiomass(x, y_ticks = 3, print_it = FALSE, ...)
            p4 <- plotM2(x, print_it = FALSE, ...)
            p5 <- plotFMort(x, print_it = FALSE, ...)
            grid::grid.newpage()
            glayout <- grid::grid.layout(3, 2) # widths and heights arguments
            vp <- grid::viewport(layout = glayout)
            grid::pushViewport(vp)
            vplayout <- function(x, y)
              grid::viewport(layout.pos.row = x, layout.pos.col = y)
            print(p1 + theme(legend.position = "none"), vp = vplayout(1, 1))
            print(p3 + theme(legend.position = "none"), vp = vplayout(1, 2))
            print(p4 + theme(legend.position = "none"), vp = vplayout(2, 1))
            print(p5 + theme(legend.position = "none"), vp = vplayout(2, 2))
            print(p2 + theme(legend.position = "right",
                             legend.key.size = unit(0.1, "cm")),
                  vp = vplayout(3, 1:2))
          }
)




##@#%@#$%^&*# Romain's plots #@%*^*@##---------------------

plotDynamics <- function(object, time_range = c(min(as.numeric(dimnames(object@n)$time)),max(as.numeric(dimnames(object@n)$time))), 
                         phenotype = TRUE, species = NULL, trait = NULL, SpIdx = NULL, print_it = T, returnData = F, save_it = F, nameSave = "Biomass.png", ylimit = c(NA,NA)){
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")) # colorful gradient
  min_value <- 1e-30
  # get the phenotype biomass through time (need at least 2 time steps for now)
  biomass <- getBiomass(object)
  time_elements <- get_time_elements(object,time_range)
  biomass <- biomass[time_elements,]
  
  # getting rid of the species that went extinct during the initialisation
  if (is.null(SpIdx))
    for (i in unique(object@params@species_params$species))
      if (sum(biomass[, i]) != 0) # & dim(object@params@species_params[object@params@species_params$species == i, ])[1] != 1)
        SpIdx = c(SpIdx, i)
  
  # sum the phenotype biomass per species
  biomassSp = NULL
  biomassTemp = biomass
  colnames(biomassTemp) = object@params@species_params$species
  for (i in SpIdx)
  {
    biomassPhen = biomassTemp[,which(colnames(biomassTemp) == i)]
    if(!is.null(dim(biomassPhen))) biomassPhen = apply(biomassPhen,1,sum)
    biomassSp = cbind(biomassSp,biomassPhen)
  }
  colnames(biomassSp) = SpIdx
  
  # apply SpIdx on biomass as well
  spSub <- object@params@species_params$ecotype[object@params@species_params$species %in% SpIdx]
  biomass <- biomass[,as.numeric(dimnames(biomass)$species) %in% spSub]
  
  plotBiom <- function(x)
  {
    Biom <- melt(x) # melt for ggplot
    colnames(Biom) = c("time","phen","value")
    
    # create a species column
    Biom$sp = sapply(Biom$phen, function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
    return(Biom)
  }
  BiomSp <- plotBiom(biomassSp)
  BiomSp <- BiomSp[BiomSp$value >= min_value,]
  
  
  if (phenotype) 
  {
    
    BiomPhen <- plotBiom(biomass)
    if(!is.null(trait))
      BiomPhen$trait <- rep(trait, each = length(unique(BiomPhen$time)))
    
    BiomPhen <- BiomPhen[BiomPhen$value >= min_value,]
    
    p <- ggplot(BiomSp) +
      geom_line(aes(x = time, y = value, colour = as.factor(sp), group = sp), size = 1.2) +
      geom_line(data = BiomPhen, aes(x = time, y = value, colour = as.factor(sp), group = phen), alpha = 0.2) +
      scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30,4,2)))) +
      scale_x_continuous(name = "Time in years") +
      labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Community biomass") 
    
    
    if (!is.null(species)) 
    {
      BiomPhen <- BiomPhen[BiomPhen$sp == species, ]
      BiomSp <- BiomSp[BiomSp$sp == species, ]
      plotTitle <- paste("Species",species)
      
      p <- ggplot(BiomPhen) +
        geom_line(aes(x = time, y = value, group = phen), alpha = .75) +
        scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30,4,2)))) +
        scale_x_continuous(name = "Time in years") +
        # scale_colour_gradientn(colours=jet.colors(9), limits = c(NA,NA))+
        theme(panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),
              legend.key = element_rect(fill = "white"))+
        ggtitle(plotTitle) 
      
      if(!is.null(trait))
      {
        plotTitle <- paste("Trait of species",species)
        
        p <- ggplot(BiomPhen) +
          geom_line(aes(x = time, y = value, colour = trait, group = trait), alpha = 1) +
          scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30,4,2)))) +
          scale_x_continuous(name = "Time in years") +
          scale_colour_gradientn(colours=jet.colors(9), limits = c(NA,NA))+
          theme(panel.background = element_rect(fill = "white", color = "black"),
                panel.grid.minor = element_line(colour = "grey92"),
                legend.key = element_rect(fill = "white"))+
          ggtitle(plotTitle) 
      }
    }
    
  } else {
    
    p <- ggplot(BiomSp) +
      geom_line(aes(x = time, y = value, colour = as.factor(sp), group = sp), size = 1.2) +
      scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(-30:4))) +
      scale_x_continuous(name = "Time in years") +
      labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Community biomass") 
  }
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(BiomSp) else if(print_it) return(p)
}


plotSS <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), min_w =min(object@params@w)/100, ylim = c(NA,NA),
                   biomass = TRUE, print_it = TRUE, species = TRUE, community = FALSE, save_it = FALSE, nameSave = "SizeSpectrum.png", returnData = FALSE, ...){
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  # min_w = 0.001
  time_elements <- get_time_elements(object,time_range)
  spec_n <- apply(object@n[time_elements,,,drop=FALSE],c(2,3), mean)
  pkt_n <- apply(object@n_pp[time_elements,,drop=FALSE],2,mean)
  alg_n <- apply(object@n_aa[time_elements,,drop=FALSE],2,mean)
  ben_n <- apply(object@n_bb[time_elements,,drop=FALSE],2,mean)
  
  y_axis_name = "Abundance"
  if (biomass){
    spec_n <- sweep(spec_n,2,object@params@w,"*")
    pkt_n <- pkt_n * object@params@w_full
    alg_n <- alg_n * object@params@w_full
    ben_n <- ben_n * object@params@w_full
    y_axis_name = "Biomass"
  }
  # Make data.frame for plot
  plot_datSP <- data.frame(value = c(spec_n), Species = dimnames(spec_n)[[1]], w = rep(object@params@w, each=nrow(object@params@species_params)), bloodline = object@params@species_params$species)
  plot_datPkt <- data.frame(value = c(pkt_n), Species = "Phytoplankton", w = object@params@w_full)
  plot_datAlg <- data.frame(value = c(alg_n), Species = "Algae", w = object@params@w_full)
  plot_datBen <- data.frame(value = c(ben_n), Species = "Benthos", w = object@params@w_full)
  
  if(community) plot_datSP <- data.frame(value = apply(spec_n, 2, sum), w = object@params@w)

  else if (species)
  {
    dimnames(spec_n)$species = object@params@species_params$species
    SpIdx = unique(object@params@species_params$species)
    spec_sp = matrix(data = NA, ncol = dim(spec_n)[2], nrow = length(SpIdx), dimnames = list(as.character(SpIdx),dimnames(spec_n)$size))
    names(dimnames(spec_sp))=list("species","size")
    
    for (i in 1:dim(spec_sp)[1])
    {
      temp = spec_n # save to manip
      temp[which(rownames(spec_n) != i), ] = 0 # make everything but the targeted species to go 0 to have correct normalisation
      temp = apply(temp, 2, sum)
      spec_sp[i, ] = temp
    }
    plot_datSP <- data.frame(value = c(spec_sp), Species = dimnames(spec_sp)[[1]], w = rep(object@params@w, each=length(SpIdx)))
  }
  # lop off 0s in background and apply min_w
  plot_datSP <- plot_datSP[(plot_datSP$value > 0) & (plot_datSP$w >= min_w),]
  plot_datPkt <- plot_datPkt[(plot_datPkt$value > 0) & (plot_datPkt$w >= min_w),]
  plot_datAlg <- plot_datAlg[(plot_datAlg$value > 0) & (plot_datAlg$w >= min_w),]
  plot_datBen <- plot_datBen[(plot_datBen$value > 0) & (plot_datBen$w >= min_w),]
  #getPalette = colorRampPalette(brewer.pal(9, "Set1"))# increase the number of colors used
  
  if(community)
  {
    p <- ggplot(plot_datSP) + 
      geom_line(aes(x=w, y = value)) + 
      geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species),alpha = 0.5, color = "blue", size = 1.5) +
      geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
      geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
      scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
      scale_y_continuous(name = "Abundance density in individuals.m^-3", limits = ylim, trans = "log10") +
      # labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Community spectrum") 
  }
  else if (species)
  {
    p <- ggplot(plot_datSP) + 
      geom_line(aes(x=w, y = value, colour = as.factor(Species), group = Species)) + 
      geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species),alpha = 0.5, color = "blue", size = 1.5) +
      geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
      geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
      scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
      scale_y_continuous(name = "Abundance density in individuals.m^-3", limits = ylim, trans = "log10") +
      labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Size spectrum")
  }
  
  else
  {
    p <- ggplot(plot_datSP) + 
      geom_line(aes(x=w, y = value, colour = as.factor(bloodline), group = Species)) + 
      geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species),alpha = 0.5, color = "blue", size = 1.5) +
      geom_line(data = plot_datAlg, aes(x = w, y = value, group = Species),alpha = 0.5, color = "green", size = 1.5) +
      geom_line(data = plot_datBen, aes(x = w, y = value, group = Species),alpha = 0.5, color = "yellow", size = 1.5) +
      scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5)))+
      scale_y_continuous(name = "Abundance density in individuals.m^-3", limits = ylim, trans = "log10") +
      labs(color='Species') +
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"),
            legend.key = element_rect(fill = "white"))+
      scale_colour_manual(values=cbPalette)+ # colorblind
      ggtitle("Size spectrum")
    
  }
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_datSP) else if(print_it) return(p)
}

plotFood <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, throughTime = F, start = 1000, every = 1000, 
                     print_it = T, returnData = F, save_it =F, nameSave = "Feeding.png"){
  
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  if (throughTime)
  {
    time_range = seq(start,max(as.numeric(dimnames(object@n)$time)),every)
    time_range = c(time_range,max(as.numeric(dimnames(object@n)$time))) # so it counts the last time step which is probably not even
    time_range = unique(time_range)
    feeding = array(data = NA, dim = c(length(unique(object@params@species_params$species)),100,length(time_range)),  
                    dimnames = list(as.character(unique(object@params@species_params$species)),object@params@w,time_range)) 
    Critfeeding = matrix(data=NA, nrow = length(time_range), ncol= 100, dimnames = list(time_range,object@params@w))
    for (i in time_range)
    {
      
      feed_time <- getFeedingLevel(object=object, time_range=i, drop=FALSE)#, ...) # get the feeding time
      feed <- apply(feed_time, c(2,3), mean) # average on the time frame
      
      Cfeed_time <- getCriticalFeedingLevel(object=object, time_range=i, drop=FALSE)#, ...) # get the critical feeding level
      Critfeed <- apply(Cfeed_time, c(2,3), mean) # average on the time frame
      Critfeed <- Critfeed[1,] # all rows the same
      
      dimnames(feed)$sp = object@params@species_params$species
      SpIdx = unique(object@params@species_params$species) # get the species names
      feed_sp = matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(feed)$w)) # prepare the new object
      names(dimnames(feed_sp))=list("species","size")
      
      for (j in SpIdx)
      {
        temp = feed # save to manip
        temp[which(rownames(feed) != j), ] = 0 # keep the ecotypes from the species only
        temp = apply(temp, 2, sum)
        temp = temp / length(which(rownames(feed)==j)) # do the mean (in 2 steps)
        feed_sp[which(rownames(feed_sp)==j), ] = temp
      }
      feeding[,,which(dimnames(feeding)[[3]] == i)] = feed_sp
      Critfeeding[which(dimnames(Critfeeding)[[1]] == i),] = Critfeed
    }
    a <- c(object@params@species_params$w_inf[SpIdx]) # to get vline of different col, need to create a data frame
    vlines <- data.frame(xint = a,grp = SpIdx)
    
    plot_dat = melt(feeding)
    colnames(plot_dat) = c("species","size","time","value")
    plot_crit = melt(Critfeeding)
    colnames(plot_crit) = c("time","size","value")
    p <- ggplot(plot_dat) + 
      geom_line(aes(x=size, y = value, colour = as.factor(species))) + 
      geom_line(data = plot_crit, aes(x = size, y = value), linetype = "dashed") +
      scale_x_log10(name = "Size") + 
      scale_y_continuous(name = "Feeding Level", lim=c(0,1))+
      geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
      facet_grid(time ~ .)+
      scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
      theme(panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"))+
      
      ggtitle("Feeding level through time")
    
    if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
    
    if (returnData) return(list(plot_dat,plot_crit)) else if(print_it) return(p)
    
  }
  
  feed_time <- getFeedingLevel(object=object, time_range=time_range, drop=FALSE) #, ...) # get the feeding time
  feed <- apply(feed_time, c(2,3), mean) # average on the time frame
  
  Cfeed_time <- getCriticalFeedingLevel(object=object, time_range=time_range, drop=FALSE)#, ...) # get the critical feeding level
  Critfeed <- apply(Cfeed_time, c(2,3), mean) # average on the time frame
  Critfeed <- Critfeed[1,] # all rows the same
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(feed)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    feed_sp = matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(feed)$w)) # prepare the new object
    names(dimnames(feed_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = feed # save to manip
      temp[which(rownames(feed) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(feed)==i)) # do the mean (in 2 steps)
      feed_sp[which(rownames(feed_sp)==i), ] = temp
    }
    feed = feed_sp
  }
  
  a <- c(object@params@species_params$w_inf[SpIdx]) # to get vline of different col, need to create a data frame
  vlines <- data.frame(xint = a,grp = SpIdx)
  
  plot_dat <- data.frame(value = c(feed), species = dimnames(feed)[[1]], size = rep(object@params@w, each=length(dimnames(feed)[[1]])))
  
  name = paste("Feeding level at time",time_range,sep=" ")
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=size, y = value, colour = as.factor(species))) + 
    geom_hline(yintercept = Critfeed[1], linetype = "dashed", color = "red") +
    geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
    scale_x_log10(name = "Size", breaks = c(1 %o% 10^(-3:5)))  + 
    scale_y_continuous(name = "Feeding Level", lim=c(0,1))+
    scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),
          legend.key = element_rect(fill = "white"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave,width = 18, height = 18,units = "cm" )
  
  if (returnData) return(list(plot_dat,Critfeed)) else if(print_it) return(p)
}

plotGrowth <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, print_it = T, returnData = F, save_it = F, ylim = c(NA,NA),
                       nameSave = "Growth.png",...){
  
  time_elements <- get_time_elements(object,time_range)
  growth_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    growth <- getEGrowth(object@params, n=n, n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                         intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
    return(growth)})
  
  #growth <- apply(growth_time, c(2,3), mean) # use this when I will have time_range on more than one time
  growth = growth_time
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(growth)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    growth_sp = matrix(data = NA, ncol = dim(growth)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(growth)$w)) # prepare the new object
    names(dimnames(growth_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = growth # save to manip
      temp[which(rownames(growth) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(growth)==i)) # do the mean (in 2 steps)
      growth_sp[which(rownames(growth_sp)==i), ] = temp
    }
    growth = growth_sp
  }
  
  name = paste("Growth level at time",time_range,sep=" ")
  plot_dat <- data.frame(value = c(growth), Species = dimnames(growth)[[1]], w = rep(object@params@w, each=length(dimnames(growth)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
    scale_y_continuous(name = "instantaneous growth", trans ="log10", limits = ylim)+
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotStarvation <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, print_it = T, returnData = F, save_it =F, 
                           nameSave = "Starvation.png"){
  
  
  cbPalette <- c("#999999","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7") #9 colors for colorblind
  
  # death_time <- getSmort(object=object, time_range=what_time, drop=FALSE)
  
  
  time_elements <- get_time_elements(object,time_range)
  death_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    starv <- getSMort(object@params, n=n, n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                      intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
    return(starv)})
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(death_time)[[1]] = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    death_sp = matrix(data = NA, ncol = dim(death_time)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(death_time)$w)) # prepare the new object
    names(dimnames(death_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = death_time # save to manip
      temp[which(rownames(death_time) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(death_time)==i)) # do the mean (in 2 steps)
      death_sp[which(rownames(death_sp)==i), ] = temp
    }
    death = death_sp
  }
  
  # a <- c(object@params@species_params$w_inf[1:9]) # to get vline of different col, need to create a data frame
  # vlines <- data.frame(xint = a,grp = c(1:9))
  
  plot_dat <- data.frame(value = c(death), species = dimnames(death)[[1]], size = rep(object@params@w, each=length(dimnames(death)[[1]])))
  
  name = paste("Starvation at time",time_range,sep=" ")
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=size, y = value, colour = as.factor(species))) + 
    #geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
    scale_x_log10(name = "Size", breaks = c(1 %o% 10^(-3:5)))  + 
    scale_y_continuous(name = "Instantaneous starvation mortality")+
    scale_colour_manual(values=cbPalette, name = "Species")+ # colorblind
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),
          legend.key = element_rect(fill = "white"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotScythe <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),print_it = TRUE, returnData = F, comments = T){
  
  # effort can be in 2 forms
  
  if(is.matrix(object@effort)) effort = object@effort[time_range,]
  else effort = object@effort[time_range]
  
  z <- getZ(object = object@params, n = object@n[time_range,,], n_pp = object@n_pp[time_range,],n_aa = object@n_aa[time_range,],n_bb = object@n_bb[time_range,],
            effort = effort, 
            intakeScalar = object@intTempScalar[,,time_range], metScalar = object@metTempScalar[,,time_range], morScalar = object@morTempScalar[,,time_range])
  dimnames(z)$prey = object@params@species_params$species
  #SpIdx = sort(unique(object@params@species_params$species)) # get the species names
  
  # need to get rid of the extinct species at that time in SpIdx
  a <- apply(object@n[time_range,,],1,sum)
  
  names(a) <- sapply(names(a), function(x) as.numeric(unlist(strsplit(as.character(x), "")))[1])
  
  d <- rowsum(a, group = names(a))
  
  if (sum(d[,1] == 0)) 
  {
    
    d <- d[-which(d[,1] == 0),]
    SpIdx <- as.numeric(names(d))
    
    
  } else {SpIdx <- as.numeric(rownames(d))}
  
  z_sp = matrix(data = NA, ncol = dim(z)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(z)$w_prey)) # prepare the new object
  names(dimnames(z_sp))=list("prey","w_prey")
  
  for (i in SpIdx)
  {
    temp = z # save to manip
    temp[which(rownames(z) != i), ] = 0 # keep the ecotypes from the species only
    temp = apply(temp, 2, sum)
    temp = temp / length(which(rownames(z)==i)) # do the mean (in 2 steps)
    z_sp[which(rownames(z_sp)==i), ] = temp
  }
  z = z_sp
  
  name = paste("Total Mortality at time",time_range,sep=" ")
  
  plot_dat <- data.frame(value = c(z), Species = dimnames(z)[[1]], w = rep(object@params@w, each=length(dimnames(z)[[1]])))
  
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
    scale_y_continuous(name = "Mortality", lim=c(0,max(plot_dat$value))) +
    theme(legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotSpawn <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, print_it = T, returnData = F, save_it = F, 
                      nameSave = "Spawn.png",...){
  
  time_elements <- get_time_elements(object,time_range)
  spawn_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    spawn <- getESpawning(object@params, n=n, n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                          intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
    return(spawn)})
  
  #spawn <- apply(spawn_time, c(2,3), mean) # use this when I will have time_range on more than one time
  spawn = spawn_time
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(spawn)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    spawn_sp = matrix(data = NA, ncol = dim(spawn)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(spawn)$w)) # prepare the new object
    names(dimnames(spawn_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = spawn # save to manip
      temp[which(rownames(spawn) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(spawn)==i)) # do the mean (in 2 steps)
      spawn_sp[which(rownames(spawn_sp)==i), ] = temp
    }
    spawn = spawn_sp
  }
  
  a <- c(object@params@species_params$w_inf[1:9]) # to get vline of different col, need to create a data frame
  vlines <- data.frame(xint = a,grp = c(1:9))
  
  name = paste("Spawn level at time",time_range,sep=" ")
  plot_dat <- data.frame(value = c(spawn), Species = dimnames(spawn)[[1]], w = rep(object@params@w, each=length(dimnames(spawn)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
    scale_y_continuous(name = "Energy allocated to spawning", trans ="log10")+
    geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotPredRate <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, ylim = c(NA,NA),
                         print_it = T, returnData = F, save_it = F, nameSave = "PredRate.png",...){
  
  time_range = max(as.numeric(dimnames(object@n)$time))
  time_elements <- get_time_elements(object,time_range)
  spawn_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    spawn <- getPredRate(object@params, n=n, n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                         intakeScalar = object@intTempScalar[,,x])
    return(spawn)})
  
  #spawn <- apply(spawn_time, c(2,3), mean) # use this when I will have time_range on more than one time
  spawn = spawn_time
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(spawn)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    spawn_sp = matrix(data = NA, ncol = dim(spawn)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(spawn)$w)) # prepare the new object
    names(dimnames(spawn_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = spawn # save to manip
      temp[which(rownames(spawn) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(spawn)==i)) # do the mean (in 2 steps)
      spawn_sp[which(rownames(spawn_sp)==i), ] = temp
    }
    spawn = spawn_sp
  }
  
  a <- c(object@params@species_params$w_inf[1:9]) # to get vline of different col, need to create a data frame
  vlines <- data.frame(xint = a,grp = c(1:9))
  
  name = paste("Predation rate at time",time_range,sep=" ")
  plot_dat <- data.frame(value = c(spawn), Species = dimnames(spawn)[[1]], w = rep(object@params@w_full, each=length(dimnames(spawn)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-7:5))) + 
    scale_y_continuous(name = "Potential death rate from predator", trans ="log10", limits = ylim)+
    geom_vline(data = vlines,aes(xintercept = xint,colour = as.factor(grp)), linetype = "dashed") + 
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}


plotCohort <- function(object, dt = 0.1, t_steps = 5, iSpecies = NULL, effort = 0, traitID = object@params@species_params$ed_int,
                       cohortSpan = c(1,dim(object@n)[1]-t_steps/dt), simStart = 0,
                       print_it = T, returnData = F, save_it = F, nameSave = paste("CohortSpecies",iSpecies,".png",sep=""))
{
  srrMatrix <- function(rdi,species_params) # srr functions handling matrixes instead of vectors
  {
    rdiSum <- apply(rdi,2,sum)
    rdiNormal = vector(mode = "numeric", length = length(rdiSum))
    names(rdiSum) <- species_params$species
    for (i in unique(species_params$species))
    {
      rdiSp = rdiSum # save to manip
      rdiSp[which(names(rdiSum) != i)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
      
      for (i in 1:length(rdiSp))
        # in case of NA
        if (is.na(rdiSp[i]) == TRUE)
          rdiSp[i] = 1e-30
      
      if (sum(rdiSp) != 0)
        rdiNormal = rdiNormal + rdiSp / sum(rdiSp)
    }
    r_maxN = species_params$r_max * rdiNormal
    
    for (i in 1:length(r_maxN))
      # do not want to divide by 0
      if (r_maxN[i] == 0)
        r_maxN[i] = species_params$r_max[i]
    
    return(t(t(rdi) * r_maxN) / t(t(rdi) + r_maxN) )
  }
  # setting up some parameters
  tic()
  sex_ratio = 0.5
  T = t_steps/dt; # number of time steps you want to follow cohort for
  

  
  no_Phen <- dim(object@n)[2] # number of phenotypes for matrix dim
  PhenName<- object@params@species_params$ecotype # names of phenotypes
  
  if(length(PhenName) != no_Phen) #conflict between biomass and species params
  {
    # getting out of hands as now the species df follows over multiple sims (bad idea?) | commented fix is here for non updated sims
    # keeping only the right portion of the df
    sp_params <- filter(object@params@species_params, extinct >= (simStart/dt) | extinct == F)
    # no_Phen <- dim(sp_params)[1] # number of phenotypes for matrix dim
    PhenName<- sp_params$ecotype # names of phenotypes
    # PhenIdx <- which(object@params@species_params$ecotype %in% PhenName) # index of surving phen in the biomass matrix and such
    # #updating it in the object (won't be saved)
    object@params@species_params <- sp_params
    # traitVec <- object@params@species_params[,traitID]
  }
  
  
  # PhenIdx <- seq(1,no_sp)
  
  # no_Phen = length(PhenIdx)
  fitness <- array(0,c(no_Phen, length(cohortSpan)), dimnames = list(PhenName,cohortSpan)) #collect the total spawn per time (start of cohort) per species
  names(dimnames(fitness)) <- list("species","cohort")
  for (t_start in cohortSpan)
  {
    cat(sprintf("Cohort number %g\n",t_start))
    
    # Initialising matrixes
    cohortW = array(0, c(no_Phen, T+1)); # row vector for following cohort weight
    cohortS = array(0, c(no_Phen, T+1)); # vector for cohort survival
    cohortR = array(0, c(no_Phen, T+1)); # vector for cohort spawning
    cohortR_sol = array(0, c(no_Phen, T+1)); # vector for cohort spawn at size
    cohortW[,1] = object@params@w[1]; # log weight initially (newborn)
    cohortS[,1] = object@n[t_start,,1]; # initial population in spectrum
    cohortInitPop <- object@n[t_start,,1] * object@params@dw[1] # initial number of individuals per species
    cohortInitPop[which(cohortInitPop == 0)] <- 1 # the phenotypes with 0 abundance get 1 so dividing by 0 doesn't end the universe
    
    for (t in seq(1,T)){ # within time period you're interested in
      # vector of the previous size bin for every phenotypes
      cohortWprev = unlist(lapply(lapply(cohortW[,t], FUN = function(x) x-object@params@w), FUN = function(x) max(which(x>= 0)))) # yolo
      # growth matrix
      growth = getEGrowth(object@params,n = object@n[t_start+t-1,,],n_pp = object@n_pp[t_start+t-1,],n_aa = object@n_aa[t_start+t-1,],n_bb = object@n_bb[t_start+t-1,], 
                          intakeScalar = object@intTempScalar[,,t_start+t-1], metScalar = object@metTempScalar[,,t_start+t-1])
      # update the new size bin with the growth / growth is how much more mass you get in the size bin so you just ad both to get the next size bin
      cohortW[,t+1] = cohortW[,t]+dt*diag(growth[,cohortWprev])
      # mortality matrix
      z = getZ(object = object@params, n = object@n[t_start+t-1,,],n_pp = object@n_pp[t_start+t-1,],n_aa = object@n_aa[t_start+t-1,],n_bb = object@n_bb[t_start+t-1,], 
               effort = effort,intakeScalar = object@intTempScalar[,,t_start+t-1], metScalar = object@metTempScalar[,,t_start+t-1], morScalar = object@morTempScalar[,,t_start+t-1] )
      # update the amount surviving the time-step
      cohortS[,t+1] = cohortS[,t]*exp(-dt*diag(z[,cohortWprev])) # where does the exponential comes from? you cannot just substract as Z is just an instantaneous mortality value?
      # need to take global n to have the right amount of resources available but need to take the right fraction at the end for the fitness, 
      # as not all the individuals reproducing are part of the cohort.
      # need to prepare n for no NAN, I just want the n of the specific cohort so I extract the right fraction
      n = object@n[t_start+t-1,,]
      # get the rdi manually to have it spread over size bins
      e_spawning <- getESpawning(object = object@params, n = n,n_pp = object@n_pp[t_start+t-1,],n_aa = object@n_aa[t_start+t-1,],n_bb = object@n_bb[t_start+t-1,], 
                                 intakeScalar = object@intTempScalar[,,t_start+t-1], metScalar = object@metTempScalar[,,t_start+t-1])
      e_spawning_pop <- apply((e_spawning*n),1,"*",object@params@dw)
      rdi <- sex_ratio*(e_spawning_pop * object@params@species_params$erepro)/object@params@w[object@params@w_min_idx] # global rdi
      rdi <- srrMatrix(rdi = rdi, species_params = object@params@species_params) # what happens with rmax?
      # get the proportion of abundance from the followed cohort. It's not t+1 as repro and death happens at the same time so we take the value from previous time step
      cohortF <- cohortS[,t]/cohortS[,1]
      cohortF[!is.finite(cohortF)] <- 0
      # update the total spawn for fitness
      cohortR[,t+1] = cohortR[,t] + dt*diag(rdi[cohortWprev,])*cohortF/cohortInitPop
      #cohortR_sol[q,t+1] = dt*rdi[cohortWprev,q] # do not sum the spawn so it is the spawn at time
      #print(cohortR[PhenIdx,1:t+1])
    }
    fitness[,which(t_start==cohortSpan)] = cohortR[,T] # fitness is the total spawn within the time period
  }
  
  rownames(fitness) <- PhenName # this is their name
  # make it a dataframe and add species and mat size for processing later
  # print(fitness)
  fitness <- as.data.frame(fitness)
  fitness$trait <- object@params@species_params[,traitID]
  fitness$species <- object@params@species_params$species
  
  fitness <- fitness[!rowSums(fitness[,-c(dim(fitness)[2]-1,dim(fitness)[2])]) == 0,] # get rid of phenotypes not appeared yet
  
  fitness_dat <- fitness # fitness is the whole data, fitness dat is the species we want to plot here
  
  if (!is.null(iSpecies)) fitness_dat <- fitness_dat[which(fitness_dat$species == iSpecies),] # select the right species if asked
  fitness_dat$species <- NULL # don't need this anymore
  fitness_dat[fitness_dat==0] <- NA # clearer plot
  
  plot_dat <- melt(fitness_dat, id = "trait")
  plot_dat$variable <- as.numeric(as.character(plot_dat$variable))
  
  colfunc <- colorRampPalette(c("black", "orange"))
  colGrad <- colfunc(length(unique(plot_dat$trait)))
  
  p <- ggplot(plot_dat) +
    geom_line(aes(x=as.numeric(variable),y=value,group = as.factor(trait), color = as.factor(trait))) +
    scale_y_continuous(trans = "log10") +
    scale_x_continuous(name = "Time") +
    scale_color_manual(name = "trait", values = colGrad)+
    #geom_vline(xintercept = 3000, linetype = "dashed") +
    theme(legend.title=element_text(),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),
          legend.justification=c(1,1),
          legend.position = "none",
          legend.key = element_rect(fill = "white"))+
    ggtitle("orange > black")
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  toc()
  if (returnData) return(fitness) else if(print_it) return(p)
}

# require plotCohort Output
plotFitness <- function(object, save_it = F, print_it = T, returnData = F, SpIdx = 1:9, ylim = c(NA,NA), iTime = NA, plotName = "fitnessTemp.png", returnPlot = T)
{
  ### data processing ###
  no_sp = length(unique(object$species))
  plotSpList <- vector("list",no_sp)
  
  counterSp = 0
  for(iSpecies in SpIdx) #c(1,6,9)) #
  {
    counterSp = counterSp +1
    myMat <- as.data.frame(object[which(object$species == iSpecies),])
    
    # then process
    myMat$species <- NULL
    myMat$group <- NULL
    myMat[myMat==0] <- NA
    
    # change slightly the w_mat of the original species so they do not overlap across sim
    #myMat$w_mat[myMat$w_mat == myMat$w_mat[1]] <- sapply(myMat$w_mat[myMat$w_mat == myMat$w_mat[1]],function(x) x*rnorm(1,1,0.00001))
    
    plot_dat <- melt(myMat, id = "trait")
    plot_dat$variable <- as.numeric(as.character(plot_dat$variable))
    
    colnames(plot_dat) <- c("trait","Time","Fitness")
    plot_dat$Species <- iSpecies
    
    # if wish to keep one time only
    # whichTime <- unique(plot_dat$Time)#[c(1,3,5,6,9,10,11)]
    # plotTimeList <- vector("list",length(whichTime))
    # iTime = 5500
    # plot_dat_time <- plot_dat[which(plot_dat$Time == iTime),]
    
    plotSpList[[iSpecies]] <- plot_dat
  }
  
  plotDat <- rbind(data.table::rbindlist(plotSpList)) #,data.table::rbindlist(plotSpListF),data.table::rbindlist(plotSpListN2),data.table::rbindlist(plotSpListF2))
  
  # plotDat <- plotDat[-which(is.na(plotDat$Fitness)),]
  
  
  # for (iSpecies in )
  # {
  #   spName <- paste("Species",iSpecies)
  #   plotList[[iSpecies]] <- ggplot(plotDat[(which(plotDat$Species == iSpecies)),]) +
  #     geom_point(aes(x=w_mat,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #     scale_y_continuous(trans = "log10", name = spName) +
  #     scale_x_continuous(name = NULL) +
  #     scale_color_gradient(name = "time", low = "black", high = "orange")+
  #     # facet_grid(.~scenario, scales = "free") +
  #     theme(legend.title=element_text(),
  #           panel.background = element_rect(fill = "white", color = "black"),
  #           # panel.grid.minor = element_line(colour = "grey92"),
  #           # legend.justification=c(1,1),
  #           legend.position = "bottom",
  #           strip.background = element_blank(),
  #           strip.text.x = element_text(c("no predation","predation")),
  #           legend.key = element_rect(fill = "white"))+
  #     ggtitle(NULL)
  #   
  #   
  # }
  plotList <- vector("list",no_sp)
  
  if(!is.na(iTime)) plotDat <- plotDat[(which(plotDat$Time == iTime)),]
  
  for(iSpecies in 1:no_sp) # improve that later
  {
    
    plotList[[iSpecies]] <- ggplot(plotDat[(which(plotDat$Species == iSpecies)),]) +
      geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
      # geom_line(aes(x=trait,y=Fitness, color = Time, group = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
      scale_y_continuous(trans = "log10", name = paste("Species",iSpecies), limits = ylim) +
      scale_x_continuous(name = NULL) +
      scale_color_gradient(name = "time", low = "black", high = "orange")+
      # facet_grid(.~scenario, scales = "free") +
      theme(legend.title=element_text(),
            panel.background = element_rect(fill = "white", color = "black"),
            # panel.grid.minor = element_line(colour = "grey92"),
            # legend.justification=c(1,1),
            legend.position = "bottom",
            strip.background = element_blank(),
            strip.text.x = element_text(c("no predation","predation")),
            legend.key = element_rect(fill = "white"))+
      ggtitle(NULL)
    
  }
  
  if (returnPlot) return(plotList)
  # p2 <- ggplot(plotDat[(which(plotDat$Species == 2)),]) +
  #   geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   scale_y_continuous(trans = "log10", name = "Species 2") +
  #   scale_x_continuous(name = NULL) +
  #   scale_color_gradient(name = "time", low = "black", high = "orange")+
  #   # facet_grid(.~scenario, scales = "free") +
  #   theme(legend.title=element_text(),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         # panel.grid.minor = element_line(colour = "grey92"),
  #         # legend.justification=c(1,1),
  #         legend.position = "bottom",
  #         strip.background = element_blank(),
  #         strip.text.x = element_text(c("no predation","predation")),
  #         legend.key = element_rect(fill = "white"))+
  #   ggtitle(NULL)
  # 
  # 
  # p3 <- ggplot(plotDat[(which(plotDat$Species == 3)),]) +
  #   geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   scale_y_continuous(trans = "log10", name = "Species 3") +
  #   scale_x_continuous(name = NULL) +
  #   scale_color_gradient(name = "time", low = "black", high = "orange")+
  #   # facet_grid(.~scenario, scales = "free") +
  #   theme(legend.title=element_text(),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         # panel.grid.minor = element_line(colour = "grey92"),
  #         # legend.justification=c(1,1),
  #         legend.position = "bottom",
  #         strip.background = element_blank(),
  #         strip.text.x = element_text(c("no predation","predation")),
  #         legend.key = element_rect(fill = "white"))+
  #   ggtitle(NULL)
  # 
  # p4 <- ggplot(plotDat[(which(plotDat$Species == 4)),]) +
  #   geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   scale_y_continuous(trans = "log10", name = "Species 4") +
  #   scale_x_continuous(name = NULL) +
  #   scale_color_gradient(name = "time", low = "black", high = "orange")+
  #   # facet_grid(.~scenario, scales = "free") +
  #   theme(legend.title=element_text(),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         # panel.grid.minor = element_line(colour = "grey92"),
  #         # legend.justification=c(1,1),
  #         legend.position = "bottom",
  #         strip.background = element_blank(),
  #         strip.text.x = element_text(c("no predation","predation")),
  #         legend.key = element_rect(fill = "white"))+
  #   ggtitle(NULL)
  # 
  # p5 <- ggplot(plotDat[(which(plotDat$Species == 5)),]) +
  #   geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   geom_line(aes(x=trait,y=Fitness, color = Time, group = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   
  #   scale_y_continuous(trans = "log10", name = "Species 5") +
  #   scale_x_continuous(name = NULL) +
  #   scale_color_gradient(name = "time", low = "black", high = "orange")+
  #   # facet_grid(.~scenario, scales = "free") +
  #   theme(legend.title=element_text(),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         # panel.grid.minor = element_line(colour = "grey92"),
  #         # legend.justification=c(1,1),
  #         legend.position = "bottom",
  #         strip.background = element_blank(),
  #         strip.text.x = element_text(c("no predation","predation")),
  #         legend.key = element_rect(fill = "white"))+
  #   ggtitle(NULL)
  # 
  # p6 <- ggplot(plotDat[(which(plotDat$Species == 6)),]) +
  #   geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   scale_y_continuous(trans = "log10", name = "Species 6") +
  #   scale_x_continuous(name = NULL) +
  #   scale_color_gradient(name = "time", low = "black", high = "orange")+
  #   # facet_grid(.~scenario, scales = "free") +
  #   theme(legend.title=element_text(),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         # panel.grid.minor = element_line(colour = "grey92"),
  #         # legend.justification=c(1,1),
  #         legend.position = "bottom",
  #         strip.background = element_blank(),
  #         strip.text.x = element_text(c("no predation","predation")),
  #         legend.key = element_rect(fill = "white"))+
  #   ggtitle(NULL)
  # 
  # p7 <- ggplot(plotDat[(which(plotDat$Species == 7)),]) +
  #   geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   scale_y_continuous(trans = "log10", name = "Species 7") +
  #   scale_x_continuous(name = NULL) +
  #   scale_color_gradient(name = "time", low = "black", high = "orange")+
  #   # facet_grid(.~scenario, scales = "free") +
  #   theme(legend.title=element_text(),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         # panel.grid.minor = element_line(colour = "grey92"),
  #         # legend.justification=c(1,1),
  #         legend.position = "bottom",
  #         strip.background = element_blank(),
  #         strip.text.x = element_text(c("no predation","predation")),
  #         legend.key = element_rect(fill = "white"))+
  #   ggtitle(NULL)
  # 
  # p8 <- ggplot(plotDat[(which(plotDat$Species == 8)),]) +
  #   geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   scale_y_continuous(trans = "log10", name = "Species 8") +
  #   scale_x_continuous(name = NULL) +
  #   scale_color_gradient(name = "time", low = "black", high = "orange")+
  #   # facet_grid(.~scenario, scales = "free") +
  #   theme(legend.title=element_text(),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         # panel.grid.minor = element_line(colour = "grey92"),
  #         # legend.justification=c(1,1),
  #         legend.position = "bottom",
  #         strip.background = element_blank(),
  #         strip.text.x = element_text(c("no predation","predation")),
  #         legend.key = element_rect(fill = "white"))+
  #   ggtitle(NULL)
  # 
  # p9 <- ggplot(plotDat[(which(plotDat$Species == 9)),]) +
  #   geom_point(aes(x=trait,y=Fitness, color = Time), size = 0.5) +#,group = Time, color = Time), size = 1) +
  #   scale_y_continuous(trans = "log10", name = "Species 9") +
  #   scale_x_continuous(name = NULL) +
  #   scale_color_gradient(name = "time", low = "black", high = "orange")+
  #   # facet_grid(.~scenario, scales = "free") +
  #   theme(legend.title=element_text(),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         # panel.grid.minor = element_line(colour = "grey92"),
  #         # legend.justification=c(1,1),
  #         legend.position = "bottom",
  #         strip.background = element_blank(),
  #         strip.text.x = element_text(c("no predation","predation")),
  #         legend.key = element_rect(fill = "white"))+
  #   ggtitle(NULL)
  
  
  g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}
  
  mylegend<-g_legend(plotList[[1]])
  
  p10 <- grid.arrange(arrangeGrob(plotList[[1]] + theme(legend.position="none"),
                                  plotList[[2]] + theme(legend.position="none"),
                                  plotList[[3]] + theme(legend.position="none"),
                                  plotList[[4]] + theme(legend.position="none"),
                                  plotList[[5]] + theme(legend.position="none"),
                                  plotList[[6]] + theme(legend.position="none"),
                                  plotList[[7]] + theme(legend.position="none"),
                                  plotList[[8]] + theme(legend.position="none"),
                                  plotList[[9]] + theme(legend.position="none"),
                                  nrow=3),
                      mylegend, nrow=2,heights=c(9.5,.5),
                      left=textGrob("Fitness", rot = 90, vjust = 1))
  
  
  if(save_it) ggsave(p10,filename = plotName,device = "png",width = 21,height = 21,units = "cm")
  
  if(returnData) return(plotDat) else if(print_it) return(p10)
  
}

plotTemperature <- function(object, temperature = seq(1,30,.2), iSpecies = 1, size = NULL, ylim = c(NA,NA),
                            f0 = object@params@f0, netEnergy = F,
                            save_it = F, print_it = T, returnData = F, 
                            plotName = "thermotoleranceTemp.png", plotTitle = "thermotolerance")
{
  if(is.null(size)) size <- # size slot closer to maturation size
      which(abs(object@params@w - object@params@species_params$w_mat[iSpecies]) == 
              min(abs(object@params@w - object@params@species_params$w_mat[iSpecies])))
  
  
  intakeScalar <- tempFun(w = object@params@w, temperature = temperature, t_ref = object@params@t_ref, t_d =  object@params@species_params$t_d[iSpecies], 
                          Ea = object@params@species_params$ea_int[iSpecies], c_a = object@params@species_params$ca_int[iSpecies], 
                          Ed = object@params@species_params$ed_int[iSpecies], c_d = object@params@species_params$cd_int[iSpecies])
  
  metabScalar <- tempFun(w = object@params@w, temperature = temperature, t_ref = object@params@t_ref , t_d = object@params@species_params$t_d[iSpecies], 
                         Ea = object@params@species_params$ea_met[iSpecies], c_a = object@params@species_params$ca_met[iSpecies], 
                         Ed = object@params@species_params$ed_met[iSpecies], c_d = object@params@species_params$cd_met[iSpecies])
  # temperatureScalar has size in rows and temperature in cols
  
  # multiplication of 2 matrices to get a 3 dim array
  b <- NULL
  for (iSp in 1:dim(object@params@intake_max)[1])
  {
    a <- object@params@intake_max[iSp,] * intakeScalar
    b <- abind(b,a,along = 3)
  }
  
  e <- sweep(f0 * b, c(1,2), object@params@species_params$alpha, "*", check.margin = FALSE)
  
  b <- NULL
  for (iSp in 1:dim(object@params@intake_max)[1])
  {
    a <- object@params@metab[iSp,] * metabScalar
    b <- abind(b,a,along = 3)
  }
  
  e <- e - b
  
  # normalise between 0 and 2 for the plot
  netE_dat <- e[size,,iSpecies]
  netE_dat <- netE_dat/max(netE_dat)*2
  
  myData <- data.frame("temperature" = temperature, "intake" = intakeScalar[size,], "metabolism" = metabScalar[size,], "netEnergy" = netE_dat)
  
  # myData <- data.frame("temperature" = temperature-273, "scalar" = temperatureScalar[20,])
  t_ref <- object@params@species_params$t_d[iSpecies] - (5 + 0.25*object@params@species_params$t_d[iSpecies] )

  p1 <- ggplot(myData)+
    geom_line(aes(x = temperature, y = intake), color = "blue")+
    geom_line(aes(x = temperature, y = metabolism), color = "red")+
    geom_vline(xintercept = t_ref, linetype = "dashed", alpha = 0.5) +
    annotate("text", x = t_ref, y = 0.1, label = "t_ref") + 
    geom_vline(xintercept = object@params@species_params$t_d[iSpecies], linetype = "dashed", alpha = 0.5,color = "red") +
    annotate("text", x = object@params@species_params$t_d[iSpecies], y = 0.1, label = "t_d") + 
    geom_hline(yintercept = 0, alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    scale_y_continuous(name = "scalar", limits = ylim)+
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(plotTitle)
  
  if(netEnergy) p1 <- p1 + geom_line(aes(x = temperature, y = netEnergy), color = "green")
  
  if(save_it) ggsave(p1, file = plotName, scale = 1.5)
  
  if(returnData) return(myData) else if(print_it) return(p1)
  
}

plotTemperaturePhen <- function(object, temperature = seq(1,30), iSpecies = 1, ylim = c(NA,NA), size = NULL,
                                save_it = F, print_it = T, returnData = F, 
                                plotName = "thermotolerancePhen.png", plotTitle = "thermotolerance")
{
  
  phenParams <- object@params@species_params[which(sim@params@species_params$species == iSpecies),]
  
  
  
  metTempScalar <- matrix(NA, nrow = dim(phenParams)[1], ncol = length(temperature), dimnames = list(phenParams$ea_met,temperature)) 
  intTempScalar <- matrix(NA, nrow = dim(phenParams)[1], ncol = length(temperature), dimnames = list(phenParams$ed_int,temperature)) 
  
  if(is.null(size)) size <- # size slot closer to maturation size
    which(abs(object@params@w - object@params@species_params$w_mat[iSpecies]) == 
            min(abs(object@params@w - object@params@species_params$w_mat[iSpecies])))
  
  for(iPhen in 1:dim(phenParams)[1]) # fill the scalars arrays
  {
    
    
    
    intTempScalar[iPhen,] <- tempFun(w = size, temperature = temperature, t_ref = object@params@t_ref , t_d = phenParams$t_d[iPhen], 
                                     Ea = phenParams$ea_int[iPhen], c_a = phenParams$ca_int[iPhen], 
                                     Ed = phenParams$ed_int[iPhen], c_d = phenParams$cd_int[iPhen])
    
    metTempScalar[iPhen,] <- tempFun(w = size, temperature = temperature, t_ref = object@params@t_ref , t_d = phenParams$t_d[iPhen], 
                                     Ea = phenParams$ea_met[iPhen], c_a = phenParams$ca_met[iPhen], 
                                     Ed = phenParams$ed_met[iPhen], c_d = phenParams$cd_met[iPhen])
    # temperatureScalar has size in rows and temperature in cols
  }
  
  plotDat <- melt(intTempScalar)
  colnames(plotDat) <- c("phenotypes","temperature","intake")
  tempDat <- melt(metTempScalar)
  plotDat$metabolism <- tempDat$value
  
  
  
  p1 <- ggplot(plotDat)+
    geom_line(aes(x = temperature, y = intake, group = phenotypes, color = phenotypes)) +
    # geom_line(aes(x = temperature, y = metabolism), color = "red")+
    geom_vline(xintercept = object@params@t_ref, linetype = "dashed", alpha = 0.5) +
    annotate("text", x = object@params@t_ref, y = 0.1, label = "t_ref") +
    geom_vline(xintercept = object@params@species_params$t_d[iSpecies], linetype = "dashed", alpha = 0.5,color = "red") +
    annotate("text", x = object@params@species_params$t_d[iSpecies], y = 0.1, label = "t_d") +
    geom_vline(xintercept = mean(object@temperature), linetype = "dashed", alpha = 0.5) +
    annotate("text", x = mean(object@temperature), y = 0.3, label = "temperate avg") +
    geom_hline(yintercept = 0, alpha = 0.5) +
    labs(color='ed_int') +
    # geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    scale_y_continuous(name = "scalar", limits = ylim)+
    theme(legend.title= element_text(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(plotTitle)
  
  
  if(save_it) ggsave(p1, file = plotName, scale = 1.5)
  
  if(returnData) return(plotDat) else if(print_it) return(p1)
}

plotPlankton <- function(object,
                         start_time = as.numeric(dimnames(sim@n)[[1]][1]),
                         end_time = as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]]),
                         y_ticks = 6, print_it = TRUE,
                         ylim = c(NA, NA),
                         total = TRUE, returnData = FALSE, ...){
  b <- getPlanktonBiom(object)
  if (start_time >= end_time) {
    stop("start_time must be less than end_time")
  }
  # Select time range
  b <- b[(as.numeric(dimnames(b)[[1]]) >= start_time) &
           (as.numeric(dimnames(b)[[1]]) <= end_time), , drop = FALSE]
  b_total <- rowSums(b)
  # Include total
  # if (total) {
  b <- cbind(b, Total = b_total)
  
  bm <- data.frame("time" = names(b_total), "value" = b_total)
  bm$time <- as.numeric(as.character(bm$time))
  
  p <- ggplot(bm) +
    geom_line(aes(x = time, y = value, group = 1)) +
    scale_y_continuous(name = "Biomass [g]") + #, breaks = log_breaks(n = y_ticks),  trans = "log10") 
    theme(panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"),
          legend.key = element_rect(fill = "white"))+
    ggtitle("Plankton biomass") 
  
  # } else {
  # names(dimnames(b)) <- c("time", "size")
  # bm <- reshape2::melt(b)
  # # Force Species column to be a factor (otherwise if numeric labels are
  # # used they may be interpreted as integer and hence continuous)
  # # bm$Species <- as.factor(bm$Species)
  # # Implement ylim and a minimal cutoff
  # min_value <- 1e-20
  # bm <- bm[bm$value >= min_value &
  #            (is.na(ylim[1]) | bm$value >= ylim[1]) &
  #            (is.na(ylim[2]) | bm$value <= ylim[1]), ]
  # # Select species
  # # spec_bm <- bm[bm$Species %in% species, ]
  # x_label <- "Year"
  # y_label <- "Biomass [g]"
  # 
  # p <- ggplot(bm) +
  #   geom_line(aes(x = time, y = value, group = size, color = as.numeric(as.character(size)))) +
  # scale_y_continuous(trans = "log10", breaks = log_breaks(n = y_ticks), name = y_label) 
  # 
  # }
  
  if(returnData) return(p) else if (print_it) return(p)
  
}

plotNetEnergy <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, ylim = c(NA,NA),
                          print_it = T, returnData = F, save_it = F, nameSave = "netEnergy.png",...){
  
  time_elements <- get_time_elements(object,time_range)
  energy_time <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    energy <- getEReproAndGrowth(object@params, n=n, n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                                 intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
    return(energy)})
  
  #growth <- apply(growth_time, c(2,3), mean) # use this when I will have time_range on more than one time
  energy = energy_time
  
  if (species) # if I want to display species instead of ecotypes
  {
    dimnames(energy)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    energy_sp = matrix(data = NA, ncol = dim(energy)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(energy)$w)) # prepare the new object
    names(dimnames(energy_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = energy # save to manip
      temp[which(rownames(energy) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(energy)==i)) # do the mean (in 2 steps)
      energy_sp[which(rownames(energy_sp)==i), ] = temp
    }
    energy = energy_sp
  }
  
  name = paste("energy level at time",time_range,sep=" ")
  plot_dat <- data.frame(value = c(energy), Species = dimnames(energy)[[1]], w = rep(object@params@w, each=length(dimnames(energy)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = value, colour = Species)) + 
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
    scale_y_continuous(name = "net energy value", trans ="log10", limits = ylim)+
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}

plotAvailableEnergy <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), species = T, ylim = c(NA,NA),
                                print_it = T, returnData = F, save_it = F,
                                nameSave = "AvailableEnergy.png",...){
  
  SpIdx <- as.numeric(sort(unique(object@params@species_params$species)))
  colfunc <- colorRampPalette(c("firebrick3","darkturquoise", "orange"))
  colGrad <- colfunc(length(SpIdx))
  names(colGrad) <- SpIdx
  colLine <- c("solid","dashed")
  names(colLine) <- c("Fish","Background")
  
  
  time_elements <- get_time_elements(object,time_range)
  prey <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    res <- getAvailPrey(object@params, n=n)
    return(res)})
  
  background <- aaply(which(time_elements), 1, function(x){
    # Necessary as we only want single time step but may only have 1 species which makes using drop impossible
    
    n <- array(object@n[x,,],dim=dim(object@n)[2:3])
    dimnames(n) <- dimnames(object@n)[2:3]
    res <- getAvailBackground(object@params, n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,])
    return(res)})
  
  if (species) # if I want to display species instead of ecotypes
  {
    ## take care of prey
    dimnames(prey)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    prey_sp = matrix(data = NA, ncol = dim(prey)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(prey)$w)) # prepare the new object
    names(dimnames(prey_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = prey # save to manip
      temp[which(rownames(prey) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(prey)==i)) # do the mean (in 2 steps)
      prey_sp[which(rownames(prey_sp)==i), ] = temp
    }
    prey = prey_sp
    
    # take care of background
    dimnames(background)$sp = object@params@species_params$species
    SpIdx = sort(unique(object@params@species_params$species)) # get the species names
    background_sp = matrix(data = NA, ncol = dim(background)[2], nrow = length(SpIdx), dimnames = list(SpIdx,dimnames(background)$w)) # prepare the new object
    names(dimnames(background_sp))=list("species","size")
    
    for (i in SpIdx)
    {
      temp = background # save to manip
      temp[which(rownames(background) != i), ] = 0 # keep the ecotypes from the species only
      temp = apply(temp, 2, sum)
      temp = temp / length(which(rownames(background)==i)) # do the mean (in 2 steps)
      background_sp[which(rownames(background_sp)==i), ] = temp
    }
    background = background_sp
  }
  
  name = paste("Available energy at time",time_range,sep=" ")
  plot_dat <- data.frame(prey = c(prey), background = c(background), species = dimnames(prey)[[1]], w = rep(object@params@w, each=length(dimnames(prey)[[1]])))
  p <- ggplot(plot_dat) + 
    geom_line(aes(x=w, y = prey, colour = species, group = species, linetype = "Fish")) +
    geom_line(aes(x=w, y = background, colour = species, group = species, linetype = "Background")) +
    geom_line(aes(x=w, y = (prey + background), colour = species, group = species), size = 1, alpha = .5) +
    scale_x_continuous(name = "Size", trans="log10", breaks = c(1 %o% 10^(-3:5))) + 
    scale_y_continuous(name = "instantaneous value", trans ="log10", limits = ylim)+
    scale_color_manual(name = "Species", values = colGrad)+ # color gradient
    scale_linetype_manual(name = "Prey", values = colLine)+
    theme(legend.title=element_text(),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), 
          legend.position="right",
          legend.key = element_rect(fill = "white"))+
    guides(color = guide_legend(nrow=length(SpIdx)),
           linetype = guide_legend(order = 2,override.aes = list(colour = "black")))+
    ggtitle(name)
  
  if(save_it) ggsave(plot = p, filename = nameSave,scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}


ratioReproduction <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)), RDI = F){
  
  time_range = max(as.numeric(dimnames(object@n)$time))
  SpIdx = sort(unique(object@params@species_params$species)) # get the species names
  time_elements <- get_time_elements(object,time_range)
  x <- which(time_elements)
  
  rdi <- getRDI(object@params, n=object@n[x,,], n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
  
  rdd <- getRDD(object@params, n=object@n[x,,], n_pp = object@n_pp[x,],n_aa = object@n_aa[x,],n_bb = object@n_bb[x,], 
                intakeScalar = object@intTempScalar[,,x], metScalar = object@metTempScalar[,,x])
  
  
  if(length(object@params@species_params$ecotype) > length(SpIdx))
  {
    
    names(rdi) <- object@params@species_params$species
    rdiSpecies = NULL
    for (iSpecies in SpIdx)
    {
      temp = rdi # save to manip
      temp[which(names(rdi) != iSpecies)] = 0 # keep the ecotypes from the species only
      temp = sum(temp)
      rdiSpecies <- cbind(rdiSpecies,temp)
    }
    rdi <- rdiSpecies
    colnames(rdi) <- SpIdx
    
    rddSpecies = NULL
    for (iSpecies in SpIdx)
    {
      temp = rdd # save to manip
      temp[which(names(rdd) != iSpecies)] = 0 # keep the ecotypes from the species only
      temp = sum(temp)
      rddSpecies <- cbind(rddSpecies,temp)
    }
    rdd <- rddSpecies
    colnames(rdd) <- SpIdx
    
  }
  
  ratio <- rdd/rdi
  
  
  if(RDI) return(rdi) else return(ratio)
  
}

#' This function is used to look at the performace curve of species.
#' It takes a basic sim for default but you can specify any parameters you want to reshape the curve (only intake available at the moment)
plotTemperatureProfile <- function(object, temperature = 1:30, iSpecies = 5, size = NULL, f0 = object@params@f0,
                                   t_d = NULL, t_ref = NULL, ed_int = NULL, ea_int = NULL, ca_int = NULL, cd_int = NULL,
                                   save_it = F, print_it = T, returnData = F, ylim = c(NA,NA), energyProfile = T)
{
  
  if(is.null(size)) size <- # size slot closer to maturation size
      which(abs(object@params@w - object@params@species_params$w_mat[iSpecies]) == 
              min(abs(object@params@w - object@params@species_params$w_mat[iSpecies])))
  
  
  if(is.null(t_d)) t_d <- object@params@species_params$t_d[iSpecies]
  if(is.null(t_ref)) t_ref <-   t_ref <- t_d - (5 + 0.25*t_d) #object@params@t_ref
  
  if(is.null(ed_int)) ed_int <- object@params@species_params$ed_int[iSpecies]
  if(is.null(ea_int)) ea_int <- object@params@species_params$ea_int[iSpecies]
  if(is.null(cd_int)) cd_int <- object@params@species_params$cd_int[iSpecies]
  if(is.null(ca_int)) ca_int <- object@params@species_params$ca_int[iSpecies]
  
  intakeScalar <- tempFun(w = size, temperature = temperature, t_ref = t_ref , t_d = t_d , 
                          Ea = ea_int, c_a = ca_int, 
                          Ed = ed_int, c_d = cd_int)
  
  metabScalar <- tempFun(w = size, temperature = temperature, t_ref = t_ref , t_d = t_d,
                         Ea = object@params@species_params$ea_met[iSpecies], c_a = object@params@species_params$ca_met[iSpecies],
                         Ed = object@params@species_params$ed_met[iSpecies], c_d = object@params@species_params$cd_met[iSpecies])
  
  
  
  a <- object@params@intake_max[iSpecies,size] * intakeScalar
  e <- object@params@f0 * a * object@params@species_params$alpha[iSpecies]
  b <- object@params@metab[iSpecies,size] * metabScalar       
  e <- e - b
  
  if(energyProfile) 
  {
    myData <- data.frame("temperature" = temperature, "energy" = t(e)) 
    ylabel <- "energy"
  }  else {
    myData <- data.frame("temperature" = temperature, "energy" = t(intakeScalar))
    ylabel <- "scalar"
  }
  
  p <- ggplot(myData)+
    geom_line(aes(x = temperature, y = energy))+
    geom_vline(xintercept = t_ref, linetype = "dashed", alpha = 0.5) +
    annotate("text", x = t_ref, y = 0.1, label = "t_ref") +
    geom_vline(xintercept = t_d, linetype = "dashed", alpha = 0.5,color = "red") +
    annotate("text", x = t_d, y = 0.1, label = "t_d") +
    scale_y_continuous(name = ylabel, limits = ylim)+
    theme(legend.title=element_blank(),
          legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(NULL)
  
  if(returnData) return(myData) else return(p)
}

plotEnergyperTemperature <- function(object,trait,ylim = c(NA,NA), iSpecies = 5, temperature = 1:30, size = NULL, traitValue = NULL,
                                     td = NULL, ed = NULL, ea = NULL, tref = NULL, save_it = F, print_it = T, returnData = F, 
                                     plotName = "traitThermoToleranceTemp.png", plotTitle = "trait thermotolerance", plotEnergy = T)
{
# use this temperature function where you need to give everything
  tempFunTest <- function(w, temperature, t_d, t_ref, Ea, c_a, Ed, c_d) # default are 0 for now as deactivation is buggy
  {
    k = 8.617332e-5 # Boltzmann constant 
    # converting to Kelvin from Celcius
    temperature <- temperature + 273 
    t_ref <- t_ref + 273
    t_d <- t_d + 273
  
    temperatureScalar <- t(sapply(w,FUN = function(x){x^(c_a*(temperature - t_ref))}) * exp(-Ea/k*(1/temperature - 1/t_ref)   )) *
      t(1/(sapply(w,FUN = function(x){x^(c_d*(temperature - t_d  ))}) *(exp(-Ed/k*(1/temperature - 1/t_d)) + 1)))*
      sapply(w,FUN = function(x){x^(c_d*(t_ref       - t_d  ))}) *(exp(-Ed/k*(1/t_ref       - 1/t_d)) + 1) 
    
    return(temperatureScalar)
  }
  
  
  if(is.null(td)) t_d <- object@params@species_params$t_d[iSpecies]
  if(is.null(ed)) ed_int <- object@params@species_params$ed_int[iSpecies]
  
  
  if(is.null(size))  size <- # size slot closer to maturation size
      which(abs(object@params@w - object@params@species_params$w_mat[iSpecies]) == 
              min(abs(object@params@w - object@params@species_params$w_mat[iSpecies])))
  
  switch(trait,
         "ed_int" = {
           if (is.null(traitValue)) ed_int <- seq(round(min(object@params@species_params$ed_int),1),round(max(object@params@species_params$ed_int),1),by = 0.1)
           else ed_int <- traitValue
           netE_dat <- matrix(NA,nrow = length(temperature) ,ncol = length(ed_int),dimnames = list("temperature" = temperature,"trait"= ed_int))
           for(var in ed_int)
           {
             if(is.null(ea)) ea_int <- object@params@species_params$ea_int[iSpecies]
             else ea_int <- ea(var)
             # else ea_int <- ea
             if(is.null(tref)) t_ref <- object@params@t_ref
             else t_ref <- tref(t_d)
             intakeScalar <- tempFunTest(w = size, temperature = temperature, t_ref = t_ref , t_d = t_d , 
                                     Ea = ea_int, c_a = object@params@species_params$ca_int[iSpecies], 
                                     Ed = var, c_d = object@params@species_params$cd_int[iSpecies])
             
             metabScalar <- tempFunTest(w = size, temperature = temperature, t_ref = t_ref , t_d = t_d,
                                    Ea = ea_int, c_a = object@params@species_params$ca_met[iSpecies],
                                    Ed = var, c_d = object@params@species_params$cd_met[iSpecies])
             
             
             
             a <- object@params@intake_max[iSpecies,size] * intakeScalar
             e <- object@params@f0 * a * object@params@species_params$alpha[iSpecies]
             b <- object@params@metab[iSpecies,size] * metabScalar       
             e <- e - b
             
             if(plotEnergy) netE_dat[,which(colnames(netE_dat) == var)] <- e
             else netE_dat[,which(colnames(netE_dat) == var)] <- intakeScalar
           }
         },
         "ea_int" = {
           if (is.null(traitValue)) ea_int <- seq(round(min(object@params@species_params$ea_int),1),round(max(object@params@species_params$ea_int),1),by = 0.1)
           else ea_int <- traitValue
           netE_dat <- matrix(NA,nrow = length(temperature) ,ncol = length(ea_int),dimnames = list("temperature" = temperature,"trait"= ea_int))
           for(var in ea_int)
           {
             if(is.null(ed)) ed_int <- object@params@species_params$ed_int[iSpecies]
             else ea_int <- ea(var)
             
             if(is.null(tref)) t_ref <- object@params@t_ref
             else t_ref <- tref(t_d)
             
             intakeScalar <- tempFunTest(w = size, temperature = temperature, t_ref = t_ref , t_d = t_d , 
                                     Ea = var, c_a = object@params@species_params$ca_int[iSpecies], 
                                     Ed = ed_int, c_d = object@params@species_params$cd_int[iSpecies])
             
             metabScalar <- tempFunTest(w = size, temperature = temperature, t_ref = t_ref , t_d = t_d,
                                    Ea = var, c_a = object@params@species_params$ca_met[iSpecies],
                                    Ed = ed_int, c_d = object@params@species_params$cd_met[iSpecies])
             
             
             
             a <- object@params@intake_max[iSpecies,size] * intakeScalar
             e <- object@params@f0 * a * object@params@species_params$alpha[iSpecies]
             b <- object@params@metab[iSpecies,size] * metabScalar       
             e <- e - b
             
             if(plotEnergy) netE_dat[,which(colnames(netE_dat) == var)] <- e
             else netE_dat[,which(colnames(netE_dat) == var)] <- intakeScalar
           }
         },
         "t_d" = {
           if (is.null(traitValue)) t_d <- seq(round(min(object@params@species_params$t_d),1),round(max(object@params@species_params$t_d),1),by = 0.2) 
           else t_d <- traitValue
           
           netE_dat <- matrix(NA,nrow = length(temperature) ,ncol = length(t_d),dimnames = list("temperature" = temperature,"trait"= t_d))
           
           for(var in t_d)
           {
             if(is.null(ea)) ea_int <- object@params@species_params$ea_int[iSpecies]
             else ea_int <- ea(ed_int)
             
             if(is.null(tref)) t_ref <- object@params@t_ref
             else t_ref <- tref(var)
             
             intakeScalar <- tempFunTest(w = size, temperature = temperature, t_ref = t_ref , t_d = var , 
                                     Ea = ea_int, c_a = object@params@species_params$ca_int[iSpecies], 
                                     Ed = ed_int, c_d = object@params@species_params$cd_int[iSpecies])
             
             metabScalar <- tempFunTest(w = size, temperature = temperature, t_ref = t_ref , t_d = var,
                                    Ea = ea_int, c_a = object@params@species_params$ca_met[iSpecies],
                                    Ed = ed_int, c_d = object@params@species_params$cd_met[iSpecies])
             
             
             
             a <- object@params@intake_max[iSpecies,size] * intakeScalar
             e <- object@params@f0 * a * object@params@species_params$alpha[iSpecies]
             b <- object@params@metab[iSpecies,size] * metabScalar       
             e <- e - b
             
             if(plotEnergy) netE_dat[,which(colnames(netE_dat) == var)] <- e
             else netE_dat[,which(colnames(netE_dat) == var)] <- intakeScalar
             
           }
         },
         {print("Trait specified is not in the list")}
  )

  
  plot_dat <- melt(netE_dat)
  
  p <- ggplot(plot_dat)+
    geom_line(aes(x = temperature, y = value, group = trait, color = trait)) +
    # geom_vline(xintercept = object@params@t_ref, linetype = "dashed", alpha = 0.5) +
    # annotate("text", x = object@params@t_ref, y = 0.1, label = "t_ref") +
    # geom_vline(xintercept = t_d, linetype = "dashed", alpha = 0.5) +
    # annotate("text", x = t_d, y = 0.1, label = "t_d") +
    scale_y_continuous(name = "net energy", limits = ylim)+
    theme(#legend.title=element_blank(),
      legend.justification=c(1,1),
      legend.key = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.minor = element_line(colour = "grey92"))+
    ggtitle(plotTitle)

  if(save_it) ggsave(p, file = plotName, scale = 1.5)
  
  if(returnData) return(plot_dat) else if(print_it) return(p)
}

plotNoPhen <- function(object, initPop = 1, comments = T, print_it = T, returnData = F, SpIdx = NULL, dt = 0.1)
{
  timeMax <- object@params@species_params$timeMax[1]
  if(is.null(SpIdx)) SpIdx <- 1:9
  exit_df <- data.frame()
  
  SumPar = data.frame(object@params@species_params$species,object@params@species_params$pop,object@params@species_params$extinct)
  
  colnames(SumPar) = c("species","pop","exit")
  
  SumPar$pop[which(SumPar$pop < (initPop/dt))] = initPop/dt # initial species pop at 1 by default, or defined by user if there was an initialisation
  
  for (i in 1:dim(SumPar)[1]) if (SumPar$exit[i] == 0) SumPar$exit[i] = timeMax # the not extinct ones get the end sim as extinction value
  
  DemList <- list()
  
  for (x in SpIdx) # for every species
  {
    exit <- NA
    DemCount = matrix(data = 0, nrow = (timeMax*dt), ncol =2) # create a matrix that count pop and extinction for each time step *dt
    DemCount[1,1] = 1 # fill the first row
    SpSumPar <- SumPar[which(SumPar$species == x),] # select the right species df
    
    if(dim(SpSumPar)[1] != 1) # if more than one phenotype
    {
      for (j in 2:dim(SpSumPar)[1]) # along the df
        DemCount[ceiling(SpSumPar$pop[j]*dt),1] = DemCount[ceiling(SpSumPar$pop[j-1]*dt),1] + 1 # total number of phenotypes at that time
      for (j in 1:dim(SpSumPar)[1]) # along the df 
        if (SpSumPar$exit[j] != timeMax) # do not take into account extinction at the last step (because its not)
          DemCount[ceiling(SpSumPar$exit[j]*dt),2] = DemCount[ceiling(SpSumPar$exit[j]*dt),2] -1 # an extinction happened at that time
    } else { DemCount[ceiling(SpSumPar$exit[1]*dt),2] = DemCount[ceiling(SpSumPar$exit[1]*dt),2] -1 } # if only one phenotype in sp
    
    # I have a matrix with holes, need to fill them
    ExCount = 0 
    for (i in 2:dim(DemCount)[1])
    {
      if (DemCount[i,1] == 0) DemCount[i,1] = DemCount[i-1,1]
      if (DemCount[i,2] != 0)
      {
        DemCount[i,2] = ExCount + abs(DemCount[i,2])
        ExCount = abs(DemCount[i,2]) } else  {
          DemCount[i,2] =  ExCount 
        }
    }
    
    pop <- DemCount[,1] - DemCount[,2] # pop - extinction
    DemCount <- as.data.frame(pop) # keep the alive number of phen
    colnames(DemCount) <- c("pop")
    DemCount$time <- as.numeric(rownames(DemCount))
    DemCount$species <- x
    DemCount$pop[DemCount$pop == 0] <- NA # if species goes extinct, put some NA (do not count in the mean)
    if (sum(is.na(DemCount$pop))) exit <- which(is.na(DemCount$pop))[1] # remember when the species go extinct
    DemList[[x]] <- DemCount
    
    if(is.finite(exit)) exit_df<- rbind(exit_df,c(x,exit))
    
  }
  myData <- do.call(rbind,lapply(DemList,function(x)x))
  # do an average across simulation
  # plot_dat <- data.frame(dataList[[1]][[1]]$time,dataList[[1]][[1]]$species,
  #                        apply(do.call(cbind,lapply(dataList, function(x) x[[1]]$pop)),1,mean, na.rm = T),
  #                        apply(do.call(cbind,lapply(dataList, function(x) x[[2]]$pop)),1,mean, na.rm = T))
  # colnames(plot_dat) <- c("time","species","popN","popF")                      

  # PLOTTING TIME
# SpIdx <- 1:9
colfunc <- colorRampPalette(c("firebrick3","darkturquoise", "orange"))
colGrad <- colfunc(length(SpIdx))
names(colGrad) <- SpIdx

p <- ggplot(myData) +
  stat_smooth(aes(x = time, y = pop, color = as.factor(species)), method = "loess", span = 0.15, se = F, size = 0.5)+
  #geom_point(data = exit_df, aes(x=extinction,y=0,color=as.factor(species))) +
  scale_x_continuous(name = "Time (yr)") +
  scale_y_continuous(name = "Number of phenotypes", limits = c(0,NA))+
  scale_color_manual(name = "Species", values = colGrad)+ # color gradient
  theme(legend.title=element_blank(),
        # legend.position=c(0.5,0.95),
        legend.justification=c(1,1),
        legend.box = "horizontal",
        legend.key = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", color = "black"),
        panel.grid.minor = element_line(colour = "grey92"))+
  # guides(color=guide_legend(override.aes=list(fill=NA), order =1),
  # linetype = guide_legend(order = 2,override.aes = list(colour = "black")))+
  ggtitle("Variation of phenotype's number throughout the simulation")

  # colfunc <- colorRampPalette(c("firebrick3","darkturquoise", "orange"))
  # colGrad <- colfunc(length(SpIdx))
  # names(colGrad) <- SpIdx
  
  # p <- ggplot(myData) +
  #   stat_smooth(aes(x = time, y = pop, color = as.factor(species)), method = "loess", span = 0.15, se = F, size = 0.5)+
  #   #geom_point(data = exit_df, aes(x=extinction,y=0,color=as.factor(species))) +
  #   scale_x_continuous(name = "Time (yr)") +
  #   scale_y_continuous(name = "Number of phenotypes")+
  #   scale_color_manual(name = "Species", values = colGrad)+ # color gradient
  #   scale_linetype_manual(name = "Fisheries", values = colLine)+
  #   theme(legend.title=element_blank(),
  #         # legend.position=c(0.5,0.95),
  #         legend.justification=c(1,1),
  #         legend.box = "horizontal",
  #         legend.key = element_rect(fill = "white"),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         panel.grid.minor = element_line(colour = "grey92"))+
  #   # guides(color=guide_legend(override.aes=list(fill=NA), order =1),
  #          # linetype = guide_legend(order = 2,override.aes = list(colour = "black")))+
  #   ggtitle("Variation of phenotype's number throughout the simulation")
  
  # if(length(exit_df) !=0)
  # {
  #   colnames(exit_df) <- c("species","extinction")
  #   p <- ggplot(plot_dat) +
  #     stat_smooth(aes(x = time, y = popN, color = as.factor(species), linetype = "un-fished"), method = "loess", span = 0.15, se = F, size = 0.5)+
  #     stat_smooth(aes(x = time, y = popF,color = as.factor(species), linetype = "fished"), method = "loess", span = 0.15, se = F, size = 0.5)+
  #     geom_point(data = exit_df, aes(x=extinction,y=0,color=as.factor(species))) +
  #     geom_vline(xintercept = 4000, linetype = "dashed") +
  #     scale_x_continuous(name = "Time (yr)") +
  #     scale_y_continuous(name = "Number of phenotypes")+
  #     scale_color_manual(name = "Species", values = colGrad)+ # color gradient
  #     scale_linetype_manual(name = "Fisheries", values = colLine)+
  #     theme(legend.title=element_blank(),
  #           legend.position=c(0.5,0.95),
  #           legend.justification=c(1,1),
  #           legend.box = "horizontal",
  #           legend.key = element_rect(fill = "white"),
  #           panel.background = element_rect(fill = "white", color = "black"),
  #           panel.grid.minor = element_line(colour = "grey92"))+
  #     guides(color=guide_legend(override.aes=list(fill=NA), order =1),
  #            linetype = guide_legend(order = 2,override.aes = list(colour = "black")))+
  #     ggtitle("Variation of phenotype's number throughout the simulation")
  #   
  # }
  
  
  # color gradient
  # colfunc <- colorRampPalette(c("green4","orange", "steelblue"))
  # colGrad <- colfunc(length(SpIdx))
  # names(colGrad) <- SpIdx
  # if (length(SpIdx) == 9) colLine <- c("solid","dashed","solid","dashed","solid","longdash","solid","longdash","solid") else colLine = rep("solid",3) # 3,6 and 8 are dashed
  # 
  # 
  # 
  # p <- ggplot(plot_dat) +
  #   geom_line(aes(x = time, y = popN-exitN, color = as.factor(species)))+
  #   geom_line(aes(x = time, y = popF-exitF,color = as.factor(species)), linetype = "dashed")+
  #   scale_x_continuous(name = "Time (yr)") +
  #   scale_y_continuous(name = "Number of phenotypes")+
  #   scale_color_manual(name = "Species", values = colGrad)+ # color gradient
  #   theme(legend.title=element_blank(),
  #         legend.position=c(0.19,0.95),
  #         legend.justification=c(1,1),
  #         legend.key = element_rect(fill = "white"),
  #         panel.background = element_rect(fill = "white", color = "black"),
  #         panel.grid.minor = element_line(colour = "grey92"))+
  #   guides(color=guide_legend(override.aes=list(fill=NA)))+
  #   ggtitle("Variation of phenotype's number throughout the simulation")
  
  if(returnData) return(myData) else return(p)
  
}

plotTrait <- function(object, SpIdx = NULL, dt = 0.1, Normalisation = F, returnData = T, traitID = NULL)
{
  # Beware this is gonna be hard
  
  # Initialisation // bunch of set I need to set up now
  if (is.null(SpIdx)) SpIdx <- 1:9 # species index
 
   Wmat <- object@params@species_params$w_mat[SpIdx] # maturation size to only select mature phenotypes
  
  # if (save_it & is.null(dir))  # if want to save but not specified where
  #   dir = paste(getwd(),"/temporary",sep = "")
  # 
  # ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir)), FALSE) #create the file if it does not exists
  
  
  # Set-up the data right, we need the right apparition, extinction, traits, ... it is done dynamically mouahahah
  SumPar = object@params@species_params #shortcut
  # TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),dt*SumPar$pop,dt*SumPar$extinct,SumPar[traitID]) # weird things happen without the as.numeric / *0.1 because dt / need to update this to have the traits as arguments
  TT = data.frame("Species" = SumPar$species, "Phenotype" =  as.numeric(SumPar$ecotype),"Apparition" = dt*SumPar$pop,
                  "Extinction" = dt*SumPar$extinct,SumPar[traitID]) # weird things happen without the as.numeric / *0.1 because dt / need to update this to have the traits as arguments
  
  # colnames(TT) = c("Species","Phenotype","Apparition","Extinction","Td","EdInt")
  # rownames(TT) = rownames(SumPar)
  TT = TT[order(TT[,1],decreasing=FALSE),]
  # TT = as.data.frame(TT) # I use TT later in the graphs (its like an artefact)
  
  TimeMax <- SumPar$timeMax[1] * dt
  for (i in 1:dim(TT)[1]) if (TT$Extinction[i] == 0) TT$Extinction[i] = TimeMax # Fill the extinction value for the non-extinct species
  
  # Weighted mean of trait value
  # 1) matrix of summed abundance of mature individuals at each time step
  biomass = object@n
  # put 0 in object@n when w < w_mat
  for (iTime in 1:dim(biomass)[1]) # for each time step
  {
    for (iPhen in 1:dim(biomass)[2]) # for each ecotypes
    {
      w_lim = SumPar$w_mat[iPhen] # get the maturation size of the ecotype
      S <- numeric(length(object@params@w))
      S[sapply(w_lim, function(x) which.min(abs(x - object@params@w)))] <- 1 # find what w bin is the closest of the maturation size
      NoW_mat = which(S == 1) # what is the col number of that size bin
      biomass[iTime,iPhen,1:NoW_mat-1] <-0 # everything under this value become 0
    }
  }
  
  abundanceM = apply(biomass, c(1,2),sum) # sum the abundance left 
  
  # 2) normalisation per species 
  colnames(abundanceM) = SumPar$species # phenotypes from same species have the same name
  abundanceNormal = matrix(0,nrow = dim(abundanceM)[1], ncol = dim(abundanceM)[2])
  
  # I am getting rid of the species which went instinct at the begining and that went extinct without having mutants (no trait variation)
  # SpIdx = NULL
  # for (i in unique(SumPar$species))
  #   if (sum(abundanceM[,i]) != 0 & dim(SumPar[SumPar$species == i,])[1] != 1) # if not extinct at the beginning and more than one ecotype (for the apply)
  #     SpIdx = c(SpIdx,i)
  
  # SpIdx is annoying:
  # if (length(SpIdx) > length(unique(SumPar$species))) SpIdx = unique(SumPar$species) # ok so I might have species not even reaching this point so I'm short cutting spidx automatcaly
  
  for (iSpecies in SpIdx)
  {
    abundanceSp = abundanceM # save to manipulate
    abundanceSp[,which(colnames(abundanceM) != iSpecies)] = 0 # make everything but the targeted species to go 0 to have correct normalisation
    abundanceSp = sweep(abundanceSp,1,apply(abundanceSp,1,sum),"/") # normalise
    abundanceSp[is.nan(abundanceSp)] <-0 # when I divide by 0 I get nan
    abundanceNormal = abundanceNormal + abundanceSp # I just need to add them up to get the final matrix
  }
  # if (comments) cat(sprintf("Abundance normalised\n"))
  
  # Now I have normalised abundance, I need to apply the trait I want to plot on them

  
  plotList <- vector(mode = "list", length = length(traitID))
  dataList <- vector(mode = "list", length = length(traitID))
  
  # TD
  
  for(iTrait in traitID)
  {
    plotStore <- list()
    dataStore <- list()
  abundanceT = sweep(abundanceNormal,2,as.matrix(SumPar[iTrait]),"*") # I use the normalised abundance and multiply by the trait value

  # Calculate mean at each time step
  TotMean = matrix(0, nrow = dim(abundanceT)[1], ncol = max(as.numeric((unique(SumPar$species)))), 
                   dimnames = list(rownames(abundanceT),as.character(seq(1,max(as.numeric(unique(SumPar$species))))))) #sort(unique(SumPar$species))[length(unique(SumPar$species))]
  names(dimnames(TotMean)) = list("time","species")
  
  for (i in SpIdx)
  {
    AMean = abundanceT[,which(colnames(abundanceT) == i)] # get the portion of the matrix I want (right species)
    if (is.null(dim(AMean)) ==F) AMean = apply(AMean,1,sum) # it is already weighted so I'm just summing to get the mean # if only one sp no need to do it
    TotMean[,i] = AMean
  }
  
  # Calculate variance and standard deviation
  # it is the sum of the difference between value and mean squared and multiplied by the weight
  
  statData = list() # list with the stats of all species as I need to do this for each species separatly
  
  for (i in SpIdx)
  {
    meanSp = TotMean[,i] # take the mean of the species
    traitSp = SumPar$w_mat[SumPar$species == i] # take the traits of the ecotypes in the species
    weightSp = abundanceNormal[,which(colnames(abundanceNormal) == i)] # take the proportion of each ecotypes
    stat = matrix(cbind(as.numeric(rownames(abundanceT)),meanSp,0,0), ncol = 5,nrow = length(meanSp), dimnames = list(NULL,c("time","mean","sd","percentMean","percentSd"))) # initialise the matrix
    for (j in 1:length(meanSp)) # for each time step
    {
      if (is.null(dim(weightSp))) {variance = sum(((traitSp-meanSp[j])^2)*weightSp[j])} else {variance = sum(((traitSp-meanSp[j])^2)*weightSp[j,])} # calculate the variance, condition if only one phen
      stat[j,3] = sqrt(variance) # calculate the standard deviation
      stat[j,4] = (meanSp[j] - SumPar$t_d[1])/SumPar$t_d[1] # normalisation of mean
      stat[j,5] = stat[j,3]/SumPar$t_d[1] # normalised sd
    }
    statData[[i]] = as.data.frame(stat) # put in the list
  }
  
  # I have the stats for every species, just need to plot now
  
  
  for (i in SpIdx)
  {
    stat = statData[[i]] # take the stats of the right species
    
    phenotype = TT[TT$Species == i,c("Apparition","Extinction",iTrait)] # recup the traits time values
    Phen = melt(phenotype,iTrait) # make the dataframe
    #name = paste("Maturation size of species ",i, sep = "")
    
    #short cut the data frame when species does not reach end of simulation
    if (sum(which(Phen$value == TimeMax)) == 0) # if there is at least one value equal to the end of the sim it means that the species survived until then
      stat = stat[-which(stat$mean == 0),] # first occurence of mean = 0, meaning dead
    
    name = paste("Species",i, sep = " ")
    
    # prepare the data for the ribbon
    g1 <- ggplot(stat)+
      geom_smooth(aes(x = time, y = percentMean-percentSd)) +
      geom_smooth(aes(x= time, y = percentMean+percentSd)) 
    
    gg1 <- ggplot_build(g1)
    dfRibbon <- data.frame(x = gg1$data[[1]]$x, ymin = gg1$data[[1]]$y, ymax = gg1$data[[2]]$y) #and extract the smooth data for the ribbon
    
    if (!Normalisation) stat$percentMean = stat$mean
    
    
    p = ggplot() +
      geom_smooth(data = stat, aes(x = time, y = percentMean)) +
      # geom_ribbon(data = dfRibbon, aes(x = x, ymin = ymin, ymax = ymax), fill = "grey", alpha = 0.4)+
      # geom_hline(yintercept = 0, linetype = "dashed") +
      theme(legend.title=element_blank(),panel.background = element_rect(fill = "white", color = "black"),
            panel.grid.minor = element_line(colour = "grey92"), legend.position="none", 
            legend.justification=c(1,1),legend.key = element_rect(fill = "white"))+ #none remove legend; things below change depending on the graph position
      scale_x_continuous(name = "Time", limits  = c(0,TimeMax)) +
      scale_y_continuous(name = name, limits = c(NA,NA)) + #, breaks = seq(window[1],window[2],0.1)) + #how can seq() be so bad at his only job
      ggtitle(paste("Trait =",iTrait))
    
    
    
    plotStore[[i]] = p
    dataStore[[i]] <- stat
    
  }
  
  plotList[[which(iTrait == traitID)]] <- plotStore
  dataList[[which(iTrait == traitID)]] <- dataStore
  }
  
  if(returnData)  return(dataList) else return(plotList)
  
}

plotTraitPhen <- function(object, trait = NULL, SpIdx = NULL, dt = 0.1, print_it = F, returnData = T, save_it = F, nameSave = "TraitPhen.png", ylimit = c(NA,NA))
{
  
  # Initialisation // bunch of set I need to set up now
  if (is.null(SpIdx)) SpIdx <- 1:9 # species index
  if(is.null(trait)) trait <- "td"
  
  # Set-up the data right, we need the right apparition, extinction, traits, ...
  SumPar = object@params@species_params #shortcut
  TT = cbind(SumPar$species,as.numeric(SumPar$ecotype),dt*SumPar$pop,dt*SumPar$extinct,SumPar$t_d,SumPar$ed_int) # weird things happen without the as.numeric / *0.1 because dt / need to update this to have the traits as arguments
  colnames(TT) = c("species","phenotype","apparition","extinction","td","ed_int")
  rownames(TT) = rownames(SumPar)
  TT = TT[order(TT[,1],decreasing=FALSE),]
  TT = as.data.frame(TT) # I use TT later in the graphs (its like an artefact)
  
  # need to change ecotype name for original species
  for (iSpecies in SpIdx)
  {
    nameList <- TT$phenotype
    target <- which(TT$phenotype == iSpecies)
  while (TT$phenotype[target] %in% nameList) TT$phenotype[target] = as.numeric(paste(TT$species[target],sample(seq(1:1e5),1),sep="")) # take 5 random digits to follow the digit species identity as a name
  }
  
  TimeMax <- SumPar$timeMax[1] * dt
  for (i in 1:dim(TT)[1]) if (TT$extinction[i] == 0) TT$extinction[i] = TimeMax # Fill the extinction value for the non-extinct species
  
  plot_dat <- NULL
  for(iSpecies in SpIdx)
  {
    species_dat <- filter(TT,species == iSpecies)  
    
    switch (trait,
            "td" = {selectTrait = species_dat$td},
            "ed_int" = {selectTrait = species_dat$ed_int},
            {print("trait not known")}
    )
    
    
    phen_dat <- data.frame("species" = rep(t(species_dat$species),2), "phen" = rep(t(species_dat$phenotype),2), 
                           "time" = rbind(t(t(species_dat$apparition)),t(t(species_dat$extinction))), "trait" = rep(t(selectTrait),2))
    plot_dat <- rbind(plot_dat,phen_dat)
  }
  
  
  p <- ggplot(plot_dat)+
    geom_line(aes(x=time,y=trait, group = phen, color = as.factor(species))) +
    scale_x_continuous(name = "Time in years", limits = c(NA,NA))+
    scale_y_continuous(name = "Trait value", limits = ylimit)+
    facet_grid(species ~.)+
    scale_color_manual(name = "Species", values = colGrad)+ # color gradient
    theme(legend.title=element_text(),panel.background = element_rect(fill = "white", color = "black"),
          panel.grid.minor = element_line(colour = "grey92"), legend.position="bottom",
          #legend.justification=c(1,1),
          legend.key = element_rect(fill = "white"))+
    guides(color = guide_legend(nrow=1)) +
    ggtitle("Trait every phenotypes")
  
  
  if(save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)
  
  if (returnData) return(plot_dat) else if(print_it) return(p)
}


plotDietCompPhen<-function(object, prey=dimnames(object@diet_comp)$prey, min_w=.001,
                           predator=dimnames(object@diet_comp)$predator, timeaverage=FALSE, print_it = F, returnData = T){
  
  # need to sum the phenotypes within species
  no_sp = length(unique(object@params@species_params$species))
  myDiet <- array(NA, dim = c(no_sp,dim(object@diet_comp)[2:4]), dimnames = list("predator" = 1:no_sp, "pred_size" = dimnames(object@diet_comp)$pred_size,
                                                                                 "prey" = dimnames(object@diet_comp)$prey,"prey_size" = dimnames(object@diet_comp)$prey_size))
  for(iSpecies in 1:no_sp)
    myDiet[iSpecies,,,] <- apply(object@diet_comp[which(dimnames(object@diet_comp)$predator == iSpecies),,,],c(2,3,4), sum)
  
  
  myDiet2 <- array(NA, dim = c(no_sp,dim(object@diet_comp)[2],no_sp+1, dim(object@diet_comp)[4]), dimnames = list("predator" = 1:no_sp, "pred_size" = dimnames(object@diet_comp)$pred_size,
                                                                                                                  "prey" = c(1:no_sp,"background"),"prey_size" = dimnames(object@diet_comp)$prey_size))
  
  for(iSpecies in 1:no_sp)
    myDiet2[,,iSpecies,] <- apply(myDiet[,,which(dimnames(myDiet)$prey == iSpecies),],c(1,2,4), sum)
  
  myDiet2[,,no_sp+1,] <- myDiet[,,"plankton",]
  
  predator=dimnames(myDiet2)$predator
  prey=dimnames(myDiet2)$prey
  out<-myDiet2
  
  # default function
  prey_nam<-prey
  pred_nam<-predator
  
  prey<-apply(out, c(1,2,3), FUN=sum) #Sum across size classess with in prey 
  tot<-apply(prey, c(1,2), FUN=sum) #Sum across prey species 
  
  prey_prop<-sweep(prey, c(1,2), tot, "/") # Get proportion of diet for each species
  
  no_pred<- length(dimnames(prey_prop)[[1]])
  no_pred_w<- length(dimnames(prey_prop)[[2]])
  no_prey<- length(dimnames(prey_prop)[[3]])
  
  #Stacked  bar chart 
  plot_dat<-expand.grid(dimnames(prey_prop)[[1]], dimnames(prey_prop)[[2]], dimnames(prey_prop)[[3]])
  colnames(plot_dat)<-c("predator","predsize","prey")
  plot_dat$predsize<-as.numeric(as.character(log10(as.numeric(as.character(plot_dat$predsize)))))
  
  plot_dat$value<- as.vector(prey_prop)
  
  
  species<-object@params@species_params$species
  wmin<-object@params@w[object@params@species_params$w_min_idx]
  wmax<-object@params@w[object@params@species_params$w_max_idx]
  
  for ( i in 1:length(species)){
    plot_dat$value[plot_dat$predator==species[i] & plot_dat$predsize < log10(wmin[i])]<- 0
    plot_dat$value[plot_dat$predator==species[i] & plot_dat$predsize > log10(wmax[i])]<- 0
  }
  
  
  plot_dat<-subset(plot_dat,predsize> log10(min_w))
  
  dsub<-plot_dat[ plot_dat$prey %in% prey_nam, ]
  dsub<-dsub[ dsub$predator %in% pred_nam, ]
  dsub[is.na(dsub)]<-0
  
  p<-  ggplot(data = dsub, aes(x = predsize, y = value, fill = prey)) + geom_area( position = 'stack')  + facet_wrap(~predator, ncol=5) + scale_color_brewer(palette="Set1") +
    scale_x_continuous(name = "log10 predator mass (g)") + scale_y_continuous(name = "Proportion of diet by mass (g)")
  
  if(returnData) return(dsub) else if (print_it) print(p)
  
}