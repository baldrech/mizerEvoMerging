# this file is here to run the model composed of different functions

# I have to pass all the parameters here to be able to build mutants again
myModel <- function(no_sp = 9, # number of species #param described in Andersen & Pedersen 2010
                    min_w_inf = 10, # minimum weight of sp
                    max_w_inf = 1e4, # maximum weight of sp
                    no_w = 100, # number of size bins community spectrum
                    min_w = 0.001, #min size bin of community spectrum/The smallest size of the species community size spectrum
                    max_w = max_w_inf * 1.1, #max size bin of both spectrum
                    min_w_pp = 1e-10, #min size bin of background size spectrum
                    # no_w_pp = round(no_w)*0.3, # number of size bins background spectrum
                    w_pp_cutoff = 1, # cut of size of the background spectrum
                    k0 = 50, # recruitment adjustment parameter
                    n = 0.75, # exponent of maximum intake (scaling of intake)
                    p = 0.75, # exponent of standard metabolism
                    q = 0.8, # exponent of search volume
                    eta = 0.25, # size at maturation relative to Mg (mass in grams ?)
                    r_pp = 4, # growth rate of resource spectrum (primary production)
                    kappa = 0.005, # ressource spectrum carrying capacity
                    lambda = 2+q-n, # exponent of the background spectrum.
                    alpha = 0.6, # assimilation efficiency
                    ks = 4, # factor for standard metabolism
                    z0pre = 2, # background mortality factor
                    h = 85, # factor of maximum intake
                    beta = 100, # preferred predator-prey weight ratio
                    sigma = 1, # width of selection function
                    gamma = NA,
                    f0 = 0.5, # average feeding level of the community/feeding level of small individuals feeding on background
                    knife_edge_size = 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp) * eta , #knife edge position
                    gear_names = rep("FishingStuff", no_sp),
                    t_max = 100,
                    dt = 0.1,
                    mu = 1,
                    extinct = TRUE, # extinction option
                    RMAX = TRUE, # enable egg density dependence
                    rm = NULL, # set up rmax by user
                    OptMutant = "M5", # mutation depends on number of eggs or not?
                    Trait = "unknow", # No default, needs to be specified by user
                    r_mult = 1e0, #rmax multiplier to try things
                    erepro = 1, # reproduction efficiency
                    no_run = 1, # number of sim in a row to do
                    effort = 0,
                    cannibalism = 1, # stay here for now but going away soon
                    initCondition = NULL, # if I want to input previous mizer object as initial condition
                    initTime = 1, # time for initialisation
                    param = NULL, # can input a param data.frame to do multispecies model
                    print_it = T, # if you want to display messages or not
                    normalFeeding = F, #if want to normalised feeding kernel
                    mAmplitude = 0.2, # width of distribution of new trait value
                    save_it = F, # do you want to save?
                    path_to_save = NULL, # where?
                    predMort = NULL, # if want to replace dynamics m2 by constant one
                    initPool = 0,
                    tau = 7, # exponent in psi function
                    interactionOne = 0.5, # to set_up interaction with one parameter
                    interaction = matrix(interactionOne,nrow=no_sp, ncol=no_sp), # default interaction matrix, controlled by the interactionOne param
                    diet_steps = 10, # for the diet thing
                    # extension
                    kappa_ben = 0.005,
                    r_bb = 2,
                    min_w_bb = 1e-10,  
                    w_bb_cutoff = 10,
                    lambda_alg = 2+q-n,
                    kappa_alg = 0.005,
                    r_aa = 2,
                    min_w_aa = 1e-10,  
                    w_aa_cutoff = 100,
                    t_ref = 10,
                    t_d = 25,
                    temperature = NA,
                    ea_met = NA,
                    ca_met = NA,
                    ea_int = NA,
                    ca_int = NA,
                    ea_mat = NA,
                    ca_mat = NA,
                    ea_mor = NA,
                    ca_mor = NA,
                    ed_met = NA,
                    cd_met = NA,
                    ed_int = NA,
                    cd_int = NA,
                    ed_mat = NA,
                    cd_mat = NA,
                    ed_mor = NA,
                    cd_mor = NA,
                    avail_PP = 1,
                    avail_BB = 0,
                    avail_AA = 0,
                    ...){
tic()
  
  # deactivate temperature functions if temperature is NA
  # if(is.na(temperature))
  # {
  #   ea_met = 0 # their default is not 0
  #   ea_int = 0
  #   ed_int = 0
  #   
  # }
  
  if (is.null(initCondition))
  {
    firstRun = 1
    s_max = no_run * t_max / dt
    # I'm deleting all the default from this function so it uses only the ones in myModel
    if (is.null(param))
      param <- set_trait_model(no_sp = no_sp, 
                       min_w_inf = min_w_inf, 
                       max_w_inf = max_w_inf,
                       no_w = no_w, 
                       min_w = min_w, 
                       max_w = max_w,
                       min_w_pp = min_w_pp,
                       w_pp_cutoff = w_pp_cutoff,
                       k0 = k0, 
                       n = n, 
                       p = p, 
                       q = q,
                       eta = eta, 
                       r_pp = r_pp, 
                       kappa = kappa,
                       lambda = lambda, 
                       alpha = alpha, 
                       ks = ks, 
                       z0pre = z0pre, 
                       h = h, 
                       beta = beta, 
                       sigma = sigma,
                       gamma = gamma,
                       f0 = f0, 
                       knife_edge_size = knife_edge_size,
                       gear_names = gear_names,
                       r_mult = r_mult,
                       cannibalism = cannibalism,
                       erepro = erepro,
                       s_max = s_max,
                       rm = rm,
                       normalFeeding = normalFeeding,
                       tau = tau,
                       interaction = interaction,
                       # extension
                       kappa_ben = kappa_ben,
                       r_bb = r_bb,
                       min_w_bb = min_w_bb,  
                       w_bb_cutoff = w_bb_cutoff,
                       lambda_alg = lambda_alg,
                       kappa_alg = kappa_alg,
                       r_aa = r_aa,
                       min_w_aa = min_w_aa,  
                       w_aa_cutoff = w_aa_cutoff,
                       t_ref = t_ref,
                       t_d = t_d,
                       ea_met = ea_met,
                       ca_met = ca_met,
                       ea_int = ea_int,
                       ca_int = ca_int,
                       ea_mat = ea_mat,
                       ca_mat = ca_mat,
                       ea_mor = ea_mor,
                       ca_mor = ca_mor,
                       ed_met = ed_met,
                       cd_met = cd_met,
                       ed_int = ed_int,
                       cd_int = cd_int,
                       ed_mat = ed_mat,
                       cd_mat = cd_mat,
                       ed_mor = ed_mor,
                       cd_mor = cd_mor,
                       avail_PP = avail_PP,
                       avail_BB = avail_BB,
                       avail_AA = avail_AA
                       ) 
    
    # Initialisation ---------------------
    # Mutant option
    M3List <- list() # So I'm creating this list to store parameters from user input and only have one thing to pass between functions
    if (OptMutant == "M3") # means that we need to know when user want mutation to appear
    {
      prompt <- "At what time do you want the mutation to occur?\n"
      M3List[[1]] <- as.integer(strsplit(readline(prompt), " ")[[1]])
      if (length(M3List[[1]]) == 0 )M3List[[1]] = 0
    }
    
    # Kick start the abundance
    if (print_it) cat(sprintf("Initialisation of the simulation, please wait.\n"))
    initBio <- project(param, t_max = initTime, extinct = FALSE, OptMutant="yo", RMAX = RMAX, diet_steps = 0) # init abundance
    #initBio <- project(param, t_max = 1, extinct = FALSE, OptMutant="yo", RMAX = T) # init abundance
    n_init <- initBio@n[dim(initBio@n)[1],,]
    n_pp_init <- initBio@n_pp[dim(initBio@n_pp)[1],]
    n_aa_init <- initBio@n_aa[dim(initBio@n_aa)[1],]
    n_bb_init <- initBio@n_bb[dim(initBio@n_bb)[1],]
    if (print_it) cat(sprintf("Initialisation completed, starting simulation.\n"))
    nameList = initBio@params@species_params$ecotype
    
    # Generate a base of phenotypes around each species
    if (initPool > 0)
    {
      for (iSpecies in sort(unique(param@species_params$species)))
      {
        # Generate phenotypes pool
        for (iPhenotype in seq(1, initPool))
        {
          mutant <- param@species_params[param@species_params$ecotype == iSpecies,] # perfect copy
          while (mutant$ecotype %in% nameList) mutant$ecotype = as.numeric(paste(mutant$species,sample(seq(1:1e5),1),sep="")) # take 5 random digits to follow the digit species identity as a name
          
          switch(Trait,
                 size = {
                   # Trait = asymptotic size
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$w_inf)
                   mutant$w_inf <- mutant$w_inf + rnorm(1, 0, sd) # change a bit the asymptotic size
                   mutant$w_mat <- mutant$w_inf * eta # calculate from the new w_inf value
                   mutant$z0 <- z0pre * as.numeric(mutant$w_inf) ^ (n - 1) # if I don't put as.numeric I lose the name z0
                 },
                 beta = {
                   # Trait = PPMR
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$beta)
                   mutant$beta <- mutant$beta + rnorm(1, 0, sd) # change a bit the PPMR
                   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)
                   mutant$gamma <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 sigma = {
                   # Trait = fedding kernel
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$sigma)
                   mutant$sigma <- mutant$sigma + rnorm(1, 0, sd) # change a bit the diet breadth
                   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)
                   mutant$gamma <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 predation = {
                   # PPMR
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$beta)
                   mutant$beta <- mutant$beta + rnorm(1, 0, sd) # change a bit the PPMR
                   # feeding kernel
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$sigma)
                   mutant$sigma <- mutant$sigma + rnorm(1, 0, sd) # change a bit the diet breadth
                   # recalculate gamma if necessary
                   # alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)
                   # mutant$gamma <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 eta = {
                   # Trait = eta
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$eta)
                   mutant$eta <- mutant$eta + rnorm(1, 0, sd) # change a bit eta
                   mutant$w_mat <- mutant$w_inf * mutant$eta # update
                 },
                 ed_int = {
                   # Trait = ed_int
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$ed_int)
                   mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
                   while(mutant$ed_int < 2.5) mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
                 },
                 t_d = {
                   # Trait = ed_int
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$t_d)
                   mutant$t_d <- mutant$t_d + rnorm(1, 0, sd)
                 },
                 temperature = {
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$ed_int)
                   mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
                   while(mutant$ed_int < 2.5) mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd) # ed_int cannot go lower than 2.5
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$t_d)
                   mutant$t_d <- mutant$t_d + rnorm(1, 0, sd)
                 },
                 
                 all = {
                   # Trait = asymptotic size
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$w_inf)
                   mutant$w_inf <- mutant$w_inf + rnorm(1, 0, sd) # change a bit the asymptotic size
                   mutant$w_mat <- mutant$w_inf * eta # calculate from the new w_inf value
                   mutant$z0 <- z0pre * as.numeric(mutant$w_inf) ^ (n - 1) # if I don't put as.numeric I lose the name z0
                   # Trait = predation
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$beta)
                   mutant$beta <- mutant$beta + rnorm(1, 0, sd) # change a bit the PPMR
                   sd = as.numeric(mAmplitude *  param@species_params[which(param@species_params$ecotype == iSpecies),]$sigma)
                   mutant$sigma <- mutant$sigma + rnorm(1, 0, sd) # change a bit the diet breadth
                   # calculate the new gamma
                   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)
                   mutant$gamma <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 {
                   print("Trait specified is not in the list")
                 })
          
          # integrate the mutant in the df 
          rownames(mutant) = mutant$ecotype
          param@species_params <- rbind(param@species_params, mutant) #include the mutant in the dataframe
          
          #updatenameList
          nameList = param@species_params$ecotype
          
          #mutant abundance
          n_mutant <- rep(0,no_w)
          n_init <- rbind(n_init,n_mutant) # this include the new mutant as last column
          rownames(n_init)[length(rownames(n_init))] <- rownames(mutant) # update the name of the mutant accordingly
          
          # need to change the interaction matrix
          interaction <- rbind(interaction,interaction[which(rownames(param@interaction) == iSpecies),])
          interaction <- cbind(interaction,interaction[,which(colnames(param@interaction) == iSpecies)])
        }
      }
    }
    #need to update some suff now that there is one more sp
    no_sp = dim(param@species_params)[1]
    # Recreate the "param" object needed for the projection
    

    
    if (print_it) cat(sprintf("Creating first MizerParams object.\n"))
    param <- MizerParams(param@species_params, no_w = no_w,  w_pp_cutoff = w_pp_cutoff, max_w = max_w, min_w_pp = min_w_pp,
                         n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda, t_ref = t_ref,
                         # normalFeeding = normalFeeding, tau = tau, 
                         interaction = interaction)
    if (print_it) cat(sprintf("Done\n"))

    # redistribute the abundance of the phenotypes more randomly
    template <- n_init[1:no_sp,]
    for (iSpecies in 1:no_sp) # for every species
    {
      #the total abundance is randomly distributed among all the phenotypes
      biomRandom<- runif(initPool+1,0,1)
      biomFrac <- biomRandom/sum(biomRandom) # that's the fraction to apply to the initial abundance to spread the biomass among phenotypes
      position = 0
      for(iPhen in which(param@species_params$species == iSpecies))
      {
        position = position + 1
        n_init[iPhen,] <- template[iSpecies,] * biomFrac[position]
      }
    }
    
  }
  else 
  {
    Nparam <- initCondition@params@species_params[initCondition@params@species_params$extinct == F,] # take the sp not extinct at start the sim
    interaction <- initCondition@params@interaction[which(rownames(initCondition@params@interaction) %in% Nparam$ecotype),which(colnames(initCondition@params@interaction) %in% Nparam$ecotype)] # take the interaciton matrix of non-extinct species
    
    for (i in unique(Nparam$species)) Nparam[which(Nparam$species == i),]$knife_edge_size <- knife_edge_size[i] # update knife edge
    
    Nparam$timeMax = no_run * t_max / dt # update the time max of the sim /// start from beginning
    param <- MizerParams(Nparam, no_w = no_w, w_pp_cutoff = w_pp_cutoff,  max_w = max_w, min_w_pp = min_w_pp,
                         n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda, t_ref = t_ref,
                         # normalFeeding = normalFeeding, tau = tau, 
                         interaction = interaction)
    spIndex = as.character(Nparam$ecotype)
    initCondition@n = initCondition@n[,spIndex,] # take the abundance of only the present species
    n_init = initCondition@n[dim(initCondition@n)[1],,] # take last time step of the abundance to make it first time step
    n_pp_init = initCondition@n_pp[dim(initCondition@n_pp)[1],] # same for plankton
    n_aa_init = initCondition@n_aa[dim(initCondition@n_aa)[1],] # same for plankton
    n_bb_init = initCondition@n_bb[dim(initCondition@n_bb)[1],] # same for plankton
    if (print_it) cat(sprintf("Starting simulation with previous run.\n"))
    no_run = no_run + max(Nparam$timeMax)/t_max*dt # update number of runs
    firstRun = max(Nparam$timeMax)/t_max*dt +1 # state at which run we're starting
    nameList = initCondition@params@species_params$ecotype # this list keep in memory all the species name (as I lose some in my ecotypes by getting rid of the extinct/ use to give ecotypes namee)
  }
  
  # need some specific stuff when running single species model
  oneSpMode = F
  if (no_sp == 1)
  {
    oneSpMode = T
    cat(sprintf("Simulation in mode: one species only\n"))
  }
  
  #Multiple run --------------------------------
  allRun <- list() # save all the runs
  interactionSave <- param@interaction # save the interaction matrix at the begining of the simulation
  
  # need to check temperature format
  if (length(temperature) == 1 && !is.na(temperature)) temperature = rep(temperature, times = t_max*no_run)

    for(j in firstRun:no_run){
    # I am ordering everything by appartition order. To keep that even if I stop and re initialise the sim, I need to change the run number
    # it means that if I do a sim after another one, the first run wont be one but the previous number of run + one
    if (print_it) cat(sprintf("run = %s\n",j))
    
    # Select the right temperature vector for the run/ must be t_max length
    temperature_vec <- temperature[seq((j-1)*t_max+1,j*t_max)]
    if (print_it && !is.na(temperature))
      {cat(sprintf("temperature for the run:\n"))
    print(temperature_vec)}
    
    # First run without mutants
    sim <- project(param, t_max = t_max, dt =dt, mu = mu, initial_n = n_init, initial_n_pp=n_pp_init, initial_n_aa=n_aa_init, 
                   initial_n_bb=n_bb_init, extinct = extinct, RMAX=RMAX, OptMutant=OptMutant, M3List = M3List, 
                   checkpoint = j, effort = effort, print_it = print_it, predMort = predMort, diet_steps = diet_steps, temperature = temperature_vec) # init first step
    
    # Post initialisation -------------------

      allData <- list() # this list will save the data output of all the projections
      counter = 1 # used to increment the list (I guess there is a better way to do that)
      
      while (length(sim) > 3 ) # ugly but if everything is done, length(sim) = 1, if sim dead, length =2, if a mutant appear, length = 5( sim,time,resident, n , npp)
      {
        n_init = sim$n # last time abundance, that will be modified and use as initiation for next projection
        
        # SAVE THE DATA FROM PREVIOUS PROJECTION
        allData[[counter]] <- sim$data
        counter = counter +1
        
        # CREATE MUTANTS
        for (i in 1: length(sim$resident))
        {
          resident = sim$resident[i] # this is the parent phenotype
          resident_params = sim$data@params@species_params[sim$data@params@species_params$ecotype == resident,]
          # if(!is.data.frame(resident_params) || length(resident_params) != 26) 
          # {print("error in resident")
          #   print(resident_params)}
          mutant <- resident_params # create a perfect copy
          mutant$pop = sim$i_stop + (j-1)*t_max/dt # what i_time does the mutant appear
          mutant$run = j
          
          mutant$ecotype =  as.numeric(unlist(strsplit(as.character(resident), "")))[1] # take the first digit of the parent name (species identity)
          if(!is.numeric(mutant$ecotype))
          {
            cat(sprintf("something is wrong with the mutant name: %i\n", mutant$ecotype))
          }
          while (mutant$ecotype %in% nameList) mutant$ecotype = as.numeric(paste(as.numeric(unlist(strsplit(as.character(resident), "")))[1],sample(seq(1:1e5),1),sep="")) # take 5 random digits to follow the digit species identity as a name
          # TRAITS
          switch(Trait,
                 size = {
                   # Trait = asymptotic size
                   sd = as.numeric(mAmplitude *  resident_params["w_inf"]) # standard deviation
                   mutant["w_inf"] <- resident_params["w_inf"] + rnorm(1, 0, sd) # change a bit the asymptotic size
                   mutant["w_mat"] <- mutant["w_inf"] * eta # calculate from the new w_inf value
                   mutant["z0"] <- z0pre * as.numeric(mutant["w_inf"]) ^ (n - 1) # if I don't put as.numeric I lose the name z0
                   #cat(sprintf("Its size mutes slightly.\n"))
                 },
                 beta = {
                   # Trait = PPMR
                   sd = as.numeric(mAmplitude *  resident_params["beta"]) # standard deviation
                   mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd) # change a bit the PPMR
                   # calculate the new gamma
                   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)
                   mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
                   #cat(sprintf("Its PPMR mutes slightly.\n"))
                 },
                 sigma = {
                   # Trait = fedding kernel
                   sd = as.numeric(mAmplitude *  resident_params["sigma"]) # standard deviation
                   mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd) # change a bit the diet breadth
                   # calculate the new gamma
                   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)
                   mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
                   #cat(sprintf("Its diet breadth mutes slightly.\n"))
                 },
                 predation = {
                   # PPMR
                   sd = as.numeric(mAmplitude *  resident_params["beta"]) # standard deviation
                   mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd) # change a bit the PPMR
                   # feeding kernel
                   sd = as.numeric(mAmplitude *  resident_params["sigma"]) # standard deviation
                   mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd) # change a bit the diet breadth
                   # recalculate gamma if necessary
                   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)
                   mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
                 },
                 eta = {
                   # Trait = eta
                   sd = as.numeric(mAmplitude *  resident_params["eta"]) # standard deviation
                   mutant["eta"] <- resident_params["eta"] + rnorm(1, 0, sd) # change a bit eta
                   if (mutant["eta"] >= 1) mutant["eta"] <- 0.95 # because yes it does happen
                   mutant["w_mat"] <- mutant["w_inf"] * mutant["eta"] # update
                   #cat(sprintf("Its w_mat is: %g\n",mutant["w_mat"]))
                 },
                 ed_int = {
                   # Trait = ed_int
                   sd = as.numeric(mAmplitude *  resident_params["ed_int"])
                   mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
                   while(mutant$ed_int < 2.5) mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd) # ed_int cannot go lower than 2.5
                 },
                 t_d = {
                   # Trait = ed_int
                   sd = as.numeric(mAmplitude *  resident_params["t_d"])
                   mutant$t_d <- mutant$t_d + rnorm(1, 0, sd)
                   cat(sprintf("Its name is %i and its trait value is %g\n", mutant$ecotype,mutant["t_d"]))
                 },
                 temperature = {
                   sd = as.numeric(mAmplitude * resident_params["ed_int"])
                   mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd)
                   while(mutant$ed_int < 2.5) mutant$ed_int <- mutant$ed_int + rnorm(1, 0, sd) # ed_int cannot go lower than 2.5
                   sd = as.numeric(mAmplitude * resident_params["t_d"])
                   mutant$t_d <- mutant$t_d + rnorm(1, 0, sd)
                 },
                 all = {
                   # Trait = asymptotic size
                   sd = as.numeric(mAmplitude *  resident_params["w_inf"]) # standard deviation
                   mutant["w_inf"] <- resident_params["w_inf"] + rnorm(1, 0, sd) # change a bit the asymptotic size
                   mutant["w_mat"] <- mutant["w_inf"] * eta # calculate from the new w_inf value
                   mutant["z0"] <- z0pre * as.numeric(mutant["w_inf"]) ^ (n - 1) # if I don't put as.numeric I lose the name z0
                   # Trait = predation
                   sd = as.numeric(mAmplitude *  resident_params["beta"]) # standard deviation
                   mutant["beta"] <- resident_params["beta"] + rnorm(1, 0, sd) # change a bit the PPMR
                   sd = as.numeric(mAmplitude *  resident_params["sigma"]) # standard deviation
                   mutant["sigma"] <- resident_params["sigma"] + rnorm(1, 0, sd) # change a bit the diet breadth
                   # calculate the new gamma
                   alpha_e <- sqrt(2 * pi) * mutant$sigma * mutant$beta ^ (lambda - 2) * exp((lambda - 2) ^ 2 * mutant$sigma ^ 2 / 2)
                   mutant["gamma"] <- h * f0 / (alpha_e * kappa * (1 - f0))
                   #cat(sprintf("Its traits mute slightly.\n"))
                 },
                 {
                   print("congrats, you managed to fuck up somewhere")
                 })
          
          # cat(sprintf("Its name is %i and its trait value is %g\n", mutant$ecotype,mutant["w_mat"]))
          
          # I need to specify the name myself as the dataframe way is not consistant and subject to errors. It will work as long as a parent has less than 1e5 mutants
          rownames(mutant) = mutant$ecotype
          sim$data@params@species_params <- rbind(sim$data@params@species_params, mutant) #include the mutant in the dataframe
          #need to update some suff now that there is one more sp
          no_sp = no_sp + 1
          w_inf <- as.numeric(unlist(sim$data@params@species_params["w_inf"])) # need to recreate the vector
          w_mat <-  as.numeric(unlist(sim$data@params@species_params["w_mat"]))

          # There are 2 interaction matrix: interaction which is the current one and does not have extinct species and interactionSave that has everything
          # Need to recreate the interaction matrix for every new mutant
          interaction <- rbind(interaction,interaction[which(rownames(sim$data@params@interaction) == resident),])
          interaction <- cbind(interaction,interaction[,which(colnames(sim$data@params@interaction) == resident)])
          rownames(interaction) <- sim$data@params@species_params$ecotype
          colnames(interaction) <- rownames(interaction)
          interactionSave <- rbind(interactionSave,interactionSave[which(rownames(sim$data@params@interaction) == resident),])
          interactionSave <- cbind(interactionSave,interactionSave[,which(colnames(sim$data@params@interaction) == resident)])
          rownames(interactionSave)[length(rownames(interactionSave))] <- mutant$ecotype
          colnames(interactionSave)[length(colnames(interactionSave))] <- mutant$ecotype

          if(sum(colnames(interaction) %in% rownames(interaction)) != dim(interaction)[1])
          {
            cat(sprintf("interaction names are wrong\n"))
            print(colnames(interaction))
            print(rownames(interaction))
          }
          
          if(sum(colnames(interaction) %in% colnames(interactionSave)) != dim(interaction)[1])
          {
            cat(sprintf("something happened\n"))
            print(colnames(interaction))
            print(colnames(interactionSave))
          }
          
          # Recreate the "param" object needed for the projection
          trait_params <- MizerParams(sim$data@params@species_params,  no_w = no_w, w_pp_cutoff = w_pp_cutoff, max_w = max_w,  min_w_pp = min_w_pp,
                                      n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda, t_ref = t_ref,
                                      # normalFeeding = normalFeeding, tau = tau, 
                                      interaction = interaction)
          ## TODO ### need to fix this name shenaniggan
          # print("output Mizerparam")
          # print(trait_params@interaction)
          # 
          # object <- trait_params
          # print(dimnames(object@psi)[[1]])
          # print(dimnames(object@intake_max)[[1]])
          # print(dimnames(object@search_vol)[[1]])
          # print(dimnames(object@metab)[[1]])
          # print(dimnames(object@ft_pred_kernel_e)[[1]])
          # print(dimnames(object@ft_pred_kernel_p)[[1]])
          # print(dimnames(object@mu_b)[[1]])
          # print(dimnames(object@selectivity)[[2]])
          # print(dimnames(object@catchability)[[2]])
          # print(dimnames(object@interaction)[[1]])
          # print(dimnames(object@interaction)[[2]])
          # print(object@species_params$ecotype)
          
          # # Use this piece of code if you want to update r_max, as it depends on the number of species. I'm not doing it though, r_max is inherited. 
          # # warning, beta is not updated here
          # alpha_p <- f0 * h * beta^(2 * n - q - 1) * exp((2 * n * (q - 1) - q^2 + 1) * sigma^2 / 2)
          # alpha_rec <- alpha_p / (alpha * h * f0 - ks)
          # # Calculating dw using Ken's code - see Ken's email 12/08/13
          # tmpA <- w_inf[1]
          # tmpB <- (log10(w_inf[length(w_inf)]) - log10(w_inf[1])) / (no_sp - 1) # Difference between logged w_infs, fine
          # dw_winf <- tmpB * tmpA *10^(tmpB*((1:no_sp)-1)) # ?
          # N0_max <- k0 * w_inf^(n*2-q-3+alpha_rec) * dw_winf  # Why * dw_winf, not / ? Ken confirms * in email
          # # No need to include (1 - psi) in growth equation because allocation to reproduction at this size = 0, so 1 - psi = 1
          # g0 <- (alpha * f0 * h * trait_params@w[1]^n - ks * trait_params@w[1]^p)
          # r_max <- N0_max * g0
          #trait_params@species_params$r_max <- r_max
          
          # abundance of the mutant
          n_mutant <- rep(0,no_w)
          n_mutant = 0.05 * n_init[dimnames(sim[[1]]@n)$sp ==resident,] # the initial abundance is 5% of the resident pop
          n_init[dimnames(sim[[1]]@n)$sp ==resident,]= n_init[dimnames(sim[[1]]@n)$sp ==resident,] - 0.05*n_init[dimnames(sim[[1]]@n)$sp ==resident,] # Witdraw the abundance of the mutant from its parent (we're not talking about eggs here but different ecotype already present)
          n_init <- rbind(n_init,n_mutant) # this include the new mutant as last column
          rownames(n_init)[length(rownames(n_init))] <- rownames(mutant) # update the name of the mutant accordingly
          sim$n <- n_init
        }
        
        # in the case of constant predation mortality
        if (!is.null(predMort)){
          predMort = matrix(data = predMort[1,], nrow = no_sp, ncol = dim(predMort)[2],byrow = T)
        }
        # let's start projecting again
        # sim <- project(trait_params, t_max = t_max, dt = dt, i_stop = sim$i_stop, initial_n = n_init, initial_n_pp=sim$n_pp, mu = mu, 
        #                extinct = extinct, RMAX=RMAX,OptMutant=OptMutant, M3List = M3List,checkpoint = j, effort = effort, 
        #                print_it = print_it, predMort = predMort )
        # print("sim specs")
        # print(class(sim))
        # print(summary(sim))
        sim <- project(trait_params, t_max = t_max, dt = dt, prevSim = sim, mu = mu,
                       extinct = extinct, RMAX=RMAX,OptMutant=OptMutant, M3List = M3List, checkpoint = j, effort = effort,
                       print_it = print_it, predMort = predMort, diet_steps = diet_steps, temperature = temperature_vec)
        
      }
      
      
      # if simulation went extinct
      ##TODO### need to update this I reckon, not sure if it works
      if (length(sim) == 2) 
      {
        allData[[counter]] <- sim[[1]] #
        for (i in 1:(length(allData)-1)) # change the time max of the sim as it's shorter now
        {
          allData[[i]]@params@species_params$timeMax = length(sim[[1]])*t_max/dt
        }
        # need to add 0 biomass to the last sim to make it last until t_max and then do final touch
        
        return(allData)
      }
      
      allData[[counter]] <- sim # one last time for the last projection
            
      
      # now allData has all the successive runs, lets stitch them
      biomass <- stitch(allData) # biomass is a list of n and background
      sim@n = biomass[[1]]
      dimnames(sim@n)$sp <- sim@params@species_params$ecotype # to be able to sort the alive from extincts
      sim@n_pp = biomass[[2]]
      sim@n_aa = biomass[[3]]
      sim@n_bb = biomass[[4]]
      
      # now I want to do more run with, as initial conditions, the final step of the previous run
      # but first I need to save it
      allRun[[j]] <- sim
      
      # then let's clean the sim of the extinct species and initial next sim
      Nparam = sim@params@species_params[sim@params@species_params$extinct == F,]
      # need to change the interaction matrix as well
      if(sum(sim@params@species_params$extinct != F)>0) interaction = interaction[-c(which(sim@params@species_params$extinct != F)),-c(which(sim@params@species_params$extinct != F))] #get rid of the lines in the interaction matrix when species are extinct
      
      # print(2)
      # print(param@w_full)
      # print(param@species_params)
      # print(Nparam$species)
      
      param <- MizerParams(Nparam, no_w = no_w, w_pp_cutoff = w_pp_cutoff, max_w = max_w,  min_w_pp = min_w_pp,
                           n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda, t_ref = t_ref,
                           # normalFeeding = normalFeeding, tau = tau,
                           interaction = interaction)
      # print(3)
      # print(param@w_full)
      # print(param@species_param)
      spIndex = as.character(Nparam$ecotype)
      n_init = sim@n[dim(sim@n)[1],spIndex,]
      n_pp_init = sim@n_pp[dim(sim@n_pp)[1],]
      n_aa_init = sim@n_aa[dim(sim@n_aa)[1],]
      n_bb_init = sim@n_bb[dim(sim@n_bb)[1],]
      

  }
  # all runs done
  
  # final param counting the extinct species
    a = NULL
    for (i in firstRun:length(allRun) ) a = rbind(a,allRun[[i]]@params@species_params) # bind the different dataframes
    a <- a[order(a$ecotype, a$extinct, decreasing=TRUE),] # weird 3 lines to get rid of duplicates and keep the ones with the extinction value
    a <- a[!duplicated(a$ecotype),]
    SummaryParams = a[order(a$pop,a$ecotype),]
    
    # At this stage, the paramDF does not remember phenotypes that appeared and went extinct in the same run, but the interaction matrix does
    # Lazy/easy way -> get rid of them in the interaction matrix as well
    
    interactionSave <- interactionSave[as.numeric(rownames(interactionSave)) %in% SummaryParams$ecotype,as.numeric(rownames(interactionSave)) %in% SummaryParams$ecotype]
    
    if (dim(interactionSave)[1] != dim(SummaryParams)[1] || dim(interactionSave)[2] != dim(SummaryParams)[1])
    {
      cat(sprintf("dim(interactionSave)[1] = %i\n",dim(interactionSave)[1]))
      cat(sprintf("dim(interactionSave)[2] = %i\n",dim(interactionSave)[2]))
      cat(sprintf("dim(SummaryParams)[1] = %i\n",dim(SummaryParams)[1]))
      cat(sprintf("rownames(interactionSave):\n"))
      print(rownames(interactionSave))
      cat(sprintf("SummaryParams:\n"))
      print(SummaryParams$ecotype)
      interactionSave <- matrix(0.5,ncol = dim(SummaryParams)[1], nrow = dim(SummaryParams)[1])
    }
    # Update all the other param from the dataframe
    FinalParam <- MizerParams(SummaryParams, no_w = no_w, w_pp_cutoff = w_pp_cutoff, max_w = max_w,  min_w_pp = min_w_pp,
                              n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda, t_ref = t_ref,
                              # normalFeeding = normalFeeding, tau = tau, 
                              interaction = interactionSave)
    # handle and save the final data
    sim = finalTouch(list(allRun,FinalParam),temperature = temperature, print_it = print_it)
    gc()
    sim = superOpt(sim)
    if (save_it)
    {
      if (is.null(path_to_save)) path_to_save = paste(getwd(),"/temporary",sep="")
      ifelse(!dir.exists(file.path(path_to_save)), dir.create(file.path(path_to_save),recursive = T), FALSE) #create the file if it does not exists
      save(sim,file = paste(path_to_save,"/run.Rdata",sep="")) #save it
    }
    toc()
    return(sim)
}
