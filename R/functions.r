# Functions for MCMC and approximate likelihood computations.
# Written by Jonas Knape

#' @useDynLib stageFreq
#' @importFrom Rcpp sourceCpp
NULL

##' The function is intended to be used in the calculation of Monte Carlo likelihoods.
##' @title Simulate gamma distributed stage durations.
##' @param mu A vector of length equal to one minus the number of stages giving the mean stage durations 
##'           for all but the last stage.
##' @param cv A vector of the same length as \code{mu} giving the coefficient of variation, which is equal to one 
##'           over the square root of the shape, for all but the last stage.
##' @param N The number of stage durations (individuals) to simulate.
##' @return A matrix of stage durations where each row gives tha stage durations for one individual.
##' @examples 
##' ## Simulate gamma stage durations for ten individuals.
##' rGammaSD(mu = c(5,2,4), cv = c(1,.5,.3), N = 10)
##' @export
rGammaSD = function(mu, cv, N) {
  sDur = matrix(nrow = N, ncol = length(mu) + 1)
  for (i in 1:(ncol(sDur) - 1)) 
    sDur[, i] = rgamma(N, shape = 1 / cv[i] ^ 2, 
                       scale = mu[i] * cv[i] ^ 2)
  sDur[, ncol(sDur)] = Inf # Assumes adults stay adults forever
  sDur
}

##' The function is intended to be used in the calculation of Monte Carlo likelihoods.
##' @title Simulate Weibull distributed stage durations.
##' @param shape A vector of length equal to one minus the number of stages giving the shape of the stage durations 
##'           for all but the last stage.
##' @param scale A vector of length equal to one minus the number of stages giving the scale of the stage durations 
##'           for all but the last stage.
##' @param N The number of stage durations (individuals) to simulate.
##' @return A matrix of stage durations where each row gives tha stage durations for one individual.
##' @examples 
##' ## Simulate Weibull stage durations for ten individuals.
##' rWeibSD(shape = c(3,1,2), scale = c(3,10,5), N = 10)
##' @export
rWeibSD = function(shape, scale, N) {
  sDur = matrix(nrow = N, ncol = length(scale) + 1)
  for (i in 1:(ncol(sDur) - 1)) 
    sDur[, i] = rweibull(N, shape = shape, 
                       scale = scale)
  sDur[, ncol(sDur)] = Inf # Assumes adults stay adults forever
  sDur
}

##' The function is intended to be used in the calculation of numeric likelihoods.
##' @title Compute gamma stage duration densities.
##' @param mu A vector of length equal to one minus the number of stages giving the mean stage durations 
##'           for all but the last stage.
##' @param cv A vector of the same length as \code{mu} giving the coefficient of variation, which is equal to one 
##'           over the square root of the shape, for all but the last stage.
##' @param nGrid The number of grid points to use.
##' @param tMax The maximum time point in the grid (i.e. the right endpoint).
##' @return A matrix of stage duration probabilities where the first column gives the grid 
##'         and the remaining columns give the stage duration probabilities.
##' @examples 
##' ## Simulate gamma stage durations for ten individuals.
##' dGammaSD(mu = c(5,2,4), cv = c(1,.5,.3), nGrid = 20, tMax = 15)
##' @export
dGammaSD = function(mu, cv, nGrid, tMax) {
  grid = seq(0, tMax, length.out = nGrid)
  delta = grid[2] - grid[1]
  nGrid = length(grid)
  dDur = matrix(nrow = nGrid, ncol = length(mu) + 2)
  dDur[, 1] = grid
  for (i in 1:(ncol(dDur) - 1)) {
    pg = pgamma(c(grid, grid[nGrid] + delta), shape = 1 / cv[i] ^ 2, 
                scale = mu[i] * cv[i] ^ 2)
    dDur[, i + 1] = pg[2:(nGrid + 1)] - pg[1:nGrid]
  }
  dDur[, ncol(dDur)] = 0 # Assumes adults stay adults forever
  dDur
}

##' The function is intended to be used in the calculation of numeric likelihoods.
##' @title Compute Weibull stage duration densities.
##' @param shape A vector of length equal to one minus the number of stages giving the shape of the stage durations 
##'           for all but the last stage.
##' @param scale A vector of length equal to one minus the number of stages giving the scale of the stage durations 
##'           for all but the last stage.
##' @param nGrid The number of grid points to use.
##' @param tMax The maximum time point in the grid (i.e. the right endpoint).
##' @return A matrix of stage duration probabilities where the first column gives the grid 
##'         and the remaining columns give the stage duration probabilities.
##' @examples 
##' ## Compute Weibull stage duration probabilities over a grid of 20 points.
##' dWeibSD(shape = c(3,1,2), scale = c(3,10,5), nGrid = 20, tMax = 15)
##' @export
dWeibSD = function(shape, scale, nGrid, tMax) {
  grid = seq(0, tMax, length.out = nGrid)
  delta = grid[2] - grid[1]
  nGrid = length(grid)
  dDur = matrix(nrow = nGrid, ncol = length(scale) + 2)
  dDur[, 1] = grid
  for (i in 1:(ncol(dDur) - 1)) {
    pg = pweibull(c(grid, grid[nGrid] + delta), shape = shape[i], 
                scale = scale[i])
    dDur[, i + 1] = pg[2:(nGrid + 1)] - pg[1:nGrid]
  }
  dDur[, ncol(dDur)] = 0 # Assumes adults stay adults forever
  dDur
}

##' Convert a data frame or matrix into an sfMat object, to be used with the likelihood functions
##' \code{\link{mnLogLik}}, \code{\link{mnLogLik2}} and \code{\link{poisLogLik}}.
#' @title Convert a data frame or matrix into an sfMat object.
#' @param data A matrix or data frame containing stage frequency data where each row gives the
##'             sample at a specific time point.
##' @param timeVar The name of the column of \code{data} holding the time variable,
##'                which defines at which time a sample is taken.
##' @param stageVars A vector of names of the columns of \code{data} holding the stage counts. 
##'                  The vector should be ordered from the first to the last stage.
##' @param N0 The name of the column giving the known initial number of individuals at each 
##'           sampling time. May be set to NULL if not known.
##' @param gamma The name of the column giving the sampling effort at each time. 
##'              Used by \code{\link{poisLogLik}}. May be set to NULL (the default).
##' @return A data matrix of class sfMat.
##' @export
stageFreqMat = function(data, timeVar, stageVars, N0 = NULL, gamma = NULL) {
  data = as.matrix(data)
  if (!(timeVar %in% colnames(data)))
    stop(paste0(timeVar, " not found in data."))
  attr(data, "time") = timeVar
  sCol = integer(length(stageVars))
  for (i in 1:length(stageVars)) {
    ind = which(colnames(data) == stageVars[i])
    if (length(ind) == 0)
      stop(paste0(stageVars[i], " not found in data."))
    sCol[i] = ind
  }
  attr(data, "sCol") = sCol
  if (!is.null(N0))
    attr(data, "N0") = N0
  if (!is.null(gamma))
    attr(data, "gamma") = gamma
  class(data) = "sfMat"
  data
}

##' Approximate log likelihood for Poisson sampling model of stage frequency data.
##' @title Approximate Poisson log likelihood for stage frequency data.
##' @param data A data matrix of class sfMat. 
##' @param rStageDur A function for simulating stage durations. Should return a matrix where the columns
##'                  give the simulated durations for each stage (in order), and the rows represent a simulated individual.
##'                  Simulated durations may be equalt to Inf if an indiviudal does not mature from a stage.
##' @param mortRate A vector of mortality rates for each stage. Should be of the same length as 
##'                 the number of columns returned by rStageDur.
##' @param lam0 A vector 
##' @param ... Further arguments passed to rStageDur.
##' @return The approximate (stochastic) log-likelihood value.
##' @export
poisLogLik = function(data, rStageDur, mortRate, lam0, ...) {
  sCol = attr(data, "sCol")
  if (!is.null(attr(data, "gamma")))
    gamma = data[, attr(data, "gamma")]
  else
    gamma = rep(1, nrow(data))
  nStage = length(sCol)
  sDur = do.call(rStageDur, list(...))
  sp = stageFreqMC(data[, attr(data, "time")], sDur, mortRate) / nrow(sDur)
  sum(dpois(data[,sCol], sp*lam0*gamma, log = TRUE))
}



##' Approximate Monte Carlo log likelihood for multinomial sampling with
##' known initial number of individuals at start of development.
##' @title Approximate Monte Carlo multinomial log likelihood for stage frequency data with known initial numbers.
##' @param data A data matrix of class sfMat. 
##' @param rStageDur A function for simulating stage durations. Should return a matrix where the columns
##'                  give the simulated durations for each stage (in order), and the rows represent a simulated individual.
##'                  Simulated durations may be equalt to Inf if an indiviudal does not mature from a stage.
##' @param mortRate A vector of mortality rates for each stage. Should be of the same length as 
##'                 the number of columns returned by rStageDur.
##' @param ... Further arguments passed to rStageDur.
##' @return The approximate (stochastic) log-likelihood value.
##' @export
mnLogLik = function(data, rStageDur, mortRate, ...) {
  sCol = attr(data, "sCol")
  nStage = length(sCol)
  sDur = do.call(rStageDur, list(...))
  sp = stageFreqMC(data[, attr(data, "time")], sDur, mortRate) / nrow(sDur)
  pd = 1-apply(sp, 1, sum)
  ll = 0
  for (i in 1:nrow(data)) {
    ll = ll + dmultinom(c(data[i, sCol], data[i, attr(data, "N0")] - sum(data[i, sCol])), prob = c(sp[i, ], pd[i]), log = TRUE)
  }
  ll
}

##' Approximate numerical log likelihood for multinomial sampling with
##' known initial number of individuals at start of development.
##' @title Approximate numerical multinomial log likelihood for stage frequency data with known initial numbers.
##' @param data A data matrix of class sfMat. 
##' @param dStageDur A function computing stage duration densities. Should return a matrix where the firs column 
##'                  give the grid over which numeric integration is to be performed and the remaining columns
##'                  give the density of the stage durations evaluated at the grid point.
##'                  Simulated durations may be equalt to Inf if an indiviudal does not mature from a stage.
##' @param mortRate A vector of mortality rates for each stage. Should be of the same length as 
##'                 the number of columns returned by rStageDur.                
##' @param ... Further arguments passed to rStageDur.
##' @return The approximate numeric log-likelihood value.
##' @export
mnLogLikN = function(data, dStageDur, mortRate,...) {
  sCol = attr(data, "sCol")
  nStage = length(sCol)
#  grid = seq(0, max(data[,attr(data, "time")]) + 1, length.out = nGrid)
  dDur = do.call(dStageDur, list(...))
  sp = stageFreqN(data[, attr(data, "time")], dDur[, 2:(nStage+1)], mortRate, grid = dDur[,1])
  pd = 1-apply(sp, 1, sum)
  pd[pd<=0] = 1e-300
  sp[sp<=0] = 1e-300
  ll = 0
  for (i in 1:nrow(data)) {
    if(min(sp[i,], pd[i])<0 | max(sp[i,], pd[i]) == 0) browser()
    ll = ll + dmultinom(c(data[i, sCol], data[i, attr(data, "N0")] - sum(data[i, sCol])), prob = c(sp[i, ], pd[i]), log = TRUE)
    if (!is.finite(ll)) browser()
  }
  ll
}

##' @title Approximate multinomial log likelihood for stage frequency data with fixed number of 
##'        individuals.
##' @param data A data matrix of class sfMat. 
##' @param rStageDur A function for simulating stage durations. Should return a matrix where the columns
##'                  give the simulated durations for each stage (in order), and the rows represent a simulated individual.
##'                  Simulated durations may be equalt to Inf if an indiviudal does not mature from a stage.
##' @param mortRate A vector of mortality rates for each stage. If NULL (the default), mortality rates are 
##'                 assumed to be zero.
##' @param ... Further arguments passed to rStageDur.
##' @return The approximate (stochastic) log-likelihood value.
##' @export
mnLogLik2 = function(data, rStageDur, mortRate = NULL, ...) {
  sCol = attr(data, "sCol")
  nStage = length(sCol)
  if (is.null(mortRate))
    mortRate = rep(0, nStage)
  sDur = do.call(rStageDur, list(...))
  sp = stageFreqMC(data[, attr(data, "time")], sDur, mortRate = mortRate) / nrow(sDur)
  ll = 0
  for (i in 1:nrow(data)) {
    ll = ll + dmultinom(data[i, sCol], prob = sp[i, ], log = TRUE)
  }
  ll
}

##' @title Compute stage probabilites by numeric integration.
##' @param sampleTimes Vector of times at which the cohort was sampled.
##' @param dDur A matrix of stage duration densities.
##' @param mortRate Mortality rates for each stage.
##' @param grid A grid over which numeric integration of stage duration probabilities is done.
##' @return A matrix of stage counts. Each row gives the number of individuals in each stage for a given sampling time.             
##' @export
stageFreqN = function(sampleTimes, dDur, mortRate, grid) {
  ## Bellows and Birley method
  nStage = ncol(dDur)
  nGrid = length(grid)
  nSamp = length(sampleTimes)
  if (grid[nGrid] < sampleTimes[nSamp])
    stop("End point of grid smaller than latest sampling time. Extend grid.")
  fout = matrix(ncol = nStage, nrow = nSamp)
  P = matrix(ncol = nStage, nrow = nGrid)
  f = matrix(ncol = nStage, nrow = nGrid)
  P[, 1] = c(1, rep(0, nGrid - 1))
  for (j in 2:nStage) {
    H = cumsum(dDur[,j - 1]) 
    f[,j-1] = convolve((1-H) * exp(-mortRate[j-1] * grid), rev(P[, j-1]), type = "open")[1:nGrid]
    fout[, j-1] = approx(grid, f[, j-1], xout = sampleTimes)$y
    P[,j] = convolve(dDur[, j - 1] * exp(-mortRate[j-1] * grid), rev(P[,j-1]), type = "open")[1:nGrid]
  }
  H = cumsum(dDur[,nStage]) 
  f[, nStage] = convolve((1-H) * exp(-mortRate[nStage] * grid), rev(P[, nStage]), type = "open")[1:nGrid]
  fout[, nStage] = approx(grid, f[, nStage], xout = sampleTimes)$y
  fout[fout<0] = 0
  fout
}



##' @title Simulate stage frequency data with known initial number of individuals for each sampling occasion.
##' @param sampleTimes Vector of times at which the cohort is sampled. Development is assumed to start at time 0.
##' @param N0 An integer vector of length equal to the length of \code{sampleTimes}, or of length one, giving the known
##'           number of individuals at the start of development. If the length is one, the value is taken as the number 
##'           initial  individuals for eac sampling time.
##' @param rStageDur A function simulating stage durations, and returning stage durations as a matrix where rows
##'                  represent individuals and columns represent stage durations for each stage. The function should
##'                  have an argument for the number of individuals to simulate. 
##' @param mortRate A vector of mortality rates for each stage.
##' @param Narg The name of the argument to the function \code{rStageDur} which determines the number of individuals
##'             that are simulated.
##' @param ... Further arguments passed to \code{rStageDur}.
##' @return An sfData object.          
##' @export
##' @examples 
##' data = simMNCohort(seq(1,50, length.out=15), N0 = 30, rGammaSD, mu = c(5, 15, 8), cv = c(.5,.5,.5), 
##'                    mortRate = c(0.01, 0.01,0.01,0.01), Narg = "N")
simMNCohort = function(sampleTimes, N0, rStageDur, mortRate, Narg, ...) {
  if (length(N0) == 1) 
    N0 = rep(N0, length(sampleTimes)) 
  nStage = length(mortRate)
  data = matrix(0, nrow = length(sampleTimes), ncol = nStage)
  rsArgs = list(N0[1], ...)
  names(rsArgs)[1] = Narg
  for (i in 1:nrow(data)) {
    rsArgs[[Narg]] = N0[i]
    sDur = do.call(rStageDur, rsArgs)
    Nt = stageFreqMC(sampleTimes[i], sDur, mortRate)
    Nt0 = stageFreqMC(sampleTimes[i], sDur, 0 * mortRate)
    nonZ = which(Nt > 0)
    data[i, nonZ] = rbinom(length(nonZ), size = Nt0[nonZ], prob = Nt[nonZ] / Nt0[nonZ])
  }
  colnames(data) = paste0("stage", 1:nStage)
  data = cbind(day = sampleTimes, N0 = N0, data)
  attr(data, "N0") = "N0"
  attr(data, "time") = "day"
  attr(data, "sCol") = 3:(2 + nStage)
  class(data) = "sfMat"
  data
}

##' The function is intended to be used with computationally expensive likelihoods, 
##' and little attention has been giving to speed up the adaptive MCMC algorithm itself.
##' The underlying algorithm is an adaptive MCMC sampler with adaptive global scaling 
##' and using the empirical covariance matrix to generate proposals.
##' 
##' @title Estimate a continuous posterior distribution via adaptive MCMC.
##' @param inits Vector of initial parameter values. The log likelihood function and prior need to have
##'              finite values for these parameters.
##' @param logLik A function giving the log likelihood value for the model. 
##'               The function should take parameter values to be estimated via an argument called pars.
##' @param logPrior A function giving the log prior density for the parameters. 
##'                 The function should take a single argument named pars. 
##'                 For parameter values that are out of range the function should return -Inf.
##' @param nIter    The size of the MCMC sample.
##' @param nThin The number of thining iterations in between each saved sample.
##' @param logScaleMCMC If true, the MCMC is run on the log scale of all parameters. Requires that all parameters are positive.
##' @param blocks May be used for blocking MCMC updates. In this case a vector of integers of 
##'               the same length as inits where each unique integer represents a block.
##'               If NULL, no blocking is applied.
##' @param ... Further arguments passed to \code{logLik}.
##' @return An \code{\link[coda]{mcmc}} object.          
##' @export
##' @examples
##' ## Simulate a cohort with four stages, sampled at 15 occasions with 30 initial individuals in 
##' ## each sample.
##' 
##' cdata = simMNCohort(seq(1,50, length.out=15), N0 = 30, rGammaSD, 
##'                     mu = c(5, 15, 8), cv = c(.5,.5,.5), mortRate = c(0.01, 0.01,0.01,0.01), 
##'                     Narg = "N")
##' 
##' ## Define a likelihood function using the function mnLogLik and rGammaSD. 
##' ## Note that it must have the argument named 'pars'.
##' 
##' exLogLik = function(pars, data, nSim) {
##'   mnLogLik(data, rStageDur = rGammaSD, mu = pars[1:3], cv = rep(pars[4], 3), 
##'         mortRate = rep(pars[5], 4),  N = nSim)
##' }
##' 
##' ## Define the log-prior function using bounded uniform priors.
##' 
##' exLogPrior = function(pars) {
##'  mu = pars[1:3]
##'  cv = pars[4]
##'  mr = pars[5]
##'  sum(log(all(c(cv < 2,mu < 30,mr < 1,cv > 0,mu > 0,mr > 0))))
##' }
##' 
##' \dontrun{
##' ## Fit model
##' inits = c(10,10,10,1,.001) 
##' names(inits) = c(paste0("mu", 1:3), "cv", "mr")
##' mcmc = adMCMC(inits, exLogLik, exLogPrior, nIter = 1000, nThin = 10, 
##'                     data = cdata, nSim = 1000)
##' plot(mcmc)
##' }
adMCMC = function (inits, logLik, logPrior, nIter, nThin, logScaleMCMC = TRUE,
                      blocks = NULL, ...) {
  # Define transform for MCMC (log scale or original scale) 
  if (logScaleMCMC) {
    transf = log
    transfInv = exp
    logdTransfInv = identity
  } else {
    transf = identity
    transfInv = identity
    logdTransfInv = function(x) return(0)  
  }
  logLikP = NA * numeric(nIter)
  logLikI = NA * numeric(nIter)
  if (!("pars" %in% names(formals(logLik)))) {
    stop("Log likelihood function needs to have an argument named 'pars'")
  }
  if (is.null(names(inits)))
    names(inits) = paste0("p", 1:length(inits))
  parsChain = new("adrwMCMC", varNames=names(inits), maxIter = nIter, thin=nThin, 
                  blocks = blocks, idVarScale = 0.1^2, varScale = .38^2)  
  logLikArgs = list(pars = inits, ...)
  parsChain$setInits(transf(inits), logLik = do.call(logLik, logLikArgs),    
                     logPrior = logPrior(inits) + sum(logdTransfInv(transf(inits))))
  if (!is.finite(parsChain$getLogLik())) {
    stop(paste("Non finite likelihood value", parsChain$getLogLik(), "at initial values.",
               "(Failure could be due to stochastic likelihood evaluation.)"))
  }
  if(!is.finite(parsChain$getLogPrior())) {
    stop(paste("Non finite prior value", parsChain$getLogPrior(), "at initial values."))
  }
  
  pts = proc.time()
  pb <- txtProgressBar(min = 0, max = nIter * nThin, style = 3)
  plot(0, 0, xlim = c(0,nIter), ylim = c(0,1), pch = 20, xlab = "iteration", ylab = "acceptance rate")
  for (i in 2:(nIter*nThin)) {
    if (is.null(blocks)) {
      propPars = parsChain$getProposal()
    } else
      propPars = parsChain$getProposal(block = sample(blocks, 1))
    propLogPrior = logPrior(transfInv(propPars)) + sum(logdTransfInv(propPars))
    if (!is.finite(propLogPrior) & !identical(propLogPrior, -Inf)) {
      browser()
      stop(paste("Non finite value", propLogPrior, 
                 "returned by prior function at proposed parameter values", transfInv(propPars)))
    }
    if (propLogPrior > - Inf) {
#      logLikArgs$pars = replace(pars, estInd, transfInv(propPars[repInd]))
      logLikArgs$pars = transfInv(propPars)
      propLL = do.call(logLik, logLikArgs)
    } else {
      propLL = -Inf
    }
    if (!is.finite(propLL) & !identical(propLL, -Inf)) {
      stop(paste("Non finite value", propLL, 
                 "returned by likelihood function at proposed parameter values", transfInv(propPars)))
    }
    laccP = (propLL - parsChain$getLogLik()) + (propLogPrior - parsChain$getLogPrior())
    if (log(runif(1, 0, 1)) < laccP) {
      parsChain$acceptProposal(TRUE, logLik = propLL, logPrior = propLogPrior) #, accP = exp(laccP))
    } else {
      parsChain$acceptProposal(FALSE, logLik = propLL, logPrior = propLogPrior) #, accP = exp(laccP))
    }  
    if (i %% nThin == 0) {
      setTxtProgressBar(pb, i)
      points(i / nThin, parsChain$getAcceptanceRate(), pch = 20)
    }
  }
  writeLines("\nProcessing time (seconds):")
  print(proc.time() - pts)
  writeLines(paste0("Total acceptance rate: ", 100*parsChain$getAcceptanceRate(), "%"))
  invisible(coda::as.mcmc(transfInv(parsChain$getChain())))
}