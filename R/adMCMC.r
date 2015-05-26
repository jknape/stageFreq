# A reference class for handling adaptive MCMC
# Written by Jonas Knape
library(coda)
setRefClass(Class = "MCMC",
            fields = c(".chain", ".nPar", ".thin", ".ichain", ".iteration", 
                       ".proposal", ".current", ".logLik", ".logPrior",".propSet",
                       ".accepted",".nAccept"),
            methods = list(
              initialize = function(varNames=c(), thin=1, maxIter=1000) {
                .iteration <<- 1
                .nPar <<- length(varNames)
                .chain <<- matrix(rep(NA, .nPar * maxIter), nr = maxIter, nc = .nPar)
                .ichain <<- NA
                colnames(.chain) <<- varNames
                .thin <<- thin
                .propSet <<- FALSE
                #.updType <<- integer(maxIter * thin)
                .accepted <<- logical(maxIter * thin)
                .nAccept <<- 0
              },
              setInits = function(inits, logLik = NULL, logPrior = NULL) {
                .current <<- matrix(inits, nrow = 1, ncol = .nPar)
                .proposal <<- .current
                .logLik <<- logLik
                .logPrior <<- logPrior
                .ichain <<- 1
                .chain[.ichain, ] <<- .current
              },
              getProposal = function() {
                if (.propSet) {
                  return(.proposal)
                } else {
                  stop("Proposal not set.")
                }
              },
              setProposal = function(proposal) {
                .proposal <<- proposal  
                .propSet <<- TRUE
              },
              getCurrent = function() {
                .current
              },
              getLogLik = function() {
                if(is.null(.logLik)) {
                  stop("logLik not set.")
                }
                .logLik
              },
              getLogPrior = function() {
                if(is.null(.logPrior)) {
                  stop("logPrior not set.")
                }
                .logPrior
              },
              getAcceptanceRate = function() {
                .nAccept/.iteration
              },
              acceptProposal = function(accepted, logLik = NULL, logPrior = NULL) {
                if (!.propSet) 
                  stop("Proposal not set, nothing to accept.")  
                if (accepted) {
                  .current <<- .proposal
                  .logLik <<- logLik
                  .logPrior <<- logPrior
                  .accepted[.iteration] <<- TRUE
                  .nAccept <<- .nAccept + 1
                }
                if (.iteration %% .thin == 0) {
                  .ichain <<- .ichain + 1
                  .chain[.ichain, ] <<- .current
                }
                .iteration <<- .iteration+1
                .propSet <<- FALSE
              },
              getIteration = function() {
                .iteration
              }
              ,
              getChain = function() {
                invisible(.chain)
              }
            )
)


setRefClass(Class = "adMCMC",
            fields = c(".adSchedule"),
            methods = list(
              initialize = function(thin = 1, maxIter = 1000, adSchedule = round(seq(ceiling(maxIter/1000), maxIter, length.out = 1000)), ...) {
                callSuper(thin = thin, maxIter = maxIter, ...)
                .adSchedule <<- sort(unique(adSchedule)) #[adSchedule > thin]
              }
            ),
            contains = "MCMC"
)


setRefClass(Class = "adrwMCMC",
            fields = c(".chVar", ".parBlocks", ".idPropProb", ".idVarScale", ".varScale", ".rescaleStep", ".updType", ".phase2"),
            methods = list(
              initialize = function(blocks = NULL, idProb = .05, idVarScale = 0.1^2, varScale = 2.38^2,...) {
                callSuper(...)
                .chVar <<- sqrt(.1) * diag(.nPar)
                if (is.null(blocks)){
                  .parBlocks <<- rep(1, .nPar)
                } else {
                  if (length(blocks) != .nPar)
                    stop("Length of blocks not equal to number of parameters")
                  .parBlocks <<- blocks
                }
                .phase2 <<- nrow(.chain)/20
                .idPropProb <<- idProb
                .idVarScale <<- idVarScale
                .varScale <<- varScale
                .updType <<- integer(length(.accepted))
              },
              getProposal = function(block = NULL) {
                if (!.propSet) {
                  if (.ichain %in% .adSchedule & .iteration %% .thin == 0) {
                    #if (abs(sum(.chVar) - sum(diag(.chVar))) > .00001)
                    #print(c(.iteration, .varScale, locAccR), digits = 4)
                    if (.ichain > .phase2) {
                     rs = min(c(.thin * diff(.adSchedule), .iteration - 1))
                     locAccR = sum(.accepted[(.iteration - rs - 1):(.iteration - 1)])/rs
                     .varScale <<- .varScale * exp( (.iteration * 1000 / (.thin * nrow(.chain))) ^ (-1/2) * (locAccR - .3))
                     .chVar <<- tryCatch(chol(var(.chain[.phase2:.ichain, , drop = FALSE])), 
                                        error = function (e) {diag(1, .nPar)})
                    } else {
                      if(.ichain > 1) {
                       sdd = apply(.chain[1:.ichain, , drop = FALSE], 2, sd)
                       if (min(sdd) > 0)
                        .chVar <<- diag(apply(.chain[1:.ichain, , drop = FALSE], 2, sd))
                      }
                    }
                    #.chVar <<- diag(.nPar)
                  }
                  #if(.ichain == 2000) 
                   # browser()
                  idUpdate = runif(1) < .idPropProb 
                  if(idUpdate)
                    .updType[.iteration] <<- 0
                  else
                    .updType[.iteration] <<- 1
                  if (is.null(block)) {
                    if (idUpdate) {
                      .proposal <<- .current + sqrt(.idVarScale / .nPar) * rnorm(.nPar)
                    } else {
                      .proposal <<- .current + sqrt(.varScale / .nPar) * rnorm(.nPar) %*% .chVar
                    }
                  } else {
                    bind = which(.parBlocks == block)
                    if (idUpdate) {
                      .proposal[bind] <<- .current[bind] +  sqrt(.idVarScale / length(bind)) * rnorm(length(bind))
                    } else { 
                      .proposal[bind] <<- .current[bind] +  sqrt(.varScale / length(bind)) * rnorm(length(bind)) %*% .chVar[bind, bind]
                    }
                  }
                  .propSet <<- TRUE
                }
                callSuper()
              }
            ),
            contains = "adMCMC"
)
