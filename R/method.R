setMethod("initialize", "RankData",
	function(.Object, ...){
        arg = list(...)
		fields = names(arg)
        # init ranking
        if("ranking" %in% fields){
            .Object@ranking = arg[["ranking"]]
        } else if ("ordering" %in% fields){
            .Object@ranking = OrderingToRanking(arg[["ordering"]])
        } else {
          stop("Either ordering or ranking matrix should be given")
        }
        # init nobj
        if("nobj" %in% fields){
            .Object@nobj = arg[["nobj"]]
        } else {
            .Object@nobj = max(.Object@ranking)
        }
        # init count
        if("count" %in% fields){
            .Object@count = arg[["count"]]
        } else {
            .Object@count = rep(1,nrow(.Object@ranking))
        }
        # init nobs
        if("nobs" %in% fields){
            .Object@nobs = arg[["nobs"]]
            stopifnot(.Object@nobs==sum(.Object@count))
        } else {
            .Object@nobs = sum(.Object@count)
        }
        # init ndistinct
        .Object@ndistinct = nrow(.Object@ranking)
        # handle topq
        if ("topq" %in% fields){  # the field is given
            if (min(arg[["topq"]]) < .Object@nobj-1){  # true topq case
                if (any(arg[["topq"]]>=.Object@nobj)){
                    warning("topq value should range between 1 and nobs-1")
                    arg[["topq"]][arg[["topq"]]>=.Object@nobj] = .Object@nobj-1
                }
                max_rank = apply(.Object@ranking,1,max)
                max_rank = as.numeric(max_rank) - 1
                if (!setequal(max_rank,arg[["topq"]])){
                    warning("The supplied top-q vector is not valid")
                    .Object@topq = unique(max_rank)
                } else {
                    .Object@topq = arg[["topq"]]
                }
                .Object@q_ind = c(1,cumsum(as.numeric(table(max_rank)[as.character(.Object@topq)]))+1)
                .Object@subobs = numeric(length(arg[["topq"]]))
                for (i in 1:length(arg[["topq"]])){
                    .Object@subobs[i] = sum(.Object@count[ .Object@q_ind[i]: (.Object@q_ind[i+1]-1) ])
                }
            } else {
                .Object@q_ind = -1
                .Object@subobs = -1
            }
        } else {
            .Object@topq = .Object@nobj - 1
            .Object@q_ind = -1
            .Object@subobs = -1
        }
        .Object
	}
)




setGeneric("SingleClusterModel",
        def=function(dat,init,ctrl,modal_ranking){standardGeneric("SingleClusterModel")}
)

setMethod(
  "SingleClusterModel",
  signature = c("RankData", "RankInit", "RankControlWtau"),
  definition = function( dat, init, ctrl, modal_ranking){
    dset = data.frame(dat@ranking, dat@count)
    nitem <- ncol(dset)-1
    
    test <- matrix(data = 0, nrow = factorial(nitem), ncol = nitem, byrow = TRUE)
    temp1 <- 1:nitem
    i <- 1
    w <- rep(1,nitem)
    modal <- 1:nitem
    
    ## generate a list of all possible rankings
    for (j in 1:(nitem^nitem-1)){
      temp1[1] <- nitem - j%%nitem
      temp2 <- j - j%%nitem
      for (k in nitem:2){
        temp1[k] <- nitem - temp2%/%(nitem^(k-1))
        temp2 <- temp2 - (nitem-temp1[k])*(nitem^(k-1))
      }
      temp2 <- 0
      for (l in 1:nitem){
        for (m in 1:nitem){
          if (temp1[l] == temp1[m] && l != m){
            temp2 <- 1
          }
        }
      }
      if (temp2 == 0){
        for (p in 1:nitem){
          test[i,p] = temp1[p]
        }
        i <- i+1
      }
    }
    
    n <- rep(0,factorial(nitem))
    for (j in 1:factorial(nitem)){
      for (k in 1:nrow(dset)){
        temp_ind <- 0
        for (l in 1:nitem){
          if (test[j,l] != dset[k,l]) {temp_ind <- temp_ind + 1}
        }
        if (temp_ind == 0) {n[j] <- dset[k,nitem+1]}
      }
    }
    test2 <- cbind(test, n)
    
    wdist <- function(x,y,w){
      d <- 0
      for (j in 1:(nitem-1)){
        for (k in (j+1):nitem){
          if ((x[j]-x[k])*(y[j]-y[k]) < 0) {d <- d+w[modal[k]]*w[modal[j]]}
        }
      }
      d
    }
    
    tempw <- rep(1,nitem)
    td <- rep(0,factorial(nitem))
    for (j in 1:factorial(nitem)){
      td[j] <- 0
      for (k in 1:factorial(nitem)){
        td[j] <- td[j] + n[k]*wdist(test[j,],test[k,],tempw)
      }
    }
    test3 <- cbind(test2, td)
    
    ## compute temporary modal ranking
    temp1 <- max(td)
    for (j in 1:factorial(nitem)){
      if (td[j] == min(td)){
        for (k in 1:nitem){
          modal[k] <- test3[j,k]
        }
      }
    }
    
    temp_modal <- 1:nitem
    mup <- 1
    while (mup == 1){
      ## loglikelihood function
      loglik_wdbm <- function(lambda){
        ed <- rep(0,factorial(nitem))
        for (j in 1:factorial(nitem)){
          ed[j] <- exp(-wdist(modal,test3[j,1:nitem],lambda))
        }
        pc <- sum(ed)
        pr <- rep(0,factorial(nitem))
        for (j in 1:factorial(nitem)){
          pr[j] <- ed[j]/pc
        }
        ll <- rep(0,factorial(nitem))
        for (j in 1:factorial(nitem)){
          ll[j] <- -log(pr[j])*test3[j,(nitem+1)]
        }
        sum(ll)
      }
      
      up <- optim(rep(1,nitem), loglik_wdbm, NULL, method = "BFGS", hessian = TRUE)
      w_up <- up$par
      for (j in 1:nitem){
        temp_modal[j] <- modal[j]
      }
      
      twd <- rep(0,factorial(nitem))
      for (j in 1:factorial(nitem)){
        twd[j] <- 0
        for (k in 1:factorial(nitem)){
          twd[j] <- twd[j] + n[k]*wdist(test[j,],test[k,],w_up)
        }
      }
      test3 <- cbind(test2, twd)
      
      temp1 <- max(td)
      for (j in 1:factorial(nitem)){
        if (twd[j] == min(twd)){
          for (k in 1:nitem){
            modal[k] <- test3[j,k]
          }
        }
      }
      mup <- 0
      # loop until tmp modal and modal is the same
      for (j in 1:nitem){
        if (temp_modal[j] != modal[j]) {mup <- 1}
      }
    }
    
    list(
      param.est = wToparam(up$par) ,log_likelihood = -1*up$value, pi0.ranking = modal
    )
    
  }
)







# single cluster model method for Weighted Kendall Distance
setMethod(
  "SingleClusterModel",
  signature = c("RankData","RankInit","RankControlWeightedKendall"),
  definition = function(dat,init,ctrl,modal_ranking) {
    param_len = max(dat@topq)
    param.coeff = CWeightGivenPi(dat@ranking,modal_ranking)
    param.coeff = matrix(param.coeff,ncol = dat@ndistinct,byrow = TRUE) %*%
      dat@count
    param.coeff = as.numeric(param.coeff)[1:param_len]
    
    if (length(dat@topq) == 1 && dat@topq == dat@nobj-1) {
        
      obj = function(param) {
        a = -1 * param %*% param.coeff - dat@nobs * LogC(c(param,rep(0,dat@nobj -
                                                                       1 - param_len)))
        as.numeric(-1 * a)
      }
      tt = t.gen(param_len)
      gradiant = function(param) {
        grad = GHC(param,tt)
        dat@nobs * grad + param.coeff
      }
      opt_res = optimx::optimx(
        par = init@param.init[[init@clu]][1:param_len],fn = obj,gr = gradiant,lower =
          rep(0,param_len),upper = rep(Inf,param_len),method = "L-BFGS-B",control =
          ctrl@optimx_control
      )
      param.est = unlist(opt_res[1:param_len])
      log_likelihood = -1 * opt_res[[param_len + 1]]
    } else {
      obj = function(param) {
        norm_vec = numeric(length(dat@topq))
        for (i in 1:length(dat@topq)) {
          j = dat@topq[i]
          norm_vec[i] = LogC(c(param[1:j],rep(0,dat@nobj-1 - j)))
        }
        a = -1 * param %*% param.coeff - dat@subobs %*% norm_vec
        as.numeric(-1 * a)
      }
      opt_res = optimx::optimx(
        par = init@param.init[[init@clu]][1:param_len],fn = obj,lower = rep(0,param_len),upper =
          rep(Inf,param_len),method = "L-BFGS-B",control = ctrl@optimx_control
      )
      param.est = unlist(opt_res[1:param_len])
      log_likelihood = -1 * obj(param.est)
    }
    param.est = c(param.est,rep(0,dat@nobj - 1 - param_len))
    list(
      param.est = param.est,w.est = paramTow(param.est),log_likelihood = log_likelihood
    )
  }
)

# remove this and uncomment above function
# setMethod("SingleClusterModel",
    # signature = c("RankData","RankInit","RankControlWeightedKendall"),
    # definition = function(dat,init,ctrl,modal_ranking){
        # param.coeff = CWeightGivenPi(dat@ranking,modal_ranking)
        # param.coeff = matrix(param.coeff,ncol = dat@ndistinct,byrow = TRUE)%*%dat@count
        # param.coeff = as.numeric(param.coeff)
        # param_len = dat@nobj-1
        # if (dat@topq>0){
            # param.coeff = param.coeff[1:dat@topq]
            # param_len = dat@topq
        # }
		# ind = c(1,6,26,86,206)
		# count_vec = numeric(4)
		# for (i in 1:4){
			# count_vec[i] = sum(dat@count[ind[i]:(ind[i+1]-1)])
		# }
        # obj = function(param){
			# norm_vec = numeric(param_len)
			# for (i in 1:param_len){
				# norm_vec[i] = LogC(c(param[1:i],rep(0,param_len-i)))
			# }
            # a = -1*param%*%param.coeff - count_vec%*%norm_vec
            # as.numeric(-1*a)
        # }

		# opt_res = optimx::optimx(par=init@param.init[[init@clu]][1:param_len],fn=obj,lower=rep(0,param_len),upper=rep(Inf,param_len),method="L-BFGS-B",control=ctrl@optimx_control)
		# param.est = unlist(opt_res[1:param_len])
		# log_likelihood=-1*opt_res[[param_len+1]]

		# param.est = c(param.est,rep(0,dat@nobj-1-param_len))
        # list(param.est=param.est,w.est=paramTow(param.est),log_likelihood=log_likelihood)
    # }
# )

# single cluster model method for Kendall distance
setMethod("SingleClusterModel",
    signature = c("RankData","RankInit","RankControlKendall"),
    definition = function(dat,init,ctrl,modal_ranking){
        param.coeff <- FindV(dat@ranking,modal_ranking)
        param.coeff <- rowSums(param.coeff)%*%dat@count
        param.coeff <- as.numeric(param.coeff)
        
        obj <- function(param){
            param*param.coeff + dat@nobs*LogC_Component(rep(param,dat@nobj-1))
        }
        
        opt_res = stats::optimize(f=obj,interval =c(0,100))
        list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
        }
)

# single cluster model method for Phi Component Model
setMethod("SingleClusterModel",
    signature = c("RankData","RankInit","RankControlPhiComponent"),
    definition = function(dat,init,ctrl,modal_ranking){
        param.coeff <- FindV(dat@ranking,modal_ranking)
        param.coeff <- t(param.coeff)%*%dat@count
        param.coeff <- as.numeric(param.coeff)
        param_len <- dat@nobj-1
        obj <- list()
        opt_res <- list()
        for (i in 1:param_len){
            obj[[i]] <- function(param){
                rhs <- exp(-param)/(1-exp(-param))-(dat@nobj+1-i)*exp(-(dat@nobj+1-i)*param)/(1-exp(-(dat@nobj+1-i)*param))
                ret <- param.coeff[i] - dat@nobs* rhs
                ret
            }
            opt_res[[i]] <- stats::uniroot(f=obj[[i]],interval=c(0,100),f.lower=-Inf)
        }
        param.est <- vapply(opt_res,function(x)x$root,numeric(1))
        log_likelihood <- -param.coeff%*%param.est - dat@nobs*LogC_Component(param.est)
        list(param.est=param.est,log_likelihood=log_likelihood)
        }
)

setGeneric("FindProb",
        def=function(dat,ctrl,modal_ranking,param){standardGeneric("FindProb")}
)


setMethod("FindProb",
        signature=c("RankData","RankControlWeightedKendall"),
        definition = function(dat,ctrl,modal_ranking,param){
            distance = param %*% matrix(CWeightGivenPi(dat@ranking,modal_ranking),ncol = dat@ndistinct,byrow = TRUE)
            if (length(dat@topq) == 1 && dat@topq == dat@nobj-1) {
                C = exp(LogC(param))
                prob = exp(-1*distance)/C
                if(dat@topq>0){
                    prob = prob*factorial(dat@nobj-dat@topq)
                }
            } else {
                cond_prob = dat@subobs/dat@nobs
                prob = exp(-1*distance)
                for (i in 1:length(dat@topq)) {
                    j = dat@topq[i]
                    norm_c = exp(LogC(c(param[1:j],rep(0,length(param) - j))) - lgamma(dat@nobj-j+1))
                    prob[dat@q_ind[i]:(dat@q_ind[i+1]-1)] = prob[dat@q_ind[i]:(dat@q_ind[i+1]-1)]/norm_c*cond_prob[i]
                }
                
            }
            prob
        }
)

# tmp version
# TODO: delete this and uncomment the function above
# setMethod("FindProb",
#         signature=c("RankData","RankControlWeightedKendall"),
#         definition = function(dat,ctrl,modal_ranking,param){
#             distance = param %*% matrix(CWeightGivenPi(dat@ranking,modal_ranking),ncol = dat@ndistinct,byrow = TRUE)
#             prob = exp(-1*distance)
# 			ind = c(1,6,26,86,206)
# 			cond_prob = c(0.3327723,0.1593631,0.1364490,0.3714156)
# 			for (i in 1:4){
# 				prob[ind[i]:(ind[i+1]-1)] = prob[ind[i]:(ind[i+1]-1)]/sum(prob[ind[i]:(ind[i+1]-1)])*cond_prob[i]
# 			}
#             prob
#         }
# )
setMethod("FindProb",
        signature=c("RankData","RankControlKendall"),
        definition = function(dat,ctrl,modal_ranking,param){
            param = param[1]
            distance = FindV(dat@ranking,modal_ranking) %*% rep(param,dat@nobj-1)
            C = exp(LogC_Component(rep(param,dat@nobj-1)))
            prob = exp(-1*distance)/C
            prob
        }
)


setMethod("FindProb",
        signature=c("RankData","RankControlPhiComponent"),
        definition = function(dat,ctrl,modal_ranking,param){
            distance = FindV(dat@ranking,modal_ranking) %*% param
            C = exp(LogC_Component(param))
            prob = exp(-1*distance)/C
            prob
        }
)


setMethod("FindProb",
          signature=c("RankData","RankControlWtau"),
          definition = function(dat,ctrl,modal_ranking,param){
            dset = data.frame(dat@ranking, dat@count)
            nitem = dat@nobj
            test <- matrix(data = 0, nrow = factorial(nitem), ncol = nitem, byrow = TRUE)
            temp1 <- 1:nitem
            i <- 1
            wdist <- function(x,y,w){
              d <- 0
              for (j in 1:(nitem-1)){
                for (k in (j+1):nitem){
                  if ((x[j]-x[k])*(y[j]-y[k]) < 0) {d <- d+w[modal_ranking[k]]*w[modal_ranking[j]]}
                }
              }
              d
            }
            
            for (j in 1:(nitem^nitem-1)){
              temp1[1] <- nitem - j%%nitem
              temp2 <- j - j%%nitem
              for (k in nitem:2){
                temp1[k] <- nitem - temp2%/%(nitem^(k-1))
                temp2 <- temp2 - (nitem-temp1[k])*(nitem^(k-1))
              }
              temp2 <- 0
              for (l in 1:nitem){
                for (m in 1:nitem){
                  if (temp1[l] == temp1[m] && l != m){
                    temp2 <- 1
                  }
                }
              }
              if (temp2 == 0){
                for (p in 1:nitem){
                  test[i,p] = temp1[p]
                }
                i <- i+1
              }
            }
            
            n <- rep(0,factorial(nitem))
            for (j in 1:factorial(nitem)){
              for (k in 1:nrow(dset)){
                temp_ind <- 0
                for (l in 1:nitem){
                  if (test[j,l] != dset[k,l]) {temp_ind <- temp_ind + 1}
                }
                if (temp_ind == 0) {n[j] <- dset[k,nitem+1]}
              }
            }
            test2 <- cbind(test, n)
            
            twd <- rep(0,factorial(nitem))
            for (j in 1:factorial(nitem)){
              twd[j] <- 0
              for (k in 1:factorial(nitem)){
                twd[j] <- twd[j] + n[k]*wdist(test[j,],test[k,],param)
              }
            }
            test3 <- cbind(test2, twd)
            ed <- rep(0,factorial(nitem))
            for (j in 1:factorial(nitem)){
              ed[j] <- exp(-wdist(modal_ranking,test3[j,1:nitem],param))
            }
            pc <- sum(ed)
            fitted <- rep(0,factorial(nitem))
            for (j in 1:factorial(nitem)){
              fitted[j] <- ed[j]/pc * sum(test3[,(nitem+1)])
            }
            prob = fitted/dat@nobs
            prob
          }
)















