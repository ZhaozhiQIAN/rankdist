
setGeneric("SingleClusterModel",
        def=function(dat,init,ctrl,modal_ranking){standardGeneric("SingleClusterModel")}
)

# single cluster model method for Weighted Kendall Distance
setMethod("SingleClusterModel",
    signature = c("RankData","RankInit","RankControlWeightedKendall"),
    definition = function(dat,init,ctrl,modal_ranking){
        param.coeff = CWeightGivenPi(dat@ranking,modal_ranking)
        param.coeff = matrix(param.coeff,ncol = dat@ndistinct,byrow = TRUE)%*%dat@count
        param.coeff = as.numeric(param.coeff)
        param_len = dat@nobj-1
        if (dat@topq>0){
            param.coeff = param.coeff[1:dat@topq]
            param_len = dat@topq
        }
        obj = function(param){
            a = -1*param%*%param.coeff - dat@nobs*LogC(param)
            as.numeric(-1*a)
        }
        tt = t.gen(param_len)
        gradiant = function(param){
            grad = GHC(param,tt)
            dat@nobs*grad + param.coeff
        }
        opt_res = optimx::optimx(par=init@param.init[[init@clu]],fn=obj,gr=gradiant,lower=rep(0,param_len),upper=rep(Inf,param_len),method="L-BFGS-B",control=ctrl@optimx_control)
        param.est = unlist(opt_res[1:param_len])
        list(param.est=param.est,w.est=paramTow(param.est),log_likelihood=-1*opt_res[[param_len+1]])
        }
)

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
        
        opt_res = optimize(f=obj,interval =c(0,100))
        list(param.est=opt_res$minimum,log_likelihood=-1*opt_res$objective)
        }
)

# single cluster model method for Phi Component Model
# in dev
#' @export
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
            opt_res[[i]] <- uniroot(f=obj[[i]],interval=c(0,100),f.lower=-Inf)
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
            C = exp(LogC(param))
            prob = exp(-1*distance)/C
            prob
        }
)

setMethod("FindProb",
        signature=c("RankData","RankControlKendall"),
        definition = function(dat,ctrl,modal_ranking,param){
            distance = FindV(dat@ranking,modal_ranking) %*% rep(param,dat@nobj-1)
            C = exp(LogC_Component(rep(param,dat@nobj-1)))
            prob = exp(-1*distance)/C
            prob
        }
)

# in dev
#' @export
setMethod("FindProb",
        signature=c("RankData","RankControlPhiComponent"),
        definition = function(dat,ctrl,modal_ranking,param){
            distance = FindV(dat@ranking,modal_ranking) %*% param
            C = exp(LogC_Component(param))
            prob = exp(-1*distance)/C
            prob
        }
)

















