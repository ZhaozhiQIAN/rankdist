#' @useDynLib rankdist
#' @importFrom Rcpp sourceCpp
NULL


#' Create Hash Value for Rank 
#'
#' Sometimes it is handier to deal with rankings into a hash value. \code{RanktoHash} returns hash values for ranks. Maximum 52 objects are supported.
#' @param r  a vector or matrix of rankings. Each row of the matrix represents a ranking.
#'   The ranking should be a integer from one to number of objects. No NA is allowed
#' @return a vector of character strings representing the hash values.
#' @seealso \code{\link{HashtoRank}} for a reverse operation.
#' @export
RanktoHash <- function(r){
    char <- c(letters,LETTERS)
    if (!is.matrix(r)){
        a <- paste(char[r],collapse="")
    } else {
        a <- apply(r,1,function(x)paste(char[x],collapse=""))
    }
    as.vector(a,mode="character")
}


#' Obtain Ranking from Hash Value
#'
#' \code{HashToRank} returns rankings from given hash values. Maximum 52 objects are supported.
#' @param h  A vector of hash values. 
#' @return a matrix of rankings if input has more than one element or a vector of rankings if input has only one element
#' @seealso \code{\link{RanktoHash}} for a reverse operation.
#' @export
HashtoRank <- function(h){
    char <- c(letters,LETTERS)
    transvec <- seq_along(char)
    names(transvec) <- char
    nobj <- nchar(h[1])
    numrank <- transvec[unlist(strsplit(h,split=""))]
    if (length(h)==1){
        as.vector(numrank,mode="numeric")
    } else {
        matrix(data=numrank,ncol=nobj,byrow=TRUE)
    }
    
}


#' Transformation between Rankings and Orderings
#' 
#' \code{OrderingToRanking} transforms between ranking representation and ordering representation.
#'
#'    Ranking representation encodes the position of objects. Ordering representation is an ordered sequence of objects.
#'    For example ranking (2 3 1 4) is equivalent to ordering (3 1 2 4), which means object 3 is first, object 1 is second, followed by object 2 and 4.
#' @param ordering  a matrix of orderings or rankings. Each row contains an obseravtion.
#' @return  a matrix of transformed rankings or orderings. Each row contains an obseravtion.
#' @export
OrderingToRanking <- function(ordering){
    if (is.matrix(ordering)){
        ranking <- t(apply(ordering,1,order))
    } else {
        ranking <- order(ordering)
    }
    ranking
}



#' Fit A Mixture of Distance-based Models
#' 
#' \code{RankDistModel} fits a mixture of ranking models based on weighted Kendall distance.
#' 
#' The procedure will estimate central rankings, the probability of each cluster and weights.
#' 
#' @param dat A \linkS4class{RankData} object.
#' @param init A \linkS4class{RankInit} object.
#' @param ctrl A \linkS4class{RankControl} object.
#' 
#' @return A list containing the following components:
#' \describe{
#' \item{\code{modal_ranking.est}}{the estimated pi0 for each cluster.}
#' \item{\code{p}}{the probability of each cluster.}
#' \item{\code{w.est}}{the estimated weights of each cluster.}
#' \item{\code{param.est}}{the param parameterization of weights of each cluster.}
#' \item{\code{SSR}}{the sum of squares of pearson residuals}
#' \item{\code{log_likelihood}}{the fitted log_likelihood}
#' \item{\code{BIC}}{the fitted Bayesian Information Criterion value}
#' \item{\code{free_params}}{the number of free parameters in the model}
#' \item{\code{expectation}}{the expected value of each observation given by the model}
#' \item{\code{iteration}}{the number of EM iteration}
#' \item{\code{model.call}}{the function call}
#' }
#' @export
RankDistanceModel <- function(dat,init,ctrl){
    # get parameters
    tt <- dat@nobj # number of objects
    n <- dat@nobs    # total observations
    distinctn <- dat@ndistinct	# total distinct observations
    clu <- init@clu
    count <- dat@count
    p_clu <- matrix(ncol = distinctn,nrow = clu)
    # setting up return values
    modal_ranking.est <- list()
    modal_ranking.est.last <- list()
    p <- rep(1/clu,clu)
    p.last <- numeric(clu)
    param <- list()
    param.last <- list()
    func.call <- match.call()
    loopind=0
    
    if (clu > 1L){  # use EM to fit multicluster model
        # further initialization
        modal_ranking.est <- init@modal_ranking.init
        param <- init@param.init
        init.clu = list()
        for (i in 1:clu){
            init.clu[[i]] <- new("RankInit",param.init=list(param[[i]]),modal_ranking.init=list(modal_ranking.est[[i]]),clu=1L)
        }
        
        while(TRUE){
            loopind <- loopind+1
            # E step
            z <- matrix(nrow = distinctn,ncol = clu)
            for (i in 1:clu){
                z[,i] <- p[i]*FindProb(dat,ctrl,modal_ranking.est[[i]],param[[i]])
            }
            sums <- rowSums(z)
            z <- z/sums
            
            # M step
            p <- t(z) %*% count
            p <- p/sum(p)
            if (any(p<0.03)){
                warning("One cluster has probability smaller than 3%; Try fewer clusters")
                return (list(p=p,modal_ranking.est=modal_ranking.est,param=param,iteration=loopind))
            }
            for ( i in 1:clu){
                dat.clu <- dat
                dat.clu@count <- z[,i] * count
                dat.clu@nobs <- sum(z[,i] * count)
                # need change 
                init.clu[[i]]@param.init <- list(param[[i]])
                init.clu[[i]]@modal_ranking.init <- list(modal_ranking.est[[i]])
                
                solve.clu <- SearchPi0(dat.clu,init.clu[[i]],ctrl)
                modal_ranking.est[[i]] <- solve.clu$pi0.ranking
                param[[i]] <- solve.clu$param.est
                p_clu[i,] <- FindProb(dat,ctrl,modal_ranking.est[[i]],param[[i]])*p[i]
            }
            log_likelihood_clu <- sum(log(colSums(p_clu))%*%dat@count)
            # break?
            if(loopind != 1){
                if (loopind < ctrl@EM_limit) {
                    cond1 <- all( p.last - p < ctrl@EM_epsilon)
                    cond2 <- all( unlist(modal_ranking.est.last) - unlist(modal_ranking.est) < ctrl@EM_epsilon)
                    cond3 <- all( unlist(param.last) - unlist(param) < ctrl@EM_epsilon)
                    if (cond1 && cond2 && cond3){
                        break
                    }else if(abs(log_likelihood_clu.last - log_likelihood_clu)<0.01){
                        break
                    }
                } else {
                    print(paste("Algorithm did not converge in",ctrl@EM_limit,"iterations"))
                    return(NULL)
                }
            }
            p.last <- p
            modal_ranking.est.last <- modal_ranking.est
            param.last <- param
            log_likelihood_clu.last <- log_likelihood_clu
        } # inf loop
        
        est_prob <- colSums(p_clu)
    } else {  # fit single cluster model
        sigle_cluster_mod <- SearchPi0(dat,init,ctrl)
        modal_ranking.est <- list(sigle_cluster_mod$pi0.ranking)
        p <- 1
        param <- list(sigle_cluster_mod$param.est)
        log_likelihood_clu.last <- sigle_cluster_mod$log_likelihood
        est_prob <- FindProb(dat,ctrl,modal_ranking.est[[1]],param[[1]])*p[1]
    }
    # finishing up with parameter estimation and summary statistics
    p <- as.numeric(p)
    v <- length(p) - 1 + sum(unlist(param)!=0)
    bic <- 2*log_likelihood_clu.last - v*log(n)
    expectation <- as.numeric(est_prob*n)
    SSR <- sum((expectation-dat@count)^2/expectation)
    ret <- list(p=p,modal_ranking.est=modal_ranking.est,w.est=lapply(param,paramTow),param.est=param,SSR=SSR,log_likelihood=log_likelihood_clu.last,free_params=v,BIC=bic,expectation=expectation,iteration=loopind,model.call=func.call)
    ret
}


#' Find Initial Values of param
#'
#' \code{MomentEst} finds the initial values of param which can be used in the subsequent optimization problems.
#' Linear model is fitted to the log odds of rankings.
#' 
#' @param dat a RankData object
#' @param size the number of samples to take in the linear model
#' @param pi0 an optional argument showing the location of central ranking. 
#' If not provided, Borda Count method is used to estimate the central ranking.
#' @return estimated param
#' @examples MomentsEst(apa_obj,40)
#' @export
MomentsEst <- function(dat,size,pi0=NULL){
    # estimating pi0
    avg_rank <- dat@count %*% dat@ranking;
    if (is.null(pi0)){
        modal_ranking <- sort(dat@ordering[1,])[order(avg_rank)]
    } else {
        modal_ranking <- pi0
    }
    nparam <- dat@nobj - 1
    # construct equation set
    prob_obs <- dat@count/dat@nobs
    # row-wise
    distance_mat <- matrix(CWeightGivenPi(dat@ranking,modal_ranking),nrow = dat@ndistinct)
    pair <- sample(1:nrow(distance_mat),2*size,replace=TRUE)
    pair_mat <- matrix(pair,nrow=2)
    logodd <- numeric(size)
    design_mat <- matrix(ncol=nparam, nrow=size)
    for (i in 1:size){
        logodd[i] <- prob_obs[pair_mat[1,i]]/prob_obs[pair_mat[2,i]]
        design_mat[i,] <- distance_mat[pair_mat[1,i],]-distance_mat[pair_mat[2,i],]
    }
    param.est <- lm(logodd~design_mat-1)
    param.est <- param.est$coefficients
    param.est[param.est<0] <- 0
    names(param.est) <- NULL
    param.est
}

#' American Psychological Association (APA) election data
#'
#' A dataset containing 5738 complete votes in APA election. There are 5 candidates in total.
#'
#' @format a RankData object
#' 
#' @references Marden, J. I. (1995). Analyzing and Modeling Rank Data (94-96). Chapman Hall, New York.
"apa_obj"

#' American Psychological Association (APA) election data (partial rankings included)
#'
#' A dataset containing 5738 complete votes and 9711 partial votes in APA election. There are 5 candidates in total.
#'
#' @format a RankData object
#' 
#' @references Marden, J. I. (1995). Analyzing and Modeling Rank Data (94-96). Chapman Hall, New York.
"apa_partial_obj"

