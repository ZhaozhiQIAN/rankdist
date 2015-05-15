#' @useDynLib rankdist
#' @importFrom Rcpp sourceCpp
NULL


#' @title RankData Class
#' @description A S4 class to represent ranking data
#' 
#' It is well understood that the ranking representation and ordering representation of ranking data can easily get confused.
#' I thus use a S4 class to store all the information about the ranking data. This can avoid unnessary confusions.
#' 
#' @slot nobj an number to store number of ranked objects.
#' @slot nobs the number of observations.
#' @slot ndistinct the number of distinct rankings.
#' @slot ordering a matrix that stores the ordering representation of distinct rankings. Each row contains one ordering.
#' @slot ranking a matrix that stores the ranking representation of distinct rankings. Each row contains one ranking.
#' @slot count the number of observations for each distinct ranking.
#' @slot topq an integer to store whether the data set is top-q ranking. If value is -1, the data is complete rankings. 
#' Otherwise, the value is number of ranked object (the q in top-q ranking).
#' @examples
#' # creating a random rank data set
#' rankmat = replicate(10,sample(1:52,52), simplify = "array")
#' ordermat = OrderingToRanking(rankmat)
#' countvec = sample(1:52,52,replace=TRUE)
#' rankdat = new("RankData",nobj=52,nobs=sum(countvec),ndistinct=10,
#'      ranking=rankmat,ordering=ordermat,count=countvec)
#' @aliases RankData RankData-class
#' @seealso \code{\link{RankInit}}, \code{\link{RankControl}}
#' @export
setClass( "RankData",
          representation = representation(
              nobj = "numeric",
              nobs = "numeric",
              ndistinct = "numeric",
              ordering = "matrix",
              ranking = "matrix",
              count = "numeric",
              topq = "numeric" # if topq >0, else = -1
          ),
          prototype = prototype(
              topq = -1L
          )
)


#' @title RankInit Class
#' @description A S4 class to store initialization information of model fitting
#' 
#' The \code{RankInit} class is used to give initial values of model fitting procedures.
#' 
#' @slot fai.init a list containing initial values of the fai parameterization of weights.
#' @slot modal_ranking.init a list containing starting points for the modal ranking search.
#' @slot clu an integer containing the number of clusters used in the model.
#' @slot p.init a numeric vector containing the initial values for cluster probabilities.
#' @examples c1init = new("RankInit",fai.init=list(rep(1,4)),
#'      modal_ranking.init=list(c(2,3,4,1,5)),clu=1L)
#' c2init = new("RankInit",fai.init=list(rep(0.1,4),rep(0.1,4)),
#'      modal_ranking.init = list(c(2,3,4,1,5),c(2,5,1,4,3)),clu=2L,p.init=c(0.5,0.5))
#' @aliases RankInit RankInit-class
#' @seealso \code{\link{RankData}}, \code{\link{RankControl}}
#' @export
setClass( "RankInit",
          representation = representation(
              fai.init = "list",
              modal_ranking.init = "list",
              clu = "integer",
              p.init = "numeric"
          ),
          prototype = prototype(
              clu = 1L,
              p.init = 1
          )
)

#' @title RankControl Class
#' @description A S4 class to store control parameters for model fitting
#' 
#' @slot EM_limit maximum number of EM iteration
#' @slot EM_epsilon convergence error for weights and cluster probabilities in EM iteration
#' @slot SearchPi0_limit maximum number of iterations in the local seach of pi0.
#' @slot SearchPi0_FUN a function object that gives a goodness of fit criterion. The default is log likelihood.
#' @slot SearchPi0_fast_traversal a logical value. If TRUE (by default), immediately traverse to the neighbour if it is better than the current pi0. Otherwise, check all neighbours and traverse to the best one.
#' @slot SearchPi0_show_message a logical value. If TRUE, the location of the current pi0 is shown.
#' @slot SearchPi0_neighbour a character string specifying which type of neighbour to use in the local search. Supported values are: "Cayley" to use neighbours in terms of Cayley distance or "Kendall" to use neighbours in terms of Kendall distance.
#' Note that Kendall neighbours are a subset of Cayley neighbours
#' @slot optimx_control a list to be passed to \code{\link[optimx]{optimx}}. The list must not contain a component \code{maximize=TRUE} since internally the negation of the likelihood function is minimized.
#' @details The control parameters that start with prefix \code{EM_} are intended for the EM iteration. The ones with prefix \code{SeachPi0} control the behaviour of searching model ranking.
#' The function object specified in the \code{SearchPi0_control} takes a list as argument. The components in the list include the following. \code{obs}: the number of observations.
#' \code{w.est}: the estimated weights. \code{log_likelihood}: the estimated log_likelihood. With this information, most of the popular information criterions can be supported and customized criterions can also be defined.
#' A larger returned value indicates a better fit. 
#' @examples # enabling messages and warnings
#' testctrl = new("RankControl",SearchPi0_show_message=TRUE, optimx_control=list(dowarn=TRUE))
#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}
#' @aliases RankControl RankControl-class
#' @export
setClass( "RankControl",
          representation = representation(
              EM_limit = "numeric",
              EM_epsilon = "numeric",
              SearchPi0_limit = "numeric",
              SearchPi0_FUN = "function",
              SearchPi0_fast_traversal = "logical",
              SearchPi0_show_message = "logical",
              SearchPi0_neighbour = "character",
              optimx_control = "list"
          ),
          prototype = prototype(
              EM_limit=500,
              EM_epsilon=1e-5,
              SearchPi0_limit=500,
              SearchPi0_FUN = function(x){x$log_likelihood},
              SearchPi0_fast_traversal=TRUE,
              SearchPi0_show_message=FALSE,
              SearchPi0_neighbour="Cayley",
              optimx_control = list(maximize=FALSE,starttests=TRUE,trace=0,dowarn=FALSE)
          )
)

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
#' \item{\code{fai.est}}{the fai parameterization of weights of each cluster.}
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
    fai <- list()
    fai.last <- list()
    func.call <- match.call()
    loopind=0
    
    if (clu > 1L){  # use EM to fit multicluster model
        # further initialization
        modal_ranking.est <- init@modal_ranking.init
        fai <- init@fai.init
        init.clu = list()
        for (i in 1:clu){
            init.clu[[i]] <- new("RankInit",fai.init=list(fai[[i]]),modal_ranking.init=list(modal_ranking.est[[i]]),clu=1L)
        }
        
        while(TRUE){
            loopind <- loopind+1
            # E step
            z <- matrix(nrow = distinctn,ncol = clu)
            for (i in 1:clu){
                z[,i] <- p[i]*FindProb(dat,modal_ranking.est[[i]],fai[[i]])
            }
            sums <- rowSums(z)
            z <- z/sums
            
            # M step
            p <- t(z) %*% count
            p <- p/sum(p)
            if (any(p<0.03)){
                warning("One cluster has probability smaller than 3%; Try fewer clusters")
                return (list(p=p,modal_ranking.est=modal_ranking.est,fai=fai,iteration=loopind))
            }
            for ( i in 1:clu){
                dat.clu <- dat
                dat.clu@count <- z[,i] * count
                dat.clu@nobs <- sum(z[,i] * count)
                # need change 
                init.clu[[i]]@fai.init <- list(fai[[i]])
                init.clu[[i]]@modal_ranking.init <- list(modal_ranking.est[[i]])
                
                solve.clu <- SearchPi0(dat.clu,init.clu[[i]],ctrl)
                modal_ranking.est[[i]] <- solve.clu$pi0.ranking
                fai[[i]] <- solve.clu$fai.est
                p_clu[i,] <- FindProb(dat,modal_ranking.est[[i]],fai[[i]])*p[i]
            }
            log_likelihood_clu <- sum(log(colSums(p_clu))%*%dat@count)
            # break?
            if(loopind != 1){
                if (loopind < ctrl@EM_limit) {
                    cond1 <- all( p.last - p < ctrl@EM_epsilon)
                    cond2 <- all( unlist(modal_ranking.est.last) - unlist(modal_ranking.est) < ctrl@EM_epsilon)
                    cond3 <- all( unlist(fai.last) - unlist(fai) < ctrl@EM_epsilon)
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
            fai.last <- fai
            log_likelihood_clu.last <- log_likelihood_clu
        } # inf loop
        
        est_prob <- colSums(p_clu)
    } else {  # fit single cluster model
        sigle_cluster_mod <- SearchPi0(dat,init,ctrl)
        modal_ranking.est <- list(sigle_cluster_mod$pi0.ranking)
        p <- 1
        fai <- list(sigle_cluster_mod$fai.est)
        log_likelihood_clu.last <- sigle_cluster_mod$log_likelihood
        est_prob <- FindProb(dat,modal_ranking.est[[1]],fai[[1]])*p[1]
    }
    # finishing up with parameter estimation and summary statistics
    p <- as.numeric(p)
    v <- length(p) - 1 + sum(unlist(fai)!=0)
    bic <- 2*log_likelihood_clu.last - v*log(n)
    est_prob <- colSums(p_clu)
    expectation <- as.numeric(est_prob*n)
    SSR <- sum((expectation-dat@count)^2/expectation)
    ret <- list(p=p,modal_ranking.est=modal_ranking.est,w.est=lapply(fai,faiTow),fai.est=fai,SSR=SSR,log_likelihood=log_likelihood_clu.last,free_params=v,BIC=bic,expectation=expectation,iteration=loopind,model.call=func.call)
    ret
}


#' Find Initial Values of Fai
#'
#' \code{MomentEst} finds the initial values of fai which can be used in the subsequent optimization problems.
#' Linear model is fitted to the log odds of rankings.
#' 
#' @param dat a RankData object
#' @param size the number of samples to take in the linear model
#' @param pi0 an optional argument showing the location of central ranking. 
#' If not provided, Borda Count method is used to estimate the central ranking.
#' @return estimated fai
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
    nfai <- dat@nobj - 1
    # construct equation set
    prob_obs <- dat@count/dat@nobs
    # row-wise
    distance_mat <- matrix(CWeightGivenPi(dat@ranking,modal_ranking),nrow = dat@ndistinct)
    pair <- sample(1:nrow(distance_mat),2*size,replace=TRUE)
    pair_mat <- matrix(pair,nrow=2)
    logodd <- numeric(size)
    design_mat <- matrix(ncol=nfai, nrow=size)
    for (i in 1:size){
        logodd[i] <- prob_obs[pair_mat[1,i]]/prob_obs[pair_mat[2,i]]
        design_mat[i,] <- distance_mat[pair_mat[1,i],]-distance_mat[pair_mat[2,i],]
    }
    fai.est <- lm(logodd~design_mat-1)
    fai.est <- fai.est$coefficients
    fai.est[fai.est<0] <- 0
    names(fai.est) <- NULL
    fai.est
}

#' American Psychological Association (APA) election data
#'
#' A dataset containing 5738 complete votes in APA election. There are 5 candidates in total.
#'
#' @format a RankData object
#' 
#' @references Marden, J. I. (1995). Analyzing and Modeling Rank Data (94-96). Chapman Hall, New York.
"apa_obj"


