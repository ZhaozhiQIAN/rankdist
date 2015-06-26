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
              topq = "numeric", # a numeric vector: topq 
			  subobs = "numeric", # number of observations for each topq
			  q_ind = "numeric"  # starting point of each topq location
          ),
          prototype = prototype(
              topq = -1,
			  subobs = 1,
			  q_ind = 1
          )
)


#' @title RankInit Class
#' @description A S4 class to store initialization information of model fitting
#' 
#' The \code{RankInit} class is used to give initial values of model fitting procedures.
#' 
#' @slot param.init a list containing initial values of the param parameterization of weights.
#' @slot modal_ranking.init a list containing starting points for the modal ranking search.
#' @slot clu an integer containing the number of clusters used in the model.
#' @slot p.init a numeric vector containing the initial values for cluster probabilities.
#' @examples c1init = new("RankInit",param.init=list(rep(1,4)),
#'      modal_ranking.init=list(c(2,3,4,1,5)),clu=1L)
#' c2init = new("RankInit",param.init=list(rep(0.1,4),rep(0.1,4)),
#'      modal_ranking.init = list(c(2,3,4,1,5),c(2,5,1,4,3)),clu=2L,p.init=c(0.5,0.5))
#' @aliases RankInit RankInit-class
#' @seealso \code{\link{RankData}}, \code{\link{RankControl}}
#' @export
setClass( "RankInit",
          representation = representation(
              param.init = "list",
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
#' @description A virtual S4 class to store control parameters for model fitting. 
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
#' @details RankControl class must be extended to reflect what distance metric should be used. Possibles extendions are \code{\link{RankControlWeightedKendall}}. The control parameters that start with prefix \code{EM_} are intended for the EM iteration. The ones with prefix \code{SeachPi0} control the behaviour of searching model ranking.
#' @section User-defined Criterion:
#' You can specify user-defined criterion to choose modal rankings. The function object SearchPi0_FUN takes a list as argument. The components in the list include the following. \code{obs}: the number of observations.
#' \code{w.est}: the estimated weights. \code{log_likelihood}: the estimated log_likelihood. With this information, most of the popular information criterions can be supported and customized criterions can also be defined.
#' A larger returned value indicates a better fit. Note that if you are fitting a mixture model the EM algorthm always tries to maximized the log likelihood. Thus the default value should be used in this case.
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
              optimx_control = "list",
              "VIRTUAL"
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


#' @title RankControlWeightedKendall Class
#' @description A S4 class to store control parameters for Weighted Kendall distance model fitting. It is derived from class \code{\link{RankControl-class}}. 
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
#' @slot param_len a numeric value that should be set to (number of object - 1). It specifies the length of parameter.
#' @details \code{RankControlWeightedKendall} is derived from virtual class \code{\link{RankControl}}. All slots in \code{\link{RankControl}} are still valid. Only one slot \code{param_len} is added. 
#' This control class tells the solver to fit a model based on Weighted Kendall distance.  
#' The control parameters that start with prefix \code{EM_} are intended for the EM iteration. The ones with prefix \code{SeachPi0} control the behaviour of searching model ranking.
#' @examples # enabling messages and warnings
#' testctrl = new("RankControlWeightedKendall",SearchPi0_show_message=TRUE, optimx_control=list(dowarn=TRUE),param_len=9)
#' @seealso \code{\link{RankData}}, \code{\link{RankInit}}, \code{\link{RankControl}}
#' @aliases RankControlWeightedKendall RankControlWeightedKendall-class
#' @export
setClass( "RankControlWeightedKendall",
        representation = representation(
            param_len = "numeric"
        ),
        contains = "RankControl"
)


#' @export
setClass( "RankControlKendall",
        representation = representation(
            param_len = "numeric"
        ),
        contains = "RankControl"
)

#' @export
setClass( "RankControlPhiComponent",
        representation = representation(
            param_len = "numeric"
        ),
        contains = "RankControl"
)












