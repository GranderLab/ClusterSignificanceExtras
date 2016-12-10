##Functions for one line variable declaration
##http://stackoverflow.com/questions/7519790/assign-multiple-new-variables-in-a-single-line-in-r

#' Grouping the left hand side
#'
#'
#'
#' @name g
#' @rdname g
#' @aliases g
#' @param ... argument 1
#' @author Jason Serviss
#' @keywords g
#' @examples
#'
#'
#' cat("no current examples")
#'
#'
NULL
#' @export

g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

#' Binary Operator
#'
#'
#'
#' @name %=%.lbunch
#' @rdname MultipleAssignmentHelper
#' @aliases %=%.lbunch
#' @param l argument 1
#' @param r argument 2
#' @param ... argument 3
#' @author Jason Serviss
#' @keywords internal
#' @examples
#'
#'
#' cat("no current examples")
#'
#'
NULL
#' @export

'%=%.lbunch' = function(l, r, ...) {
    Envir = as.environment(-1)
    
    if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
    
    if (length(l) > length(r))  {
        warning("LHS has more args than RHS. RHS will be repeated.")
        r <- extendToMatch(r, l)
    }
    
    for (II in 1:length(l)) {
        do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
    }
}

#' Used if LHS is larger than RHS
#'
#'
#'
#' @name extendToMatch
#' @rdname extendToMatch
#' @aliases extendToMatch
#' @param l argument 1
#' @param r argument 2
#' @author Jason Serviss
#' @keywords extendToMatch
#' @examples
#'
#'
#' cat("no current examples")
#'
#'
NULL
#' @export

extendToMatch <- function(source, destin) {
    s <- length(source)
    d <- length(destin)
    
    # Assume that destin is a length when it is a single number and source is not
    if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
    
    dif <- d - s
    if (dif > 0) {
        source <- rep(source, ceiling(d/s))[1:d]
    }
    return (source)
}

#' Generic form
#'
#'
#'
#' @name %=%
#' @rdname MultipleAssignment
#' @aliases %=%
#' @param ... argument 1
#' @author Jason Serviss
#' @keywords internal
#' @examples
#'
#'
#' cat("no current examples")
#'
#'
NULL
#' @export

'%=%' = function(l, r, ...) UseMethod('%=%')


##automatic handling of ellipsis
.ellipsis <- function(...) {
    if(length(list(...)) > 0) {
        e=parent.frame()
        arguments <- list(...)
        for(i in 1:length(arguments)){
            name <- names(arguments[i])
            assign(name, arguments[[i]], envir=e)
        }
    }
}


#' Parallelize ClusterSignificance permutation
#'
#' Only works, for now, with 2 groups.
#'
#' @name parallelCS
#' @rdname parallelCS
#' @aliases parallelCS
#' @param mat matrix to ClusterSignificance
#' @param groups groups to ClusterSignificance
#' @param iter iter to ClusterSignificance
#' @param projMethod projMethod to ClusterSignificance
#' @param cores the number of cores to run parallel
#' @author Jason Serviss
#' @keywords parallelCS
#' @examples
#'
#'
#' cat("no current examples")
#'
#'
NULL
#' @export
#' @importFrom doMC registerDoMC
#' @importFrom foreach foreach %dopar% registerDoSEQ

parallelCS <- function(
    mat,
    groups,
    iter,
    projmethod,
    cores = 2,
    user.permutations = NULL
) {
    registerDoMC(cores)
    ceiling(iterations <- iter / cores)
    completed <- 0
    count <- 0
    
    while(TRUE) {
        
        loopOutput <- foreach(i = 1:cores, .combine = c) %dopar% {
            output <- permute(
                mat = mat,
                groups = groups,
                iter = iterations,
                projmethod = projmethod,
                user.permutations=user.permutations
            )
        }
        
        if(count == 0) {
            output <- loopOutput
        } else {
            output <- c(output, loopOutput)
        }
        
        completed <- length(getData(output, "scores.vec")[[1]])
        if(completed >= iter) {break}
        
        iterations <- ceiling((iterations - completed) / cores)
        message("adding iterations")
        count <- count + 1
    }
    output <- .trimIterations(output, iter)
    
    registerDoSEQ()
    return(output)
}

# if length(scores.vec) > iter
.trimIterations <- function(output, iter) {
    scores.vec <- getData(output, "scores.vec")[[1]]
    
    if(length(scores.vec) == iter) {
        return(output)
    } else {
        message("trimming iterations")
        output <- methods::new(
            "PermutationResults",
            scores.real=getData(output, "scores.real"),
            scores.vec=list(scores.vec[1:iter])
        )
        return(output)
    }
}



