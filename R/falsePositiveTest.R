
#' False positive test
#'
#' This function performs the false positive test. This test is designed to
#' determine the false positive rate of the Pcp and Mlp methods. To do this
#' we utilize the hundred dataset which was derived using the \code{\link{zeroOrHundred}}
#' function with argument \code{overlap}=100. The dataset includes matrices where
#' the groups overlap is 100% and thus we anticipate that no significant
#' seperations should be found with Pcp or Mlp. For each matrice, the
#' permute function is run, using both the Pcp and Mlp methods afterwhich,
#' the resulting pvalue is recorded.
#'
#'
#' @name falsePositiveTest
#' @rdname falsePositiveTest
#' @aliases falsePositiveTest
#' @author Jason Serviss
#' @keywords falsePositiveTest
#' @param interval The interval of points amounts that should be tested. i.e. 1 tests points
#' amounts \code{1:100}, 10 tests \code{seq(1, 100, 10)}, etc.
#' @param iterations The number of iterations for each call to \code{\link[ClusterSignificance]{permute}}.
#' @param cores The number of cores the test should be parallalized with.
#' @param save Should the output be saved?
#' @param verbose Should the function be verbose?
#' @param outPath The location the file should be saved to if \code{save=TRUE}.
#' @examples
#' falsePositiveTest(interval=50, iterations=1, cores=1, save=FALSE)
NULL
#' @export


falsePositiveTest <- function(
    interval = 1,
    iterations = 10000,
    cores = 4,
    save=TRUE,
    verbose = TRUE,
    outPath = 'data'
) {
    output <- data.frame()
    tests <- seq(5, 100, 1)
    tests <- subset(tests, tests %% interval == 0)
    
    for(jj in 1:length(tests)) {
        
        if(verbose == TRUE) {print(paste("points: ", tests[jj], sep=""))}
        
        test <- paste(tests[jj], "points", sep=".")
        data(hundred)
        pointsArray <- hundred[[test]]
        
        pointsAmount <- gsub(
            "([0-9]*).*",
            "\\1",
            names(hundred[test])
        )
        
        repitition <- gsub(
            "([0-9]*).*",
            "\\1",
            dimnames(pointsArray)[[3]]
        )
        
        for(yy in 1:dim(pointsArray)[3]) {
            if(verbose == TRUE) {print(paste("reps: ", yy, sep=""))}
            
            mat <- list(pointsArray[,,yy])
            groups <- makeGroups(mat, paste("grp", 1:2, sep=""))
            
            if(verbose == TRUE) {print("Pcp")}
            
            permMat <- permMatrix(mat[[1]], groups, iterations)
            
            permutePcp <- parallelCS(
                mat = mat[[1]],
                groups = groups,
                iter = iterations,
                projmethod="pcp",
                cores = cores,
                user.permutations = permMat
            )
            
            if(verbose == TRUE) {print("Mlp")}

            permuteMlp <- parallelCS(
                mat = mat[[1]],
                groups = groups,
                iter = iterations,
                projmethod="mlp",
                cores = cores,
                user.permutations = permMat
            )
            
            pPcp <- pvalue(permutePcp)
            pMlp <- pvalue(permuteMlp)
            
            tmp <- data.frame(
                points = pointsAmount,
                rep = repitition[yy],
                pValuePcp = pPcp,
                pValueMlp = pMlp,
                Pcp.scores.real = getData(permutePcp, "scores.real")[[1]],
                Pcp.scores.vec = I(list(getData(permutePcp, "scores.vec")[[1]])),
                Mlp.scores.real = getData(permuteMlp, "scores.real")[[1]],
                Mlp.scores.vec = I(list(getData(permuteMlp, "scores.vec")[[1]]))
            )
            
            output <- rbind(output, tmp)
        }
    }
    rownames(output) <- 1:nrow(output)
    
    if(save == TRUE) {
        fpTestResults <- output
        save(
            fpTestResults,
            file=paste(outPath, "fpTestResults.rda", sep="/"),
            compress="bzip2"
        )
    }
    return(output)
}

