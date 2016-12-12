
#' Sensitivity test
#'
#' This function performs the sensitivity test. This test is designed to
#' determine the true positive rate of the Pcp and Mlp methods. To do this
#' we utilize the zero dataset which was derived using the \code{\link{zeroOrHundred}}
#' function with argument \code{overlap}=0. The dataset includes matrices where
#' the groups overlap is 0% and thus we anticipate that significant
#' seperations should be found with Pcp or Mlp for all matrices. For each
#' matrice, the permute function is run, using both the Pcp and Mlp methods
#' afterwhich, the resulting pvalue is recorded.
#'
#'
#' @name sensitivityTest
#' @rdname sensitivityTest
#' @aliases sensitivityTest
#' @author Jason Serviss
#' @keywords sensitivityTest
#' @param interval The interval of points amounts that should be tested. i.e. 1 tests points
#' amounts \code{1:100}, 10 tests \code{seq(1, 100, 10)}, etc.
#' @param iterations The number of iterations for each call to \code{\link[ClusterSignificance]{permute}}.
#' @param cores The number of cores the test should be parallalized with.
#' @param save Should the output be saved?
#' @param verbose Should the function be verbose?
#' @param outPath The location the file should be saved to if \code{save=TRUE}.
#' @examples
#' sensitivityTest(interval=50, iterations=1, cores=1, save=FALSE)
NULL
#' @export

sensitivityTest <- function(
    interval = 1,
    iterations = 10000,
    cores = 4,
    save=TRUE,
    verbose=TRUE,
    outPath='data'
) {
    output <- data.frame()
    tests <- seq(5, 100, 1)
    tests <- subset(tests, tests %% interval == 0)
    
    for(jj in 1:length(tests)) {
        
        if(verbose == TRUE) {print(paste("points: ", tests[jj], sep=""))}

        test <- paste(tests[jj], "points", sep=".")
        data(zero)
        pointsArray <- zero[[test]]
        
        pointsAmount <- gsub(
            "([0-9]*).*",
            "\\1",
            names(zero[test])
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
            
            pPcp <- pvalue(permutePcp)
            
            tmp <- data.frame(
                points = pointsAmount,
                rep = repitition[yy],
                pValuePcp = pPcp
            )
            
            output <- rbind(output, tmp)
        }
    }
    rownames(output) <- 1:nrow(output)
    
    if(save == TRUE) {
        sensitivityTestResults <- output
        save(
            sensitivityTestResults,
            file=paste(outPath, "sensitivityTest.rda", sep="/"),
            compress="bzip2"
        )
    }
    return(output)
}
