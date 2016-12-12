
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


#' Specificity test
#'
#' This function performs the specificity test. This test is designed to
#' determine the true negative rate of the Pcp method. To do this
#' we utilize the hundred dataset which was derived using the \code{\link{zeroOrHundred}}
#' function with argument \code{overlap}=100. The dataset includes matrices where
#' the group's range overlap is 100% and thus we anticipate that no significant
#' seperations should be found with Pcp. For each matrice, the
#' permute function is run, using the Pcp method afterwhich,
#' the resulting pvalue is recorded.
#'
#'
#' @name specificityTest
#' @rdname specificityTest
#' @aliases specificityTest
#' @author Jason Serviss
#' @keywords specificityTest
#' @param interval The interval of points amounts that should be tested. i.e. 1 tests points
#' amounts \code{1:100}, 10 tests \code{seq(1, 100, 10)}, etc.
#' @param iterations The number of iterations for each call to \code{\link[ClusterSignificance]{permute}}.
#' @param cores The number of cores the test should be parallalized with.
#' @param save Should the output be saved?
#' @param verbose Should the function be verbose?
#' @param outPath The location the file should be saved to if \code{save=TRUE}.
#' @examples
#' specificityTest(interval=50, iterations=1, cores=1, save=FALSE)
NULL
#' @export


specificityTest <- function(
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
        specificityTestResults <- output
        save(
        specificityTestResults,
        file=paste(outPath, "specificityTest.rda", sep="/"),
        compress="bzip2"
        )
    }
    return(output)
}

#' Sensitivity and Specificity Figure
#'
#' This function plots the specificity and sensivity testing results.
#'
#'
#' @name sensSpecPlot
#' @rdname sensSpecPlot
#' @aliases sensSpecPlot
#' @author Jason Serviss
#' @keywords sensSpecPlot
#' @param sens Sensivitity testing results data.frame.
#' @param spec Specificity testing results data.frame.
#' @examples
#'
#' sensSpecPlot()
#'
#'
NULL
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggthemes theme_few scale_colour_economist
#' @importFrom plyr ddply summarize

sensSpecPlot <- function(sens=NULL, spec=NULL) {
    
    if(is.null(spec) | is.null(sens)) {
        sensitivityTestResults <- sensitivityTestResults
        specificityTestResults <- specificityTestResults
    } else {
        sensitivityTestResults <- sens
        specificityTestResults <- spec
    }
    
    sensitivityTestResults$testType <- "Sensitivity test"
    specificityTestResults$testType <- "Specificity test"
    
    sensitivityTestResults$Pcp.BHcorrectedP <- ddply(
        sensitivityTestResults,
        "points",
        summarize,
        fdr=p.adjust(pValuePcp, method="fdr")
    )$fdr
    
    specificityTestResults$Pcp.BHcorrectedP <- ddply(
        specificityTestResults,
        "points",
        summarize,
        fdr=p.adjust(pValuePcp, method="fdr")
    )$fdr
    
    data <- rbind(
        sensitivityTestResults,
        specificityTestResults
    )[ ,c("points", "rep", "Pcp.BHcorrectedP", "testType")]
    
    m <- melt(data, id.vars=c("points", "rep", "testType"))
    m$points <- as.numeric(as.character(m$points))
    
    m <- subset(m, variable == "Pcp.BHcorrectedP")
    p <- ggplot(m, aes(x=as.factor(points), y=-log10(value), colour=testType))+
        geom_boxplot(outlier.size = 1)+
        geom_hline(yintercept=-log10(0.05), linetype=2, colour="darkgrey")+
        theme_few()+
        theme(
            legend.position="top",
            legend.title=element_blank(),
            legend.text=element_text(size=15),
            axis.title=element_text(size=17),
            axis.text=element_text(size=15),
            plot.title=element_text(
                hjust=0.5,
                family="ArialMT",
                face="plain",
                size=24
            ),
            strip.text.x = element_text(size = 17)
        )+
        labs(
            x="Points per group",
            y="-log10(adjusted p-Value)",
            title="Specificity/Sensitivity test results"
        )+
        scale_colour_economist()+
        scale_x_discrete(breaks = round(seq(min(m$points), max(m$points), by = 10),1))
    
    p
    return(p)
}
