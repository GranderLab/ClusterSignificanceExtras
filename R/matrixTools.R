
#' Make Matrix
#'
#' This function is specifically designed for the generation of 2 groups of points within a matrix.
#'
#' Makes a matrix containing data for two groups of points.
#' Specifying a specific number of points, dimensions, and range
#' of points for each group is suported. The resulting matrix will have group
#' 1 occupying the top half of the matrix and group 2 the bottom half.
#'
#'
#' @name makeMatrix
#' @rdname makeMatrix
#' @aliases makeMatrix
#' @param x A numerical vector of length 2 specifying the range of values possible for group 1.
#' @param y A numerical vector of length 2 specifying the range of values possible for group 2.
#' @param dim The number of dimensions (columns) to generate.
#' @param points The number of data points in each group.
#' @param by Indicates the frequency of ponts which are subsequently randomly picked according
#' to the \emph{points} argument.
#' @param seed Sets the seed to make results reproducible.
#' @param replace Passed to the \code{\link[base]{sample}} function.
#' @param ... Arguments to pass on.
#' @author Jason Serviss
#' @keywords makeMatrix
#' @examples
#' mat <- makeMatrix(x=c(0,1), y=c(0,1), dim=2, points=10)
NULL
#' @export

makeMatrix <- function(
    x,
    y,
    dim,
    points,
    by = 0.001,
    seed = 3,
    replace = FALSE,
    ...
) {
    set.seed(seed)
    
    #below, a call to sapply facilitates a iteration for each dim
    #for each iteration sample is called twice, one for x and one for y
    #results from the sample calls are concatenated
    #calls to sample randomly pick numbers between the values for x and y
    #the number of points picked by sample is governed by the points argument
    #finally, results from sample calls are formated into a matrix with the expected
    # ClusterSignificance structure.
    mat <- matrix(
        c(
            sapply(1:dim, function(xi)
                c(
                    sample(
                        seq(
                            x[1],
                            x[2],
                            by
                        ),
                        replace=replace,
                        size = points
                    ),
                    sample(
                        seq(
                            y[1],
                            y[2],
                            by
                        ),
                        replace=replace,
                        size = points
                    )
                )
            )
        ), ncol = dim
    )
    output <- list(mat, x, y) #matrix and arguments output
    class(output) <- "makeMatrix" #class assignment to output
    return(output)
}

#' Make Groups
#'
#' This function makes the \code{groups} input variable to ClusterSignificance from a
#' \code{\link{makeMatrix}} generated matrix.
#'
#' This function is specifically designed for the generation of the groups argument corresponding
#' to the output from the \code{\link{makeMatrix}} function. It always results in 2 groups with
#' each group being assigned an equal amount of points in the matrix.
#'
#' @name makeGroups
#' @rdname makeGroups
#' @aliases makeGroups
#' @param mat The matrix for which the groups should be generated.
#' @param names A character vector indicating the names of the groups.
#' @param ... Arguments to pass on.
#' @author Jason Serviss
#' @keywords makeGroups
#' @examples
#' mat <- makeMatrix(x=c(0,1), y=c(0,1), dim=2, points=10)
#' groups <- makeGroups(mat, names=c("grp1", "grp2"))
NULL
#' @export

makeGroups <- function(
    mat,
    names,
    ...
) {
    #input checks
    if(length(mat[[1]]) == 1) {
        stop("the matrix you submitted is not properly formated")
    }
    
    if(length(names) < 2) {
        stop("the names you submitted are < length 2")
    }
    
    mat <- mat[[1]]
    
    #below a call to sapply provides iteration over the names argument
    #the names are evenly divided over the rows of the matrix
    #this is done so that the first name gets the first division of rows
    #the second name gets the second division of rows, et.
    groups <- c(
        sapply(
            1:length(names),
            function(ix, iy)
                rep(
                    iy[ix],
                    nrow(mat)/length(names)
                ),
                iy = names
        )
    )
    return(groups)
}

#' Overlap Specific Matrices
#'
#' Generates matrices where the groups have a specific percentage of overlap with each other.
#'
#' The function generates overlap specific matrices. It can do this either for a specific
#' overlap or many at once. It is also capable of generating replicates for a single overlap value.
#' All matrices must pass multiple tests before they are accepted with which tests they must pass being
#' dependant on the desired overlap. This function has only been tested with the \code{dim}, \code{points},
#' and \code{by} values set to their defaults as it utilizes the \code{\link{overlapProbability}} output and
#' the \code{\link{dynamicXY}} function.
#'
#'
#' @name overlapMatrices
#' @rdname overlapMatrices
#' @aliases overlapMatrices
#' @param reps The desired number of overlap specific matrices to generate for each overlap.
#' @param overlaps A numerical vector containing all of the overlaps to be generated.
#' @param dim The number of dimensions (columns) to generate. Passed to the makeMatrix function.
#' @param points The number of data points in each group. Passed to the makeMatrix function.
#' @param by Passed to the \code{\link{makeMatrix}} function.
#' @param verbose Should the function be verbose?
#' @param save Should output be saved?
#' @author Jason Serviss
#' @keywords overlapMatrices
#' @examples
#' matrices <- overlapMatrices(reps=2, overlaps=0, dim=2, points=10, save=FALSE, verbose=FALSE)
NULL
#' @export

overlapMatrices <- function(
    reps = 10,
    overlaps = seq(0,100, by=10),
    dim = 2,
    points = 102,
    by = 0.9,
    verbose = TRUE,
    save = TRUE,
    ...
) {
    
    #setup a empty array to hold the data
    multiplyMatrix <- array(
        NA,
        dim=c(
            points*2,
            dim,
            reps * length(overlaps)
        )
    )
    
    dimnames(multiplyMatrix)[[3]] <- paste(
        sort(
            rep(
                overlaps,
                reps
            )
        ),
        1:reps,
        sep="."
    )
    
    #run a iteration for each of the desired overlaps
    for(y in overlaps) {
        
        if(verbose == TRUE) {message(y)}
        
        ##make otherMatrices variable for this y
        ##the otherMatrices variable keeps track of previously generated matrices
        ##it is fed into the matrix tests which checks to make sure
        ##that no matrices in different reprtitions are identical.
        otherMatrices <- array(
            0,
            dim=c(
                points*2,
                dim,
                reps
            )
        )
        
        #run one iteration for each repetition
        for(i in 1:reps) {
            
            #make the mat and group input
            mat <- makeMatrix(
                x=c(0,0),
                y=c(0,0),
                dim=dim,
                points=points,
                replace=TRUE
            )
            
            groups <- makeGroups(
                mat,
                paste(
                    "grp",
                    1:dim,
                    sep=""
                )
            )
            
            #specify which tests should be performed on the matrix and perform them
            #if the desired overlap is 100, additional tests are added
            #otherwise the "standard" tests are used
            #note that if the matrix does NOT pass one or more of the tests a new
            #matrix is generated within the testing function and re-tested and that
            #this is repeated until the matrix passes all tests, then it is returned here
            if(y==100){
                tests <- "all"
                interspersionThreshold <- 50
            } else {
                tests <- c(
                    "checkIdenticalDims",
                    "checkIdenticalMatrix",
                    "overlapSpecificMatrix",
                    "checkIdenticalPoints"
                )
                interspersionThreshold <- ""
            }
            
            mat2 <- matrixTests(
                tests = tests,
                mat = mat,
                groups = groups,
                dim = dim,
                points = points,
                by = by,
                otherMatrices = otherMatrices,
                desiredOverlap = y,
                interspersionThreshold = interspersionThreshold,
                idPointsThreshold = 3,
                verbose=FALSE
            )
            
            
            ##collect matrixes for this overlap in an array
            #as mentioned above, this tracking is used to make sure that
            #matrices in different repetitions are not identical
            otherMatrices[,,i] <- mat2[[1]][[1]]
            
            ##cache the generated matrix for subsequent return
            multiplyMatrix[,,paste(y,i,sep=".")] <- mat2[[1]][[1]]
        }
    }
    
    #save the results if so desired
    if(save == TRUE) {
        overlapMatrixData <- multiplyMatrix
        save(
            overlapMatrixData,
            file="data/overlapMatrixData.rda",
            compress="bzip2"
        )
    }
    
    return(list(multiplyMatrix))
}


#' Dynamic X and Y setter
#'
#' This function allows the dynamic adjustment of X and Y values in the \code{\link{makeMatrix}}
#' function dependant on the desired group overlap. It uses a previously generated table from the
#' \code{\link{overlapProbability}} function. The use of this function allows much faster generation of
#' matrices with a specific group overlap than if the computer were to need to randomly try
#' X and Y values until the overlap was generated. Y values are currently always (0,101) and
#' X2 values are always x1 + 100. Is appropriate only to use when the points argument to
#' \code{\link{makeMatrix}} is 102 and the by argument in 0.9. The returned X and Y values have
#' been optimized to have the highest probability to generate the same overlap in both dimensions
#' of the matrix i.e. when the \code{\link{makeMatrix}} argument \code{dim} is set to 2.
#'
#'
#' @name dynamicXY
#' @rdname dynamicXY
#' @aliases dynamicXY
#' @param desiredOverlap Overlap for which X and Y values should be generated.
#' @author Jason Serviss
#' @keywords dynamicXY
#' @examples
#' XY <- dynamicXY(overlapTable, desiredOverlap = 0)
#' mat <- makeMatrix(XY[[1]], XY[[2]], dim=2, points=102)[[1]]
#' groups <- makeGroups(mat, names=c("grp1", "grp2"))
#' calculateOverlap(mat[ ,1], groups) #overlap dimension 1
#' calculateOverlap(mat[ ,2], groups) #overlap dimension 2
NULL
#' @export

dynamicXY <- function(
    desiredOverlap
){
    
    #below, y is set to the values used for the calculation of
    #the overlapTable
    yOut <- c(0,101)
    
    #below, we search the overlapTable for the corresponding
    #x1 value according to the desired overlap sepcified in the
    #desiredOverlap argument
    #x2 always equals x1+100
    
    ##if handles the situation where a specific x is available
    ##in the table for the desired overlap and that overlap is
    ##not 100
    
    ##the first else if handles a 100 percent overlap
    ##the 2nd else if handles a 0 percent overlap

    ##else handles the case where the overlap is not found
    if(
        desiredOverlap != 100 &
        desiredOverlap != 0
    ) {
        if(
            nrow(subset(
                overlapTable,
                overlap == desiredOverlap
            )) == 1
        ){
            x1 <- as.numeric(
                subset(
                    overlapTable,
                    overlap == desiredOverlap,
                    select="x1"
                )[1:1,]
            )
            
            x2 <- x1 + 100
            xOut <- c(x1, x2)
        }
        
    } else if(desiredOverlap == 100) {
        
        xOut <- yOut
        
    } else if(desiredOverlap == 0) {
        
        xOut <- c(101, 101+100)
        
    } else {
        
        stop(
            paste(
                "Could not find entry for ",
                desiredOverlap,
                " in the look-up table",
                sep=""
            )
        )
        
    }
    
    return(list(xOut, yOut))
}

#' Add data point
#'
#' Adds one point to each group in the matrix.
#'
#' This function is designed to add one point in each dimension per group to the current matrix.
#' The points that are added are done so with a call to \code{\link{makeMatrix}} where the \code{x}
#' and \code{y} arguments are identical to the origional matrix. The function is primairly designed
#' to work with the \code{\link{zeroOrHundred}} function.
#'
#' @name addOnePoint
#' @rdname addOnePoint
#' @aliases addOnePoint
#' @param mat from makeMatrix
#' @param groups from makeGroups
#' @param seed passed to set.seed
#' @author Jason Serviss
#' @keywords addOnePoint
#' @examples
#' mat <- makeMatrix(x=c(0,100), y=c(0, 100), dim=2, points=10)
#' groups <- makeGroups(mat, names=c("grp1", "gpr2"))
#' matPlus <- addOnePoint(mat, groups)
NULL
#' @export

addOnePoint <- function(
    mat,
    groups,
    seed = 3
){
    
    #input check
    if(length(mat[[1]]) == 1){
        stop("the matrix you submitted is not properly formated")
    }
    
    #extract the matris and the x and y values used to generate it
    m <- mat[[1]]
    x <- mat[[2]]
    y <- mat[[3]]
    
    #etract the values corresponding to both groups from the matrix
    t <- ifelse(
        groups == unique(groups)[1],
        TRUE,
        FALSE
    )
    
    m1 <- m[t, ]
    m2 <- m[!t, ]
    
    #make the points to add to the input matrix
    #utilize the same x and y values as the input matrix
    newPoints <- makeMatrix(
        x=x,
        y=y,
        dim=dim(m)[2],
        points=1,
        seed=seed
    )[[1]]
    
    #add the new points to the input matrix
    n1 <- rbind(m1, newPoints[1,])
    n2 <- rbind(m2, newPoints[2,])
    mat <- rbind(n1, n2)
    
    #retrun the results
    out <- list(mat, x, y)
    class(out) <- "makeMatrix"
    return(out)
}

#' Zero or Hundred Matrix
#'
#' Creates the \emph{zero} or \emph{hundred} dataset.
#'
#' Creates one of two datasets; one where the 2 groups have 0 percent overlap with each other or
#' one where the two groups have 100 percent overlap with each other. 0 percent, in this context,
#' is defined as no datapoints from one group occupy the same space as the points of the other group.
#' 100 percent overlap is defined as one group completley containing the other group within the space
#' that it occupies. For ClusterSignificance testing this function was run with the default paramaters
#' for overlaps 0 and 100.
#'
#' @name zeroOrHundred
#' @rdname zeroOrHundred
#' @aliases zeroOrHundred
#' @param overlap The desired numerical overlap to be generated. Can be either 0 or 100.
#' @param dim The number of desired spatial dimensions for each groups points to occupy.
#' @param reps The number of times a dataset for a specific number of points should be generated.
#' @param p The interval of points per group for which datasets should be generated.
#' @param verbose Logical indicating if the program should report status updates.
#' @param save Logical indicating if the output should be saved to the data folder of the current directory.
#'
#' @author Jason Serviss
#' @keywords zeroOrHundred
#' @examples
#' zeroOrHundred(overlap=0, reps=2, p=5:6, save=FALSE)
NULL
#' @export

zeroOrHundred <- function(
    overlap,
    dim = 2,
    reps = 10,
    p = 5:100,
    verbose = TRUE,
    checksVerbose = FALSE,
    save=TRUE
) {
    
    #get the x and y arguments for makeMatrix depending on the overlap argument
    start <- dynamicXY(
        desiredOverlap = overlap
    )
    x <- start[[1]]
    y <- start[[2]]
    
    output <- list()
    previousPoints <- array(NA, dim=c((min(p)*2), dim, reps))
    seed <- 3
    
    #one iteration for each points amount
    for(i in p) {
        if(verbose == TRUE){print(paste("points: ", i, sep=""))}
        
        #for the first points amount (i=5) "seed matrices"
        #are made using the overlapMatrices function.
        #For each additional, points amount, instead of
        #creating a totally new matrix, values are added to
        #the seed matrices.
        if(i == min(p)){
            
            otherMatrices <- overlapMatrices(
                reps=reps,
                overlaps=overlap,
                dim=2,
                points=i,
                save=FALSE,
                verbose=FALSE
            )[[1]]
            
            dimnames(otherMatrices) <- NULL
            
        } else {
            
            #get the matrices corresponding to the previous points amount value so they
            #can be added to.
            previousPoints <- output[[i-5]]
            otherMatrices <- array(0, dim=c(i*2, dim, reps))
            
            #run one iteration per points amount replicate
            for(j in 1:reps){
                if(verbose == TRUE) {print(paste("rep: ", j, sep=""))}
                
                #extract the matrix for the corresponding replicate with the previous points amount
                #generate corresponding groups variable
                mat <- previousPoints[,,j]
                groups1 <- makeGroups(list(mat), names=c("grp1", "grp2"))
                mat <- list(mat, x, y)
                
                ##add points to the matrix corresponding to the
                #previous points amount so that it contains
                #the number of points corresponding to the
                #current points amount.
                #run checks, if the checks are ok, continue, otherwise,
                #generate different points and re-check
                while(TRUE) {
                    
                    mat1 <- addOnePoint(mat, groups1, seed)
                    groups2 <- makeGroups(mat1, names=c("grp1", "grp2"))
                    
                    if(overlap == 100){
                        tests <- "all"
                        interspersionThreshold <- 50
                    } else {
                        tests <- c(
                            "checkIdenticalDims",
                            "checkIdenticalMatrix",
                            "overlapSpecificMatrix",
                            "checkIdenticalPoints"
                        )
                    }
                    
                    checks <- matrixTests(
                        tests = tests,
                        mat = mat1,
                        groups = groups2,
                        dim = dim,
                        points = i,
                        by = 0.9,
                        seed = seed,
                        otherMatrices = otherMatrices,
                        desiredOverlap = overlap,
                        interspersionThreshold = interspersionThreshold,
                        idPointsThreshold = 1,
                        verbose = ifelse(checksVerbose == TRUE, TRUE, FALSE)
                    )[[9]]
                    
                    if(checks == 0) {break}
                    seed <- seed + 1
                }
                #add this repitition to the otherMatrices variable
                otherMatrices[,,j] <- mat1[[1]]
            }
        }
        #name and add all the matrices corresponding to the
        #current points amount to the output variable
        dimnames(otherMatrices)[[3]] <- paste(1:reps, "repitition", sep=".")
        output[[i-4]] <- otherMatrices
    }
    
    #give the output variable meaningful names and save
    names(output) <- paste(p, "points", sep=".")
    
    if(save == TRUE & overlap == 0) {
        zero <- output
        save(
            zero,
            file="data/zero.rda",
            compress="bzip2"
        )
    }
    
    if(save == TRUE & overlap == 100) {
        hundred <- output
        save(
            hundred,
            file="data/hundred.rda",
            compress="bzip2"
        )
    }
    return(output)
}

#' Overlap Probability Function
#'
#' This function generates a table that contains the probability of generating a specific overlap between groups
#' when starting with a specific x value.
#'
#' The function works by examining probability of every x value from 0-100 by 0.9 to generate a matrix with a specific group
#' overlap. x, in this case, is x[1] supplied to \code{\link{makeMatrix}} and x[2] = x[1] + 100. A matrix is generated \code{reps} times for each x
#' using the \code{\link{makeMatrix}} function. The overlap between groups is then calculated. When this process completes, a lookup table is
#' generated. Each x is listed in the lookup table together with the overlap generated in dimension 1 and 2. It was found that
#' the same overlap is observed in each dimension and that the probability for each overlap is 100 percent (i.e. same overlap observed
#' in all trials). The lookup table is output as a data frame, primiarly, for use by the \code{\link{overlapMatrices}} function.
#' Note that the 100 percent overlap for both dims is missing from the table and is currently handled by the \code{\link{dymanicXY}} function.
#' As well, the lookup table does include the 0 percent overlap but this is overridden by the \code{\link{dynamicXY}} function.
#'
#' @name overlapProbability
#' @rdname overlapProbability
#' @aliases overlapProbability
#' @author Jason Serviss
#' @keywords overlapProbability
#' @examples
#' overlapProbability(save=FALSE, reps=2)
NULL
#' @export

overlapProbability <- function(
    points = 102,
    reps = 30,
    x1 = seq(0, 101, by=1),
    by = 0.9,
    save = TRUE,
    verbose = TRUE
) {
    #below is an example of "left-side variable assignment"
    #i.e. pointsOverlaps is assigned data.frame() and count is assigned 1.
    #syntactic sugar
    g(pointsOverlaps, count) %=% list(data.frame(), 1)
    
    #the first loop was included with the original thought of having a
    #overlap table with many different points amounts, although, after
    #the realization that not all overlaps can be generated with any points
    #amount, the function was always run with points = 102, due to the fact
    #that all overlaps can be generated with this points amount.
    for(c in points) {
        #set y values and x2 values to static values, only x1 varies
        g(x2, y) %=% list(x1 + 101, c(0, 101))
        
        if(verbose == TRUE) {message(c)}
        
        #iterate over the values of x1
        #these values represent all values from y1 to y2.
        #thus x acts like a sliding window over y from, in theory,
        #100 percent overlap of the y range to 0 percent overlap of
        #the y range
        for(i in 1:length(x1)) {
            
            #make an array to keep track of matrices which have already
            #been tested so they are not tested again
            otherMatrices <- array(
                rep(0, c*2*2*reps),
                dim=c(c*2, 2, reps)
            )
            
            #iterate over each desired replicate of this x1 value
            for(z in 1:reps) {
                
                #create the mat and groups variable
                cmag <- .createMatAndGroups(x1, x2, i, y, c, by)
                mat <- cmag[[1]]
                groups <- cmag[[2]]
                
                #run the matrix tests on the matrix
                #note that if the matrix does not pass one or more of the
                #tests, it is automatically re-generated by the matrixTests
                #function until it passes all tests, at which time, it is
                #returned here.
                mat <- matrixTests(
                    tests=c(
                        "checkIdenticalDims",
                        "checkIdenticalMatrix",
                        "checkIdenticalPoints"
                    ),
                    mat,
                    groups,
                    dim = 2,
                    points = c,
                    by = by,
                    otherMatrices = otherMatrices,
                    idPointsThreshold = 3,
                    verbose = ifelse(verbose == TRUE, TRUE, FALSE)
                )[[1]]
                
                #add the matrix to the tracking array, explained above
                otherMatrices[,,z] <- mat[[1]]
                
                #calculate the group overlap percentage for the current matrix
                overlap <- calculateOverlap(
                    mat,
                    groups
                )
                
                #the .savePointsOverlaps function reformats the data associated
                #with the current matrix for subsequent adding to the output variable
                add <- .savePointsOverlaps(
                    i,
                    c,
                    z,
                    x1,
                    overlap,
                    mat
                )
                
                #bind the current matrix and its associated data to the output variable
                pointsOverlaps <- rbind(
                    pointsOverlaps,
                    add
                )
                
                count <- count + 1
            }
        }
    }
    
    #below we examine, for each x1, which overlap is most commonly generated
    pointsOverlaps$overlap <- round(pointsOverlaps$overlap)
    
    freq <- plyr::ddply(
       pointsOverlaps,
       c("x1", "overlap"),
       plyr::summarize,
       freq=length(x1)
    )
    
    #below we extract the x1 with the highest probability for giving a specific overlap
    overlapTable <- plyr::ddply(
       freq,
       "overlap",
       .fun = function(x)
           plyr::summarise(
               x,
               x[which(x$freq == max(x$freq)), 'x1']
           )[1:1,]
    )
    
    #save the output os so desired
    colnames(overlapTable) <- c("overlap", "x1")
    
    if(save == TRUE) {
        save(
            overlapTable,
            file = "data/overlapTable.rda",
            compress="bzip2"
        )
    }
    
    return(overlapTable)
}

#creates the matrix and group variables
.createMatAndGroups <- function(x1, x2, i, y, c, by) {
    mat <- makeMatrix(
        x=c(x1[i], x2[i]),
        y=y,
        dim=2,
        points = c,
        by = by
    )
    
    groups <- makeGroups(
        mat,
        names=c("grp1", "grp2")
    )
    return(list(mat, groups))
}

#reformats the data associated with the current matrix
.savePointsOverlaps <- function(
    i,
    c,
    z,
    x1,
    overlap,
    mat
){
    add <- data.frame(
        points = c,
        reps = z,
        x1 = x1[i],
        overlap = overlap,
        mat = I(list(mat)),
        stringsAsFactors = FALSE
    )
    return(add)
}

#' Make permutation matrix
#'
#' The .permMatrix function from ClusterSignificance.
#'
#'
#' @name permMatrix
#' @rdname permMatrix
#' @aliases permMatrix
#' @param mat from makeMatrix
#' @param groups from makeGroups
#' @param iterations number of iterations.
#' @author Jason Serviss
#' @keywords permMatrix
#' @examples
#'
#' mat <- makeMatrix(x=c(0,100), y=c(0, 100), dim=2, points=10)
#' groups <- makeGroups(mat, names=c("grp1", "gpr2"))
#' matPlus <- addOnePoint(mat, groups)
#'
#'
NULL
#' @export

## the permutation matrix function from ClusterSignificance
## takes a matrix as input, permutes the data x times, and outputs the result.

permMatrix <- function(
    mat,
    groups,
    iterations
){
    
    uniq.groups <- combn(
        unique(groups),
        2
    )
    
    permats <- lapply(
        1:ncol(uniq.groups),
        function(y)
            lapply(
                1:iterations,
                function(x)
                    sapply(
                        1:ncol(mat[groups %in% uniq.groups[, y], ]),
                        function(i)
                            mat[groups %in% uniq.groups[, y], i] <-
                                sample(
                                    mat[groups %in% uniq.groups[, y], i],
                                    replace=FALSE
                                )
                    )
            )
    )
    
    return(permats)
}

############################################################################################################
#                                                                                                          #
#                                               Checks                                                     #
#                                                                                                          #
############################################################################################################


#' Matrix Tests
#'
#' A wrapper for all matrix checks. The function also handles "regeneration", i.e.
#' if a matrix does \strong{not} pass the checks that it should, it is regenerated
#' until it does. Which tests should be run is controled by the \code{tests} arguement.
#'
#' @name matrixTests
#' @rdname matrixTests
#' @aliases matrixTests
#' @author Jason Serviss
#' @keywords matrixTests
#' @param tests Specifies which tests should be run.
#' @param mat The input matrix.
#' @param groups The groups corresponding to the input matrix.
#' @param dim The dimensions (colums) of the input matrix.
#' @param points The number of points per group in the input matrix.
#' @param by The \code{by} variable to the \code{\link{makeMatrix}} function.
#' @param seed The \code{seed} variable to the \code{\link{makeMatrix}} function.
#' @param otherMatrices Input to the \code{\link{checkIdenticalMatrix}} function.
#' @param desiredOverlap The desired group overlap.
#' @param interspersionThreshold Input to the \code{\link{interspersion}} function.
#' @param idPointsThreshold Input to the \code{\link{checkIdenticalDims}} function.
#' @param verbose Logical if the function should be verbose. Mostly for testing.
NULL
#' @export

matrixTests <- function(
    tests = "all",
    mat,
    groups,
    dim = ncol(mat[[1]]),
    points = nrow(mat)/2,
    by = 0.9,
    seed = 3,
    otherMatrices = "",
    desiredOverlap = "",
    interspersionThreshold = "",
    idPointsThreshold = "",
    verbose = FALSE
) {
    
    if(is.numeric(desiredOverlap)) {
        XY <- dynamicXY(desiredOverlap)
        x <- XY[[1]]
        y <- XY[[2]]
    } else {
        x <- mat[[2]]
        y <- mat[[3]]
    }
    
    iterations = 0
    
    #the while loop below runs the desired tests and checks that the input matrix passes all tests
    #at 2 points, the input matrix is checked to see if it has passed the tests thus far, if not
    #a new matrix is generated (with the .regenerate funciton) and the tests are re-run. Once the
    #matrix passes all tests it is returned.
    while(TRUE) {
        testOutput <- list()
        
        if(verbose == TRUE) {
            print(paste("XY:", x, y, sep=" "))
        }
        
        #run identical dimensions test
        if("checkIdenticalDims" %in% tests | "all" %in% tests) {
            testOutput <- c(testOutput, checkIdenticalDims(mat, groups))
            if(verbose == TRUE) {
                print(paste("checkIdenticalDims: ", checkIdenticalDims(mat, groups), sep=""))
            }
        }
        
        #run identical matrix test
        if("checkIdenticalMatrix" %in% tests | "all" %in% tests) {
            testOutput <- c(testOutput, checkIdenticalMatrix(mat, otherMatrices))
            if(verbose == TRUE) {
                print(paste("checkIdenticalMatrix: ", checkIdenticalMatrix(mat, otherMatrices), sep=""))
            }
        }
        
        #checks if the matrix has passed all tests thus far and, if not, regenerates the matrix
        if(any(testOutput == FALSE)) {
            new <- .regenerate(x, y, dim, points, by, seed, verbose)
            mat <- new[[1]]
            groups <- new[[2]]
            seed <- new[[3]]
            iterations <- iterations + 1
            next
        }
        
        #check the group overlap in the current matrix
        if("overlapSpecificMatrix" %in% tests | "all" %in% tests) {
            overlap1 <- calculateOverlap(mat[[1]][,1], groups)
            overlap2 <- calculateOverlap(mat[[1]][,2], groups)

            testOutput <- c(testOutput, overlap1 == desiredOverlap & overlap2 == desiredOverlap)
            
            if(verbose == TRUE) {
                print(paste("overlaps:", overlap1, overlap2, sep=" "))
            }
        }
        
        #run the interspersion test
        if("interspersion" %in% tests | "all" %in% tests) {
            intValues <- numeric()
            for(i in 1:dim) {
                int <- interspersion(mat, groups, i)
                intValues <- c(intValues, int)
            }
            testOutput <- c(testOutput, all(intValues > interspersionThreshold))
            if(verbose == TRUE) {
                print(paste("interspersion: ", intValues, sep=""))
            }
        }
        
        #run the check identical points test
        if("checkIdenticalPoints" %in% tests | "all" %in% tests) {
            idPoints <- logical()
            for(q in 1:dim){
                ident <- checkIdenticalPoints(mat, groups, q, idPointsThreshold)
                idPoints <- c(idPoints, ident)
            }
            testOutput <- c(testOutput, all(idPoints == TRUE))
            if(verbose == TRUE) {
                print(paste("checkIdenticalPoints: ", ident, sep=""))
            }
        }
        
        #if the matrix passes all tests break the while loop
        if(all(testOutput == TRUE)) {
            if(verbose == TRUE){print("complete")}
            break
        }
        
        if(verbose == TRUE) {
            print(paste("testOutput", testOutput, sep=": "))
        }
        
        ##if the matrix didn't pass all tests, generate a new matrix
        new <- .regenerate(x, y, dim, points, by, seed, verbose)
        mat <- new[[1]]
        groups <- new[[2]]
        seed <- new[[3]]
        iterations <- iterations + 1
    }
    
    return(list(mat, groups, otherMatrices, desiredOverlap, dim, points, by, seed, iterations))
}

#handles matrix regeneration
.regenerate <- function(x, y, dim, points, by, seed, verbose) {
    if(verbose == TRUE){print("re-generating")}
    seed <- seed + 1
    
    mat <- makeMatrix(
        x = x,
        y = y,
        dim = dim,
        points = points,
        by = by,
        seed = seed
    )
    
    groups <- makeGroups(
        mat,
        names = paste("grp", 1:dim, sep="")
    )
    
    return(list(mat, groups, seed))
}

#' Calculate Overlap
#'
#' Calculate the overlap between two groups.
#'
#' This function accepts a matrix containing data points for 2 groups where group 1
#' occupies the top half of the matrix and group 2, the bottom half. It calculates and
#' returns the overlap between the 2 groups. \emph{Overlap} in this case, is calculated dependant
#' on the class of the mat argument. If the class of the mat argument is numeric and is one dimension
#' form a \code{\link{makeMatrix}} matrix, then overlap refers to
#' \eqn{100 * ((grp 1 points in grp2 range + grp2 points in grp1 range) / total number of points)}.
#' If, on the other hand, the class of mat is \emph{makeMatrix} then the percentage of intersecting area
#' between the two groups is calculated.
#'
#' @name calculateOverlap
#' @rdname calculateOverlap
#' @aliases calculateOverlap
#' @param mat The numerical vector.
#' @param groups Defines the location of the groups within the vector.
#' @param ... Arguments to pass on.
#' @author Jason Serviss
#' @keywords calculateOverlap
#' @references \url{https://en.wikipedia.org/wiki/Convex_hull}
#' @examples
#' mat <- makeMatrix(x=c(0,1), y=c(1.1,2.1), dim=2, points=20)
#' groups <- makeGroups(mat, c("g1", "g2"))
#' calculateOverlap(mat, groups)
NULL
#' @export
#' @importFrom gpclib area.poly intersect
#' @importClassesFrom gpclib gpc.poly

calculateOverlap <- function(
    mat,
    groups,
    ...
) {
    #input checks
    if(class(mat) != "makeMatrix" & class(mat) != "numeric"){
        stop("the matrix you submitted is not properly formated")
    }
    
    if(class(mat) == "numeric") {
        spl <- split(
            mat,
            groups
        )
        
        #the if statment handles 100 percent overlap, i.e. one group's
        #range is totally contained with the other groups range.
        
        #else if handles 0 percent overlap, i.e. the range on both
        #groups is not overlapping
        
        #finally, else handles all other scenerios
        
        if(
            all(spl[[2]] > min(spl[[1]]) & spl[[2]] < max(spl[[1]])) |
            all(spl[[1]] > min(spl[[2]]) & spl[[1]] < max(spl[[2]]))
        ){
            percent_overlap <- 100
            return(percent_overlap)
            
        } else if(
            all(spl[[2]] > max(spl[[1]])) |
            all(spl[[1]] > max(spl[[2]]))
        ) {
            percent_overlap <- 0
            return(percent_overlap)
            
        } else {
            #calculate the number of points from each group within the other
            #groups range as a logical and then as a number
            grp1 <- spl[[1]] > min(spl[[2]]) & spl[[1]] < max(spl[[2]])
            grp2 <- spl[[2]] > min(spl[[1]]) & spl[[2]] < max(spl[[1]])
            total <- c(grp1, grp2)
            count <- length(total[total == TRUE])
            
            #calculate the percentage of overlap
            percent_overlap <- round(
                (count / length(mat)) * 100
            )
            
            return(percent_overlap)
        }

    }
    
    #intersecting are calculation
    #this first uses "convex hulls" to determine a shape for each group
    #after the shape is determined, the intersecting area of the 2 shapes is calculated
    if(class(mat) == "makeMatrix") {
        #make a variable holding the matrix values for each group
        mat <- mat[[1]]
        spl <- split(mat, groups)
        g1 <- matrix(spl[[1]], ncol=2)
        g2 <- matrix(spl[[2]], ncol=2)
        
        #find the convex hull for each group
        hull1 <- chull(g1)
        hull2 <- chull(g2)
        
        #use the convex hull to make a polygon representing each group
        p1 <- as(g1[hull1, ], "gpc.poly")
        p2 <- as(g2[hull2, ], "gpc.poly")
        
        #calculate the area of each polygon
        area1 <- area.poly(p1)
        area2 <- area.poly(p2)
        
        #find the intersecting area of each polygon and calculate the percentage of overlap
        areaIntersection <- area.poly(intersect(p1,p2))
        percentOverlap <- round(((areaIntersection * 2) / (area1 + area2)) * 100)
    }
    
    return(percentOverlap)
}

#' Calculate interspersion
#'
#' Calculates the interspersion between 2 groups. For example, if
#' the the groups are alternating when the points are ordered, the
#' interspersion will be 100%. Alternativley, if one groups points
#' lay within the other groups minimum point and its minimum - 1
#' point, the result will be 0%. See the examples.
#'
#'
#' @name interspersion
#' @rdname interspersion
#' @aliases interspersion
#' @param mat The mattrix to be checked. Usually a \code{\link{makeMatrix}} matrix.
#' @param groups A character vector indicating the groups location within the matrix.
#' Can be generated with the makeGroups function.
#' Can only accept 2 groups.
#' @param dim Which dimention of the matrix to check.
#' @author Jason Serviss
#' @keywords interspersion
#' @examples
#' mat <- makeMatrix(x=c(0,1), y=c(0,1), dim=2, points=10, seed=15)
#' groups <- makeGroups(mat, names=c("g1", "g2"))
#' interspersion(mat, groups, 1)
#'
#' mat <- makeMatrix(x=c(0,1), y=c(0.9,0.95), dim=2, points=10, seed=15)
#' groups <- makeGroups(mat, names=c("g1", "g2"))
#' interspersion(mat, groups, 1)
NULL
#' @export

interspersion <- function(
    mat,
    groups,
    dim
) {
    if(length(mat[[1]]) == 1){
        stop("the matrix you submitted is not properly formated")
    }
    
    mat <- mat[[1]]
    vector <- mat[, dim]
    
    l <- length(groups) / 2
    
    ##max score; this represents 100% interspersion
    ##indexes of grp1 and grp2 alternate every element
    names <- rep(
        c("1", "2"),
        l
    )
    
    max <- sum(
        which(names == "1")
    ) * sum(
        which(names == "2")
    )
    
    ##min score; this represents 0% interspersion
    ##indexes of group2 are within one group1 interval
    names <- c(
        "1",
        rep("2", l),
        rep("1", l - 1)
    )
    
    min <- sum(
        which(names == "1")
    ) * sum(
        which(names == "2")
    )
    
    ##scale; here i divide range between the min and max possible interspersion into
    ##100 elements. By performing the same calculation with the real data and then seeing
    ##where it falls on the scale, I can determine the "percentage of interspersion".
    s <- seq(
        min,
        max,
        by=(max-min)/100
    )
    
    ##real score; calculates the real score
    names(vector) <- groups
    g <- unique(groups)
    vector <- sort(vector)
    
    real <- sum(
        which(names(vector) == g[1])
    ) * sum(
        which(names(vector) == g[2])
    )
    
    ##score; finds where the real score falls on the scale of the theoretical scale
    score <- which.min(
        abs(s - real)
    ) - 1
    
    return(score)
}

#' Check Identical Dimensions
#'
#' This function takes an input matrix and checks if any
#' single dimension for a single group is equal to any other
#' single dimension for a single group including itself.
#' For example, with 2 groups with 2 dimensions each, it checks
#' if x dimension for group 1 is exactly equal to dimension y for
#' group 1 or dimension x or y for group 2.
#' \strong{Returns TRUE if none are not identical and FALSE if they are}.
#'
#'
#' @name checkIdenticalDims
#' @rdname checkIdenticalDims
#' @aliases checkIdenticalDims
#' @param mat \code{\link{makeMatrix}} function output.
#' @param groups \code{\link{makeGroups}} function output.
#' @param ident A boolean value. Do not alter.
#' @author Jason Serviss
#' @keywords checkIdenticalDims
#' @examples
#' mat <- makeMatrix(x=c(0,1), y=c(0,1), dim=2, points=10)
#' groups <- makeGroups(mat, names=c("grp1", "grp2"))
#' checkIdenticalDims(mat, groups)
#'
#' mat <- list(matrix(rep(1,40), ncol=2))
#' groups <- makeGroups(mat, c("grp1", "grp2"))
#' checkIdenticalDims(mat, groups)
NULL
#' @export

checkIdenticalDims <- function(
    mat,
    groups,
    ident = FALSE
) {
    
    if(length(mat[[1]]) == 1){
        warning("the matrix you submitted is not properly formated")
    }
    
    mat <- mat[[1]]
    
    #subset points in mat corresponding to group 1 and 2 into seperate matrices
    #reformat them into a dataframe with 1 coloum per dimension per group
    g <- ifelse(groups == unique(groups)[1], TRUE, FALSE)
    m1 <- mat[g,]
    m2 <- mat[!g,]
    dat <- data.frame(m1, m2)
    
    #find all possible combinations of 2 of groups and dimensions
    c <- combn(dat, 2, simplify=FALSE)
    
    #check if any of the dimensions for any of the groups are equal to any
    #other dimension of any other groups
    bools <- lapply(
        c,
        function(x)
            sapply(
                as.list(x)[-1],
                FUN = setequal,
                as.list(x)[[1]]
            )
    )
    
    #check that no dimensions are equal to any other dimension
    ident <- all(bools == FALSE)
    
    return(ident)
}

#' Check Identical Matrices
#'
#' This function takes an input matrix (mat), typically from \code{\link{makeMarix}},
#' and an array contaning other matrices (otherMatrices) of the same dimensionality.
#' The function checks if the matrix is identical to any of the other matrices submitted.
#' This was designed for running replicate experiments when each replicate should be unique
#' and is used with the \code{\link{overlapMatrices}}, \code{\link{overlapProbability}}, and
#' \code{\link{zeroOrHundred}} functions.
#'
#' @name checkIdenticalMatrix
#' @rdname checkIdenticalMatrix
#' @aliases checkIdenticalMatrix
#' @param mat The matrix to be checked.
#' @param otherMatrices an array of the other matrices to be checked against.
#' @param ident A boolean value. Do not alter.
#' @author Jason Serviss
#' @keywords checkIdenticalDims
#' @examples
#' mat <- makeMatrix(x=c(0,1), y=c(0,1), dim=2, points=10, seed=11)
#' otherMatrices <- array(NA, dim=c(20, 2, 1))
#' otherMatrices[,,1] <- mat[[1]]
#' checkIdenticalMatrix(mat, otherMatrices)
#'
#' otherMatrices <- array(NA, dim=c(20, 2, 1))
#' otherMatrices[,,1] <- makeMatrix(x=c(0,1), y=c(0,1), dim=2, points=10, seed=50)[[1]]
#' checkIdenticalMatrix(mat, otherMatrices)
NULL
#' @export

checkIdenticalMatrix <- function(
    mat,
    otherMatrices,
    ident = FALSE
) {
    
    if(length(mat[[1]]) == 1){
        warning("the matrix you submitted is not properly formated")
    }
    
    mat <- mat[[1]]
    
    #checks if mat is equal to any matrix in the otherMatrices variable
    bools <- sapply(
        1:dim(otherMatrices)[3],
        function(x)
            all((otherMatrices[,,x] == mat) == TRUE)
    )
    
    ident <- all(bools == FALSE)
    
    return(ident)
}

#' Check Identical Points
#'
#' This function takes an input matrix and checks if, for each group, the frequency
#' of identical data points is greater than what is specified by the \emph{identicalPointsThreshold}
#' argumet. If the number of identical points are less than the identicalPointsThreshold
#' the function returns TRUE and, otherwise, returns FALSE.
#'
#'
#' @name checkIdenticalPoints
#' @rdname checkIdenticalPoints
#' @aliases checkIdenticalPoints
#' @param mat a matrix from makeMatrix function.
#' @param groups a character vector from makeGroups function.
#' @param dim Which dimension of the matrix to check.
#' @param identicalPointsThreshold Sets the threshold of identical points.
#' @author Jason Serviss
#' @keywords checkIdenticalPoints
#' @examples
#' mat <- makeMatrix(x=c(0,1), y=c(0,1), dim=2, points=10, seed=11)
#' groups <- makeGroups(mat, names=c("grp1", "grp2"))
#' checkIdenticalPoints(mat, groups, 1, 2)
#'
#' mat <- list(rbind(mat[[1]], mat[[1]][20,], mat[[1]][20,]))
#' groups <- makeGroups(mat, c("g1", "g2"))
#' checkIdenticalPoints(mat, groups, 1, 2)
NULL
#' @export

checkIdenticalPoints <- function(
    mat,
    groups,
    dim,
    identicalPointsThreshold
) {
    
    #input checks
    if(length(mat[[1]]) == 1){
        warning("the matrix you submitted is not properly formated")
    }
    
    if(!is.numeric(identicalPointsThreshold)){
        warning("the identicalPointsThreshold you submitted is not numeric")
    }
    
    #split the matrix into invididual matrices for each group
    mat <- mat[[1]]
    m <- mat[, dim]
    spl <- split(m, groups)
    
    #calculate the frequency for points in each group
    #as well, check if any of the requencies are greater than identicalPointsThreshold
    bools <- logical()
    for(j in 1:length(spl)) {
        freq <- as.data.frame(table(spl[[j]]))$Freq
        bools <- c(bools, any(freq > identicalPointsThreshold))
    }
    
    output <- all(bools == FALSE)
    return(output)
}

############################################################################################################
#                                                                                                          #
#                                               Plotting                                                   #
#                                                                                                          #
############################################################################################################

#' Plot Overlaps
#'
#' Plots the group overlaps in 2 dimensions for a given matrix and groups.
#'
#'
#' @name visualizeOverlaps
#' @rdname visualizeOverlaps
#' @aliases visualizeOverlaps
#' @param mat The matrix to be checked.
#' @param groups A character vector indicating the groups location within the matrix. Can be generated with the makeGroups function.
#' @author Jason Serviss
#' @keywords visualizeOverlaps
#' @examples
#'
#' x1 <- 0
#' x2 <- x1+1
#' mat <- makeMatrix(x=c(x1,x2), y=c(0,1), dim=2, points=10)[[1]]
#' groups <- makeGroups(mat, names=c("grp1", "grp2"))
#'
#' visualizeOverlaps(mat, groups)
#'
#'
NULL
#' @export
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom ggthemes theme_few scale_colour_ptol
#' @importFrom plyr ddply summarize

visualizeOverlaps <- function(
    mat,
    groups,
    cex=1,
    alpha=0.6,
    title=NULL,
    subtitle=NULL,
    plotType="lines"
){
    
    if(length(mat[[1]]) == 1){
        stop("the matrix you submitted is not properly formated")
    }
    
    mat <- mat[[1]]
    mat <- data.frame(
        dim1 = mat[,1],
        dim2 = mat[,2],
        groups = groups
    )
    
    m <- melt(mat)
    
    m2 <- ddply(
        m,
        c(
            "groups",
            "variable"
        ),
        plyr::summarize,
        median = median(value),
        min = min(value),
        max = max(value)
    )
    
    if(plotType == "lines") {
        p <- ggplot(
            m2,
            aes(
                x = groups,
                y = median,
                colour = groups,
                ymin = min,
                ymax = max
            )
        )+
        geom_pointrange()+
        facet_grid(.~variable)+
        theme_few()+
        scale_colour_ptol()+
        theme(
            legend.position="none",
            legend.title=element_blank(),
            legend.text=element_text(size=15),
            axis.title=element_text(size=17),
            axis.text=element_text(size=15),
            plot.title=element_text(
                hjust=0.5,
                family="ArialMT",
                face="bold",
                size=24
            )
        )
        p
    } else {
        p <- ggplot(
            mat,
            aes(
                x=dim1,
                y=dim2,
                colour=groups
            )
        )+
        geom_point(size=cex, alpha=alpha)+
        theme_few()+
        scale_colour_ptol()+
        theme(
            legend.position="top",
            legend.title=element_text(size=17),
            legend.text=element_text(size=15),
            axis.title=element_text(size=17),
            axis.text=element_text(size=15),
            plot.title=element_text(
                hjust=0.5,
                family="ArialMT",
                face="plain",
                size=24
            ),
            plot.subtitle=element_text(
                hjust=0.5,
                family="ArialMT",
                face="italic",
                size=17
            )
        )+
        labs(
            title=title,
            subtitle=subtitle
        )
        p
    }
    return(p)
}






