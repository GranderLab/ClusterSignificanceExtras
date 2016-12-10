context("Matrix Tools")

test_that("check that the makeMatrix function returns the correct output", {
  
  ####TEST1####
  ##prepare normal input data
  x1 <- 1.1
  x2 <- x1 + 1
  y1 <- 0
  y2 <- y1 + 1
  dim <- 2
  points <- 10

  #setup expected data
  min.g1 <- 1.1
  min.g2 <- 0
  max.g1 <- 2.1
  max.g2 <- 1
  
  ##run function
  matOut <- makeMatrix(x=c(x1, x2), y=c(y1, y2), dim=dim, points=points)
  mat <- matOut[[1]]
  x <- matOut[[2]]
  y <- matOut[[3]]
  
  groups <- makeGroups(matOut, names=c("grp1", "grp2"))
  split <- split(mat, groups)
  
  ##test
  expect_equal(x1, x[1])
  expect_equal(x2, x[2])
  expect_equal(y1, y[1])
  expect_equal(y2, y[2])
  expect_true(min.g1 <= min(split[[1]]))
  expect_true(min.g2 <= min(split[[2]]))
  expect_true(max.g1 >= max(split[[1]]))
  expect_true(max.g2 >= max(split[[2]]))
})

test_that("check that the makeGroups function returns the correct output", {
    
    ####TEST1####
    ##prepare normal input data
    mat <- list(matrix(rep(0, 8), ncol=2))
    
    #setup expected data
    exp <- c(rep("grp1", 2), rep("grp2", 2))
    
    ##run function
    output <- makeGroups(mat, names=c("grp1", "grp2"))
    
    ##test
    expect_equal(exp, output)
    
    ####TEST2 bad matrix format####
    ##prepare normal input data
    mat <- matrix(rep(0, 8), ncol=2)
    
    ##test
    expect_error(makeGroups(mat, names=c("grp1", "grp2")))
    
    ####TEST3 short names####
    ##prepare normal input data
    mat <- list(matrix(rep(0, 8), ncol=2))
    
    ##test
    expect_error(makeGroups(mat, names=c("grp1")))
})

test_that("check that the dynamicXY function returns the correct output", {
    
    ####TEST1 100####
    ##prepare normal input data
    desiredOverlap <- 100
    
    #setup expected data
    exp <- list(c(1, 100), c(0, 101))
    
    ##run function
    output <- dynamicXY(desiredOverlap)
    
    ##test
    expect_equal(exp, output)
    
    ####TEST2 0####
    ##prepare normal input data
    desiredOverlap <- 0
    
    #setup expected data
    exp <- list(c(101, 201), c(0, 101))
    
    ##run function
    output <- dynamicXY(desiredOverlap)
    
    ##test
    expect_equal(exp, output)

    ####TEST3 message####
    ##prepare normal input data
    desiredOverlap <- 1000
    
    ##test
    expect_error(dynamicXY(desiredOverlap))
})

test_that("check that the addOnePoint function returns the correct output", {
    
    ####TEST1####
    ##prepare normal input data
    mat <- matrix(rep(0,8), ncol=2)
    x <- c(0,100)
    y <- c(0,100)
    mat <- list(mat, x, y)
    groups <- c(rep("1", 2), rep("2", 2))
    seed <- 3
    
    ##run function
    output <- addOnePoint(mat, groups, seed)[[1]]
    
    ##test
    expect_equal(nrow(output), 6)
    expect_equal(ncol(output), 2)

})

test_that("check that the calculateOverlap function returns the correct output", {
    
    ####TEST1 0% overlap####
    ##prepare normal input data
    g1 <- seq(0.0, 1.0, 0.1)
    g2 <- seq(1.0, 2.0, 0.1)
    mat <- list(matrix(c(g1, g2), ncol=1))
    groups <- c(rep("grp1", length(g1)), rep("grp2", length(g2)))
    
    #setup expected data
    expected <- 0
    
    ##run function
    output <- calculateOverlap(mat, groups, 1)
    
    ##test
    expect_equal(output, expected)
    
    ####TEST2 100% overlap####
    ##prepare normal input data
    g1 <- seq(0.0, 10, 1.0)
    g2 <- seq(0.1, 1.1, 0.1)
    mat <- list(matrix(c(g1, g2), ncol=1))
    groups <- c(rep("grp1", length(g1)), rep("grp2", length(g2)))
    
    #setup expected data
    expected <- 100
    
    ##run function
    output <- calculateOverlap(mat, groups, 1)
    
    ##test
    expect_equal(output, expected)
    
    ####TEST3 50% overlap####
    ##prepare normal input data
    g1 <- seq(0.0, 1.0, 1.0)
    g2 <- seq(0.5, 1.5, 1.0)
    mat <- list(matrix(c(g2, g1), ncol=1))
    groups <- c(rep("grp1", length(g1)), rep("grp2", length(g2)))
    
    #setup expected data
    expected <- 50
    
    ##run function
    output <- calculateOverlap(mat, groups, 1)
    
    ##test
    expect_equal(output, expected)
    
    ####TEST4 100% overlap group switch####
    ##prepare normal input data
    g1 <- seq(0.0, 10, 1.0)
    g2 <- seq(0.1, 1.1, 0.1)
    mat <- list(matrix(c(g2, g1), ncol=1))
    groups <- c(rep("grp1", length(g1)), rep("grp2", length(g2)))
    
    #setup expected data
    expected <- 100
    
    ##run function
    output <- calculateOverlap(mat, groups, 1)
    
    ##test
    expect_equal(output, expected)
    
    ####TEST5 0% overlap group switch####
    ##prepare normal input data
    g1 <- seq(0.0, 1.0, 0.1)
    g2 <- seq(1.0, 2.0, 0.1)
    mat <- list(matrix(c(g2, g1), ncol=1))
    groups <- c(rep("grp1", length(g1)), rep("grp2", length(g2)))
    
    #setup expected data
    expected <- 0
    
    ##run function
    output <- calculateOverlap(mat, groups, 1)
    
    ##test
    expect_equal(output, expected)
    
    ####TEST6 wrong matrix format####
    ##prepare normal input data
    g1 <- seq(0.0, 1.0, 0.1)
    g2 <- seq(1.0, 2.0, 0.1)
    mat <- matrix(c(g2, g1), ncol=1)
    groups <- c(rep("grp1", length(g1)), rep("grp2", length(g2)))
    
    expect_error(calculateOverlap(mat, groups, 1))
})

test_that("check that the checkIdenticalDims function returns the correct output", {
    
    ####TEST1 all identical; FALSE####
    ##prepare normal input data
    dim1 <- seq(1, 10, 1)
    dim2 <- seq(1, 10, 1)
    vec <- rep(c(dim1, dim2), 2)
    mat <- matrix(vec, ncol=2)
    groups <- c(rep("1", nrow(mat)/2), rep("2", nrow(mat)/2))
    
    ##run function
    output <- checkIdenticalDims(list(mat), groups)
    
    ##test
    expect_false(output)
    
    ####TEST2 identical within groups; FALSE####
    ##prepare normal input data
    dim1 <- seq(1, 10, 1)
    dim2 <- seq(2, 11, 1)
    vec <- rep(c(dim1, dim2), 2)
    mat <- matrix(vec, ncol=2)
    groups <- c(rep("1", nrow(mat)/2), rep("2", nrow(mat)/2))
    
    ##run function
    output <- checkIdenticalDims(list(mat), groups)
    
    ##test
    expect_false(output)
    
    ####TEST3 one dim TRUE one dim; FALSE####
    ##prepare normal input data
    dim1a <- seq(1, 10, 1)
    dim2a <- seq(1, 10, 1)
    dim1b <- seq(1, 10, 1)
    dim2b <- seq(2, 11, 1)
    vec <- c(dim1a, dim2a, dim1b, dim2b)
    mat <- matrix(vec, ncol=2)
    groups <- c(rep("1", nrow(mat)/2), rep("2", nrow(mat)/2))
    
    ##run function
    output <- checkIdenticalDims(list(mat), groups)
    
    ##test
    expect_false(output)
    
    ####TEST4 reverse one; FALSE####
    ##prepare normal input data
    dim1a <- seq(1, 10, 1)
    dim2a <- seq(10, 1, -1)
    dim1b <- seq(3, 12, 1)
    dim2b <- seq(2, 11, 1)
    vec <- c(dim1a, dim2a, dim1b, dim2b)
    mat <- matrix(vec, ncol=2)
    groups <- c(rep("1", nrow(mat)/2), rep("2", nrow(mat)/2))
    
    ##run function
    output <- checkIdenticalDims(list(mat), groups)
    
    ##test
    expect_false(output)
    
    ####TEST5 none identical TRUE####
    ##prepare normal input data
    dim1a <- seq(1, 10, 1)
    dim2a <- seq(2, 11, 1)
    dim1b <- seq(3, 12, 1)
    dim2b <- seq(4, 13, 1)
    vec <- c(dim1a, dim2a, dim1b, dim2b)
    mat <- matrix(vec, ncol=2)
    groups <- c(rep("1", nrow(mat)/2), rep("2", nrow(mat)/2))
    
    ##run function
    output <- checkIdenticalDims(list(mat), groups)
    
    ##test
    expect_true(output)
})

test_that("check that the checkIdenticalMatrix function returns the correct output", {
    
    ####TEST1 FALSE####
    ##prepare normal input data
    dim1 <- seq(1, 10, 1)
    dim2 <- seq(1, 10, 1)
    vec <- rep(c(dim1, dim2), 2)
    mat <- matrix(vec, ncol=2)
    
    otherMatrices <- array(c(rep(0, 40)), dim=c(nrow(mat), ncol(mat), 1))
    otherMatrices[,,1] <- mat
    
    ##run function
    output <- checkIdenticalMatrix(list(mat), otherMatrices)
    
    ##test
    expect_false(output)
    
    ####TEST2 TRUE####
    ##prepare normal input data
    dim1 <- seq(1, 10, 1)
    dim2 <- seq(2, 11, 1)
    vec <- rep(c(dim1, dim2), 2)
    mat <- matrix(vec, ncol=2)
    
    otherMatrices <- array(c(rep(0, 40)), dim=c(nrow(mat), ncol(mat), 1))

    ##run function
    output <- checkIdenticalMatrix(list(mat), otherMatrices)
    
    ##test
    expect_true(output)
    
    ####TEST3 multiple others FALSE####
    ##prepare normal input data
    dim1 <- seq(1, 10, 1)
    dim2 <- seq(1, 10, 1)
    vec <- rep(c(dim1, dim2), 2)
    mat <- matrix(vec, ncol=2)
    
    otherMatrices <- array(c(rep(0, 40)), dim=c(nrow(mat), ncol(mat), 2))
    otherMatrices[,,1] <- mat
    otherMatrices[,,2] <- mat

    ##run function
    output <- checkIdenticalMatrix(list(mat), otherMatrices)
    
    ##test
    expect_false(output)
    
    ####TEST4 multiple others TRUE####
    ##prepare normal input data
    dim1 <- seq(1, 10, 1)
    dim2 <- seq(2, 11, 1)
    vec <- rep(c(dim1, dim2), 2)
    mat <- matrix(vec, ncol=2)
    
    otherMatrices <- array(c(rep(0, 80)), dim=c(nrow(mat), ncol(mat), 2))
    
    ##run function
    output <- checkIdenticalMatrix(list(mat), otherMatrices)
    
    ##test
    expect_true(output)
    
    ####TEST5 only one FALSE####
    ##prepare normal input data
    dim1 <- seq(1, 10, 1)
    dim2 <- seq(1, 10, 1)
    vec <- rep(c(dim1, dim2), 2)
    mat <- matrix(vec, ncol=2)
    
    otherMatrices <- array(c(rep(0, 40)), dim=c(nrow(mat), ncol(mat), 2))
    otherMatrices[,,1] <- mat
    
    ##run function
    output <- checkIdenticalMatrix(list(mat), otherMatrices)
    
    ##test
    expect_false(output)
    
    ####TEST6 send inproper mat format####
    ##prepare normal input data
    dim1 <- seq(1, 10, 1)
    dim2 <- seq(1, 10, 1)
    vec <- rep(c(dim1, dim2), 2)
    mat <- matrix(vec, ncol=2)
    
    otherMatrices <- array(c(rep(0, 40)), dim=c(nrow(mat), ncol(mat), 2))
    otherMatrices[,,1] <- mat
    
    ##test
    expect_warning(checkIdenticalMatrix(mat, otherMatrices))
})


test_that("check that the overlapProbability function returns the correct output", {
    
    ####TEST1####
    #setup input
    x1 <- seq(0, 101, by=1)
    points <- 5
    reps <- 2
    by <- 0.9
    
    #setup expected data
    out.nrow <- 2 * 102
    out.ncol <- 6
    x1.s <- sort(rep(seq(0, 101, by=1) , 2))
    overlaps <- c(0, 20, 30, 40, 50, 60, 70, 80, 100) ##all possible overlaps with 5 points

    ##run function
    output <- overlapProbability(points=points, reps=reps, x1=x1, by=by, save=FALSE)
    mat <- output[1:1, 'mat'][[1]]
    groups <- makeGroups(mat, names=c("1", "2"))
    spl <- split(mat[[1]], groups)
    
    ##test
    expect_equal(nrow(output), out.nrow)
    expect_equal(ncol(output), out.ncol)
    expect_equal(output$x1, x1.s)
    expect_true(all(output$points == 5))
    expect_true(all(output$overlapPC1 %in% overlaps))
    expect_true(all(output$overlapPC2 %in% overlaps))
    expect_false(identical(output[1:1, 'mat'], output[2:2, 'mat']))
    expect_identical(unique(output$reps), 1:reps)
    expect_identical(unique(output$points), points)
    expect_identical(unique(output$x1), seq(0, 101, by=1))
})


test_that("check that the interspersion function returns the correct output", {
    
    ####TEST1 0 percent####
    ##prepare normal input data
    mat <- list(matrix(c(seq(1,10,1), 0.1, seq(11,19,1)), ncol=1))
    groups <- c(rep("1", 10), rep("2", 10))
    
    ##run function
    output <- interspersion(mat, groups, 1)
    
    ##test
    expect_equal(output, 0)
    
    ####TEST2 100 percent####
    ##prepare normal input data
    mat <- list(matrix(seq(0,19, by=1), ncol=1))
    groups <- rep(c("1", "2"), 10)
    
    ##run function
    output <- interspersion(mat, groups, 1)
    
    ##test
    expect_equal(output, 100)
    
    ####TEST3 wrong matrix format####
    ##prepare normal input data
    mat <- matrix(seq(0,19, by=1), ncol=1)
    groups <- rep(c("1", "2"), 10)
    
    ##test
    expect_error(interspersion(mat, groups, 1))
})

test_that("check that the matrixTests function returns the correct output", {
    
    ####TEST1 0 percent overlap####
    ##prepare normal input data
    g1 <- c(seq(1,10,1), 10, 10)
    g2 <- c(seq(1,10,1), 10, 10)
    vec <- rep(c(g1, g2), 2)
    
    mat <- matrix(vec, ncol=2)
    mat <- list(mat, c(1,10), c(11,20))
    class(mat) <- "makeMatrix"
    
    tests <- c(
        "checkIdenticalDims",
        "checkIdenticalMatrix",
        "overlapSpecificMatrix",
        "checkIdenticalPoints"
    )
    
    groups <- c(rep("1", 12), rep("2", 12))
    dim <- 2
    points <- 12
    by=1
    seed=3
    otherMatrices <- array(vec, dim=c(24, 2, 1))
    desiredOverlap <- 0
    interspersionThreshold <- ""
    idPointsThreshold <- 2
    verbose <- FALSE
    
    ##run function
    output <- matrixTests(
        tests,
        mat,
        groups,
        dim,
        points,
        by,
        seed,
        otherMatrices,
        desiredOverlap,
        interspersionThreshold,
        idPointsThreshold,
        verbose
    )[[1]]
    
    ##test
    expect_equal(calculateOverlap(output, groups, 1), desiredOverlap)
    expect_equal(calculateOverlap(output, groups, 2), desiredOverlap)
    expect_true(checkIdenticalDims(output, groups))
    expect_true(checkIdenticalMatrix(output, otherMatrices))
    expect_true(checkIdenticalPoints(output, groups, 1, idPointsThreshold))
    expect_true(checkIdenticalPoints(output, groups, 2, idPointsThreshold))
    
    ####TEST2 100 percent overlap####
    ##prepare normal input data
    g1 <- c(seq(1,10,1), 10, 10)
    g2 <- c(seq(11,20,1), 11, 11)
    vec <- rep(c(g1, g2), 2)
    
    mat <- matrix(vec, ncol=2)
    mat <- list(mat, c(1,10), c(11,20))
    class(mat) <- "makeMatrix"
    
    tests <- c(
    "checkIdenticalDims",
    "checkIdenticalMatrix",
    "overlapSpecificMatrix",
    "interspersion",
    "checkIdenticalPoints"
    )
    
    groups <- c(rep("1", 12), rep("2", 12))
    dim <- 2
    points <- 12
    by=1
    seed=3
    otherMatrices <- array(vec, dim=c(24, 2, 1))
    desiredOverlap <- 100
    interspersionThreshold <- 50
    idPointsThreshold <- 2
    verbose <- FALSE
    
    ##run function
    output <- matrixTests(
        tests,
        mat,
        groups,
        dim,
        points,
        by,
        seed,
        otherMatrices,
        desiredOverlap,
        interspersionThreshold,
        idPointsThreshold,
        verbose
    )[[1]]
    
    ##test
    expect_equal(calculateOverlap(output, groups, 1), desiredOverlap)
    expect_equal(calculateOverlap(output, groups, 2), desiredOverlap)
    expect_true(interspersion(output, groups, 1) > interspersionThreshold)
    expect_true(interspersion(output, groups, 2) > interspersionThreshold)
    expect_true(checkIdenticalDims(output, groups))
    expect_true(checkIdenticalMatrix(output, otherMatrices))
    expect_true(checkIdenticalPoints(output, groups, 1, idPointsThreshold))
    expect_true(checkIdenticalPoints(output, groups, 2, idPointsThreshold))
})

test_that("check that the checkIdenticalPoints function returns the correct output", {
    
    ####TEST1 TRUE####
    ##prepare normal input data
    vec <- seq(1,10,1)
    mat <- matrix(rep(vec, 4), ncol=2)
    mat <- list(mat)
    groups <- makeGroups(mat, names=c("1", "2"))
    dim <- 1
    identicalPointsThreshold <- 1
    
    ##run function
    output <- checkIdenticalPoints(mat, groups, dim, identicalPointsThreshold)
    
    ##test
    expect_true(output)
    
    ####TEST2 FALSE####
    ##prepare normal input data
    vec <- rep(seq(1,5,1),2)
    mat <- matrix(rep(vec, 4), ncol=2)
    mat <- list(mat)
    groups <- makeGroups(mat, names=c("1", "2"))
    dim <- 1
    identicalPointsThreshold <- 1
    
    ##run function
    output <- checkIdenticalPoints(mat, groups, dim, identicalPointsThreshold)
    
    ##test
    expect_false(output)
    
    ####TEST3 identicalPointsThreshold warning####
    ##prepare normal input data
    vec <- rep(seq(1,5,1),2)
    mat <- matrix(rep(vec, 4), ncol=2)
    mat <- list(mat)
    groups <- makeGroups(mat, names=c("1", "2"))
    dim <- 1
    identicalPointsThreshold <- ""
    
    ##test
    expect_warning(checkIdenticalPoints(mat, groups, dim, identicalPointsThreshold))
})


























