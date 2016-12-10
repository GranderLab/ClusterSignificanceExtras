#' Output from the falseNegativeTest function with default arguments.
#'
#'
#' @title False Negative Test Results
#' @docType data
#' @name fnTestResults
#' @format data.frame
#' \describe{
#' \item{points}{The number of points per group.}
#' \item{rep}{The repetition number for a specific number of points.}
#' \item{pValuePcp}{The resulting p-value from ClusterSignificance using the Pcp method.}
#' \item{pValueMlp}{The resulting p-value from ClusterSignificance using the Mlp method.}
#' \item{Pcp.scores.real}{The resulting scores.real value from ClusterSignificance using the Pcp method.}
#' \item{Pcp.scores.vec}{The resulting scores.vec values from ClusterSignificance using the Pcp method.}
#' \item{Mlp.scores.real}{The resulting scores.real value from ClusterSignificance using the Mlp method.}
#' \item{Mlp.scores.vec}{The resulting scores.vec values from ClusterSignificance using the Mlp method.}
#' }
#' @usage fnTestResults
#' @return data.frame
#' @examples
#' fnTestResults
#'

NULL