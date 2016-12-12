#' tSNE results from the hematological malagnancies anaysis.
#'
#'
#' @title tSNE results dataset.
#' @docType data
#' @name hemCancData
#' @format data.frame
#' \describe{
#' \item{X1}{The first tSNE dimension.}
#' \item{X2}{The second tSNE dimension.}
#' \item{X3}{The third tSNE dimension.}
#' \item{geo_accession}{The sample GEO accession number.}
#' \item{characteristics_ch1.1}{The samples cancer type.}
#' }
#' @usage hemCancData
#' @return data.frame
#' @examples
#' hemCancData
NULL

#' Output from the zeroOrHundred function with overlap argument = 100 and other default arguments.
#'
#'
#' @title Hundred overlap dataset.
#' @docType data
#' @name hundred
#' @format list of arrays
#' \describe{
#' \item{list names}{The number of points per group.}
#' \item{list elements}{An array containing all repititions of the ClusterSignificance input matrix for the specific points number indicated in the list element name}
#' \item{array dimension 1}{Rows of the input matrix to ClusterSignificance.}
#' \item{array dimension 2}{Columns of the input matrix to ClusterSignificance.}
#' \item{array dimension 3}{The repititions for each points per group amounts.}
#' }
#' @usage hundred
#' @return list of arrays
#' @examples
#' hundred
#'

NULL

#' A character vector of lncRNA ensembl IDs included in the hematological cancer analysis.
#'
#'
#' @title lncRNAs dataset.
#' @docType data
#' @name lncGenes
#' @format character vector
#' @usage lncGenes
#' @return character vector
#' @examples
#' lncGenes
NULL

#' ClusterSignificance permute output from the hematological cancer dataset.
#'
#'
#' @title Hematological malagnancies ClusterSignificance result.
#' @docType data
#' @name pe
#' @format PermutationResults
#' @usage pe
#' @return PermutationResults
#' @examples
#' pe
NULL

#' Output from the sensitivityTest function with default arguments.
#'
#'
#' @title Sensitivity test results dataset.
#' @docType data
#' @name sensitivityTestResults
#' @format data.frame
#' \describe{
#' \item{points}{The number of points per group.}
#' \item{rep}{The repitition.}
#' \item{pValuePcp}{The pValue from ClusterSignificance with 10^4 iterations.}
#' }
#' @usage sensitivityTestResults
#' @return data.frame
#' @examples
#' sensitivityTestResults
NULL

#' Output from the specificityTest function with default arguments.
#'
#'
#' @title Specificity test results dataset.
#' @docType data
#' @name specificityTestResults
#' @format data.frame
#' \describe{
#' \item{points}{The number of points per group.}
#' \item{rep}{The repitition.}
#' \item{pValuePcp}{The pValue from ClusterSignificance with 10^4 iterations.}
#' }
#' @usage specificityTestResults
#' @return data.frame
#' @examples
#' specificityTestResults
NULL

#' Output from the zeroOrHundred function with overlap argument = 0 and other default arguments.
#'
#'
#' @title Zero overlap dataset.
#' @docType data
#' @name zero
#' @format list of arrays
#' \describe{
#' \item{list names}{The number of points per group.}
#' \item{list elements}{An array containing all repititions of the ClusterSignificance input matrix for the specific points number indicated in the list element name}
#' \item{array dimension 1}{Rows of the input matrix to ClusterSignificance.}
#' \item{array dimension 2}{Columns of the input matrix to ClusterSignificance.}
#' \item{array dimension 3}{The repititions for each points per group amounts.}
#' }
#' @usage zero
#' @return list of arrays
#' @examples
#' zero
#'

NULL