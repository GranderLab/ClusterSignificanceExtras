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

    data <- rbind(sensitivityTestResults, specificityTestResults)[ ,c("points", "rep", "Pcp.BHcorrectedP", "testType")]
    
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