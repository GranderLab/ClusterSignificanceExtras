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
#' @importFrom plyr ddply

sensSpecPlot <- function(method=NULL) {
    fnTestResults$testType <- "Sensitivity test"
    fpTestResults$testType <- "Specificity test"
    
    fpTestResults$Pcp.BHcorrectedP <- ddply(fpTestResults, "points", summarize, fdr=p.adjust(pValuePcp, method="fdr"))$fdr
    fnTestResults$Pcp.BHcorrectedP <- ddply(fnTestResults, "points", summarize, fdr=p.adjust(pValuePcp, method="fdr"))$fdr

    fpTestResults$Mlp.BHcorrectedP <- ddply(fpTestResults, "points", summarize, fdr=p.adjust(pValueMlp, method="fdr"))$fdr
    fnTestResults$Mlp.BHcorrectedP <- ddply(fnTestResults, "points", summarize, fdr=p.adjust(pValueMlp, method="fdr"))$fdr

    data <- rbind(fnTestResults, fpTestResults)[ ,c("points", "rep", "Pcp.BHcorrectedP", "Mlp.BHcorrectedP", "testType")]
    
    m <- melt(data, id.vars=c("points", "rep", "testType"))
    m$points <- as.numeric(as.character(m$points))
    
    if(method == "Mlp") {
        m <- subset(m, variable == "Mlp.BHcorrectedP")
        p <- ggplot(m, aes(x=as.factor(points), y=-log10(value), colour=testType))+
            geom_boxplot(outlier.size = 1)+
            geom_hline(yintercept=-log10(0.05), linetype=2, colour="darkgrey")+
            labs(
                x="Points per group",
                y="-log10(adjusted p-Value)",
                title="Specificity/Sensitivity test results"
            )+
            theme_few()+
            theme(
                legend.position="top",
                legend.title=element_blank(),
                legend.text=element_text(size=15),
                axis.title=element_text(size=17),
                axis.text=element_text(size=15),
                axis.text.x=element_text(angle=90),
                plot.title=element_text(
                    hjust=0.5,
                    family="ArialMT",
                    face="bold",
                    size=24,
                    margin=margin(b=15)
                ),
                strip.text.x = element_text(size = 17)
            )+
            scale_colour_economist()+
            scale_x_discrete(breaks = round(seq(min(m$points), max(m$points), by = 10),1))
        
        p
        
    } else if(method == "Pcp") {
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
    } else {
        p <- ggplot(m, aes(x=as.factor(points), y=-log10(value), colour=testType))+
            geom_boxplot(outlier.size = 1)+
            geom_hline(yintercept=-log10(0.05), linetype=2, colour="darkgrey")+
            facet_grid(
                .~variable,
                labeller = as_labeller(
                    c(
                        `Pcp.BHcorrectedP` = 'Pcp method',
                        `Mlp.BHcorrectedP` = 'Mlp method'
                    )
                )
            )+
            labs(
                x="Points per group",
                y="-log10(adjusted p-Value)",
                title="Specificity/Sensitivity test results"
            )+
            theme_few()+
            theme(
                legend.position="top",
                legend.title=element_blank(),
                legend.text=element_text(size=15),
                axis.title=element_text(size=17),
                axis.text=element_text(size=15),
                axis.text.x=element_text(angle=90),
                plot.title=element_text(
                    hjust=0.5,
                    family="ArialMT",
                    face="bold",
                    size=24,
                    margin=margin(b=15)
                ),
                strip.text.x = element_text(size = 17)
            )+
            scale_colour_economist()+
            scale_x_discrete(breaks = round(seq(min(m$points), max(m$points), by = 10),1))
        
        p
    }
    return(p)
}