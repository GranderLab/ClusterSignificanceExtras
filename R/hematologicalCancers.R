#' Hematological Cancers
#'
#' Downloads, processes, clusters, and runs ClusterSignificance on
#' the hematological malagnancies dataset.
#'
#' @name hematologicalCancers
#' @rdname hematologicalCancers
#' @aliases hematologicalCancers
#' @param iterations Number of permutations.
#' @param ... Arguments to pass on.
#' @author Jason Serviss
#' @keywords hematologicalCancers
#' @examples
#' hematologicalCancers()
NULL
#' @export
#' @import GEOquery
#' @import biomaRt
#' @importFrom Rtsne Rtsne

hematologicalCancers <- function(iterations=10000, ...) {
    
    tmp <- .downloadAndProcessData()
    exp <- tmp[[1]]
    pheno <- tmp[[2]]
    
    exp <- .annotateBiotypes(exp)
    nc <- .filterLongNonCoding(exp)
    lncGenes <- nc$ensembl_gene_id
    pheno <- .renamePhenos(pheno)
    plot <- .runTsne(nc, pheno)
    
    tmp <- .runPcp(plot)
    mat <- tmp[[1]]
    groups <- tmp[[2]]
    prj <- tmp[[3]]
    group.color <- tmp[[4]]
    
    #plot, group.color needed for t-sne plots
    #prj needed for projection plots
    
    cl <- classify(prj)
    
    #cl needed for classification plots
    
    pe <- permute(mat, groups, iter=iterations, projmethod="pcp")
    pValues <- as.data.frame(pvalue(pe))
    colnames(pValues) <- "pValue"
    
    #pValues needed for pValue table
    
    return(list(plot, group.color, prj, cl, pe, pValues, mat, groups, nc, lncGenes))
    
}


.downloadAndProcessData <- function() {
    gset <- getGEO("GSE13159", GSEMatrix =TRUE)[[1]]
    
    #extract expression values, remove incomplete cases, and get phenotype data
    fvarLabels(gset) <- make.names(fvarLabels(gset))
    ex <- as.data.frame(exprs(gset))
    exp <- ex[complete.cases(ex),]
    pheno <- pData(phenoData(gset))
    
    return(list(exp, pheno))
}

.annotateBiotypes <- function(exp) {
    
    ensembl = biomaRt::useMart(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset="hsapiens_gene_ensembl",
        host = "jul2015.archive.ensembl.org"
    )
    
    affyids = rownames(exp)
    
    IDs <- biomaRt::getBM(
        attributes=c(
            'affy_hg_u133_plus_2',
            'ensembl_gene_id',
            'external_gene_name',
            'chromosome_name',
            'start_position',
            'end_position',
            'gene_biotype'
        ),
        filters = 'affy_hg_u133_plus_2',
        values = affyids,
        mart = ensembl
    )
    
    colnames(IDs) = c(
        "probe_id",
        "ensembl_gene_id",
        "external_gene_name",
        "chromosome_name",
        "start_position",
        "end_position",
        "gene_biotype"
    )
    
    #remove probes that are duplicated if any
    IDs = subset(
        IDs,
        duplicated(IDs$probe_id) == FALSE &
        duplicated(IDs$probe_id, fromLast=TRUE) == FALSE
    )
    
    #merge annotation with expression data and extract lncRNA biotypes
    exp <- merge(exp, IDs, by.x=0, by.y="probe_id")
    
    return(exp)
}

.filterLongNonCoding <- function(exp) {
    #specify long non-coding
    ncGenes <- c(
        "3prime_overlapping_ncrna",
        "antisense",
        "IG_V_pseudogene",
        "lincRNA",
        "misc_RNA",
        "polymorphic_pseudogene",
        "processed_pseudogene",
        "processed_transcript",
        "pseudogene",
        "sense_intronic",
        "sense_overlapping",
        "TR_V_pseudogene",
        "transcribed_processed_pseudogene",
        "transcribed_unitary_pseudogene",
        "transcribed_unprocessed_pseudogene",
        "unitary_pseudogene",
        "unprocessed_pseudogene"
    )
    
    #extract long non-coding
    nc <- exp[exp$gene_biotype %in% ncGenes, ]
    rownames(nc) <- nc$Row.names
    nc$Row.names <- NULL
    return(nc)
}

.renamePhenos <- function(pheno) {
    #change names of genetic subtypes to something more readable
    pheno$characteristics_ch1.1 <- gsub(
        "leukemia class: (.*)",
        "\\1",
        pheno$characteristics_ch1.1
    )
    
    #concatenate genetic subtypes into their respective cancer types
    #this is done due to the fact that, if the individual genetic subtypes were retained, ClusterSignificance would need to make 153 comparisons instead of 21.
    #153 comparisons is way too many for the publication and will provide a overwhelming amount of results, potentially being hard to interpret.
    
    Ball <- c("ALL with hyperdiploid karyotype", "ALL with t(1;19)", "ALL with t(12;21)", "c-ALL/Pre-B-ALL with t(9;22)", "c-ALL/Pre-B-ALL without t(9;22)", "mature B-ALL with t(8;14)", "Pro-B-ALL with t(11q23)/MLL")
    aml <- c("AML complex aberrant karyotype", "AML with inv(16)/t(16;16)", "AML with normal karyotype + other abnormalities", "AML with t(11q23)/MLL", "AML with t(15;17)", "AML with t(8;21)")
    
    
    pheno$characteristics_ch1.1 <- ifelse(
        pheno$characteristics_ch1.1 %in% Ball,
        "B-ALL",
        ifelse(
            pheno$characteristics_ch1.1 %in% aml,
            "AML",
            ifelse(
                pheno$characteristics_ch1.1 == "T-ALL",
                "T-ALL",
                ifelse(
                    pheno$characteristics_ch1.1 == "Non-leukemia and healthy bone marrow",
                    "Non-leukemia and healthy bone marrow",
                    ifelse(
                        pheno$characteristics_ch1.1 == "CML",
                        "CML",
                        ifelse(
                            pheno$characteristics_ch1.1 == "CLL",
                            "CLL",
                            "MDS"
                        )
                    )
                )
            )
        )
    )
    
    pheno$characteristics_ch1.1 <- ifelse(
        pheno$characteristics_ch1.1 == "Non-leukemia and healthy bone marrow",
        "normal",
        pheno$characteristics_ch1.1
    )
    
    return(pheno)
}

.runTsne <- function(nc, pheno) {
    tsne <- Rtsne(t(
        nc[ ,!colnames(nc) %in% c(
            "Row.names",
            "probe_id",
            "ensembl_gene_id",
            "external_gene_name",
            "chromosome_name",
            "start_position",
            "end_position",
            "gene_biotype"
        )]
    ), dims=3)
    
    #put tSNE results and phenotype data in a data.frame
    plot <- data.frame(tsne$Y, pheno)
    return(plot)
}

.runPcp <- function(plot) {
    mat <- as.matrix(plot[ ,c("X1", "X2", "X3")])
    groups <- plot$characteristics_ch1.1
    prj <- pcp(mat, groups)
    group.color <- getData(prj, "group.color")
    return(list(mat, groups, prj, group.color))
}

#' tSNE plots
#'
#' Performs the tSNE plots for the hematological malagnancies dataset.
#'
#' @name tsnePlots
#' @rdname tsnePlots
#' @aliases tsnePlots
#' @param ... Arguments to pass on.
#' @author Jason Serviss
#' @keywords tsnePlots
#' @examples
#' #tsnePlots()
NULL
#' @export
#' @import scatterplot3d

tsnePlots <- function(plot, dim) {
    mat <- as.matrix(plot[ ,c("X1", "X2", "X3")])
    groups <- plot$characteristics_ch1.1
    prj <- pcp(mat, groups)
    group.color <- getData(prj, "group.color")

    if(dim == "X") {
        .dimX(plot, group.color, groups)
    } else if(dim == "Y") {
        .dimY(plot, group.color, groups)
    } else if(dim == "Z") {
        .dimZ(plot, group.color, groups)
    } else {
        stop("The specified dim could not be found.")
    }
}

.dimX <- function(plot, group.color, groups) {
    p <- scatterplot3d(
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X1"],
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X2"],
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X3"],
        color=rgb(
            group.color["red",   colnames(group.color)[1]],
            group.color["green", colnames(group.color)[1]],
            group.color["blue",  colnames(group.color)[1]],
            150,
            maxColorValue=255
        ),
        pch=16,
        box=FALSE,
        xlab="X",
        ylab="Y",
        zlab="Z",
        xlim=c(min(plot$X1), max(plot$X1)),
        ylim=c(min(plot$X2), max(plot$X2)),
        zlim=c(min(plot$X3), max(plot$X3))
    
    )
    
    invisible(
        sapply(2:length(unique(groups)), function(x)
            p$points3d(
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X1"],
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X2"],
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X3"],
                col=rgb(
                    group.color["red",   colnames(group.color)[x]],
                    group.color["green", colnames(group.color)[x]],
                    group.color["blue",  colnames(group.color)[x]],
                    150,
                    maxColorValue=255
                ),
                pch=16,
                cex = 1
            )
        )
    )
}

.dimY <- function(plot, group.color, groups) {
    p <- scatterplot3d(
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X2"],
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X1"],
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X3"],
        color=rgb(
            group.color["red",   colnames(group.color)[1]],
            group.color["green", colnames(group.color)[1]],
            group.color["blue",  colnames(group.color)[1]],
            150,
            maxColorValue=255
        ),
        pch=16,
        box=FALSE,
        xlab="Y",
        ylab="X",
        zlab="Z",
        xlim=c(min(plot$X2), max(plot$X2)),
        ylim=c(min(plot$X1), max(plot$X1)),
        zlim=c(min(plot$X3), max(plot$X3))
    )
    
    invisible(
        sapply(2:length(unique(groups)), function(x)
            p$points3d(
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X2"],
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X1"],
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X3"],
                col=rgb(
                    group.color["red",   colnames(group.color)[x]],
                    group.color["green", colnames(group.color)[x]],
                    group.color["blue",  colnames(group.color)[x]],
                    150,
                    maxColorValue=255
                ),
                pch=16,
                cex = 1
            )
        )
    )
}

.dimZ <- function(plot, group.color, groups) {
    p <- scatterplot3d(
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X3"],
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X2"],
        plot[plot$characteristics_ch1.1 == colnames(group.color)[1], "X1"],
        color=rgb(
            group.color["red",   colnames(group.color)[1]],
            group.color["green", colnames(group.color)[1]],
            group.color["blue",  colnames(group.color)[1]],
            150,
            maxColorValue=255
        ),
        pch=16,
        box=FALSE,
        xlab="Z",
        ylab="Y",
        zlab="X",
        xlim=c(min(plot$X3), max(plot$X3)),
        ylim=c(min(plot$X2), max(plot$X2)),
        zlim=c(min(plot$X1), max(plot$X1))
    )
    
    invisible(
        sapply(2:length(unique(groups)), function(x)
            p$points3d(
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X3"],
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X2"],
                plot[plot$characteristics_ch1.1 == colnames(group.color)[x], "X1"],
                col=rgb(
                    group.color["red",   colnames(group.color)[x]],
                    group.color["green", colnames(group.color)[x]],
                    group.color["blue",  colnames(group.color)[x]],
                    150,
                    maxColorValue=255
                ),
                pch=16,
                cex = 1
            )
        )
    )
}


#' normal vs. MDS
#'
#' Plots the clustering of the normal and MDS cancer types.
#'
#' @name normalMDS
#' @rdname normalMDS
#' @aliases normalMDS
#' @param ... Arguments to pass on.
#' @author Jason Serviss
#' @keywords normalMDS
#' @examples
#' #normalMDS()
NULL
#' @export

normalMDS <- function(plot, view) {
    mat <- as.matrix(plot[ ,c("X1", "X2", "X3")])
    groups <- plot$characteristics_ch1.1
    prj <- pcp(mat, groups)
    group.color <- getData(prj, "group.color")
    inv <- subset(plot, characteristics_ch1.1 == "normal" | characteristics_ch1.1 == "MDS")
    
    if(view == 1) {
        .view1(inv, group.color, groups)
    } else if(view == 2) {
        .view2(inv, group.color, groups)
    } else if(view == 3) {
        .view3(inv, group.color, groups)
    } else {
        stop("The view you requested cannot be found.")
    }
}

.view1 <- function(inv, group.color, groups) {
    p <- scatterplot3d(
        inv[inv$characteristics_ch1.1 == "normal", "X1"],
        inv[inv$characteristics_ch1.1 == "normal", "X2"],
        inv[inv$characteristics_ch1.1 == "normal", "X3"],
        color=rgb(
            group.color["red",   "normal"],
            group.color["green", "normal"],
            group.color["blue",  "normal"],
            150,
            maxColorValue=255
        ),
        pch=16,
        box=FALSE,
        xlab="X",
        ylab="Y",
        zlab="Z",
        xlim=c(min(inv$X1), max(inv$X1)),
        ylim=c(min(inv$X2), max(inv$X2)),
        zlim=c(min(inv$X3), max(inv$X3))
    )
    
    invisible(
        sapply(2:length(unique(groups)), function(x)
            p$points3d(
                inv[inv$characteristics_ch1.1 == "MDS", "X1"],
                inv[inv$characteristics_ch1.1 == "MDS", "X2"],
                inv[inv$characteristics_ch1.1 == "MDS", "X3"],
                col=rgb(
                    group.color["red",   "MDS"],
                    group.color["green", "MDS"],
                    group.color["blue",  "MDS"],
                    150,
                    maxColorValue=255
                ),
                pch=16,
                cex = 1
            )
        )
    )
}

.view2 <- function(inv, group.color, groups) {
    p <- scatterplot3d(
        inv[inv$characteristics_ch1.1 == "normal", "X2"],
        inv[inv$characteristics_ch1.1 == "normal", "X1"],
        inv[inv$characteristics_ch1.1 == "normal", "X3"],
        color=rgb(
            group.color["red",   "normal"],
            group.color["green", "normal"],
            group.color["blue",  "normal"],
            150,
            maxColorValue=255
        ),
        pch=16,
        box=FALSE,
        xlab="X",
        ylab="Y",
        zlab="Z",
        xlim=c(min(inv$X2), max(inv$X2)),
        ylim=c(min(inv$X1), max(inv$X1)),
        zlim=c(min(inv$X3), max(inv$X3))
    )
    
    invisible(
        sapply(2:length(unique(groups)), function(x)
            p$points3d(
                inv[inv$characteristics_ch1.1 == "MDS", "X2"],
                inv[inv$characteristics_ch1.1 == "MDS", "X1"],
                inv[inv$characteristics_ch1.1 == "MDS", "X3"],
                col=rgb(
                    group.color["red",   "MDS"],
                    group.color["green", "MDS"],
                    group.color["blue",  "MDS"],
                    150,
                    maxColorValue=255
                ),
                pch=16,
                cex = 1
            )
        )
    )
}

.view3 <- function(inv, group.color, groups) {
    p <- scatterplot3d(
        inv[inv$characteristics_ch1.1 == "normal", "X3"],
        inv[inv$characteristics_ch1.1 == "normal", "X2"],
        inv[inv$characteristics_ch1.1 == "normal", "X1"],
        color=rgb(
            group.color["red",   "normal"],
            group.color["green", "normal"],
            group.color["blue",  "normal"],
            150,
            maxColorValue=255
        ),
        pch=16,
        box=FALSE,
        xlab="X",
        ylab="Y",
        zlab="Z",
        xlim=c(min(inv$X3), max(inv$X3)),
        ylim=c(min(inv$X2), max(inv$X2)),
        zlim=c(min(inv$X1), max(inv$X1))
    )
    
    invisible(
        sapply(2:length(unique(groups)), function(x)
            p$points3d(
                inv[inv$characteristics_ch1.1 == "MDS", "X3"],
                inv[inv$characteristics_ch1.1 == "MDS", "X2"],
                inv[inv$characteristics_ch1.1 == "MDS", "X1"],
                col=rgb(
                    group.color["red",   "MDS"],
                    group.color["green", "MDS"],
                    group.color["blue",  "MDS"],
                    150,
                    maxColorValue=255
                ),
                pch=16,
                cex = 1
            )
        )
    )
}