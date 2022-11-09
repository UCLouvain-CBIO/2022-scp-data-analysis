
# BiocManager::install(c("scpdata", "scp", "tidyverse", "patchwork",
#                        "SCP.replication", "reticulate", "zellkonverter",
#                        "scuttle"))
library("scpdata")
library("scp")
library("tidyverse")
library("patchwork")
library("SCP.replication")
library("biomaRt")
library("scater")
library("scran")
library("igraph")
library("viridis")
library("cluster")
library("RColorBrewer")

####---- Define computational workflows for data processing ----####

## We wrap the data processing workflows in meta-functions that take a
## QFeatures object as input as well as arguments that point to the different
## variable names in the QFeatures object

## The metafunctions were constructed based on the replication results available
## from the SCP.replication website: https://uclouvain-cbio.github.io/SCP.replication/

runSceptre <- function(x, colvarSampleType, scPattern, rowvarProtein,
                       colvarMSbatch, colvarChannel) {
    require("reticulate")
    require("zellkonverter")
    require("scuttle")
    spt <- import("sceptre")
    if (!identical("proteins", names(x)))
        stop("'x' should contain a single assay called 'proteins'")
    ## Keep only single cells
    x <- subsetByColData(x, grepl(scPattern, colData(x)[, colvarSampleType]))
    ## Batch correction 
    cat("Batch correction\n\n")
    sce <- getWithColData(x, "proteins")
    sce$`File ID` <- colData(sce)[, colvarMSbatch]
    sce$Channel <- colData(sce)[, colvarChannel]
    adata <- SCE2AnnData(sce, X_name = 1)
    spt$normalize(adata)
    sce <- AnnData2SCE(adata)
    colData(sce) <- NULL
    x <- addAssay(x, sce, name = "proteins_norm")
    x <- addAssayLink(x, from = "proteins", to = "proteins_norm",
                      varFrom = rowvarProtein, varTo = rowvarProtein)
    ## Sample QC
    cat(paste0("Sample QC\nNumber cells before QC: ",
               ncol(x[["proteins_norm"]]),"\n"))
    qc <- perCellQCMetrics(x[["proteins_norm"]], assay.type = 1)
    qc$outlier <- isOutlier(qc$sum, nmads = 2, log = TRUE, type = "both")
    mads <- attr(qc$outlier,"thresholds")
    colData(x)[, colnames(qc)] <- qc[rownames(colData(x)), ]
    x <- subsetByColData(x, !x$outlier & x$detected > 700)
    cat(paste0("Number cells after QC: ",
               ncol(x[["proteins_norm"]]),"\n\n"))
    ## Feature QC
    cat(paste0("Protein QC\nNumber proteins before QC: ",
               nrow(x[["proteins_norm"]]),"\n"))
    sce <- x[["proteins_norm"]]
    qc <- perFeatureQCMetrics(sce, assay.type = 1)
    qc$ncells <- qc$detected / 100 * ncol(sce)
    rowData(x[["proteins_norm"]])[, colnames(qc)] <- qc
    x <- filterFeatures(x, ~ ncells >= 3)
    cat(paste0("Number proteins after QC: ",
               nrow(x[["proteins_norm"]]),"\n\n"))
    ## Normalization
    cat("Normalization\n\n")
    sf <- colSums(assay(x, "proteins_norm"))
    sf <- sf / median(sf)
    x <- sweep(x, MARGIN = 2, STATS = sf, FUN = "/",
               i = "proteins_norm", name = "proteins_medshift")
    ## Log-transformation
    cat("Log-transformation\n\n")
    x <- logTransform(x, base = 2, pc = 1, i = "proteins_medshift", 
                      name = "proteins_log")
    ## Imputation
    cat("Imputation\n\n")
    sce <- getWithColData(x, "proteins_log")
    adata <- SCE2AnnData(sce, X_name = 1)
    spt$impute(adata)
    sce <- AnnData2SCE(adata)
    colData(sce) <- NULL
    x <- addAssay(x, sce, name = "proteins_imput")
    x <- addAssayLinkOneToOne(x, from = "proteins_log", to = "proteins_imput")
    ## Normalization
    cat("Normalization\n\n")
    means <- rowMeans(assay(x, "proteins_imput"))
    x <- sweep(x, MARGIN = 1, STATS = means, FUN = "-",
               i = "proteins_imput", name = "proteins_centered")
    stds <- rowSds(assay(x, "proteins_centered"))
    sweep(x, MARGIN = 1, STATS = stds, FUN = "/",
          i = "proteins_centered", name = "proteins_scaled")
}

runSCoPE2 <- function(x, colvarSampleType, colvarMSbatch, colvarScType, 
                      rowvarProtein, rowvarPeptide, rowvarPEP, rowvarPIF = NULL,
                      carrierPattern, scPattern, ncPattern, refPattern) {
    ## Feature QC
    cat(paste0("Feature QC\n", 
               "Number PSMs before QC: ", sum(dims(x)[1, ]),"\n"))
    x <- pep2qvalue(x, i = names(x), PEP = rowvarPEP, rowDataName = "qvalue_psm")
    x <- pep2qvalue(x, i = names(x), groupBy = rowvarProtein,
                    PEP = rowvarPEP, rowDataName = "qvalue_protein")
    x <- computeSCR(x, i = names(x), colvar = colvarSampleType,
                    carrierPattern = carrierPattern, 
                    samplePattern = paste0(scPattern, "|", ncPattern),
                    rowDataName = "MeanSCR")
    x <- filterFeatures(x,
                        ~ qvalue_psm < 0.01 & qvalue_protein < 0.01 &
                            !is.na(MeanSCR) & !is.infinite(MeanSCR) & MeanSCR < 0.1)
    if (!is.null(rowvarPIF)) {
        filt <- as.formula(paste0("~ !is.na(", rowvarPIF, ") & ", rowvarPIF, "> 0.8"))
        x <- filterFeatures(x, filt)
    }
    cat(paste0("Number PSMs after QC: ", sum(dims(x)[1, ]),"\n\n"))
    ## Batch correction by reference
    cat(paste0("Batch correction by reference\n\n"))
    x <- divideByReference(x, i = names(x), colvar = colvarSampleType,
                           samplePattern = ".", refPattern = refPattern)
    ## Aggregate PSMs
    cat(paste0("Aggregate PSMs\n\n"))
    remove.duplicates <- function(x)
        apply(x, 2, function(xx) xx[which(!is.na(xx))[1]] )
    peptideAssays <- paste0("peptides_", names(x))
    x <- aggregateFeaturesOverAssays(x, i = names(x), fcol = rowvarPeptide,
                                     name = peptideAssays, fun = remove.duplicates)
    ## Clean data 
    cat(paste0("Clean and join data\n\n"))
    x <- infIsNA(x, i = peptideAssays)
    x <- zeroIsNA(x, i = peptideAssays)
    x <- joinAssays(x, i = peptideAssays, name = "peptides")
    ## Sample QC
    cat(paste0("Sample QC\nNumber cells before QC: ",
               sum(grepl(scPattern, colData(x)[, colvarSampleType])),"\n"))
    x <- medianCVperCell(x, i = peptideAssays, groupBy = rowvarProtein,
                         nobs = 6, na.rm = TRUE, colDataName = "MedianCV",
                         norm = "SCoPE2")
    x <- subsetByColData(x, x$MedianCV < 0.365 & 
                             grepl(scPattern, colData(x)[, colvarSampleType]))
    cat(paste0("Number cells after QC: ",
               sum(grepl(scPattern, colData(x)[, colvarSampleType])),"\n\n"))
    ## Normalization
    cat("Normalization\n\n")
    x <- normalizeSCP(x, i = "peptides", method = "div.median",
                      name = "peptides_norm1")
    x <- sweep(x, i = "peptides_norm1", MARGIN = 1, FUN = "/",
               STATS = rowMeans(assay(x[["peptides_norm1"]]),
                                na.rm = TRUE),
               name = "peptides_norm2")
    ## Missing data filtering
    cat(paste0("Missing data filtering\n",
               "Number peptides before filter: ", nrow(x[["peptides_norm2"]]),"\n"))
    x <- filterNA(x, i = "peptides_norm2", pNA = 0.99)
    cat(paste0("Number peptides after filter: ", nrow(x[["peptides_norm2"]]),"\n\n"))
    ## Log-transformation
    cat("Log-transformation\n\n")
    x <- logTransform(x, base = 2, i = "peptides_norm2",
                      name = "peptides_log")
    ## Peptide aggregation
    cat("Peptide aggregation\n\n")
    x <- aggregateFeatures(x, i = "peptides_log", name = "proteins",
                           fcol = rowvarProtein,
                           fun = matrixStats::colMedians, na.rm = TRUE)
    ## Normalization
    cat("Normalization\n\n")
    x <- normalizeSCP(x, i = "proteins", method = "center.median",
                      name = "proteins_norm1")
    x <- sweep(x, i = "proteins_norm1", MARGIN = 1, FUN = "-",
               STATS = rowMeans(assay(x[["proteins_norm1"]]),
                                na.rm = TRUE),
               name = "proteins_norm2")
    ## Imputation
    cat("Imputation\n\n")
    # x <- impute(x, i = "proteins_norm2", method = "knn",
    #             name = "proteins_impd",
    #             k = 3, rowmax = 1, colmax= 1,
    #             maxp = Inf, rng.seed = 1234)
    x <- imputeKnnSCoPE2(x, i = "proteins_norm2", 
                         name = "proteins_impd", k = 3)
    ## Batch correction
    cat("Batch correction\n\n")
    sce <- getWithColData(x, "proteins_impd")
    batch <- as.character(colData(sce)[, colvarMSbatch])
    formul <- as.formula(paste0("~", colvarScType))
    ## Remove confounded batches
    tab <- table(batch = batch, interest = colData(x)[, colvarScType])
    sel <- rownames(tab)[rowSums(tab != 0) > 1]
    sce <- sce[, colData(sce)[, colvarMSbatch] %in% sel]
    model <- model.matrix(formul, data = colData(sce))
    batch <- as.character(colData(sce)[, colvarMSbatch])
    # assay(sce) <- sva::ComBat(dat = assay(sce), batch = batch, mod = model)
    assay(sce) <- ComBatv3.34(dat = assay(sce), batch = batch, mod = model)
    colData(sce) <- NULL
    x <- addAssay(x, y = sce, name = "proteins_batchC")
    addAssayLinkOneToOne(x, from = "proteins_impd", to = "proteins_batchC")
}

####---- Process specht data using SCoPE2 workflow ----####

## First, we must prepare the data for the workflow. Some steps are dataset
## specific and were not included in the metafunctions.

## Retrieve the specht2019v3 dataset
specht <- specht2019v3()
specht <- selectRowData(specht, c("dart_PEP", "protein", "peptide", "PIF",
                                  "dart_qval"))

## Remove the peptide and protein data that was provided by the authors
specht <- specht[, , -(178:179)]
## Solve the peptide to protein mapping ambiguity 
ppMap <- data.frame(rbindRowData(specht,
                                 i = names(specht))) %>%
    group_by(peptide) %>%
    ## The majority vote happens here
    mutate(protein = names(sort(table(protein),
                                decreasing = TRUE))[1]) %>%
    select(peptide, protein) %>%
    filter(!duplicated(peptide, protein))
consensus <- lapply(names(specht), function(i) {
    ind <- match(rowData(specht[[i]])$peptide, ppMap$peptide)
    DataFrame(protein = ppMap$protein[ind])
})
names(consensus) <- names(specht)
rowData(specht) <- consensus
## Remove failed runs
nPSMs <- dims(specht)[1, ]
specht <- specht[, , nPSMs > 500]
## Remove contaminants
specht <- filterFeatures(specht, ~ !grepl("REV|CON", protein))
## Run the data processing workflow
spechtProcessedScope2 <-
    runSCoPE2(specht,
              colvarSampleType = "SampleType",
              colvarMSbatch = "Set",
              rowvarPEP = "dart_PEP",
              rowvarProtein = "protein",
              rowvarPeptide = "peptide",
              rowvarPIF = "PIF",
              colvarScType = "SampleType",
              carrierPattern = "Carrier",
              scPattern = "Mono|Macro",
              ncPattern = "Blank",
              refPattern = "Reference")
save(spechtProcessedScope2, file = "data/spechtProcessedScope2.Rda")

####---- Process specht data using sceptre workflow ----####

specht2 <- specht
## Apply 1% FDR filter, the SCeptre workflow assumes the FDR filtering is already 
## performed
specht2 <- filterFeatures(specht2, ~ dart_qval < 0.01)
## We need to construct the protein data
proteinAssays <- paste0(names(specht2), "_proteins")
specht2 <- aggregateFeaturesOverAssays(specht2, i = names(specht2),
                                       fcol = "protein",
                                       name = proteinAssays,
                                       fun = matrixStats::colMedians,
                                       na.rm = TRUE)
specht2 <- joinAssays(specht2, i = proteinAssays, name = "proteins")
spechtproteins <- specht2[, , "proteins"]
## Run the data processing workflow
spechtProcessedSceptre <-
    runSceptre(spechtproteins,
               scPattern = "Mono|Macro",
               colvarSampleType = "SampleType",
               colvarMSbatch = "Set",
               colvarChannel = "Channel")
save(spechtProcessedSceptre, file = "data/spechtProcessedSceptre.Rda")

####---- Process schoof data using SCoPE2 workflow ----####

schoofall <- schoof2021()
schoofall <- selectRowData(schoofall, c("isContaminant", "Percolator.PEP", "Master.Protein.Accessions", "Annotated.Sequence",
                                  "Isolation.Interference.in.Percent"))

schoofpsms <- schoofall[, , 1:192]
schoofproteins <- schoofall[, , "proteins"]
## Add PIF information
rds <- rowData(schoofpsms)
pifs <- lapply(rds, function(rd) {
    rd$PIF <- 1 - rd$Isolation.Interference.in.Percent / 100
    rd[, "PIF", drop = FALSE]
})
rowData(schoofpsms) <- pifs
## Remove failed runs
nPSMs <- dims(schoofpsms)[1, ]
schoofpsms <- schoofpsms[, , nPSMs > 3000]
## Remove contaminants
schoofpsms <-
    filterFeatures(schoofpsms,
                   ~ !isContaminant & Master.Protein.Accessions != "")
## Rename protein names to keep protein groups
rds <- rowData(schoofpsms)
prots <- lapply(rds, function(rd) {
    rd$Master.Protein.Accessions <- 
        gsub(";", "_or_", rd$Master.Protein.Accessions)
    rd[, "Master.Protein.Accessions", drop = FALSE]
})
rowData(schoofpsms) <- prots
## Run the data processing workflow
schoofProcessedScope2 <-
    runSCoPE2(schoofpsms,
              colvarSampleType = "SampleType",
              colvarMSbatch = "File.ID",
              colvarScType = "Population",
              rowvarPEP = "Percolator.PEP",
              rowvarProtein = "Master.Protein.Accessions",
              rowvarPeptide = "Annotated.Sequence",
              rowvarPIF = "PIF",
              carrierPattern = "booster",
              scPattern = "^sc$",
              ncPattern = "neg",
              refPattern = "norm")
save(schoofProcessedScope2, file = "data/schoofProcessedScope2.Rda")

####---- Process schoof data using sceptre workflow ----####

## Remove failed runs
sel <- names(nPSMs)[nPSMs > 3000]
schoofproteins <- subsetByColData(schoofproteins, 
                                  schoofproteins$File.ID %in% sel) 
## Remove contaminants
schoofproteins <- filterFeatures(schoofproteins, ~ !isContaminant)
## Run the data processing workflow
schoofProcessedSceptre <- 
    runSceptre(schoofproteins, scPattern = "sc",
               colvarSampleType = "SampleType", 
               colvarMSbatch = "File.ID",
               colvarChannel = "Channel")
save(schoofProcessedSceptre, file = "data/schoofProcessedSceptre.Rda")

####---- Prepare results for comparison ----####

## First we extract the last assay and save it in a new file (avoids the need to
## recompute all the processing every time we want to make changes to the figure)

load("data/spechtProcessedSceptre.Rda")
protsSpechtSceptre <- getWithColData(spechtProcessedSceptre, "proteins_scaled")
protsSpechtSceptre$CellType <- protsSpechtSceptre$SampleType
protsSpechtSceptre$MSbatch <- protsSpechtSceptre$Set
metadata(protsSpechtSceptre) <- list(dataset = "specht2021", workflow = "SCeptre")
save(protsSpechtSceptre, file = "data/protsSpechtSceptre.Rda")

load("data/spechtProcessedScope2.Rda")
protsSpechtScope2 <- getWithColData(spechtProcessedScope2, "proteins_batchC")
protsSpechtScope2$CellType <- protsSpechtScope2$SampleType
protsSpechtScope2$MSbatch <- protsSpechtScope2$Set
metadata(protsSpechtScope2) <- list(dataset = "specht2021", workflow = "SCoPE2")
save(protsSpechtScope2, file = "data/protsSpechtScope2.Rda")

load("data/schoofProcessedSceptre.Rda")
protsSchoofSceptre <- getWithColData(schoofProcessedSceptre, "proteins_scaled")
protsSchoofSceptre$CellType <- protsSchoofSceptre$Population
protsSchoofSceptre$MSbatch <- protsSchoofSceptre$File.ID
metadata(protsSchoofSceptre) <- list(dataset = "schoof2021", workflow = "SCeptre")
save(protsSchoofSceptre, file = "data/protsSchoofSceptre.Rda")

load("data/schoofProcessedScope2.Rda")
protsSchoofScope2 <- getWithColData(schoofProcessedScope2, "proteins_batchC")
## The Schoof dataset requires to convert uniprot ID to protein symbol
mart <- useMart("ensembl")
human <- useDataset("hsapiens_gene_ensembl", mart)
upid <- unique(unlist(strsplit(rownames(protsSchoofScope2), "_or_ ")))
symbol <- getBM(attributes = c("uniprot_gn_symbol", "uniprot_gn_id"), 
                filters = "uniprot_gn_id",
                values = upid,
                mart = human)
rnames <- rownames(protsSchoofScope2)
for(i in seq_along(rnames)) {
    idsvec <- strsplit(rnames[[i]], "_or_ ")[[1]]
    for (j in seq_along(idsvec)) {
        s <- symbol$uniprot_gn_symbol[match(idsvec[[j]], symbol$uniprot_gn_id)]
        if (!is.na(s)) idsvec[[j]] <- s
    }
    rnames[[i]] <- paste0(idsvec, collapse = "; ")
}
rownames(protsSchoofScope2) <- rnames
protsSchoofScope2$CellType <- protsSchoofScope2$Population
protsSchoofScope2$MSbatch <- protsSchoofScope2$File.ID
metadata(protsSchoofScope2) <- list(dataset = "schoof2021", workflow = "SCoPE2")
save(protsSchoofScope2, file = "data/protsSchoofScope2.Rda")

## Load data (to avoid rerunning the above)
load("data/protsSpechtSceptre.Rda")
load("data/protsSpechtScope2.Rda")
load("data/protsSchoofSceptre.Rda")
load("data/protsSchoofScope2.Rda")

## Combine all datasets in a single object
prots <- list(SpechtSceptre = protsSpechtSceptre,
              SpechtScope2 = protsSpechtScope2,
              SchoofSceptre = protsSchoofSceptre,
              SchoofScope2 = protsSchoofScope2)

####---- Compare dimension reduciton ----####


## Compare PCA
prots <- lapply(prots, runPCA, exprs_values = 1)
pcaPlots <- lapply(prots, plotPCA, colour_by = "CellType")

## Compare t-SNE
## (Not shown in the article)
prots <- lapply(prots, function(x) {
    runTSNE(x, exprs_values = 1,
            ## parametrized using guidelines from Kobak et Berens 2019
            perplexity = ncol(x)/100, eta = ncol(x)/12)
})
tsnePlots <- lapply(prots, plotTSNE, colour_by = "CellType")

####---- Perform clustering ----####

prots <- lapply(prots, function(x){
    ## Cell types are known a priori
    x$SupervisedLabel <- as.integer(as.factor(x$CellType)) 
    ## Build neighbor graph
    metadata(x)$SNNgraph <- buildSNNGraph(scale(assay(x)),
                                          transposed = FALSE,
                                          d = Inf, ## Use all variables
                                          k = 15, ## number of nearest neighbors
                                          type = "jaccard") ## "rank", "number"
    ## Perform Louvain cluserting
    metadata(x)$LouvainCluster <- cluster_louvain(metadata(x)$SNNgraph)
    x$UnsupervisedLabel <- metadata(x)$LouvainCluster$membership
    x
})

####---- Compare unsupervised clustering ----####

## Create contingency tables to compare clusterings
## This function takes two SingleCellExperiment objects, x and y, and compares
## the cluster partitions stored under "UnsupervisedLabel" in the colData
compareClusters <- function(x, y, wfName1, wfName2) {
    ## Find the set of common cells across x and y
    sharedCells <- intersect(colnames(x), colnames(y))
    ## Get the contigency table
    clTable <- table(x = colData(x)[sharedCells, "UnsupervisedLabel"],
                     y = colData(y)[sharedCells, "UnsupervisedLabel"])
    clTable <- as.data.frame(clTable)
    ## Plot
    ggplot(clTable) +
        aes(x = x,
            y = y,
            fill = Freq) +
        geom_tile() +
        xlab(paste("Clusters from", wfName1)) +
        ylab(paste("Clusters from", wfName2)) +
        scale_fill_viridis(discrete = FALSE, name = "Number of cells") +
        theme_minimal()
}
## Create plots for both datasets
clustPlot <- list()
clustPlot[["schoof2021"]] <- compareClusters(prots$SchoofSceptre, prots$SchoofScope2,
                                             "SCeptre", "SCoPE2")
clustPlot[["specht2021"]] <- compareClusters(prots$SpechtSceptre, prots$SpechtScope2,
                                             "SCeptre", "SCoPE2")

####---- Compare supervised clustering ----####

## We consider the known cell types as clusters

## Compare silhouette widths
prots <- lapply(prots, function(x){
    ## Get the neighbourhood network from the computed clustering
    graphRes <- metadata(x)$SNNgraph
    ## Compute the silhouette width using the known cell types
    x$SilhouetteWidth <- silhouette(x$SupervisedLabel, dist = 1 - similarity(graphRes, method = "jaccard"))[, "sil_width"]
    x
})
## Store the silhouette widths for both workflows on both datasets in a single 
## table
swDf <- do.call(rbind, lapply(prots, function(x, sw) {
    data.frame(colData(x)[, c("SilhouetteWidth", "CellType", "MSbatch", "Channel")],
               dataset = metadata(x)$dataset,
               workflow = metadata(x)$workflow)
}))
## For each dataset, plot silhouette widths per workflow and per cell type
swPlots <- list()
for (i in c("schoof2021", "specht2021")) {
    swPlots[[i]] <- filter(swDf, dataset == i) %>% 
        ggplot() +
        aes(x = CellType,
            y = SilhouetteWidth,
            fill = CellType) +
        geom_violin() +
        geom_point(data = . %>% group_by(workflow, CellType) %>% 
                       summarise(SilhouetteWidth = median(SilhouetteWidth))) +
        xlab("") + ylab("Silhouette width") +
        facet_wrap(~ workflow) +
        theme_minimal()
}

## Compute correlation within cell types
prots <- lapply(prots, function(x) {
    celltypes <- unique(x$CellType)
    withinCor <- lapply(celltypes, function (ct) {
        ctind <- which(x$CellType == ct)
        cors <- cor(assay(x)[, ctind], use = "pairwise.complete.obs")
        cors[lower.tri(cors, diag = TRUE)] <- NA ## avoid duplicate points
        cors <- as.vector(cors)
        cors[!is.na(cors)]
    })
    names(withinCor) <- celltypes
    metadata(x)$WithinCellTypeCorrelations <- withinCor
    x
})
## Collect the correlations in a single table
corDf <- do.call(rbind, lapply(prots, function(x) {
    md <- metadata(x)
    cors <- md$WithinCellTypeCorrelations
    do.call(rbind, lapply(names(cors), function(ct) {
        data.frame(Correlation = cors[[ct]],
                   CellType =  ct,
                   workflow = md$workflow,
                   dataset = md$dataset)
    }))
}))

## For each dataset, plot the correlations per workflow and per cell type
corPlots <- list()
for (i in c("schoof2021", "specht2021")) {
    corPlots[[i]] <- filter(corDf, dataset == i) %>% 
        ggplot() +
        aes(x = CellType,
            y = Correlation,
            fill = CellType) +
        geom_violin() +
        geom_point(data = . %>% group_by(workflow, CellType) %>% 
                       summarise(Correlation = median(Correlation))) +
        xlab("") + ylab("Within cell type correlation") +
        facet_wrap(~ workflow) +
        theme_minimal()
}

####---- Compare differential abundance analysis ----####

daDf <- do.call(rbind, lapply(prots, function(x){
    wf <- metadata(x)$workflow
    ds <- metadata(x)$dataset
    cat("Differental abundance for the", ds, "data processed by", wf, "\n")
    ctypes <- unique(x$CellType)
    df <- do.call(rbind, lapply(seq_len(nrow(x)), function(i){
        xx <- assay(x)[i, ]
        group1 <- xx[x$CellType == ctypes[1]]
        group2 <- xx[x$CellType == ctypes[2]]
        res <- t.test(group1, group2)
        data.frame(pval = res$p.value,
                   logFC = unname(diff(res$estimate)),
                   workflow = wf,
                   dataset = ds)
    }))
    df$padj <- p.adjust(df$pval, "BH")
    df
}))

## For each dataset, plot the volcano plots per workflow
daPlots <- list()
for (i in c("schoof2021", "specht2021")) {
    daPlots[[i]] <- filter(daDf, dataset == i) %>% 
        ggplot() +
        aes(x = logFC,
            y = -log10(padj)) +
        geom_point() +
        xlab("Log fold change") + ylab("-log10 adjusted p-value") +
        geom_vline(xintercept = c(-1, 1)) + 
        geom_hline(yintercept = -log10(0.05)) + 
        facet_wrap(~ workflow) +
        theme_minimal()
}
i <- "specht2021"; daPlots[[i]] + ggtitle(paste("Dataset:", i))

####---- Create final plots ----####

## Colors
cols <- brewer.pal(8, name = "Set2")

## Figure 2
(plotSchoof2021 <- 
        (swPlots$schoof2021 +
             corPlots$schoof2021 &
             theme(axis.text.x = element_text(angle = 40, hjust = 1),
                   legend.position = "none",
                   plot.tag = element_text(vjust = -8),
                   strip.text = element_text(face = "bold", size = 12)) &
             scale_fill_manual(values = cols)) +
        (clustPlot$schoof2021 + 
             ggtitle("Unsupervised clustering") +
             theme(plot.title = element_text(face = "bold", size = 12, 
                                             hjust = 0.5, vjust = -7),
                   plot.tag = element_text(vjust = -8),
                   axis.title.x = element_text(hjust = 0, vjust = 12))) +
        (pcaPlots$SchoofSceptre + ggtitle("SCeptre") +
             pcaPlots$SchoofScope2 + ggtitle("SCoPE2")  +
             plot_layout(guides = "collect") &
             scale_color_manual(values = cols, name = "Cell type") &
             theme(plot.title = element_text(face = "bold", size = 12, 
                                             hjust = 0.5))) + 
        plot_layout(heights = c(0.4, 0.6), design = "
                112233
                444444") +
        plot_annotation(tag_levels = "A",
                        title = "Dataset: schoof2021",
                        theme = theme(
                            plot.title = element_text(size = 20, vjust = -5,
                                                      hjust = 0.5, face = "bold"))) &
        theme(plot.tag = element_text(face = "bold", size = 18)))
ggsave(filename = "figs/SCoPE2_vs_sceptre_schoof2021.pdf", plotSchoof2021,
       width = 9, height = 7)

## Figure 3
(plotSpecht2021 <- 
        (swPlots$specht2021 +
             corPlots$specht2021 &
             theme(axis.text.x = element_text(angle = 30, hjust = 1),
                   legend.position = "none", 
                   plot.tag = element_text(vjust = -8),
                   strip.text = element_text(face = "bold", size = 12)) &
             scale_fill_manual(values = cols)) + 
        (clustPlot$specht2021 + 
             ggtitle("Unsupervised clustering") +
             theme(plot.title = element_text(face = "bold", size = 12, 
                                             hjust = 0.5, vjust = -7),
                   plot.tag = element_text(vjust = -8),
                   axis.title.x = element_text(hjust = 0, vjust = 8))) +
        (pcaPlots$SpechtSceptre + ggtitle("SCeptre") +
             pcaPlots$SpechtScope2 + ggtitle("SCoPE2") +
             plot_layout(guides = "collect") &
             scale_color_manual(values = cols, name = "Cell type") &
             theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)))  +
        plot_layout(heights = c(0.4, 0.6), design = "
                112233
                444444") +
        plot_annotation(tag_levels = "A",
                        title = "Dataset: specht2021",
                        theme = theme(plot.title = element_text(size = 20, vjust = -5,
                                                                hjust = 0.5, face = "bold"))) &
        theme(plot.tag = element_text(face = "bold", size = 18)))
ggsave(filename = "figs/SCoPE2_vs_sceptre_specht2021.pdf", plotSpecht2021,
       width = 9, height = 7)
