library("scpdata")
library("scp")
library("tidyverse")
library("scater")
library("sva")
library("patchwork")
library("scuttle")
library("SCP.replication")

####---- williams2020_tmt ----####

## Replicate the TMT data processing by Williams et al. 2020

## Load the data
williams <- williams2020_tmt()
## Keep data of insterest
williams <- subsetByColData(williams, williams$Amount == "1cell")
## Feature quality control: remove contaminants
williams <- filterFeatures(williams, 
                           ~ Reverse != "+" &
                               Potential.contaminant != "+")
## Missing data cleaning
williams <- zeroIsNA(williams, i = "proteins_corrected")
## Log-transformation
williams <- logTransform(williams, 
                         i = "proteins_corrected",
                         name = "proteins_log", 
                         base = 2)
## Normalization by median centering
williams <- normalizeSCP(williams,
                         i = "proteins_log",
                         name = "proteins_norm",
                         method = "center.median")
## Batch correction by median centering
sceNorm <- getWithColData(williams, "proteins_norm")
for (batch in sceNorm$Batch) {
    ind <- which(sceNorm$Batch == batch) 
    rowMeds <- rowMedians(assay(sceNorm[, ind]), na.rm = TRUE)
    assay(sceNorm[, ind]) <- sweep(assay(sceNorm[, ind]), FUN = "-", 
                                   STATS = rowMeds, MARGIN = 1)
}
williams <- addAssay(williams, sceNorm, name = "proteins_bc")
williams <- addAssayLinkOneToOne(williams, from = "proteins_norm", 
                                 to = "proteins_bc")
## Feature quality control: remove highly missing proteins
williams <- filterNA(williams, i = "proteins_bc", pNA = 0.3)
## Feature quality control: keep proteins with at least 2 unique peptides
prots <- rowData(williams[["peptides_corrected"]])$Leading.razor.protein
prots <- unique(prots[duplicated(prots)])
prots <- prots[prots %in% rownames(williams)[["proteins_bc"]]]
williams <- williams[prots, , ]

####---- brunner2022 ----####

## Replicate the data processing by Brunner et al. 2022

## Load and format the data (not yet in scpdata)
df <- read.delim("~/PhD/.localdata/SCP/brunner2022/20210919_DIANN_SingleCellOutput.pg_matrix.tsv")
# df <- read.delim("~/PhD/.localdata/SCP/brunner2022/20210919_DIANN_SingleCellOutput.tsv")
coln <- colnames(df)[grepl("Andrea", colnames(df))]
colnames(df) <- sub("^D.*_(\\d*).d", "\\1", colnames(df))
prots <- readSingleCellExperiment(df, ecol = grep("^\\d", colnames(df)), fnames = "Protein.Names")
prots <- prots[rownames(prots) != ""]
rownames(prots) <- make.unique(rownames(prots))
cd <- DataFrame(MsBatch = sub("^D.*RawData.(.*)$", "\\1", coln))
annots <- strsplit(cd$MsBatch, "_")
annots <- lapply(annots, function(x) {
    if (length(x) == 12) x <- x[-4]
    x
})
annots <- do.call(rbind, annots)
annots <- DataFrame(annots)
colnames(annots) <- c("Date", "MsInstrument", "Purification", "User",
                      "SampleAnnotation", "SampleType", "CellCycleStage", ".1", 
                      "PlatePosition", "CellNumber", "RunID")
annots$RunID <- sub(".d", "", annots$RunID)
rownames(annots) <- annots$RunID
brunner2022 <- QFeatures(ExperimentList(proteins = prots), 
                         colData = annots)
## Cell quality control
brunner2022 <- impute(brunner2022, i = "proteins", method = "zero", 
                      name = "proteins_zero")
cellqc <- perCellQCMetrics(brunner2022[["proteins_zero"]], assay.type = 1)
colData(brunner2022)[, colnames(cellqc)] <- cellqc[rownames(colData(brunner2022)), ]
brunner2022 <- subsetByColData(brunner2022, brunner2022$detected >= 600)
## Feature quality control
protqc <- perFeatureQCMetrics(brunner2022[["proteins_zero"]], assay.type = 1)
rowData(brunner2022) <- List(proteins = protqc, proteins_zero = protqc)
brunner2022 <- filterFeatures(brunner2022, ~ detected > 0.15)
## Log-transformation
brunner2022 <- logTransform(brunner2022, i = "proteins", 
                            name = "proteins_log", 
                            base = exp(1), pc = 1)
## Imputation by downshift

imputeMatrixWithdownShiftedNormal <- function(x, scale = 0.3, shift = 1.8) {
    m <- mean(x, na.rm = TRUE)
    std <- sd(x, na.rm = TRUE)
    nmis <- sum(is.na(x))
    repl <- rnorm(nmis, mean = m - shift * std, sd = scale * std)
    x[is.na(x)] <- repl
    x
}
brunner2022 <- impute(brunner2022, 
                      i = "proteins_log", 
                      name = "proteins_dsnImpd",
                      FUN = imputeMatrixWithdownShiftedNormal)
brunner2022$CellCycleStage <- recode(brunner2022$CellCycleStage,
                                     TB = "G1-S",
                                     NB = "G2-M",
                                     UB = "other")

####---- Plot MS drift ----####

specht <- specht2019v3()
peps <- getWithColData(specht, "peptides")
peps <- peps[, grepl("Mono|Macro", peps$SampleType) &
                 peps$lcbatch == "LCB7"]
feat <- "_KLLEGEESR_3" 
## Peptide with MS drift without cell type effect: "_HVLVTLGEK_3"
ordering <- gsub("^(\\d{6}).*16plex(.*)_Set_(\\d*).*$", "\\1-\\2-\\3", peps$Set)
ordering <- gsub("--", "-1-", ordering)
ordering <- sapply(ordering, function(x){
    x <- strsplit(x, "-")[[1]]   
    if (nchar(x[[3]]) == 1) x[[3]] <- paste0("0", x[[3]])
    paste(x[c(1, 3, 4, 2)], collapse = "")
})
(msdrift <- data.frame(value = assay(peps)[feat, ],
                       order = as.numeric(as.factor(ordering)),
                       colData(peps)) %>% 
        ggplot() +
        aes(x = order, 
            y = value,
            colour = SampleType) +
        geom_point() +
        geom_smooth(colour = "grey40", method = loess, 
                    formula = y ~ x, se = FALSE) +
        theme_minimal() +
        ggtitle("MS drift") +
        xlab("Batch index") +
        ylab("Log2 peptide intensity"))

####---- Confounder effect: cell type correlated with tag ----####

sce <- getWithColData(williams, "proteins_bc")
sce <- runPCA(sce, exprs_values = 1)
levels(sce$Channel) <- paste0("TMT", c("131C", "130N", "128C",
                                       "128N", "131N", "129C",
                                       "129N", "130C"))
pvar <- attr((reducedDim(sce, "PCA")), "percentVar")
ccols <- c("green4", "purple3", "blue3", "yellow3", "orange3",
           "darkseagreen3", "olivedrab2", "dodgerblue")
(confounderTag <- plotPCA(sce, colour_by = "Channel", shape_by = "SampleType") +
        ggtitle("Effect of TMT label") +
        scale_color_manual(values = ccols,
                           name = "TMT channel") +
        theme_minimal() +
        xlab(paste0("PC1 (", round(pvar[[1]]), "%)")) +
        ylab(paste0("PC2 (", round(pvar[[2]]), "%)")) +
        guides(shape = guide_legend(title = "Cell line")))

####---- Confounder effect: cell type correlated with date ----####

sce <- getWithColData(brunner2022, "proteins_dsnImpd")
sce <- sce[, sce$CellCycleStage != "other"]
sce <- runPCA(sce, ncomponents = 50,
              ntop = Inf,
              scale = TRUE,
              exprs_values = 1,
              name = "PCA")
sce$Date <- format(as.Date(sce$Date, format("%Y%m%d")), "%d %b %Y")
pvar <- attr((reducedDim(sce, "PCA")), "percentVar")
(confounderDate <- plotPCA(sce, shape_by = "CellCycleStage", 
                           colour_by = "Date") +
        xlab(paste0("PC1 (", round(pvar[[1]]), "%)")) + 
        ylab(paste0("PC2 (", round(pvar[[2]]), "%)")) + 
        ggtitle("Effect of date") +
        guides(shape = guide_legend(title = "Cell cycle")) +
        theme_minimal())

####---- Residual batch effect ----####

load("scripts/data/protsSpechtSceptre.Rda") ## generated in make_figure2and3.R
set.seed(1234)
protsSpechtSceptre <- runPCA(protsSpechtSceptre, exprs_values = 1)
l <- length(unique(protsSpechtSceptre$Set))
labels <- breaks <- unique(protsSpechtSceptre$Set)[seq(1, l, length.out = 11)]
labels <- gsub("_.*_", "_", labels)
labels[seq(2, length(labels), 2)] <- "..."
pvar <- attr((reducedDim(protsSpechtSceptre, "PCA")), "percentVar")
(residBC <-
        plotPCA(protsSpechtSceptre, colour_by = "Set") +
        ggtitle("Residual batch effect") +
        xlab(paste0("PC1 (", round(pvar[[1]]), "%)")) +
        ylab(paste0("PC2 (", round(pvar[[2]]), "%)")) +
        scale_color_discrete(name = "MS acquisition",
                             breaks = breaks, labels = labels) +
        scale_y_reverse() +
        theme_minimal())

####---- Make figure 4 ----####

(pl <- msdrift +
     confounderDate +
     confounderTag +
     residBC +
     plot_annotation(tag_levels = "A") &
     theme(plot.tag = element_text(size = 20, face = "bold"),
           legend.position = "right"))
ggsave(pl, filename = "figs/Batch_effects.pdf",
       width = 9, height = 7)
