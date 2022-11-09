library(tidyverse)
library(patchwork)
library(RColorBrewer)

datasets <- c("zhu2018", "budnik2018", "dou2019", "zhu2019", 
              "cong2020", "tsai2020", "williams2020_lfq",
              "williams2020_tmt", "liang2020", 
              "specht2021", "cong2021", "schoof2021", "woo2021",
              "brunner2022", "leduc2022", "woo2022", "webber2022",
              "derks2022")
datasets <- factor(datasets, levels = datasets)

steps <- c("Feature quality control", "Sample quality control",
           "Log-transform", "Normalization", "Batch correction",
           "Imputation", "Aggregation", "Identification update", "Other")
steps <- factor(steps, levels = steps)

levels <- c("preprocessing", "PSM/precursor", "peptide", "protein", NA)
levels <- factor(levels, levels = levels)

## Custom colors
colors <- c("purple", "dodgerblue2",  "gold2", "#447821ff",
            "#88bd66ff", "red3", "darkorange", "#d7aff0", "grey40")
## Color-blind friendly
colors <- c(RColorBrewer::brewer.pal(8, "Dark2"), "#6A3D9A")
colors <- colors[c(2, 4, 5, 6, 7, 3, 1, 9, 8)]

names(colors) <- steps
coln <- c("dataset", "level", "order", "step", "method")

## Step specific methods
methods <- list(c("Contaminants and decoys", "Empty signal", "FDR filtering",
                  "Missing data", "Unique peptides per protein",
                  "Spectral purity", "Sample to carrier ratio"),
                c("Poorly documented", "Failed runs", "Missing data",
                  "Median CV", "MAD sum intensity", "Grubb's test"),
                c("log2", "log1p"),
                c("Poorly documented", "Sample median", "Feature mean", 
                  "Feature median", "Total", "Quantile", "Standardization"),
                c("Reference", "Median", "ComBat", "limma"),
                c("Poorly documented", "KNN", "zero", "Normal distribution"),
                c("Poorly documented", "Random selection", "Summed", "Top N",
                  "Median", "iBAQ", "maxLFQ"),
                c("Match between runs", "Percolator rescoring", "DART-ID rescoring"),
                c("Isotope impurity correction", "Quantification thresholding",
                  "S/N computation"))
names(methods) <- steps

####---- zhu2018 ----####

zhu2018 <- data.frame(rbind(
    c("zhu2018", NA, NA, "Identification update", "Match between runs"),
    c("zhu2018", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("zhu2018", "preprocessing", 2, "Aggregation", "iBAQ"),
    c("zhu2018", "protein", 1, "Feature quality control", "Contaminants and decoys")
))
colnames(zhu2018) <- coln

####---- budnik2018 ----####

budnik2018 <- data.frame(rbind(
    c("budnik2018", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("budnik2018", "peptide", 1, "Feature quality control", "Contaminants and decoys"),
    c("budnik2018", "peptide", 2, "Normalization", "Sample median"),
    c("budnik2018", "peptide", 3, "Normalization", "Feature mean"),
    c("budnik2018", "peptide", 4, "Aggregation", "Median"),
    c("budnik2018", "protein", 1, "Imputation", "KNN"),
    c("budnik2018", "protein", 2, "Batch correction", "Reference")
))
colnames(budnik2018) <- coln


####---- dou2019 ----####

dou2019 <- data.frame(rbind(
    c("dou2019", "PSM/precursor", 1, "Feature quality control", "Contaminants and decoys"),
    c("dou2019", "PSM/precursor", 2, "Feature quality control", "FDR filtering"),
    c("dou2019", "PSM/precursor", 3, "Feature quality control", "Unique peptides per protein"),
    c("dou2019", "PSM/precursor", 4, "Other", "Isotope impurity correction"),
    c("dou2019", "PSM/precursor", 5, "Aggregation", "Poorly documented"),
    c("dou2019", "protein", 1, "Log-transform", "log2"),
    c("dou2019", "protein", 2, "Sample quality control", "Poorly documented"),
    c("dou2019", "protein", 3, "Normalization", "Sample median"),
    c("dou2019", "protein", 4, "Feature quality control", "Missing data"),
    c("dou2019", "protein", 5, "Batch correction", "ComBat"),
    c("dou2019", "protein", 6, "Imputation", "Poorly documented")
))
colnames(dou2019) <- coln

####---- zhu2019 ----####

zhu2019 <- data.frame(rbind(
    c("zhu2019", NA, NA, "Identification update", "Match between runs"),
    c("zhu2019", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("zhu2019", "preprocessing", 2, "Aggregation", "iBAQ"),
    c("zhu2019", "protein", 1, "Feature quality control", "Contaminants and decoys"),
    c("zhu2019", "protein", 2, "Feature quality control", "Missing data"),
    c("zhu2019", "protein", 3, "Normalization", "Poorly documented"),
    c("zhu2019", "protein", 4, "Log-transform", "log2"),
    c("zhu2019", "protein", 5, "Imputation", "zero"),
    c("zhu2019", "protein", 6, "Batch correction", "ComBat")
))
colnames(zhu2019) <- coln

####---- cong2020 ----####

cong2020 <- data.frame(rbind(
    c("cong2020", NA, NA, "Identification update", "Match between runs"),
    c("cong2020", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("cong2020", "preprocessing", 2, "Aggregation", "maxLFQ"),
    c("cong2020", "protein", 1, "Feature quality control", "Contaminants and decoys")
))
colnames(cong2020) <- coln

####---- tsai2020 ----####

tsai2020 <- data.frame(rbind(
    c("tsai2020", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("tsai2020", "preprocessing", 2, "Aggregation", "Summed"),
    c("tsai2020", "protein", 1, "Log-transform", "log2"),
    c("tsai2020", "protein", 2, "Batch correction", "Reference"),
    c("tsai2020", "protein", 3, "Normalization", "Sample median"),
    c("tsai2020", "protein", 4, "Normalization", "Quantile")
))
colnames(tsai2020) <- coln

####---- williams2020_lfq ----####

williams2020_lfq <- data.frame(rbind(
    c("williams2020_lfq", NA, NA, "Identification update", "Match between runs"),
    c("williams2020_lfq", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("williams2020_lfq", "preprocessing", 2, "Aggregation", "iBAQ"),
    c("williams2020_lfq", "protein", 1, "Feature quality control", "Contaminants and decoys"),
    c("williams2020_lfq", "protein", 2, "Feature quality control", "Missing data"),
    c("williams2020_lfq", "protein", 3, "Log-transform", "log2")
))
colnames(williams2020_lfq) <- coln

####---- williams2020_tmt ----####

williams2020_tmt <- data.frame(rbind(
    c("williams2020_tmt", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("williams2020_tmt", "preprocessing", 2, "Aggregation", "iBAQ"),
    c("williams2020_tmt", "preprocessing", 3, "Batch correction", "Reference"),
    c("williams2020_tmt", "protein", 1, "Feature quality control", "Contaminants and decoys"),
    c("williams2020_tmt", "protein", 2, "Log-transform", "log2"),
    c("williams2020_tmt", "protein", 3, "Normalization", "Sample median"),
    c("williams2020_tmt", "protein", 4, "Batch correction", "Median"),
    c("williams2020_tmt", "protein", 5, "Feature quality control", "Missing data"),
    c("williams2020_tmt", "protein", 6, "Feature quality control", "Unique peptides per protein")
))
colnames(williams2020_tmt) <- coln

####---- specht2021 ----####

specht2021 <- data.frame(rbind(
    c("specht2021", NA, NA, "Identification update", "DART-ID rescoring"),
    c("specht2021", "preprocessing", 1, "Other", "Isotope impurity correction"),
    c("specht2021", "PSM/precursor", 1, "Sample quality control", "Failed runs"),
    c("specht2021", "PSM/precursor", 2, "Feature quality control", "Contaminants and decoys"),
    c("specht2021", "PSM/precursor", 3, "Feature quality control", "Spectral purity"),
    c("specht2021", "PSM/precursor", 4, "Feature quality control", "FDR filtering"),
    c("specht2021", "PSM/precursor", 5, "Feature quality control", "Sample to carrier ratio"),
    c("specht2021", "PSM/precursor", 6, "Batch correction", "Reference"),
    c("specht2021", "PSM/precursor", 7, "Aggregation", "Random selection"),
    c("specht2021", "peptide", 1, "Sample quality control", "Median CV"),
    c("specht2021", "peptide", 2, "Normalization", "Sample median"),
    c("specht2021", "peptide", 3, "Normalization", "Feature mean"),
    c("specht2021", "peptide", 4, "Feature quality control", "Missing data"),
    c("specht2021", "peptide", 5, "Log-transform", "log2"),
    c("specht2021", "peptide", 6, "Aggregation", "Median"),
    c("specht2021", "protein", 1, "Normalization", "Sample median"),
    c("specht2021", "protein", 2, "Normalization", "Feature mean"),
    c("specht2021", "protein", 3, "Imputation", "KNN"),
    c("specht2021", "protein", 4, "Batch correction", "ComBat"),
    c("specht2021", "protein", 5, "Normalization", "Sample median"),
    c("specht2021", "protein", 6, "Normalization", "Feature mean")
))
colnames(specht2021) <- coln

####---- liang2020 ----####

## cf https://github.com/PayneLab/Lymphocytes1/blob/master/B_cells_versus_T_cells_FP.ipynb

liang2020 <- data.frame(rbind(
    c("liang2020", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("liang2020", "preprocessing", 2, "Aggregation", "Summed"),
    c("liang2020", "protein", 1, "Feature quality control", "Contaminants and decoys"),
    c("liang2020", "protein", 2, "Log-transform", "log2"),
    c("liang2020", "protein", 3, "Normalization", "Sample median"),
    c("liang2020", "protein", 4, "Feature quality control", "Missing data") ## min 3 cells per group
))
colnames(liang2020) <- coln

####---- schoof2021 ----####

schoof2021 <- data.frame(rbind(
    c("schoof2021", NA, NA, "Identification update", "Percolator rescoring"),
    c("schoof2021", "preprocessing", 1, "Other", "Isotope impurity correction"),
    c("schoof2021", "preprocessing", 2, "Other", "S/N computation"),
    c("schoof2021", "preprocessing", 3, "Feature quality control", "FDR filtering"),
    c("schoof2021", "preprocessing", 4, "Aggregation", "Poorly documented"),
    c("schoof2021", "protein", 1, "Feature quality control", "Contaminants and decoys"),
    c("schoof2021", "protein", 2, "Sample quality control", "Failed runs"),
    c("schoof2021", "protein", 3, "Batch correction", "Median"),
    c("schoof2021", "protein", 4, "Other", "Quantification thresholding"),
    c("schoof2021", "protein", 5, "Sample quality control", "Missing data"),
    c("schoof2021", "protein", 6, "Sample quality control", "MAD sum intensity"),
    c("schoof2021", "protein", 7, "Feature quality control", "Missing data"),
    c("schoof2021", "protein", 8, "Normalization", "Sample median"),
    c("schoof2021", "protein", 9, "Log-transform", "log2"),
    c("schoof2021", "protein", 10, "Imputation", "KNN"),
    c("schoof2021", "protein", 11, "Normalization", "Standardization")
))
colnames(schoof2021) <- coln


####---- cong2021 ----####

cong2021 <- data.frame(rbind(
    c("cong2021", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("cong2021", "preprocessing", 2, "Normalization", "Poorly documented"),
    c("cong2021", "preprocessing", 3, "Aggregation", "Poorly documented"),
    c("cong2021", "protein", 1, "Log-transform", "log2"),
    c("cong2021", "protein", 2, "Feature quality control", "Missing data"),
    c("cong2021", "protein", 3, "Imputation", "Normal distribution")
))
colnames(cong2021) <- coln

####---- woo2021 ----####

woo2021 <- data.frame(rbind(
    c("woo2021", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("woo2021", "preprocessing", 2, "Other", "Isotope impurity correction"),
    c("woo2021", "preprocessing", 3, "Aggregation", "Summed"),
    c("woo2021", "protein", 1, "Feature quality control", "Contaminants and decoys"),
    c("woo2021", "protein", 2, "Log-transform", "log2"),
    c("woo2021", "protein", 3, "Feature quality control", "Missing data"),
    c("woo2021", "protein", 4, "Imputation", "Normal distribution"),
    c("woo2021", "protein", 5, "Normalization", "Quantile"),
    c("woo2021", "protein", 6, "Batch correction", "ComBat")
))
colnames(woo2021) <- coln

####---- brunner2022 ----####

brunner2022 <- data.frame(rbind(
    c("brunner2022", NA, NA, "Identification update", "Match between runs"),
    c("brunner2022", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("brunner2022", "protein", 1, "Sample quality control", "Missing data"),
    c("brunner2022", "protein", 2, "Feature quality control", "Missing data"),
    c("brunner2022", "protein", 3, "Log-transform", "log1p"),
    c("brunner2022", "protein", 4, "Imputation", "Normal distribution")
))
colnames(brunner2022) <- coln

####---- leduc2022 ----####

leduc2022 <- data.frame(rbind(
    c("leduc2022", NA, NA, "Identification update", "DART-ID rescoring"),
    c("leduc2022", "preprocessing", 1, "Other", "Isotope impurity correction"),
    c("leduc2022", "PSM/precursor", 1, "Feature quality control", "Contaminants and decoys"),
    c("leduc2022", "PSM/precursor", 2, "Feature quality control", "Spectral purity"),
    c("leduc2022", "PSM/precursor", 3, "Feature quality control", "FDR filtering"),
    c("leduc2022", "PSM/precursor", 4, "Feature quality control", "Sample to carrier ratio"),
    c("leduc2022", "PSM/precursor", 5, "Feature quality control", "Empty signal"),
    c("leduc2022", "PSM/precursor", 6, "Batch correction", "Reference"),
    c("leduc2022", "PSM/precursor", 7, "Aggregation", "Random selection"),
    c("leduc2022", "peptide", 1, "Sample quality control", "Median CV"),
    c("leduc2022", "peptide", 2, "Normalization", "Sample median"),
    c("leduc2022", "peptide", 3, "Normalization", "Feature median"),
    c("leduc2022", "peptide", 4, "Feature quality control", "Missing data"),
    c("leduc2022", "peptide", 5, "Log-transform", "log2"),
    c("leduc2022", "peptide", 6, "Aggregation", "Median"),
    c("leduc2022", "protein", 1, "Normalization", "Sample median"),
    c("leduc2022", "protein", 2, "Normalization", "Feature median"),
    c("leduc2022", "protein", 3, "Imputation", "KNN"),
    c("leduc2022", "protein", 4, "Batch correction", "limma"),
    c("leduc2022", "protein", 5, "Normalization", "Sample median"),
    c("leduc2022", "protein", 6, "Normalization", "Feature median")
))
colnames(leduc2022) <- coln

####---- woo2022 ----####

woo2022 <- data.frame(rbind(
    c("woo2022", "preprocessing", 1, "Feature quality control", "FDR filtering"),
    c("woo2022", "preprocessing", 2, "Aggregation", "iBAQ"),
    c("woo2022", "protein", 1, "Log-transform", "log2"),
    c("woo2022", "protein", 2, "Feature quality control", "Missing data"),
    c("woo2022", "protein", 3, "Normalization", "Sample median"),
    c("woo2022", "protein", 4, "Batch correction", "ComBat"),
    c("woo2022", "protein", 5, "Imputation", "Normal distribution")
))
colnames(woo2022) <- coln

####---- webber2022 ----####

webber2022 <- data.frame(rbind(
    c("webber2022", "preprocessing", 1, "Aggregation", "Top N"),
    c("webber2022", "protein", 1, "Sample quality control", "Grubb's test"),
    c("webber2022", "protein", 2, "Normalization", "Total"),
    c("webber2022", "protein", 3, "Imputation", "Normal distribution")
))
colnames(webber2022) <- coln

####---- derks2022 ----####

derks2022 <- data.frame(rbind(
    c("derks2022", "PSM/precursor", 1, "Other", "Isotope impurity correction"),
    c("derks2022", "PSM/precursor", 2, "Feature quality control", "FDR filtering"),
    c("derks2022", "PSM/precursor", 3, "Feature quality control", "Missing data"),
    c("derks2022", "PSM/precursor", 4, "Batch correction", "Median"),
    c("derks2022", "PSM/precursor", 5, "Normalization", "Sample median"),
    c("derks2022", "PSM/precursor", 6, "Normalization", "Feature mean"),
    c("derks2022", "PSM/precursor", 7, "Aggregation", "Median"),
    c("derks2022", "protein", 1, "Batch correction", "Median"),
    c("derks2022", "protein", 2, "Normalization", "Sample median"),
    c("derks2022", "protein", 3, "Normalization", "Feature mean"),
    c("derks2022", "protein", 4, "Log-transform", "log2"),
    c("derks2022", "protein", 5, "Feature quality control", "Missing data"),
    c("derks2022", "protein", 6, "Sample quality control", "Missing data"),
    c("derks2022", "protein", 7, "Imputation", "KNN"),
    c("derks2022", "protein", 8, "Batch correction", "ComBat"),
    c("derks2022", "protein", 9, "Normalization", "Sample median"),
    c("derks2022", "protein", 10, "Normalization", "Feature mean")
))
colnames(derks2022) <- coln

####---- Combine all tables ----####

wfTable <- rbind(zhu2018, budnik2018, dou2019, zhu2019, cong2020, 
                 williams2020_lfq, williams2020_tmt, specht2021, 
                 liang2020, schoof2021, tsai2020, cong2021, woo2021, 
                 brunner2022, leduc2022, woo2022, webber2022, derks2022)
## Check validity of table
stopifnot(all(wfTable$dataset %in% datasets))
stopifnot(all(wfTable$level %in% levels))
stopifnot(all(wfTable$step %in% steps))
wrongEntries <- lapply(steps, function(step){
    sel <- wfTable$step == step & !wfTable$method %in% methods[[step]]
    wfTable[sel, ]  
})
wrongEntries <- do.call(rbind, wrongEntries)
stopifnot(nrow(wrongEntries) == 0)

wfTable$order <- as.numeric(wfTable$order)
## Store the method's shape
wfTable$shape <- apply(wfTable, 1, function(line) {
    which(methods[[line[["step"]]]] == line[["method"]])
})
wfTable$shape <- as.factor(wfTable$shape)
## Store the step's color
wfTable$color <- sapply(wfTable$step, function(s) colors[[s]])

####---- Utility functions  ----####

plotStep <- function(step) {
    table <- wfTable
    availMethods <- methods[[step]]
    ds <- datasets
    ## Subset table
    table <- table[table$step == step, ]
    ## Control factor level ordering
    table$method <- factor(table$method, levels = availMethods)
    table$dataset <- factor(table$dataset, levels = ds)
    ## Plot
    ggplot(table) +
        aes(x = method,
            y = dataset,
            shape = shape) +
        geom_point(size = 3, color = table$color, fill = table$color) +
        ggtitle(step) +
        scale_y_discrete(drop = FALSE)
}

plotWorklfow <- function() {
    table <- wfTable[wfTable$step != "Identification update", ]
    ds <- datasets
    ## Control factor level ordering
    table$dataset <- factor(table$dataset, levels = ds)
    table$level <- factor(table$level, levels = levels)
    table$order <- factor(table$order, levels = 1:max(table$order))
    ## Plot
    ggplot(table) +
        aes(x = order,
            y = dataset,
            shape = shape) +
        facet_grid(~ level, scales = "free", space = "free") +
        geom_point(size = 3, color = table$color, fill = table$color) +
        ggtitle("Workflows") +
        scale_y_discrete(drop = FALSE)
}


####---- Identification update ----####

(pl_id <- plotStep("Identification update") +
     xlab("") + ylab("") +
     theme_minimal() +
     theme(legend.position = "none",
           plot.title = element_text(hjust = 0.5),
           axis.text.x = element_text(angle = 35, hjust = 1)))
# ggsave(pl_id, file = "figs/Compare_identification_update.svg", 
#        height = 3, width = 3)


####---- Quantitative data processing workflows ----####

(pl_quant <- plotStep("Sample quality control") +
     plotStep("Feature quality control") +
     plotStep("Imputation") +
     plotStep("Log-transform") +
     plotStep("Aggregation") +
     plotStep("Normalization") +
     plotStep("Batch correction") +
     plotStep("Other") +
     plotWorklfow() +
     plot_layout(design = "
    00001111122233
    44445555566677
    88888888888888
                ") &
     scale_shape_manual(values = c(21:25, 3:4)) &
     xlab("") & ylab("") &
     theme_minimal() &
     theme(legend.position = "none",
           plot.title = element_text(hjust = 0.5),
           axis.text.x = element_text(angle = 35, hjust = 1)))

addAnnotations <- function() {
    library(grid)
    annotPar <- gpar(fontface = "bold", fontsize = 18)
    grid.text("A", x = 0.05, y = 0.98, gp = annotPar)
    grid.text("B", x = 0.35, y = 0.98, gp = annotPar)
    grid.text("C", x = 0.62, y = 0.98, gp = annotPar)
    grid.text("D", x = 0.85, y = 0.98, gp = annotPar)
    grid.text("E", x = 0.05, y = 0.635, gp = annotPar)
    grid.text("F", x = 0.35, y = 0.635, gp = annotPar)
    grid.text("G", x = 0.62, y = 0.635, gp = annotPar)
    grid.text("H", x = 0.85, y = 0.635, gp = annotPar)
    grid.text("I", x = 0.05, y = 0.29, gp = annotPar)
}

## Save plot
pdf("figs/Workflows_overview.pdf",
    height = 13, width = 10)
pl_quant
addAnnotations()
dev.off()
