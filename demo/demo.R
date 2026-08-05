library(SingleR)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

#### Set-up ####

# demo bulk TPM matrix
bulk_mat <- readRDS("demo/bulk_demo.rds")

# demo single cell Seurat object
sc_obj <- readRDS("demo/sc_demo.rds")

# functions
source("ssp_bkr.R")
source("singler_subtype.R")
source("penalized_module_score.R")

# make output file folder
dir.create("demo/demo_outputs")


#### Single sample classifier for bulk RNA-seq ####

classifier_genes <- predict.class.bkr(mode="mat")
writeLines(classifier_genes$gg, con="demo/demo_outputs/demo_bulk_classifier_genes.txt")

pred <- predict.class.bkr(bulk_mat)
write.csv(cbind(sample=colnames(bulk_mat), pred), file="demo/demo_outputs/demo_bulk_classifier_output.csv")


#### Correlation-based classifier for scRNA-seq ####

mat <- sc_obj@assays$RNA@data
out <- singler_subtype(mat, reference_file="references/singler_reference.RData")

write.csv(out$data, file="demo/demo_outputs/demo_sc_classifier_data_output.csv")
print(head(out$singler_pred))

sc_obj <- AddMetaData(sc_obj, out$data)


#### Penalized gene set scoring for scRNA-seq ####

# Use the built-in cell cycle gene set list to test 
print(cc.genes)

sc_obj <- add_penalized_module_score(sc_obj, genelists=cc.genes)
names(sc_obj@meta.data)[match(names(cc.genes), names(sc_obj@meta.data))] <- paste0(names(cc.genes), "_penalized")

write.csv(sc_obj[[paste0(names(cc.genes), "_penalized")]], file="demo/demo_outputs/demo_sc_penalized_score_output.csv")

sc_obj <- AddModuleScore(sc_obj, features=cc.genes, name="cc.genes")
names(sc_obj@meta.data)[match(paste0("cc.genes", 1:2), names(sc_obj@meta.data))] <- paste0(names(cc.genes), "_original")


#### sc Visualization ####

pdf("demo/demo_outputs/demo_sc_plots.pdf")
DimPlot(sc_obj, group.by="subtype", cols=c(Classical1="#1C6CAB", Basal1="#814C42", Classical2="#A4C0E5", Basal2="#FF7311", Mixed="grey40")) + 
  ggtitle("singler_subtype") + 
  theme(aspect.ratio=1, axis.ticks=element_blank(), axis.text=element_blank())
FeaturePlot(sc_obj, features=c(paste0(names(cc.genes), "_penalized"), paste0(names(cc.genes), "_original"))) *
  theme(aspect.ratio=1, axis.ticks=element_blank(), axis.text=element_blank())
dev.off()

