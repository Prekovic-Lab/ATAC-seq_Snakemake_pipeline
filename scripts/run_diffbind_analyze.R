library(DiffBind)
library(dplyr)
library(rtracklayer)
library(BiocParallel)
register(SerialParam()) 

# 1. Load data and params
object.counted <- readRDS(snakemake@input[["rds"]])
col_name       <- snakemake@params[["column"]]
target_cond    <- snakemake@params[["target"]]
reference_cond <- snakemake@params[["reference"]]

# 2. Normalize and Setup Contrast dynamically
object.normalize <- dba.normalize(object.counted)
object.contrast  <- dba.contrast(object.normalize, 
                                 contrast=c(col_name, target_cond, reference_cond))

# 3. Analyze
object.analysed <- dba.analyze(object.contrast, method=DBA_DESEQ2, bParallel=FALSE)

# 4. Save PCA Plot (specifically for this subset of samples)
pdf(snakemake@output[["pca"]])
# We dynamically pass the column string to attributes so DiffBind colors by it
dba.plotPCA(object.analysed, attributes=col_name) 
dev.off()

# 5. Extract Results (FDR < 0.05 is default in DiffBind)
res <- dba.report(object.analysed, contrast=1)
res_df <- as.data.frame(res)

# Save CSV
write.csv(res_df, file=snakemake@output[["diff"]], row.names=FALSE)

# Save Volcano Plot
pdf(snakemake@output[["volc"]])
dba.plotVolcano(object.analysed, contrast=1)
dev.off()

# 6. Save BED files for downstream analysis (e.g., Homer, MEME)
# Fold > 0 means enriched in Target. Fold < 0 means enriched in Reference.
bed_up <- res_df %>% 
  filter(FDR < 0.05 & Fold > 0) %>% 
  select(seqnames, start, end)

bed_down <- res_df %>% 
  filter(FDR < 0.05 & Fold < 0) %>% 
  select(seqnames, start, end)

write.table(bed_up, file=snakemake@output[["bed_up"]], sep="\t", quote=F, row.names=F, col.names=F)
write.table(bed_down, file=snakemake@output[["bed_down"]], sep="\t", quote=F, row.names=F, col.names=F)

message("Analysis complete and all outputs saved.")