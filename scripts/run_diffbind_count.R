library(DiffBind)
library(BiocParallel)
register(SerialParam()) 

# 1. Fetch dynamic params
sample_sheet_path <- snakemake@input[["samples"]]
outdir            <- snakemake@params[["outdir"]]   # NEW: Get the absolute path
col_name          <- snakemake@params[["column"]]
target_cond       <- snakemake@params[["target"]]
reference_cond    <- snakemake@params[["reference"]]

# 2. Read the raw TSV 
samples_raw <- read.delim(sample_sheet_path, sep="\t", stringsAsFactors=FALSE)

if (!(col_name %in% colnames(samples_raw))) {
  stop(sprintf("Error: Column '%s' not found in %s!", col_name, sample_sheet_path))
}

# 3. Subset to only the target and reference conditions
subset_samples <- samples_raw[samples_raw[[col_name]] %in% c(target_cond, reference_cond), ]

message(sprintf("Running counting based on %s: %s vs %s", col_name, target_cond, reference_cond))
message("Samples included: ", paste(subset_samples$name, collapse=", "))

# 4. TRANSLATE TO DIFFBIND FORMAT (Using Absolute Paths)
dba_sheet <- data.frame(
  SampleID   = subset_samples$name,
  Condition  = subset_samples$condition,
  Replicate  = subset_samples$replicate,
  PeakCaller = subset_samples$peak_type,
  
  # Point DiffBind to the absolute path of the files
  bamReads   = file.path(outdir, "4.filtered", paste0(subset_samples$name, ".final.bam")),
  Peaks      = file.path(outdir, "6.peak_calls", paste0(subset_samples$name, "_peaks.narrowPeak")),
  
  stringsAsFactors = FALSE
)

# 5. Create DBA object
object <- dba(sampleSheet=dba_sheet)

# 6. Count reads 
object.counted <- dba.count(object, bUseSummarizeOverlaps=TRUE, bParallel=FALSE)

# 7. Save RDS
saveRDS(object.counted, file=snakemake@output[["rds"]])
message("Counting complete.")