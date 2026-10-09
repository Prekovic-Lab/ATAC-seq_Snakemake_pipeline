# Redirect all output and errors to the Snakemake log file
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

cat("Loading required libraries...\n")
suppressPackageStartupMessages({
    library(DiffBind)
    library(BiocParallel)
    library(dplyr)
    library(readr)
})

# Setup parallel processing based on Snakemake threads
threads <- snakemake@threads
register(MulticoreParam(workers = threads))
cat("Using", threads, "threads for processing.\n")

# 1. Read the sample metadata
cat("Reading sample sheet...\n")
sample_meta <- read_tsv(snakemake@params[["sample_sheet"]], show_col_types = FALSE)

# 2. Build the DiffBind Sample Sheet dataframe
# DiffBind requires specific column names: SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller
cat("Constructing DiffBind sample sheet...\n")

# Extract file paths from snakemake
bam_files <- snakemake@input[["bams"]]
peak_files <- snakemake@input[["peaks"]]

# Create a mapping data frame for the files
file_map <- data.frame(
    name = basename(bam_files) %>% gsub(".final.bam", "", .),
    bamReads = bam_files,
    Peaks = peak_files,
    stringsAsFactors = FALSE
)

# Merge metadata with file paths
dba_sheet <- sample_meta %>%
    inner_join(file_map, by = "name") %>%
    rename(
        SampleID = name,
        Condition = condition,
        Replicate = replicate
    ) %>%
    mutate(
        PeakCaller = "narrow" # Specify that we are using narrowPeak format
    ) %>%
    select(SampleID, Condition, Replicate, bamReads, Peaks, PeakCaller)

cat("Sample sheet created with", nrow(dba_sheet), "samples.\n")
print(dba_sheet)

# 3. Initialize DiffBind object
cat("\nInitializing DiffBind object...\n")
db_obj <- dba(sampleSheet = dba_sheet)

# 4. Count reads (This is the most time-consuming step)
cat("Counting reads in consensus peaks (this may take a while)...\n")
db_obj <- dba.count(db_obj, bParallel = TRUE)

# 5. Normalize
cat("Normalizing data...\n")
db_obj <- dba.normalize(db_obj)

# 6. Global PCA Plot
cat("Generating Global PCA plot...\n")
pdf(snakemake@output[["pca"]], width = 8, height = 6)
dba.plotPCA(db_obj, DBA_CONDITION, label=DBA_ID)
dev.off()

# 7. Setup Contrasts
cat("Setting up design and contrasts...\n")
# DiffBind 3+ modern way of setting up designs
db_obj <- dba.contrast(db_obj, design="~Condition")

# Run the analysis for all contrasts based on the design
db_obj <- dba.analyze(db_obj)

# Save the master DBA object so you don't have to rerun this 4-hour job to change a plot!
cat("Saving master DiffBind RDS object...\n")
saveRDS(db_obj, file = snakemake@output[["dba_obj"]])

# 8. Loop through requested contrasts from Snakemake config
contrasts_list <- snakemake@params[["contrasts"]]
outdir <- snakemake@params[["outdir"]]

for (contrast_name in names(contrasts_list)) {
    cat(sprintf("\nProcessing contrast: %s...\n", contrast_name))
    
    target <- contrasts_list[[contrast_name]][["target"]]
    reference <- contrasts_list[[contrast_name]][["reference"]]
    
    # In DiffBind >=3.0, you can extract specific contrasts using the contrast parameter in dba.report
    # format: c("Factor", "Numerator/Target", "Denominator/Reference")
    res <- dba.report(db_obj, contrast = c("Condition", target, reference))
    
    # Define output file names
    csv_file <- file.path(outdir, paste0(contrast_name, "_DiffPeaks.csv"))
    volcano_file <- file.path(outdir, paste0(contrast_name, "_Volcano.pdf"))
    
    # Save CSV
    if (!is.null(res) && length(res) > 0) {
        cat(sprintf("Found %d significant differential peaks.\n", sum(res$FDR < 0.05)))
        # Convert GRanges to dataframe and write
        res_df <- as.data.frame(res)
        write_csv(res_df, csv_file)
        
        # Save Volcano Plot
        pdf(volcano_file, width = 6, height = 6)
        dba.plotVolcano(db_obj, contrast = c("Condition", target, reference))
        title(main = contrast_name)
        dev.off()
    } else {
        cat("No significant peaks found for this contrast. Creating empty files.\n")
        write_csv(data.frame(Message="No significant peaks found"), csv_file)
        pdf(volcano_file)
        plot.new()
        text(0.5, 0.5, "No differential peaks to plot")
        dev.off()
    }
}

cat("\nDiffBind analysis completed successfully!\n")

# Close sink
sink(type="message")
sink()