library(DESeq2)

###### read in data and put in DESeq dataset #####
htseq_dir <- "/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/output/HTSeq/GENCODE/release_25/"
sampleFiles <- list.files(htseq_dir)
# change below
sampleCondition <- ifelse(grepl('t0',sampleFiles),'untreated','treated')
sampleName <- gsub(".txt","",sampleFiles)
# change below
timepoint <- str_extract(sampleName,"t[0-9]+")
# change below
batch <- str_extract(sampleName, "batch[0-9]")
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = sampleCondition, stringsAsFactors=FALSE,
                          timepoint = timepoint,
                          batch=batch)
sampleTable$condition <- as.factor(sampleTable$condition)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = htseq_dir,
                                  design= ~ condition)

sample_names <- gsub('(.*_t[0-9].*?)_.*','\\1',colnames(dds))
dds_short_names <- dds
colnames(dds_short_names) <- sample_names

included <- c("batch2_TCC-03-1", "batch2_TCC-04-1", "batch3_TCC-01-1", "batch3_TCC-01-2", "batch3_TCC-01-3",
              "batch3_TCC-04-2", "batch3_TCC-05-1", "batch3_TCC-06-1", "batch3_TCC-07-1", "batch3_TCC-10-1",
              "batch3_TCC-11-1", "batch3_TCC-12-1", "batch3_TCC-13-1", "batch4_TCC-02-1", "batch4_TCC-03-2",
              "batch4_TCC-14-1", "batch4_TCC-15-1", "batch4_TCC-16-1", "batch4_TCC-17-1", "batch4_TCC-18-2",
              "batch4_TCC-19-1", "batch5_TCC-04-3")
dds_filtered <- dds_short_names[,gsub("(.*)_t.*",'\\1', colnames(dds_short_names)) %in% included]
vsd_filtered_no_estimations <- varianceStabilizingTransformation(dds_filtered)
dds_filtered <- estimateSizeFactors(dds_filtered)
dds_filtered <- estimateDispersions(dds_filtered)
vsd_filtered <- varianceStabilizingTransformation(dds_filtered)