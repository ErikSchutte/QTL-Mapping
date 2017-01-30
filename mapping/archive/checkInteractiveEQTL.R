covs <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/covariates_leiden_samples.txt", header=T,row.names = 1)
gt1 <- read.table("~/Dropbox/Erik Schutte Internship 2016/Data/filtered-data/snps.t.tsv",header=T, sep=" ")
gt1 <- cbind.data.frame(gt1, gt1, gt1, gt1)

load("~/Dropbox/Erik Schutte Internship 2016/Data/eQTL-data/gsTcell_vst_88samples.Rdata")

load("~/Dropbox/Erik Schutte Internship 2016/Results/Cis/GENCODE/gsTcell_all_timepoints_interactive-cis-eQTLs_0.05.RData")

load("~/Dropbox/Erik Schutte Internship 2016/Results/Trans/GENCODE/gsTcell_all_timepoints_interactive-trans-eQTLs_0.05.RData")
vst <- vst[,-grep("TCC-17", colnames(vst))]
colnames(vst) <- gsub("batch[0-9]_T", "T", colnames(vst))
colnames(vst) <- gsub("_t[0-9]+", "", colnames(vst))

condition = unlist(covs[3,])
age = unlist(covs[1,])
gender= unlist(covs[2,])


gt <- unlist(gt1["rs61907765",])
expr <- unlist(vst["ENSG00000200763.1",])

m1 <- lm( expr ~ age + gender + condition + gt + condition * gt )
m0 <- lm( expr ~ age + gender + condition + gt )
p.I <- anova(m1,m0)[2,6]
