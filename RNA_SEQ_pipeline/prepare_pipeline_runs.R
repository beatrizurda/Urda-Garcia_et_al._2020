##########################################################################
# Prepare runs to run the pipeline for all selected diseases in the Mare
# 
#
#
#
##########################################################################

######## AT THE DISEASE LEVEL #######-----------------------------------
setwd("Analysis/")

# Create a file with the diseases to analyze: cdisease_to_analyze.txt
diseases <- list.files(path="GREIN_SE_good_quality/")
diseases <- gsub("_grein_se.rds","",diseases)
write.table(diseases,file="RNA_SEQ_pipeline/cdisease_to_analyze.txt",sep="\t",row.names=F)

# Import list of diseases to be run
diseases <- read.csv2('RNA_SEQ_pipeline/cdisease_to_analyze.txt',header=TRUE, sep="\t",stringsAsFactors = F)
diseases <- diseases$x

# Generate the running lines and save them in a file
runs <- paste("Rscript ","run_rnaseq_pipeline_for_disease.R ",diseases,"_grein_se.rds",sep="") 
write.table(runs,file="RNA_SEQ_pipeline/new_run_rnaseq_pipeline_dl.txt",sep="\n",row.names=F, quote = F, col.names = FALSE)

######## AT THE ICD9 CODE LEVEL #######-----------------------------------
setwd("~/Desktop")
icds <- list.files("SE_icd9_good_quality/")

# Generate the running lines and save them in a file
runs <- paste("Rscript ","run_rnaseq_pipeline_for_disease.R ",icds," icd9",sep="") 
write.table(runs,file="new_run_rnaseq_pipeline_icd9_level.txt",sep="\n",row.names=F, quote = F, col.names = FALSE)
