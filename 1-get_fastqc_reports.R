# CCDL ALSF 2018
# C. Savonen 
# 
# Get sequence quality reports from the FASTQC runs in a summary
library(optparse)
library(fastqcr)
option_list <- list( 
make_option(opt_str = c("-d", "--dir"), type = "character", default = NULL, 
            help = "Directory of the fastqc reports",
            metavar = "character"))

opt <- parse_args(OptionParser(option_list=option_list))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Aggregate the reports
qc <- fastqcr::qc_aggregate(qc.dir = opt$dir)
setwd(opt$dir)

# Write full report results to a csv file
write.csv(qc_stats(qc), file = "../fastqc_quality_report_full.csv")

# Filter out samples that have failed the quality tests
qc %>%
  dplyr::select(sample, module, status) %>%    
  dplyr::filter(status %in% c("WARN", "FAIL")) %>%
  dplyr::arrange(sample)

write.csv(qc_stats(qc), file = "../fastqc_quality_report_filtered.csv")
