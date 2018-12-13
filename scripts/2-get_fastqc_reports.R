# CCDL ALSF 2018
# C. Savonen
#
# Purpose: Get a summary file of individual sequence quality reports
# from the FASTQC files.
#
# Options:
# "-d" - directory/path of where the fastqc reports have been placed.
# "-o" - directory/path where the output summary of fastqc will be placed
#
# Example bash usage:
#
# Rscript scripts/2-get_fastqc_reports.R \
# -d darmanis_data/fastqc_reports \
# -o results
#
# Get sequence quality reports from the FASTQC runs in a summary
library(optparse)
library(fastqcr)

# Get options using optparse
option_list <- list( 
  make_option(opt_str = c("-d", "--dir"), type = "character",
              default = NULL,
              help = "Directory of the fastqc reports",
              metavar = "character"),
  make_option(opt_str = c("-o", "--output"), type = "character", 
              default = getwd(),
              help = "Directory where results should be placed",
              metavar = "character"))

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Aggregate the reports
qc <- fastqcr::qc_aggregate(qc.dir = opt$dir)

# Write full report results to a csv file
write.csv(qc_stats(qc), file = file.path(opt$output,
                                        "fastqc_quality_report_full.csv"))

# Filter out samples that have failed the quality tests
qc_filtered <- data.frame(qc) %>%
  dplyr::select(sample, module, status) %>%    
  dplyr::filter(status %in% "PASS") %>%
  dplyr::arrange(sample)

# Write filtered results to a csv file
write.csv(qc_filtered, file = file.path(opt$output,
                                        "fastqc_quality_report_filtered.csv"))
