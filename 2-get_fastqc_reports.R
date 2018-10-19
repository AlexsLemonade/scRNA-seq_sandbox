# CCDL ALSF 2018
# C. Savonen 
# 
# Get sequence quality reports from the FASTQC run
if(!("fastqcr" %in% installed.packages())){
  install.packages('fastqcr')
}
# Magrittr pipe
`%>%` <- dplyr::`%>%`

# Aggregate the reports
qc <- fastqcr::qc_aggregate(qc.dir = "fastqc_reports")

write.csv(qc_stats(qc), file = "fastqc_quality_report_full.csv")

qc %>%
  dplyr::select(sample, module, status) %>%    
  dplyr::filter(status %in% c("WARN", "FAIL")) %>%
  dplyr::arrange(sample)

write.csv(qc_stats(qc), file = "fastqc_quality_report_filtered.csv")
