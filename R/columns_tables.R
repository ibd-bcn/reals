library("tidyverse")
library("readxl")
library("openxlsx")

bd <- read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = c("n.a.", ""))
upa <- read_xls("data/M13-740_abbvie_database_300519.xls", na = c("n.a.", ""))

dff <- read_xlsx("processed/AU_markers.xlsx")
upa2 <- mutate(upa, PCR_UPA_paper = if_else(BarCode %in% dff$Sample, "yes", "no"),
               SubjectID = as.character(SubjectID),
               cd = tolower(gsub("Container", "", ContainerName))) %>%
  left_join(dff[, c("Sample", "Patient", "General_location", "remission")],
            by = c("BarCode" = "Sample", "SubjectID" = "Patient", "cd" = "General_location")) %>%
  select(-cd) %>%
  distinct()
bd2 <- mutate(bd, PCR_UPA_paper = if_else(
  gsub("^0| reseq", "", tolower(Sample_id)) %in% dff$Sample,
                                          "yes", "no"),
  s = gsub("^0| reseq", "", tolower(Sample_id))) %>%
  left_join(dff[, c("Sample", "Patient", "remission", "General_location")],
            by = c("Pacient_id" = "Patient",
                   "sample_location" = "General_location",
                   "s" = "Sample")) %>%
  select(-s) %>%
  distinct()

write.xlsx(upa2, "processed/M13-740_abbvie_database_210619.xlsx", keepNA = FALSE,
           firstRow = TRUE, colNames = TRUE)
write.xlsx(bd2, "processed/bd_BCN_tnf_biopsies_210619.xlsx", keepNA = FALSE,
           firstRow = TRUE, colNames = TRUE)

bd3 <- select(bd2, -c(grep("\\.", colnames(bd2), value = TRUE), "fastq_file")) %>%
  filter(PCR_UPA_paper == "yes") %>%
  mutate(k = duplicated(gsub("^0| reseq", "", tolower(Sample_id))),
         k2 = if_else(seq_len(n()) %in% (which(k)-1), FALSE, TRUE)
         ) %>%
  filter(k2) %>%
  select(-PCR_UPA_paper, -k, -k2)
k <- sapply(bd3, function(x){all(is.na(x))})
write.xlsx(bd3[, !k], "processed/bd_BCN_PCR_UPA_paper.xlsx", keepNA = FALSE,
           firstRow = TRUE, colNames = TRUE, withFilter = TRUE)
