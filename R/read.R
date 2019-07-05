library("readxl")
library("dplyr")
library("stringr")
library("ggplot2")
library("tidyr")
library("writexl")
library("ggpubr")
library("patchwork")
library("lubridate")
library("broom")
library("purrr")
library("ggforce")
library("ggrepel")

bd <- read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = c("n.a.", ""))
bd_upa <- read_xls("data/M13-740_abbvie_database_300519.xls", na = c("n.a.", "NA"))
reals <- read_xls("data/SOX6 CHI3L1 PDGFD PTGDR2_20190516_115511_Results_Export.xls",
                  skip = 15, na = c("", "Undetermined"))
reals2 <- read_xls("data/THY HTR3E RETNLB COL3A1_20190524_165141_Results_Export.xls",
                   skip = 15, na = c("", "Undetermined"))
reals3 <- read_xls("data/TBX21 IL17A IFN GZMH_20190531_125712_Results_Export.xls",
                   skip = 15, na = c("", "Undetermined"))
reals4 <- read_xls("data/AQP7 DERL OSM S100A8_20190618_105927_Results_Export.xls",
                   skip = 15, na = c("", "Undetermined"))
reals5 <- read_xls("data/AQP8_20190704_165805_Results_Export.xls",
                   skip = 0, na = c("", "Undetermined"))
upa_codes <- read_xlsx("data/RT UPA DISEASE BLOC 1 2 3 500NG 20UL.xlsx")
control_codes <- read_xlsx("data/RT CONTROLS 500NG 20 UL.xlsx")
antiTNF_codes <- read_xlsx("data/RT biopsies BCN bloc A B C D 500 ng en 20 ul.xlsx")
w0_remitters_upa <- read_xlsx("processed/missing_response based on ulcers UPA.xlsx") %>%
  mutate(remitter = tolower(remitter)) %>%
  rename(Biopsy_Location = "Biopsy_Location week 0") %>%
  select(Patient, remitter, Biopsy_Location)

# Clean the original data ####
reals <- merge(reals, reals2, all = TRUE)
reals <- merge(reals, reals3, all = TRUE)
reals <- merge(reals, reals4, all = TRUE)
reals <- merge(reals, reals5, all = TRUE)

# Correct ids or remove duplicated samples ####
# Do you remember some samples that were in more than one sample?
# Now we remove a plate
reals[reals$`Sample Name` %in% c("102w14c", "73w14c")  &
        grepl("290419bcn 102 a 108", reals$`Experiment Name`), 1:5]
reals <- filter(reals,
                !(`Sample Name` %in% c("102w14c", "73w14c")  & grepl("290419bcn 102 a 108", `Experiment Name`)))

reals[reals$`Experiment Name` == "290419 bcn samples plat 3.eds" &
        reals$`Sample Name` == "44w0c", "Sample Name"] <- "41w14"

reals$`Sample Name`[grep("UP11 ", reals$`Sample Name`)] <- "UP11 AB4216302"
reals$`Sample Name`[grep("UP52 ", reals$`Sample Name`)] <- "UP52 Y46412103"
reals$`Sample Name`[grepl("^48C 42W", reals$`Sample Name`, ignore.case = TRUE)] <- "48C 52w14c"

# Y45570803 is duplicated! but it should according to what we know of the samples...

# Look for the beta actina of the same well as the genes.
nested_well <- reals %>%
  filter(`Target Name` != "BETA ACTINA") %>%
  group_by(`Experiment Name`, `Target Name`, `Sample Name`) %>%
  nest(Well)

nested_well$data <- lapply(wells, unlist)

# Look for the wells how many have values
xy <- vector("list", length = nrow(nested_well))
for (i in seq_len(nrow(nested_well))) {
  x <- nested_well[i, ]
  wells_value <- reals %>%
    filter(`Sample Name` == x[["Sample Name"]],
           `Experiment Name` == x[["Experiment Name"]],
           Well %in% x[["data"]][[1]]
    )
  xy[[i]] <- wells_value %>%
    group_by(`Target Name`) %>%
    summarise(nas = sum(is.na(`Cт`)))

}

# Look how many missing values of Beta actina per target sample, experiment
ba <- sapply(xy, function(x) {
  x[x$`Target Name` == "BETA ACTINA", "nas", drop = TRUE]})
table(ba)
# Three samples that don't have right the BETA ACTINA
# For genes AQP7, CHI3L1, IFN GAMMA

# Look how many missing values of gene per target sample, experiment
gene <- sapply(xy, function(x){
  x[x$`Target Name` != "BETA ACTINA", "nas", drop = TRUE]})
table(gene)
# Three samples that don't have right the BETA ACTINA
# For genes AQP7, CHI3L1, IFN GAMMA

# Relations between genes and beta actina
table(gene, ba)

# Missing samples: Probably bad pipetting on the wet lab
nested_well[gene == 3 & ba == 3, 1:3]


reals %>%
  filter(`Sample Name` %in% c("41C 41W14I", "C15 ILI", "UP3 AA5292203"),
         `Experiment Name`  %in% c("040619 PLAT 4 BCN.eds",
                                   "070519  bcn controls ili PLAT 1 sox6 chi pdgfd ptg2.eds",
                                   "280519 UPA PLAT 1.eds"),
         `Target Name` %in% c("AQP7", "CHI3L1", "IFN GAMMA", "BETA ACTINA")) %>%
  arrange(`Experiment Name`, Well, `Sample Name`, `Target Name`) %>%
  select(Well, `Sample Name`, `Target Name`, `Cт`)

# Add values to those with beta actine but without gene values
bigger_ct <- 40
# Three missing values on the gene but values for beta actina
cond <- gene == 3 & ba != 3
df <- as.data.frame(nested_well[cond, 1:3])
for (row in seq_len(nrow(df))) {
  l <- reals$"Experiment Name" == df[row, "Experiment Name"] &
    reals$"Target Name" == df[row, "Target Name"] &
    reals$"Sample Name" == df[row, "Sample Name"]
  reals[l, "Cт"] <- bigger_ct
}

# Verify that the samples with missing values have high values of the ones
# they have
cond <- gene == 2
df <- as.data.frame(nested_well[cond, 1:3])
o <- sapply(seq_len(nrow(df)), function(row) {
  l <- reals$"Experiment Name" == df[row, "Experiment Name"] &
    reals$"Target Name" == df[row, "Target Name"] &
    reals$"Sample Name" == df[row, "Sample Name"]
  ct <- reals[l, "Cт"]
  c(mean = mean(ct, na.rm = TRUE),
    max = max(ct, na.rm = TRUE))
})
o <- t(o)
colnames(o) <- c("mean", "max")
hist(o[, 2], xlab = "Max Ct",
     main = c("Distribution of max Ct of the genes with missing values",
     "\nBy genes having only 1 value."))


# Check if some samples are in more than one plate ####
multiple_exp <- reals %>%
  group_by(`Sample Name`) %>%
  summarise(n = n_distinct(`Experiment Name`))
multiple_exp %>%
  group_by(n) %>%
  count()


well_failed <- reals %>%
  group_by(`Experiment Name`, `Well`) %>%
  summarise(failed = sum(is.na(`Cт`)), targets = list(`Target Name`)) %>%
  ungroup()

n <- well_failed %>%
  group_by(failed) %>%
  count() %>%
  filter(failed == 3)
stopifnot(nrow(n) == 0)

wrong_genes <- well_failed %>%
  filter(failed == 2) %>%
  inner_join(reals) %>%
  select(`Sample Name`, `Target Name`) %>%
  filter(`Target Name` != "BETA ACTINA") %>%
  distinct(.keep_all = TRUE)

failed_reference <- well_failed %>%
  filter(failed == 2) %>%
  inner_join(reals) %>%
  select(`Sample Name`, `Target Name`) %>%
  group_by(`Sample Name`) %>%
  summarise(reference = all(`Target Name` == "BETA ACTINA")) %>%
  filter(reference)

stopifnot(nrow(failed_reference) == 0)

#  Superseeded by gene and ba check above.
# failed_betas <- well_failed %>%
#   filter(failed == 1) %>%
#   inner_join(reals) %>%
#   filter(`Target Name` != "BETA ACTINA",
#          !is.na(`Cт`)) %>%
#   select(`Sample Name`, `Target Name`) %>%
#   distinct()
# stopifnot(nrow(failed_betas) == 0)

# Analyse if the targets are right or wrong ####
failed <- reals %>%
  group_by(`Experiment Name`, `Target Name`, `Sample Name`) %>%
  summarise(n = sum(is.na(`Cт`))) %>%
  ungroup() %>%
  select(-`Experiment Name`) %>%
  filter(`Target Name` != "BETA ACTINA",
         n > 1)

sd_Ct <- group_by(reals, `Target Name`, `Sample Name`) %>%
  summarise(sd = sd(`Cт`, na.rm = TRUE), mean = mean(`Cт`, na.rm = TRUE),
            sem = sd/sqrt(n()))


erroneous <- filter(sd_Ct,
                    `sd` > 0.5 | mean > 25 & !is.na(sd),
                    `Target Name` == "BETA ACTINA")

if (nrow(erroneous) > 1) {
  write_xlsx(erroneous, "processed/big_sd.xlsx")

  pdf("Figures/histograms.pdf")
  hist(reals$Cт[reals$`Target Name` == "BETA ACTINA"], main = "Beta Actina", xlab = "Ct")
  hist(reals$Cт[reals$`Target Name` == "CHI3L1"], main = "CHI3L1", xlab = "Ct")
  hist(reals$Cт[reals$`Target Name` == "PDGFD"], main = "PDGFD", xlab = "Ct")
  hist(reals$Cт[reals$`Target Name` == "SOX6"], main = "SOX6", xlab = "Ct")
  dev.off()
}

write.csv(reals, "processed/reals.csv", row.names = FALSE)

# Check the results ####
multiple_exp <- reals %>%
  group_by(`Sample Name`) %>%
  summarise(n = n_distinct(`Experiment Name`))
m <-multiple_exp %>%
  group_by(n) %>%
  count()
# Normal same samples have been tested for several genes in different plates

preclean <- reals %>%
  filter(`Target Name` != "BETA ACTINA") %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(`ΔCт` = unique(`ΔCт Mean`)) %>%
  ungroup() %>%
  filter(!is.na(`ΔCт`))

antiTNF <- preclean %>%
  filter(grepl("[0-9]w", `Sample Name`, ignore.case = TRUE)) %>%
  separate(`Sample Name`, c("Id", "Sample Name"), fill = "left") %>%
  mutate(`Sample Name` = tolower(`Sample Name`),
         Pacient_id = as.character(as.numeric(gsub("w.*", "", `Sample Name`))),
         Time = str_extract(`Sample Name`, "w[0-9]*"),
         Sample_id = paste0(Pacient_id, "-", Time),
         AU = 2^(-`ΔCт`)*1000) %>%
  select(-`ΔCт`)

# To verify
incorrect <- antiTNF %>%
  group_by(Id, `Sample Name`) %>%
  count() %>%
  ungroup() %>%
  filter(!is.na(Id)) %>%  # Exclude some missing with NAs
  group_by(`Sample Name`) %>%
  count() %>%
  filter(n != 1)
stopifnot(nrow(incorrect) == 0)

UPA <- preclean %>%
  filter(grepl("^UP", `Sample Name`)) %>%
  mutate(`Sample Name` = gsub("^UP ", "UP", `Sample Name`),
         AU = 2^(-`ΔCт`)*1000) %>%
  separate(`Sample Name`, c("Id", "Pacient_id")) %>%
  select(-`ΔCт`)

# To verify
incorrect <- UPA %>%
  group_by(Id, Pacient_id) %>%
  count() %>%
  ungroup() %>%
  group_by(Id) %>%
  count() %>%
  filter(n != 1) %>%
  pull(Id)
# Should be 0 or null
if (length(incorrect) != 0){
  filter(upa_codes, `RT TUBE` %in% incorrect) %>%
    mutate(sample = paste(`RT TUBE`, `nº muestra`)) %>%
    pull(sample)
}
stopifnot(length(incorrect) == 0)


controls <- preclean %>%
  mutate(`Sample Name` = gsub(" CONT", "CONT", `Sample Name`)) %>%
  filter(grepl("(^C)|([0-9]+CONT C)", `Sample Name`)) %>%
  mutate(AU = 2^(-`ΔCт`)*1000,
         Location = case_when(grepl("[0-9]+ C|S", `Sample Name`) ~ "colon",
                              grepl(" ILI$", `Sample Name`) ~ "ileum")) %>%
  separate(`Sample Name`, c("Id", "Patient_id", "Loc"), remove = FALSE,
           fill = "left") %>% # We fill from the left because some don't have two IDs
  select(-`ΔCт`)

# To verify
controls %>%
  group_by(Id, `Sample Name`) %>%
  count() %>%
  ungroup() %>%
  filter(!is.na(Id)) %>%
  group_by(Id) %>%
  count() %>%
  filter(n != 1)
# 0

# Clean the "databases" ####
bd2 <- bd %>%
  mutate(Sample_id = tolower(Sample_id)) %>%
  filter(str_detect(Sample_id, "-w"), IBD == "CD") %>%
  select(Sample_id, CD_endoscopic_remission, CD_endoscopic_response,
         sample_location, Ulcers, biopsied_segment) %>%
  mutate(Sample_id = gsub(" reseq|rep", "", Sample_id),
         Pacient_id = as.character(as.numeric(gsub("-w.*", "", Sample_id))),
         Time = str_extract(Sample_id, "w[0-9]*"),
         Sample_id = paste0(Pacient_id, "-", Time)) %>%
  filter(Time %in% c("w0", "w14")) %>%
  distinct() %>%
  rename(Location = sample_location)

upa_codes <- filter(upa_codes, !is.na(`PLACA MICRONIC`))

# Merge them together ####

remission <- function(x) {
  if (any(x$Time == "w14")) {
    x <- unique(x)
    if (x$Ulcers[x$Time == "w14"] == "yes") {
      return("no")
    } else {
      return("yes")
    }
  } else {
    return("missing")
  }
}

## antiTNF ####

df <- merge(bd2, antiTNF, all.x = FALSE, all.y = TRUE,
            by = c("Sample_id", "Pacient_id", "Time")) %>%
  filter(!(Time == "w0" & Ulcers == "no"))

response <- df %>%
  group_by(Pacient_id, biopsied_segment) %>%
  nest(Ulcers, Time, .key = AnyUlcers) %>%
  mutate(remission = map(AnyUlcers, remission)) %>%
  unnest(remission, .drop = TRUE)

response$remission[response$remission == "missing"] <- "no"
df <- inner_join(df, response)

## UPA ####
duplic <- group_by(bd_upa, BarCode) %>%
  summarise(BarCodes = n()) %>%
  filter(BarCodes != 1) %>%
  pull(BarCode)
# Because we have duplicates we remove some data (first filtering for what do we care
bd_upa <- bd_upa %>%
  select(Week, pSES.CD, BarCode, ContainerName, SubjectID, ulcers, Biopsy_Location) %>%
  distinct()

# Prepare the data for filtering and remissions
df_upa <- merge(UPA, upa_codes,  by.x = "Id", by.y = "RT TUBE",
                all.x = TRUE, all.y = FALSE) %>%
  left_join(bd_upa, by = c("nº muestra" = "BarCode")) %>%
  mutate(ContainerName = gsub("Biopsy ", "", ContainerName),
         Location = ContainerName,
         Week = paste0("w", Week),
         pSES.CD = as.numeric(pSES.CD),
         remission_old = case_when(
           Week == "w0" ~ "w0",
           Week != "w0" & pSES.CD < 5 ~ "w14 remiters",
           Week != "w0" & pSES.CD >= 5 ~ "w14 non-remiters",
           TRUE ~ Week),
         SubjectID = as.character(SubjectID)) %>%
  filter(!(Week == "w0" & ulcers == "no")) %>%
  rename(Ulcers = ulcers, Time = Week)

response <- df_upa %>%
  group_by(SubjectID, Biopsy_Location) %>%
  nest(Ulcers, Time, .key = "AnyUlcers") %>%
  mutate(remission = map(AnyUlcers, remission), SubjectID = as.character(SubjectID)) %>%
  unnest(remission, .drop = TRUE)

response_all <- merge(response, w0_remitters_upa,
                      by.x = c("SubjectID", "Biopsy_Location"),
                      by.y = c("Patient", "Biopsy_Location"),
                      all.x = TRUE, all.y = FALSE) %>%
  mutate(remission = if_else(remission == "missing", remitter, remission)) %>%
  select(-remitter)


mi <- response_all %>%
  filter(is.na(Biopsy_Location)) %>%
  pull(SubjectID) %>%
  unique()

# Don't know where mi_loc came from
# filter(bd_upa, SubjectID %in% mi_loc, Week < 50) %>%
#   arrange(SubjectID, Biopsy_Location, Week) %>%
#   select(SubjectID, Biopsy_Location, Week, ContainerName, ulcers, BarCode,
#          pSES.CD) %>%
#   group_by(SubjectID, Biopsy_Location) %>%
#   filter(n_distinct(Week) == 1) %>%
#   ungroup()
# Not required as in week 0 it doesn't have ulcers

filter(df_upa, is.na(Biopsy_Location)) %>%
  select(Biopsy_Location, SubjectID, Location) %>%
  distinct()

df_upa <- df_upa %>%
  left_join(distinct(response_all)) %>%
  mutate(remission = if_else(is.na(remission), "yes", remission),
         Time = if_else(Time == "w52", "w14", Time)) %>%
  distinct()


coln <- c("Sample", "Patient", "Time", "remission", "Target", "AU", "Location")
antiTNF_df <- select(df, Sample_id, Pacient_id, Time, remission,
                     "Target Name", AU, biopsied_segment)
colnames(antiTNF_df) <- coln
upa_df <- select(df_upa, Pacient_id, SubjectID, Time, remission,
                 "Target Name", AU, Biopsy_Location)
colnames(upa_df) <- coln

## Controls ####
controls$Time <- "C"
controls$Remission <- "yes"
controls$ulcers <- NA
controls_df <- controls[, c("Sample Name", "Patient_id", "Time",
                            "Remission", "Target Name", "AU", "Location")]
colnames(controls_df) <- coln

## All ####
dff <- rbind(cbind(upa_df, "Study" = "UPA"),
             cbind(antiTNF_df, "Study" = "TNF"),
             cbind(controls_df, "Study" = "C")) %>%
  mutate(Location = tolower(Location)) %>%
  filter(remission != "not interesting") %>%
  # Remove a sample of a gene because there was some problem with the control
  filter(!(Target == "OSM" & Patient == "48" & Study == "TNF" &
             Location == "ileum" & Time == "w14")) %>%
  mutate(General_location = if_else(Location == "ileum" & !is.na(Location),
                                    "ileum", "colon")) %>%
  distinct()

# Strange values
sam <- filter(dff, Study == "UPA",
              AU > 6,
              Target == "SOX6",
              Location == "colon") %>%
  pull("Sample")

stopifnot(nrow(filter(dff, is.na(remission))) == 0)

a <- reals %>%
  filter(grepl("^UP", `Sample Name`)) %>%
  mutate(`Sample Name` = gsub("^UP ", "UP", `Sample Name`)) %>%
  separate(`Sample Name`, c("Id", "Pacient_id")) %>%
  filter(Pacient_id %in% sam) %>%
  arrange(desc(Pacient_id), desc(`Target Name`))


b <- reals %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(undetected = sum(is.na(Cт)), total = n()) %>%
  filter(undetected != 0)

write_xlsx(arrange(b, desc(`Sample Name`), desc(`Target Name`)),
           "processed/undetected_samples.xlsx")
write_xlsx(as.data.frame(a[, 3:7]),
           "processed/anormal_high_markers_colon.xlsx")

# Excluding the samples that at w0 didn't have ulcers!!
write_xlsx(dff, "processed/AU_markers.xlsx")


# Looking if all genes have the same number of samples
o <- dff %>%
  mutate(group = if_else(Time == "w0", "w0", remission)) %>%
  group_by(General_location, Target, Study, group) %>%
  summarise(nSample = n_distinct(Sample)) %>%
  arrange(General_location,  Target, Study, group) %>%
  group_by(Study, group, General_location) %>%
  summarise(n_distinct(nSample), min(nSample), max(nSample),
            diff = max(nSample) - min(nSample))
o
# No they don't one group has even 4 different number of samples...
table(o[, 4, drop = TRUE],
      o[, "diff", drop = TRUE])
# Most differences are of just 1 or in few instances

dff %>%
  mutate(group = if_else(Time == "w0", "w0", remission)) %>%
  filter(Study == "UPA", group == "w0", General_location == "ileum") %>%
  group_by(Target) %>%
  count()

# Some other things ####
bd_upa2 <- bd_upa %>%
  filter(BarCode %in% UPA$Pacient_id)


tb <- table(bd_upa2$Week, bd_upa2$SubjectID, bd_upa2$Biopsy_Location)
dmtb <- dimnames(tb)
check_ulcers <- sapply(dmtb[[3]], function(x) {
  s <- colSums(tb[, , x]) # See if they have samples for that location
  tb[1, , x] == 0 & s != 0
  })
k <- apply(check_ulcers, 1, any)
check_patients_upa <- names(k)[k]

check_ulcers <- colnames(tb)[tb[1, ] == 0]
  group_by(SubjectID) %>%
  summarise(nt = n_distinct(Week)) %>%
  group_by(nt) %>%
  count()

bd2b <- bd2 %>%
  filter(Pacient_id %in% antiTNF$Pacient_id)


tb <- table(bd2b$Time, bd2b$Pacient_id, bd2b$biopsied_segment)
dmtb <- dimnames(tb)
check_ulcers <- sapply(dmtb[[3]], function(x) {
  s <- colSums(tb[, , x]) # See if they have samples for that location
  tb[1, , x] == 0 & s != 0
})
k <- apply(check_ulcers, 1, any)
check_ulcers_tnf <- names(k)[k]

filter(bd, Pacient_id %in% check_ulcers_tnf) %>%
  select(Sample_id, Ulcers, biopsied_segment)
  # Patient 70 on w0 (not done the reals) has ulcers!
  # Patient 103 ulcers ileum on w0 by score on p CEDIS ileum
  # Patient 132 ulcers ileum on w0 by score on p CEDIS ileum
  # Patient 136 ulcers ileum ?? partial CEDIS is not conclusive (ileum is lower)
  # Patient 70 on w14 SESCD global 20 i global CDEIS 0 Red alert: Changed on the database
