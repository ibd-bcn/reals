## Real times Azu (UPA) review 2021 ##

library("xlsx")
library("gdata", warn.conflicts = FALSE)
library("ggplot2", warn.conflicts = FALSE)
library("ggpubr")
library("tidyr", quietly = TRUE, warn.conflicts = FALSE)
library("ggthemes", warn.conflicts = FALSE)
library("Rmisc", quietly = TRUE, warn.conflicts = FALSE, mask.ok = "mutate")
library("dplyr", quietly = TRUE, warn.conflicts = FALSE)
library("gridExtra", quietly = TRUE, warn.conflicts = FALSE)
library("grid")


# Dades
new_files <- list.files(pattern = "*.xls", path = "data/reals2021", full.names = TRUE)
dades0 <- lapply(new_files, readxl::read_xls, skip = 7, n_max = 192)
for (i in seq_along(dades0)) {
  dades0[[i]] <- cbind(Experiment = basename(new_files[i]), dades0[[i]])
}

merger <- function(x, y) {
    merge(x, y, all = TRUE, sort = FALSE)
}
dades <- Reduce(merger, dades0)
dades <- dades %>%
  filter(nchar(Well) <= 3, !is.na(Well), # Filter the end info of the tables
         !is.na(`Sample Name`))  # remove blanks

# Read other sample's ID ####
RT_BCN <- readxl::read_xlsx("data/Samples.xlsx", sheet = 2)
RT_BCN <- RT_BCN[!is.na(RT_BCN$`RT TUBE`), ]
colnames(RT_BCN)[2] <- "SAMPLE"
# Some samples are not available even if there was a RT tube
# (expected as they don't appear on the template)
RT_BCN[!(RT_BCN$`RT TUBE` %in% dades$Sample_name) & RT_BCN$`RT TENEMOS` != "NO", ]

RT_UPA <- readxl::read_xlsx("data/Samples.xlsx", sheet = 4)
RT_UPA <- RT_UPA[!is.na(RT_UPA$`RT TUBE`), ]
colnames(RT_UPA)[2] <- "SAMPLE"
RT_CTRL <- readxl::read_xlsx("data/Samples.xlsx", sheet = 6)
RT_CTRL <- RT_CTRL[!is.na(RT_CTRL$`RT TUBE`), ]
colnames(RT_CTRL)[2] <- "SAMPLE"
RT_CTRL$`RT TUBE` <- gsub(pattern = " CONT", replacement = "", x = RT_CTRL$`RT TUBE`)
RT_CTRL$`RT TUBE` <- paste0("Sample ", RT_CTRL$`RT TUBE`)
# Merge all samples together
RT <- rbind(RT_BCN[, 1:2], RT_UPA[, 1:2], RT_CTRL[, 1:2])

# Use easier colnames
colnames(dades)[colnames(dades) == "Cт"] <- "CT"
colnames(dades)[colnames(dades) == "Sample Name"] <- "Sample_name"
colnames(dades)[colnames(dades) == "Target Name"] <- "Target.Name"
dades$CT[grep(pattern = "Und", dades$CT)] <- NA
dades$CT <- as.numeric(dades$CT)

# Samples Ids ####
# Rename a duplicated sample due to piepetting
# Just 2 wells for HDC gene
stopifnot(sum(endsWith(dades$Sample_name, "BCN")) == 4)
dades$Sample_name[endsWith(dades$Sample_name, "BCN")] <- "107C"

dades2 <- merge(dades, RT,
                by.y = "RT TUBE", by.x = "Sample_name",
                all.x = TRUE, sort = FALSE)
dades2$SAMPLE <- gsub(pattern = "^0", replacement = "", x = dades2$SAMPLE)

# QC check ####
# Samples missing both genes in the same well
bad_well <- dades2 %>%
  filter(is.na(CT)) %>%
  select(Experiment, Well, Sample_name, Target.Name) %>%
  group_by(Experiment, Well) %>%
  summarize(n = n(), sample = unique(Sample_name)) %>%
  filter(n == 2) %>%
  select(-n) %>%
  filter(!is.na(sample)) # Blanks
bad_well
stopifnot(nrow(bad_well) == 3)
# Samples missing a gene from the two wells
bad_genes <- dades2 %>%
  filter(is.na(CT)) %>%
  select(Experiment, Well, Sample_name, Target.Name) %>%
  group_by(Sample_name, Target.Name) %>%
  summarize(n = n(), wells = list(Well), Exp = unique(Experiment)) %>%
  filter(n != 1, !is.na(Sample_name)) %>%
  select(-n)
bad_genes
stopifnot(nrow(bad_genes) == 8)

# Fill the missing values of ADCYAP with 40 as it is a lowly expressed gene
dades2$CT[dades2$Target.Name == "ADCYAP1" &
           is.na(dades2$CT) &
           dades2$Sample_name %in% tail(bad_genes$Sample_name, 6)] <- 40

# Check beta actina ####
nested_well <- dades2 %>%
  filter(Target.Name != "BETA ACTINA") %>%
  group_by(Experiment, Target.Name, Sample_name) %>%
  nest(data = Well)

nested_well$data <- lapply(nested_well$data, unlist)

# Look for the wells how many have values
xy <- vector("list", length = nrow(nested_well))
for (i in seq_len(nrow(nested_well))) {
  x <- nested_well[i, ]
  wells_value <- dades2 %>%
    filter(`Sample_name` == x[["Sample_name"]],
           `Experiment` == x[["Experiment"]],
           Well %in% x[["data"]][[1]]
    )
  xy[[i]] <- wells_value %>%
    group_by(`Target.Name`) %>%
    summarise(nas = sum(is.na(`CT`)))
}

# Look how many missing values of Beta actina per target sample, experiment
ba <- sapply(xy, function(x) {
  x[x$`Target.Name` == "BETA ACTINA", "nas", drop = TRUE]})
table(lengths(ba))
table(unlist(ba))

# Look how many missing values of gene per target sample, experiment
gene <- sapply(xy, function(x){
  x[x$`Target.Name` != "BETA ACTINA", "nas", drop = TRUE]})
table(gene)

# Calculate AU ####
## Barcelona ####
antiTNF <- dades2 %>%
  filter(endsWith(Sample_name, "C")) %>%
  arrange(Experiment, SAMPLE, Sample_name) %>%
  group_by(Experiment, SAMPLE, Sample_name) %>%
  summarize(
    Gene = unique(Target.Name[Target.Name != "BETA ACTINA"]),
    base_expr = mean(CT[Target.Name == "BETA ACTINA"], na.rm = TRUE),
    base_n = sum(!is.na(CT[Target.Name == "BETA ACTINA"])),
    target_expr = mean(CT[Target.Name != "BETA ACTINA"], na.rm = TRUE),
    target_n = sum(!is.na(CT[Target.Name != "BETA ACTINA"])),
    dCT = base_expr - target_expr,
    AU = 2^dCT * 1000) %>%
  ungroup()


# Check that we only have a missing sample
stopifnot(all(antiTNF$base_n <= 2))
stopifnot(all(antiTNF$target_n <= 2))
stopifnot(sum(is.na(antiTNF$AU)) == 1)
stopifnot(sum(is.na(antiTNF$AU[antiTNF$Sample_name == "107C"])) == 1)
antiTNF <- antiTNF %>%
  select(-base_n, -target_n) %>%
  filter(!is.na(AU))

## Controls ####
# Find which wells have the data for the same genes of samples of experiments
ctrl_groups <- dades2 %>%
  filter(startsWith(Sample_name, "Sample")) %>%
  arrange(Experiment, SAMPLE, Sample_name, Well) %>%
  group_by(Experiment, SAMPLE, Sample_name) %>%
  filter(Target.Name != "BETA ACTINA") %>%
  mutate(group = Target.Name) %>%
  select(Well, group)


ctrl <- dades2 %>%
  filter(startsWith(Sample_name, "Sample")) %>%
  arrange(Experiment, SAMPLE, Sample_name, Well) %>%
  full_join(ctrl_groups) %>%
  group_by(Experiment, SAMPLE, Sample_name, group) %>%
  summarize(
    Gene = unique(Target.Name[Target.Name != "BETA ACTINA"]),
    base_expr = mean(CT[Target.Name == "BETA ACTINA"], na.rm = TRUE),
    base_n = sum(!is.na(CT[Target.Name == "BETA ACTINA"])),
    target_expr = mean(CT[Target.Name != "BETA ACTINA"], na.rm = TRUE),
    target_n = sum(!is.na(CT[Target.Name != "BETA ACTINA"])),
    dCT = base_expr - target_expr,
    AU = 2^dCT * 1000) %>%
  ungroup()

# Check that we used only duplicates
stopifnot(all(ctrl$group == ctrl$Gene))
stopifnot(all(ctrl$base_n == 2))
stopifnot(all(ctrl$target_n == 2))
ctrl <- select(ctrl, -group, -base_n, -target_n)

## UPA ####

# Merge metadata ####
# To remember on this datasets there is only one segment per patient
# Read files with the information
bd <- readxl::read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = c("n.a.", ""))
bd_upa <- readxl::read_xls("data/M13-740_abbvie_database_300519.xls", na = c("n.a.", "NA"))

# Function to know if there were Ulcers
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


## Barcelona ####
bd2 <- bd %>%
  mutate(Sample_id = tolower(Sample_id)) %>%
  filter(stringr::str_detect(Sample_id, "-w"), IBD == "CD") %>%
  select(Sample_id, CD_endoscopic_remission, CD_endoscopic_response,
         sample_location, Ulcers, biopsied_segment) %>%
  mutate(Sample_id = gsub(" reseq|rep", "", Sample_id),
         Pacient_id = as.character(as.numeric(gsub("-w.*", "", Sample_id))),
         Time = stringr::str_extract(Sample_id, "w[0-9]*"),
         Sample_id = paste0(Pacient_id, "-", Time)) %>%
  filter(Time %in% c("w0", "w14")) %>%
  distinct() %>%
  rename(Location = sample_location)

antiTNF <- antiTNF %>%
  mutate(
    SAMPLE = tolower(SAMPLE),
    Pacient_id = as.character(as.numeric(gsub("-w.*", "", SAMPLE))),
    Time = stringr::str_extract(SAMPLE, "w[0-9]*"),
    Sample_id = paste0(Pacient_id, "-", Time),
  )

df <- merge(bd2, antiTNF, all.x = FALSE, all.y = TRUE,
            by = c("Sample_id", "Pacient_id", "Time")) %>%
  filter(!(Time == "w0" & Ulcers == "no"))

# Check if they are in remission
response <- df %>%
  group_by(Pacient_id, biopsied_segment) %>%
  nest(AnyUlcers = c(Ulcers, Time)) %>%
  mutate(remission = purrr::map_chr(AnyUlcers, remission)) %>%
  ungroup() %>%
  select(-AnyUlcers)
# Those missing had as no ulcers as the code on read.R says and
# I recall from exchange with Azucena Salas
response$remission[response$remission == "missing"] <- "no"
df <- left_join(df, response)

## Controls ####

ctrl$Time <- "C"
ctrl$Remission <- "yes"
ctrl$ulcers <- NA


## UPA ####

upa_codes <- readxl::read_xlsx("data/RT UPA DISEASE BLOC 1 2 3 500NG 20UL.xlsx",
                               range = cell_cols(1:7)) %>%
  filter(!is.na(`PLACA MICRONIC`))

bd_upa <- bd_upa %>%
  select(Week, pSES.CD, BarCode, ContainerName, SubjectID,
         ulcers, Biopsy_Location) %>%
  distinct()

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
  mutate(SubjectID = as.character(SubjectID)) %>%
  group_by(SubjectID, Biopsy_Location) %>%
  nest(AnyUlcers = c(Ulcers, Time)) %>%
  mutate(remission = purrr::map_chr(AnyUlcers, remission)) %>%
  unnest(AnyUlcers)

response_all <- merge(response, w0_remitters_upa,
                      by.x = c("SubjectID", "Biopsy_Location"),
                      by.y = c("Patient", "Biopsy_Location"),
                      all.x = TRUE, all.y = FALSE) %>%
  mutate(remission = unlist(remission),
         remission = if_else(remission == "missing", remitter, remission)) %>%
  select(-remitter)

df_upa <- df_upa %>%
  left_join(distinct(response_all)) %>%
  mutate(remission = if_else(is.na(remission), "yes", remission),
         Time = if_else(Time == "w52", "w14", Time)) %>%
  distinct()


#  Write
write.xlsx(
  x = dCT, file = "processed/mostres_reals21.xlsx",
  sheetName = "1", col.names = TRUE, row.names = FALSE
)

# Merge all data ####
coln <- c("Sample", "Patient", "Time", "remission", "Target", "AU", "Location")
# Reorder columns and remove some of them.

# plots ######

## colon
taula_colon <- taula[taula$ContainerName == "Biopsy Colon", ]
taula_colon <- taula_colon[taula_colon$pSES.CD == "0", ]

plot <- ggplot(data = taula_colon, aes(x = week, y = AU, fill = Biopsy_Location)) +
  geom_boxplot(
    outlier.colour = "white",
    fill = "white"
  ) +
  geom_point(shape = 21, aes(size = pSES.CD), position = position_jitterdodge(jitter.width = 0.1)) +
  labs(title = "KLRC2 - Colon Samples") +
  theme_classic()

plotcolon <- plot + stat_compare_means(
  method = "wilcox.test",
  comparisons = list(c("0", "14")),
  label = "p.format",
  tip.length = 0.01
)
plotcolon

## ileum
taula_ileum <- taula[taula$ContainerName == "Biopsy Ileum", ]
taula_ileum <- taula_ileum[taula_ileum$pSES.CD == "0", ]

plot <- ggplot(data = taula_ileum, aes(x = week, y = AU, fill = pSES.CD)) +
  geom_boxplot(
    outlier.colour = "white",
    fill = "white"
  ) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.98)) +
  scale_fill_manual(
    values = c("#ffebee", "#ff8a80", "#d32f2f", "#c62828", "#b71c1c", "#ffcdd2", "#ef9a9a", "#e57373", "#ef5350", "#ff1744", "#f44336", "#e53935", "#d50000", "gray"),
    labels = levels(taula$pSES.CD),
    name = "Partial SES-CD"
  ) +
  labs(title = "KLRC2 - Ileum Samples") +
  theme_classic()

plotileum <- plot + stat_compare_means(
  method = "wilcox.test",
  comparisons = list(c("0", "14")),
  label = "p.format",
  tip.length = 0.01
)
plotileum


# pdf #
pdf("processed/Reals21_pSESCD_0.pdf")
plotcolon
plotileum

dev.off()
