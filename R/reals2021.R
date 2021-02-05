## Real times Azu (UPA) review 2021 ##
library("xlsx")
library("gdata", warn.conflicts = FALSE)
library("ggplot2", warn.conflicts = FALSE)
library("ggpubr")
library("tidyr", quietly = TRUE, warn.conflicts = FALSE)
library("ggthemes", warn.conflicts = FALSE)
library("Rmisc", quietly = TRUE, warn.conflicts = FALSE, mask.ok = "mutate")
library("gridExtra", quietly = TRUE, warn.conflicts = FALSE)
library("grid")
library("purrr")
library("broom")
library("readxl")
library("dplyr", quietly = TRUE, warn.conflicts = FALSE)

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
stopifnot(nrow(RT_BCN[!(RT_BCN$`RT TUBE` %in% dades$Sample_name) & RT_BCN$`RT TENEMOS` != "NO", ]) == 100)

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
# Use consistent sample name
upa_logic <- endsWith(dades$Sample_name, "UP") |
  endsWith(dades$Sample_name, "UPA") |
  endsWith(dades$Sample_name, "UPA REPETICIO TPSAB1") # Watch out later
samples_upa <- dades$Sample_name[upa_logic]
dades$Sample_name[upa_logic] <- gsub(
  pattern = "Sample\\s*(\\d+)\\s*UPA?.*",
  replacement = "UP\\1",
  x = samples_upa)


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
stopifnot(nrow(bad_well) == 5) # 3 From BCN and 2 from UPA
# Samples missing a gene from the two wells
bad_genes <- dades2 %>%
  filter(is.na(CT)) %>%
  select(Experiment, Well, Sample_name, Target.Name) %>%
  group_by(Sample_name, Target.Name) %>%
  summarize(n = n(), wells = list(Well), Exp = unique(Experiment)) %>%
  filter(n != 1, !is.na(Sample_name)) %>%
  select(-n)
bad_genes
stopifnot(nrow(bad_genes) == 8 + 3) # 8 from BCN and 3 from UPA

to_be_correct_genes <- bad_genes %>%
  filter(Target.Name == "ADCYAP1") %>%
  pull(Sample_name)
# Fill the missing values of ADCYAP with 40 as it is a lowly expressed gene
dades2$CT[dades2$Target.Name == "ADCYAP1" &
           is.na(dades2$CT) &
           dades2$Sample_name %in% to_be_correct_genes] <- 40

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
stopifnot(sum(is.na(antiTNF$AU)) <= 1)
stopifnot(sum(is.na(antiTNF$AU[antiTNF$Sample_name == "107C"])) == 1)
antiTNF <- antiTNF %>%
  filter(base_n >= 1) %>%
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
  select(Well, group) %>%
  ungroup()


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
# Find which wells have the data for the same genes of samples of experiments
upa_groups <- dades2 %>%
  filter(startsWith(Sample_name, "UP")) %>%
  arrange(Experiment, SAMPLE, Sample_name, Well) %>%
  group_by(Experiment, SAMPLE, Sample_name) %>%
  filter(Target.Name != "BETA ACTINA") %>%
  mutate(group = Target.Name) %>%
  select(Well, group) %>%
  ungroup()

upa <- dades2 %>%
  filter(startsWith(Sample_name, "UP")) %>%
  arrange(Experiment, SAMPLE, Sample_name, Well) %>%
  full_join(upa_groups) %>%
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
stopifnot(all(upa$group == upa$Gene))
stopifnot(all(upa$base_n <= 2))
stopifnot(all(upa$target_n <= 2))
stopifnot(sum(is.na(upa$AU)) == 0)
upa <- upa %>%
  filter(base_n >= 1) %>%
  select(-group, -base_n, -target_n)

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

df_bcn <- merge(bd2, antiTNF, all.x = FALSE, all.y = TRUE,
            by = c("Sample_id", "Pacient_id", "Time")) %>%
  filter(!(Time == "w0" & Ulcers == "no"))

# Check if they are in remission
response <- df_bcn %>%
  group_by(Pacient_id, biopsied_segment) %>%
  nest(AnyUlcers = c(Ulcers, Time)) %>%
  mutate(remission = purrr::map_chr(AnyUlcers, remission)) %>%
  ungroup() %>%
  select(-AnyUlcers)
# Those missing had as no ulcers as the code on read.R says and
# I recall from exchange with Azucena Salas
response$remission[response$remission == "missing"] <- "no"
df_bcn <- left_join(df_bcn, response) %>%
  select(-Experiment, -SAMPLE, -base_expr, -target_expr,
         -CD_endoscopic_remission, -CD_endoscopic_response,
         -Sample_name) %>%
  mutate(Study = "BCN",
         remitter = case_when(
           Time == "w0" ~ "w0",
           Time == "w14" & remission == "yes" ~ "w14 R",
           Time == "w14" & remission == "no" ~ "w14 NR",
           TRUE ~ "strange"
         )) %>%
  ungroup()

## Controls ####

ctrl$Time <- "C"
ctrl$Study <- "C"
ctrl$remission <- "yes"
ctrl$Ulcers <- NA
ctrl$remitter <- NA
ctrl$Location <- ifelse(endsWith(ctrl$SAMPLE, "ILI"), "ileum", "colon")
ctrl <- merge(ctrl, RT_CTRL[, 1:2])
ctrl <- select(ctrl, -base_expr, target_expr, -dCT,
               -Sample_name, -Experiment)

## UPA ####
upa_codes <- read_xlsx("data/RT UPA DISEASE BLOC 1 2 3 500NG 20UL.xlsx",
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
df_upa <- merge(upa, upa_codes,  by.x = "Sample_name", by.y = "RT TUBE",
                all.x = TRUE, all.y = FALSE) %>%
  left_join(bd_upa, by = c("nº muestra" = "BarCode")) %>%
  mutate(ContainerName = gsub("Biopsy ", "", ContainerName),
         Location = ContainerName,
         Week = paste0("w", Week),
         pSES.CD = as.numeric(pSES.CD),
         remission_old = case_when(
           Week == "w0" ~ "w0",
           Week != "w0" & pSES.CD < 5 ~ "w14 R",
           Week != "w0" & pSES.CD >= 5 ~ "w14 NR",
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

w0_remitters_upa <- read_xlsx("processed/missing_response based on ulcers UPA.xlsx") %>%
  mutate(remitter = tolower(remitter)) %>%
  rename(Biopsy_Location = "Biopsy_Location week 0") %>%
  select(Patient, remitter, Biopsy_Location)
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
  distinct() %>%
  ungroup() %>%
  select(-Experiment) %>%
  mutate(Study = "UPA")

#  Write
write.xlsx(
  x = df_upa, file = "processed/mostres_reals21_upa.xlsx",
  sheetName = "1", col.names = TRUE, row.names = FALSE
)

# Merge all data ####
coln <- c("Sample", "Patient", "Time", "remission", "Target", "AU", "Location", "Study")
# Reorder columns and remove some of them.

df_upa2 <- df_upa %>%
  select(Sample_name, SAMPLE, Time, remission, Gene, AU, Location, Study) %>%
  mutate(remission = case_when(
    Time == "w0" ~ "w0",
    Time != "w0" & remission == "no"~ "w14 NR",
    Time != "w0" & remission == "yes"~ "w14 R",
    TRUE ~ "what"
  ))
df_bcn2 <- select(df_bcn,
                 Sample_id, Pacient_id, Time, remitter, Gene, AU, Location, Study)
ctrl2 <- ctrl %>%
  mutate(Pacient_id = gsub(pattern = "^(C\\d+)-.+", replacement =  "\\1",
                           x = SAMPLE)) %>%
  select(SAMPLE, Pacient_id, Time, remission, Gene, AU, Location, Study)
colnames(df_upa2) <- coln
colnames(df_bcn2) <- coln
colnames(ctrl2) <- coln
dff <- rbind(df_upa2, ctrl2, df_bcn2)
dff$Location <- tolower(dff$Location)

write.csv(dff, "processed/mostres_reals21_all.csv", row.names = FALSE)
write.xlsx(
  x = dff, file = "processed/mostres_reals21_all.xlsx",
  sheetName = "1", col.names = TRUE, row.names = FALSE
)

ctrl_bcn <- merge(ctrl, df_bcn, all = TRUE) # Set colnames
write.xlsx(ctrl, file = "processed/new_genes_ctrl.xlsx", row.names = FALSE)

# Stats ####
## study vs control ####
colnames(ctrl)[1] <- c("Sample_id")
colnames(ctrl)[3] <- c("Pacient_id")

length_total <- length(unique(dff$Target)) *
  length(unique(dff$Location)) *
  (length(unique(dff$Study)) - 1) # Remove the control group
l <- vector("list", length =  length_total)
i <- 1
for (gene in unique(dff$Target)) {
  for (site in unique(dff$Location)) {
    for (study in c("BCN", "UPA")) {
      d <- filter(dff,
                  Target == gene,
                  Location == site,
                  Time %in% c("C", "w0"),
                  Study %in% c(study, "C"))
      l[[i]] <- tidy(wilcox.test(AU ~ Study, data = d)) %>%
        mutate(Study = study,
               Location = site,
               Time = "w0",
               Target = gene)
      i <- i + 1
    }
  }
}

# Corregir fdr per numero de Targets (genes)
vsC <- do.call(rbind, l) %>%
  group_by(Location, Study) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr")) %>%
  arrange(Location, Target, -p.value) %>%
  select(Location, Target, Study, p.value, fdr, method, alternative) %>%
  mutate(test = paste(Study, "w0 vs C")) %>%
  as.data.frame()

write.xlsx(vsC, file = "processed/new_genes_vs_C.xlsx", row.names = FALSE)

## antiTNF vs UPA ####
# Falten
bcnVSupa <- dff %>%
  filter(Study != "C",
         remission == "w0") %>%
  nest_by(Target, Location, remission) %>%
  summarize(broom::tidy(wilcox.test(AU ~ Study, data = data))) %>%
  group_by(Location) %>%
  mutate(fdr = p.adjust(p.value, method = "fdr"),
         test = "BCN w0 vs UPA w0") %>%
  ungroup() %>%
  as.data.frame()
write.xlsx(bcnVSupa, file = "processed/new_genes_w0_BCN_vs_UPA.xlsx",
           row.names = FALSE)

## w0 vs w14 remitters ####
w0vsw14r <- dff %>%
  filter(Study != "C",
         remission != "w14 NR") %>%
  nest_by(Target, Location, Study) %>%
  summarize(broom::tidy(wilcox.test(AU ~ remission, data = data))) %>%
  ungroup() %>%
  group_by(Location, Study) %>%
  mutate(fdr = p.adjust(p.value, n = n(), method = "fdr"),
         test = "w0 vs w12/16 R") %>%
  ungroup() %>%
  select(-statistic) %>%
  arrange(Location, Target, Study, -p.value)

## w0 vs w14 non-remitters ####
w0vsw14nr <- dff %>%
  filter(Study != "C",
         remission != "w14 R") %>%
  nest_by(Target, Location, Study) %>%
  summarize(broom::tidy(wilcox.test(AU ~ remission, data = data))) %>%
  ungroup() %>%
  group_by(Location, Study) %>%
  mutate(fdr = p.adjust(p.value, n = n(), method = "fdr"),
         test = "w0 vs w12/16 NR") %>%
  ungroup() %>%
  select(-statistic) %>%
  arrange(Location, Target, Study, -p.value)

## Merge w0 vs w14 ####
w0vsw14 <- rbind(w0vsw14nr, w0vsw14r) %>%
  as.data.frame()
write.xlsx(w0vsw14, "processed/new_genes_statistics_vs_w0.xlsx",
           row.names = FALSE)

## controls SEM ####
ctrl_SEM <- ctrl %>%
  group_by(Location, Gene) %>%
  summarize(mean_se(AU)) %>%
  select(-y) %>%
  ungroup() %>%
  as.data.frame()
write.xlsx(x = ctrl_SEM, file = "processed/new_genes_controls_statistics_SEM.xlsx",
           row.names = FALSE)

df_bcn <- as.data.frame(df_bcn)
write.xlsx(df_bcn, file = "processed/new_genes_AU_antitnf.xlsx")

# plots #####

## Figure 3A ####
dff %>%
  filter(Location == "colon",
         Time %in% c("C", "w0")) %>%
  mutate(Study = forcats::fct_relevel(Study, c("C","UPA", "BCN"))) %>%
  ggplot() +
  geom_boxplot(aes(Study, AU)) +
  facet_wrap(~Target, scales = "free_y") +
  theme_minimal() +
  labs(title = "colon", x = element_blank())
ggsave("Figures/reals21_new_genes_studies_colon.png")

## Figure 3B ####
dff %>%
  filter(Location == "ileum",
         Time %in% c("C", "w0")) %>%
  mutate(Study = forcats::fct_relevel(Study, c("C", "UPA", "aTNF"))) %>%
  ggplot() +
  geom_boxplot(aes(Study, AU)) +
  facet_wrap(~Target, scales = "free_y") +
  theme_minimal() +
  labs(title = "ileum", x = element_blank())
ggsave("Figures/reals21_new_genes_studies_ileum.png")

## Figure 4A ####
dff %>%
  filter(Location == "colon",
         Study == "BCN") %>%
  mutate(remission = gsub("w14", "week 12/16", x = remission)) %>%
  mutate(remission = gsub("remitter", "R", x = remission)) %>%
  mutate(remission = gsub("non-", "N", x = remission)) %>%
  mutate(remission = forcats::fct_relevel(as.factor(remission),
                                          c("w0", "week 12/16 R", "week 12/16 NR"))) %>%
  ggplot() +
  geom_boxplot(aes(remission, AU)) +
  geom_hline(data = ctrl_SEM[ctrl_SEM$Location == "colon", ],
             aes(yintercept = ymin),
             linetype = 5) +
  geom_hline(data = ctrl_SEM[ctrl_SEM$Location == "colon", ],
             aes(yintercept = ymax),
             linetype = 5) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "colonic antiTNF", x = element_blank())
ggsave("Figures/reals21_new_genes_colon_antiTNF.png")

## Figure 4B ####
dff %>%
  filter(Location == "colon",
         Study == "UPA") %>%
  mutate(remission = gsub("w14", "week 12/16", x = remission)) %>%
  mutate(remission = gsub("remitter", "R", x = remission)) %>%
  mutate(remission = gsub("non-", "N", x = remission)) %>%
  mutate(remission = forcats::fct_relevel(as.factor(remission),
                                          c("w0", "week 12/16 R", "week 12/16 NR"))) %>%
  ggplot() +
  geom_boxplot(aes(remission, AU)) +
  geom_hline(data = ctrl_SEM[ctrl_SEM$Location == "colon", ],
             aes(yintercept = ymin),
             linetype = 5) +
  geom_hline(data = ctrl_SEM[ctrl_SEM$Location == "colon", ],
             aes(yintercept = ymax),
             linetype = 5) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_minimal() +
  labs(title = "colonic UPA", x = element_blank())
ggsave("Figures/reals21_new_genes_colon_UPA.png")
