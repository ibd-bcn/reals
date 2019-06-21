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

# Correct ids or remove duplicated samples ####
# Do you remember some samples that were in more than one sample?
# Now we remove a plate
reals <- filter(reals,
                !(`Sample Name` %in% c("102w14c", "73w14c")  & grepl("290419bcn 102 a 108", `Experiment Name`)))

reals[reals$`Experiment Name` == "290419 bcn samples plat 3.eds" &
        reals$`Sample Name` == "44w0c", "Sample Name"] <- "41w14"

reals$`Sample Name`[grep("UP11 ", reals$`Sample Name`)] <- "UP11 AB4216302"
reals$`Sample Name`[grep("UP52 ", reals$`Sample Name`)] <- "UP52 Y46412103"
reals$`Sample Name`[grepl("^48C 42W", reals$`Sample Name`, ignore.case = TRUE)] <- "48C 52w14c"

# Y45570803 is duplicated! but it should according to what we know of the samples...

under_expressed <- function(x, target, s, value) {
  keep <- grepl(s, x$`Sample Name`)
  x$Cт[x$`Target Name` == target & keep] <- value
  m <- mean(x$Cт[x$`Target Name` == "BETA ACTINA" & keep])
  x$`Cт Mean`[x$`Target Name` == target & keep] <- value-m
  x
}

bigger_ct <- 45
reals <- reals %>%
  under_expressed("RTNLB", "^UP112 ", bigger_ct) %>%
  under_expressed("IL17A", "^19 CONT ", bigger_ct) %>%
  under_expressed("RTNLB", "^UP84 ", bigger_ct) %>%
  under_expressed("HTR3E", "^1C ", bigger_ct) %>%
  under_expressed("HTR3E", "^39C ", bigger_ct) %>%
  under_expressed("HTR3E", "^67C ", bigger_ct) %>%
  under_expressed("HTR3E", "^UP103 ", bigger_ct)

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

failed_betas <- well_failed %>%
  filter(failed == 1) %>%
  inner_join(reals) %>%
  filter(`Target Name` != "BETA ACTINA",
         !is.na(`Cт`)) %>%
  select(`Sample Name`, `Target Name`) %>%
  distinct()
stopifnot(nrow(failed_betas) == 0)

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

## antiTNF ####


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

df <- merge(bd2, antiTNF, all.x = FALSE, all.y = TRUE,
            by = c("Sample_id", "Pacient_id", "Time")) %>%
  filter(!(Time == "w0" & Ulcers == "no"))

response <- df %>%
  group_by(Pacient_id, biopsied_segment) %>%
  nest(Ulcers, Time, .key = "AnyUlcers") %>%
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
# If we rerun the code above we won't see any sample

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

filter(bd_upa, SubjectID %in% mi_loc, Week < 50) %>%
  arrange(SubjectID, Biopsy_Location, Week) %>%
  select(SubjectID, Biopsy_Location, Week, ContainerName, ulcers, BarCode,
         pSES.CD) %>%
  group_by(SubjectID, Biopsy_Location) %>%
  filter(n_distinct(Week) == 1) %>%
  ungroup()
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

controls$Time <- "C"
controls$Remission <- "yes"
controls$ulcers <- NA
controls_df <- controls[, c("Sample Name", "Patient_id", "Time",
                            "Remission", "Target Name", "AU", "Location")]
colnames(controls_df) <- coln

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

## Compare against controls ####
l <- vector("list", length =  length(unique(dff$Target))*2*2*2)
i <- 1
for (gene in unique(dff$Target)) {
  for (site in unique(dff$Location)) {
    for (study in c("UPA", "TNF")) {
      for (rem in c("yes", "no")){
        d <- filter(dff,
                    Target == gene,
                    Location == site,
                    Study %in% c(study, "C"),
                    remission %in% c(rem, "yes"))
        l[[i]] <- tidy(wilcox.test(AU ~ remission, data = d)) %>%
          mutate(Study = study,
                 Location = site,
                 Target = gene,
                 remission = rem)

        i <- i +1
      }

    }
  }
}

# Corregir fdr per numero de Targets (genes)
vsC <- do.call(rbind, l) %>%
  group_by(Location, Study, remission) %>%
  nest() %>%
  mutate(map(data, ~mutate(., fdr = p.adjust(p.value, n = n(), method = "fdr")))) %>%
  unnest() %>%
  arrange(Location, Target, remission, -p.value) %>%
  select(Location, Target, remission, Study, p.value, fdr, method, alternative)

g <- vsC %>% group_by(Location, remission, Study) %>%
  summarise(n = n_distinct(Target)) %>%
  ungroup() %>%
  distinct(`n`) %>%
  pull()
stopifnot(g == 16)
write_xlsx(vsC, "processed/compare_with_controls.xlsx")

## Compare between studies ####
l <- vector("list", length =  length(unique(dff$Target))*2*2)
i <- 1
for (gene in unique(dff$Target)) {
  for (site in unique(dff$Location)) {
    for (rem in c("w0", "w14 remiters")){
      d <- filter(dff,
                  Target == gene,
                  Location == site,
                  remission == rem,
                  Study != "C")
      l[[i]] <- tidy(wilcox.test(AU ~ Study, data = d)) %>%
        mutate(Location = site,
               Target = gene,
               remission = rem)

      i <- i +1
    }

  }
}

# Corregir fdr per numero de Targets (genes)
between_studies <- do.call(rbind, l) %>%
  group_by(Location, remission) %>%
  nest() %>%
  mutate(map(data, ~mutate(., fdr = p.adjust(p.value, n = n(), method = "fdr")))) %>%
  unnest() %>%
  select(Location, Target, remission, p.value, fdr, method, alternative) %>%
  arrange(Location, Target, remission, -p.value)

write_xlsx(between_studies, "processed/compare_between_studies.xlsx")



## Compare whitin studies ####
l <- vector("list", length =  length(unique(dff$Target))*2*2)
i <- 1
for (gene in unique(dff$Target)) {
  for (site in unique(dff$Location)) {
    for (study in c("UPA", "TNF")){
      d <- filter(dff,
                  Target == gene,
                  Location == site,
                  remission  %in% c("w0", "w14 remiters"),
                  Study == study)
      l[[i]] <- tidy(wilcox.test(AU ~ remission, data = d)) %>%
        mutate(Location = site,
               Target = gene,
               Study = study)

      i <- i +1
    }

  }
}

# Corregir fdr per numero de Targets (genes)
within_studies <- do.call(rbind, l) %>%
  group_by(Study, Location) %>%
  nest() %>%
  mutate(map(data, ~mutate(., fdr = p.adjust(p.value, n = n(), method = "fdr")))) %>%
  unnest() %>%
  select(Location, Target, Study, p.value, fdr, method, alternative) %>%
  arrange(Location, Target, Study, -p.value)

write_xlsx(within_studies, "processed/compare_whitin_studies.xlsx")

# Plots ####

today <- format(Sys.time(), "%Y%m%d")

preplot <- dff %>%
  group_by(remission, Target, Location, Study) %>%
  summarise(meanAU = mean(AU), sem = sd(AU)/sqrt(n())) %>%
  mutate(ymax = meanAU + sem, ymin = meanAU - sem)


# Calculate the line between the studies
d <- preplot %>%
  filter(Study != "C", !grepl("non-remiters", remission)) %>%
  group_by(Target, remission, Location) %>%
  mutate(
    d = max(meanAU) - min(meanAU),
    distU = (max(ymin) - min(ymax))/2,
    orig = max(ymin) - distU/10, # Make them not overlap witht the error bars.
    final = min(ymax) + distU/10,
    center = distU + final)

# Keep just one every two lines
empty <- seq(from = 1, to = nrow(d), by = 2)
d <- d[empty, ]

dw <- merge(between_studies, d) %>%
  filter(fdr < 0.05)

ws <- within_studies %>%
  mutate(remission = "w14 remiters")

db <- preplot %>%
  filter(Study != "C", !grepl("non-remiters", remission)) %>%
  merge(ws) %>%
  unique() %>%
  filter(fdr < 0.05)



### W0 ####

# Compare studies by gene at w0 indepdendently of remission

preplot <- dff %>%
  group_by(Target, Study, remission, Time, General_location) %>%
  summarise(meanAU = mean(AU), sem = sd(AU)/sqrt(n())) %>%
  mutate(ymax = meanAU + sem, ymin = meanAU - sem,
         label = paste(Study, `remission`))

loc <- "colon"

pd <- position_dodge2(width = 0.5)
preplot %>%
  filter(Study != "C", Time == "w0", General_location == loc) %>%
  ggplot(aes(Study, meanAU)) +
  geom_point(aes(col = Study, group = Study, shape = remission),
             position = pd) +
  geom_errorbar(aes(col = Study, ymin = ymin, ymax = ymax), width = 0.2,
                position = pd)  +
  facet_wrap(~Target, scales = "free_y") +
  geom_hline(data = filter(preplot, Study == "C", General_location == loc),
             aes(yintercept = ymin), linetype = "dotted") +
  geom_hline(data = filter(preplot, Study == "C", General_location == loc),
             aes(yintercept = ymax), linetype = "dotted") +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = element_blank(), title = paste("Basal expression:", loc),
       subtitle = "Comparison between UPA and TNF")

preplot <- dff %>%
  group_by(Target, Study, General_location, remission, Time) %>%
  summarise(meanAU = mean(AU), sem = sd(AU)/sqrt(n())) %>%
  mutate(ymax = meanAU + sem, ymin = meanAU - sem)

st <- "UPA"

preplot %>%
  filter(Study != "C", General_location == loc, Study == st) %>%
  mutate(label = paste(Study, remission)) %>%
  ggplot(aes(Time, meanAU)) +
  # geom_point(aes(col = Study, shape = remission), position = pd) +
  geom_pointrange(aes(col = Study, ymin = ymin, ymax = ymax, shape = remission), position = pd) +
  geom_line(aes(group = label, col = Study), position = pd) +
  facet_wrap(~Target, scales = "free_y", ncol = 4) +
  geom_hline(data = filter(preplot, Study == "C", General_location == loc),
             aes(yintercept = ymin), linetype = "dotted") +
  geom_hline(data = filter(preplot, Study == "C", General_location == loc),
             aes(yintercept = ymax), linetype = "dotted") +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = element_blank(), title = "Expression of marker genes",
       subtitle = paste(loc, "at baseline and at week 14"))


## Other things ####



out <- dff %>%
  filter(Study != "C") %>%
  select(-Target, -AU) %>%
  distinct(.keep_all = TRUE) %>%
  group_by(Patient, Location, Study) %>%
  summarise(Times = n_distinct(Time)) %>%
  filter(Times != 2) %>%
  pull("Patient")
  # group_by(n) %>%
  # count()

out2 <- dff %>%
  filter(Patient %in% out,
         Time == "w0") %>%
  filter(Study != "C") %>%
  select(-Target, -AU) %>%
  distinct(.keep_all = TRUE) %>%
  select(Sample, Patient, Study, Location) %>%
  arrange(Study, Patient, Sample, Location)
out2 %>%
  write_xlsx("processed/missing_response.xlsx")

df_upa %>%
  mutate(SubjectID = as.character(SubjectID), Location = tolower(Location)) %>%
  inner_join(out2, by = c("SubjectID" = "Patient", "Location" = "Location")) %>%
  select(SubjectID, pSES.CD, Location, Biopsy_Location, Week) %>%
  distinct() %>%
  arrange(desc(SubjectID), Week) %>%
  write_xlsx("processed/missing_response_upa.xlsx")


a <- bd %>%
  distinct(Pacient_id, week, biopsied_segment, Ulcers) %>%
  filter(!(week == "0" & Ulcers == "no")) %>%
  group_by(Pacient_id, biopsied_segment) %>%
  nest(Ulcers, week, .key = "AnyUlcers") %>%
  mutate(remission = map(AnyUlcers, adf))




preplot %>%
  filter(Study != "C", Target %in% c("AQP7", "DERL3", "OSM", "S100A8")) %>%
  ggplot(aes(remission, meanAU)) +
  geom_point(aes(col = Study, group = Study)) +
  expand_limits(y = 0) +
  geom_errorbar(aes(col = Study, group = Study, ymin = ymin, ymax = ymax), width = 0.2) +
  geom_line(aes(col = Study, group = Study)) +
  geom_hline(data = filter(preplot, Study == "C",
                           Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
             aes(yintercept = ymin), linetype = "dotted") +
  geom_hline(data = filter(preplot, Study == "C", Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
             aes(yintercept = ymax), linetype = "dotted") +
  geom_segment(data = filter(dw,  Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
               aes(x = remission, y = orig, xend = remission, yend = final)) +
  geom_text(data = filter(dw,  Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
            aes(x = remission, y = center), label = "\t*", size = 5) + # Dirty trick to dodge the symbol
  geom_text(data = filter(db,  Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
            aes(x = remission, y = meanAU), label = "\t+", size = 3) + # Dirty trick to dodge the symbol
  facet_wrap(Target ~ Location, scales = "free_y", ncol = 4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  geom_segment(aes(xend = remission, yend = meanAU, group = remission)) +
  ylab("AU (mean\u00B1SEM)")


# tests ####
within_studies_baseline <- dff %>%
  filter(Study != "C", Time == "w0") %>%
  group_by(General_location, Target, Study) %>%
  nest(AU, remission) %>%
  mutate(model = map(data, ~ tidy(wilcox.test(AU ~ remission, data = .)))) %>%
  unnest(model) %>%
  group_by(Study, General_location) %>%
  nest(.key = "data3") %>%
  mutate(map(data3, ~mutate(., fdr = p.adjust(p.value, n = n(), method = "fdr")))) %>%
  unnest() %>%
  select(General_location, Target, Study, p.value, fdr, method, alternative) %>%
  arrange(General_location, Target, Study, desc(p.value))

write_xlsx(within_studies_baseline, "processed/within_studies_baselines.xlsx")

within_studies_time <- dff %>%
  filter(Study != "C") %>%
  group_by(General_location, Target, Study, remission) %>%
  nest(AU, Time) %>%
  mutate(model = map(data, ~ tidy(wilcox.test(AU ~ Time, data = .)))) %>%
  unnest(model) %>%
  group_by(Study, General_location) %>%
  nest(.key = "data3") %>%
  mutate(map(data3, ~mutate(., fdr = p.adjust(p.value, n = n(), method = "fdr")))) %>%
  unnest() %>%
  select(General_location, remission, Target, Study, p.value, fdr, method, alternative) %>%
  arrange(General_location, remission, Target, Study, desc(p.value))
write_xlsx(within_studies_time, "processed/within_studies_time.xlsx")
