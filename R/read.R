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
bd_upa <- read_xls("data/M13-740_abbvie_database_210918.xls")
reals <- read_xls("data/SOX6 CHI3L1 PDGFD PTGDR2_20190516_115511_Results_Export.xls",
    skip = 15, na = c("", "Undetermined"))
reals2 <- read_xls("data/THY HTR3E RETNLB COL3A1_20190524_165141_Results_Export.xls",
    skip = 15, na = c("", "Undetermined"))
upa_codes <- read_xlsx("data/RT UPA DISEASE BLOC 1 2 3 500NG 20UL.xlsx")
control_codes <- read_xlsx("data/RT CONTROLS 500NG 20 UL.xlsx")
antiTNF_codes <- read_xlsx("data/RT biopsies BCN bloc A B C D 500 ng en 20 ul.xlsx")

# Clean the original data ####
reals <- merge(reals, reals2, all = TRUE)

# Check if some samples are in more than one plate
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

well_failed %>%
    group_by(failed) %>%
    count()

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

sd_Ct <- group_by(reals, `Target Name`, `Sample Name`) %>%
    summarise(sd = sd(`Cт`), mean = mean(`Cт`))


erroneous <- filter(sd_Ct, `sd` > 0.5 | mean > 25 & !is.na(sd), `Target Name` == "BETA ACTINA")
if (nrow(erroneous) > 1) {
  write_xlsx(erroneous, "processed/big_sd.xlsx")

  pdf("Figures/histograms.pdf")
  hist(reals$Cт[reals$`Target Name` == "BETA ACTINA"], main = "Beta Actina", xlab = "Ct")
  hist(reals$Cт[reals$`Target Name` == "CHI3L1"], main = "CHI3L1", xlab = "Ct")
  hist(reals$Cт[reals$`Target Name` == "PDGFD"], main = "PDGFD", xlab = "Ct")
  hist(reals$Cт[reals$`Target Name` == "SOX6"], main = "SOX6", xlab = "Ct")
  dev.off()
}

# Correct ids or remove duplicated samples ####
# Do you remember some samples that were in more than one sample?
# Now we remove a plate
reals <- filter(reals,
            !(`Sample Name` %in% c("102w14c", "73w14c")  & grepl("290419bcn 102 a 108", `Experiment Name`)))

reals[reals$`Experiment Name` == "290419 bcn samples plat 3.eds" &
    reals$`Sample Name` == "44w0c", "Sample Name"] <- "41w14"

reals$`Sample Name`[grep("UP11 ", reals$`Sample Name`)] <- "UP11 AB4216302"
reals$`Sample Name`[grep("UP52 ", reals$`Sample Name`)] <- "UP52 Y46412103"
reals$`Sample Name`[grepl("42W", reals$`Sample Name`, ignore.case = TRUE)] <- "48C 52w14c"

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
         Location = case_when(grepl("i$", `Sample Name`) ~ "ileum",
                              grepl("c$", `Sample Name`) ~ "colon"),
         Sample_id = paste0(Pacient_id, "-", Time),
         AU = 2^(-`ΔCт`)*1000) %>%
  select(-`ΔCт`)

# looking at the table
antiTNF$Location[antiTNF$`Sample Name` == "04w14"] <- "colon"

# To verify
antiTNF %>%
  group_by(Id, `Sample Name`) %>%
  count() %>%
  ungroup() %>%
  group_by(Id) %>%
  count() %>%
  filter(n != 1)

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
stopifnot(length(incorrect) == 0)
filter(upa_codes, `RT TUBE` %in% incorrect) %>%
  mutate(sample = paste(`RT TUBE`, `nº muestra`)) %>%
  pull(sample)

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
  group_by(Id) %>%
  count() %>%
  filter(n != 1)

# Clean the "databases" ####
bd2 <- bd %>%
  select(Sample_id, CD_endoscopic_remission, CD_endoscopic_response, IBD,
         sample_location) %>%
  mutate(Sample_id = tolower(Sample_id)) %>%
  filter(str_detect(Sample_id, "-w")) %>%
  mutate(Sample_id = gsub(" reseq|rep", "", Sample_id),
         Pacient_id = as.character(as.numeric(gsub("-w.*", "", Sample_id))),
         Time = str_extract(Sample_id, "w[0-9]*"),
         Sample_id = paste0(Pacient_id, "-", Time)) %>%
  filter(Time %in% c("w0", "w14"))

upa_codes <- filter(upa_codes, !is.na(`PLACA MICRONIC`))

# Merge them together ####

## antiTNF ####
df <- merge(bd2, antiTNF, all.x = FALSE, all.y = TRUE) %>%
  mutate(remission = case_when(Time == "w0" ~ "w0",
                               Time == "w14" & CD_endoscopic_remission == "yes" ~ "w14 remiters",
                               Time == "w14" & CD_endoscopic_remission == "no" ~ "w14 non-remiters"))

# Some errors when writting the name on the reals!
filter(df, Location != sample_location) %>%
  group_by(Sample_id) %>%
  select(Location, sample_location, `Sample Name`)
# Should be null!!
# Samples 24-w14, 46-w14, and 73-w0 are know to be faulty

## UPA ####
duplic <- group_by(bd_upa, BarCode) %>%
        summarise(BarCodes = n()) %>%
        filter(BarCodes != 1) %>%
        pull(BarCode)
# Because we have duplicates we remove some data (first filtering for what do we care
bd_upa <- distinct(select(bd_upa, Week, pSES.CD, BarCode, ContainerName, SubjectID))

df_upa <- merge(UPA, upa_codes,  by.x = "Id", by.y = "RT TUBE")
stopifnot(nrow(df_upa) == nrow(UPA))
# right_join(bd_upa, df_upa  by = c("BarCode" = "nº muestra")
df_upa <- right_join(bd_upa, df_upa,
            by = c("BarCode" = "nº muestra")) %>%
  mutate(ContainerName = gsub("Biopsy ", "", ContainerName),
         Location = ContainerName,
         Week = paste0("w", Week),
         pSES.CD = as.numeric(pSES.CD),
         remission = case_when(
           Week == "w0" ~ "w0",
           Week != "w0" & pSES.CD < 5 ~ "w14 remiters",
           Week != "w0" & pSES.CD >= 5 ~ "w14 non-remiters",
           TRUE ~ Week))
stopifnot(nrow(df_upa) == nrow(UPA))

coln <- c("Sample", "Patient", "Time", "remission", "Target", "AU", "Location")
antiTNF_df <- select(df, Sample_id, Pacient_id, Time, remission,
                     "Target Name", AU, sample_location)
colnames(antiTNF_df) <- coln
upa_df <- select(df_upa, BarCode, SubjectID, Week, remission,
                     "Target Name", AU, Location)
colnames(upa_df) <- coln

controls$Time <- "C"
controls$Remission <- "C"
controls_df <- controls[, c("Sample Name", "Patient_id", "Time",
                            "Remission", "Target Name", "AU", "Location")]
colnames(controls_df) <- coln

dff <- rbind(cbind(upa_df, "Study" = "UPA"),
             cbind(antiTNF_df, "Study" = "TNF"),
             cbind(controls_df, "Study" = "C")) %>%
  mutate(Location = tolower(Location))

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
write_xlsx(dff, "processed/AU_markers.xlsx")

## Compare against controls ####
l <- vector("list", length =  length(unique(dff$Target))*2*2*2)
i <- 1
for (gene in unique(dff$Target)) {
  for (site in unique(dff$Location)) {
    for (study in c("UPA", "TNF")) {
      for (rem in c("w0", "w14 non-remiters")){
        d <- filter(dff,
                    Target == gene,
                    Location == site,
                    Study %in% c(study, "C"),
                    remission %in% c(rem, "C"))
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

vsC %>% group_by(Location, remission, Study) %>% summarise(n = n_distinct(Target))
filter(vsC, p.value == fdr)
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

between_studies %>% group_by(Location, remission) %>% summarise(n = n_distinct(Target))
filter(between_studies, p.value == fdr)
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
filter(within_studies, p.value == fdr)

ggplot(within_studies) +
  geom_abline(slope = 1, intercept = 0, col = "grey") +
  geom_point(aes(p.value, fdr, col = Target)) +
  theme_bw()
write_xlsx(within_studies, "processed/compare_whitin_studies.xlsx")

# Plots ####

today <- format(Sys.time(), "%Y%m%d")
my_comparisons <- list(c("C C", "w0 TNF"),
                       c("C C", "w0 UPA"),
                       c("w0 UPA", "w0 TNF"),
                       c("w14 remiters UPA", "w0 UPA"),
                       c("w14 remiters TNF", "w0 TNF"))

aTNF = "#a8ddb5"
aUPA = "#43a2ca"
cols = c("C" = "grey", "TNF" = aTNF, "UPA" = aUPA)

preplot <- dff %>%
  group_by(remission, Target, Location, Study) %>%
  summarise(meanAU = mean(AU), sem = sd(AU)/sqrt(n())) %>%
  mutate(ymax = meanAU + sem, ymin = meanAU - sem)


# Calculate the lin between the studies (problems when they cross??)
d <- preplot %>%
  filter(Study != "C", !grepl("non-remiters", remission)) %>%
  group_by(Target, remission, Location) %>%
  mutate(distU = (max(ymin) - min(ymax))/2,
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

preplot %>%
  filter(Study != "C", !grepl("non-remiters", remission)) %>%
  ggplot(aes(remission, meanAU)) +
  geom_point(aes(col = Study, group = Study)) +
  expand_limits(y = 0) +
  geom_errorbar(aes(col = Study, group = Study, ymin = ymin, ymax = ymax), width = 0.2) +
  geom_line(aes(col = Study, group = Study)) +
  geom_hline(data = filter(preplot, Study == "C"),
             aes(yintercept = ymin), linetype = "dotted") +
  geom_hline(data = filter(preplot, Study == "C"),
             aes(yintercept = ymax), linetype = "dotted") +
  geom_segment(data = dw, aes(x = remission, y = orig, xend = remission, yend = final)) +
  geom_text(data = dw, aes(x = remission, y = center), label = "\t*", size = 5) + # Dirty trick to dodge the symbol
  geom_text(data = db, aes(x = remission, y = meanAU), label = "\t+", size = 3) + # Dirty trick to dodge the symbol
  facet_grid(Target ~ Location, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  geom_segment(aes(xend = remission, yend = meanAU, group = remission)) +
  ylab("AU (mean\u00B1SEM)")

