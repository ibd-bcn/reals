library("xlsx")
library("readxl")
library("dplyr")
library("stringr")
library("ggplot2")
library("tidyr")
library("writexl")
library("ggpubr")
library("patchwork")

bd <- read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = c("n.a.", ""))
bd_upa <- read_xls("data/M13-740_abbvie_database_210918.xls")
reals <- read_xls("data/SOX6 CHI3L1 PDGFD PTGDR2_20190516_115511_Results_Export.xls", 
    skip = 15, na = c("", "Undetermined"))
upa_codes <- read_xlsx("data/RT UPA DISEASE BLOC 1 2 3 500NG 20UL.xlsx")

# Clean the original data ####

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
        distinct()

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

write_xlsx(erroneous, "processed/big_sd.xlsx")

pdf("Figures/histograms.pdf")
hist(reals$Cт[reals$`Target Name` == "BETA ACTINA"], main = "Beta Actina", xlab = "Ct")
hist(reals$Cт[reals$`Target Name` == "CHI3L1"], main = "CHI3L1", xlab = "Ct")
hist(reals$Cт[reals$`Target Name` == "PDGFD"], main = "PDGFD", xlab = "Ct")
hist(reals$Cт[reals$`Target Name` == "SOX6"], main = "SOX6", xlab = "Ct")
dev.off()


# Do you remember some samples that were in more than one sample?
# Now we remove a plate
reals <- filter(reals, 
            !(`Sample Name` == "73w14c" & grepl("290419bcn 102 a 108", `Experiment Name`)), 
            !(`Sample Name` == "10w14"  & grepl("290419bcn 102 a 108", `Experiment Name`)))

reals[reals$`Experiment Name` == "290419 bcn samples plat 3.eds" & 
    reals$`Sample Name` == "44w0c", "Sample Name"] <- "41w14"

# Check the results
multiple_exp <- reals %>%
    group_by(`Sample Name`) %>%
    summarise(n = n_distinct(`Experiment Name`)) 
m <-multiple_exp %>%
    group_by(n) %>%
    count()

stopifnot(m[m$n == 2, "nn"] == 1)


preclean <- reals %>%
  filter(`Target Name` != "BETA ACTINA") %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(`ΔCт` = unique(`ΔCт Mean`)) %>%
  ungroup() %>%
  filter(!is.na(`ΔCт`))

antiTNF <- preclean %>%
  filter(grepl("^[0-9]", `Sample Name`)) %>%
  mutate(`Sample Name` = tolower(`Sample Name`),
         Pacient_id = gsub("w.*", "", `Sample Name`),
         Time = str_extract(`Sample Name`, "w[0-9]*"),
         Location = case_when(grepl("i$", `Sample Name`) ~ "ileum",
                              grepl("c$", `Sample Name`) ~ "colon"),
         Sample_id = paste0(Pacient_id, "-", Time),
         AU = 2^(-`ΔCт`)*1000) %>%
  select(-`ΔCт`)

UPA <- preclean %>%
  filter(grepl("^UP", `Sample Name`)) %>%
  mutate(`Sample Name` = gsub("^UP ", "UP", `Sample Name`),         
         AU = 2^(-`ΔCт`)*1000) %>%
  separate(`Sample Name`, c("Id", "Pacient_id")) %>%
  select(-`ΔCт`)
controls <- preclean %>%
  filter(grepl("^C", `Sample Name`)) %>%
  mutate(AU = 2^(-`ΔCт`)*1000,
         Location = case_when(grepl(" C|S", `Sample Name`) ~ "colon",
                              grepl(" ILI$", `Sample Name`) ~ "ileum"),
         Patient_id = gsub(" .+$", "", `Sample Name`)) %>%
  select(-`ΔCт`)


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

# Plots ####
i_antiTNF <- filter(df, sample_location == "ileum") %>%
  ggplot(aes(remission, AU, col = remission,
             group = remission)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  facet_wrap(~`Target Name`, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(caption = "ileum")

today <- format(Sys.time(), "%Y%m%d")
ggsave(paste0("Figures/", today, "_antiTNF_ileum.png"), plot = i_antiTNF)

c_antiTNF <- filter(df, sample_location == "colon") %>%
  ggplot(aes(remission, AU, col = remission,
             group = remission)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  facet_wrap(~`Target Name`, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(caption = "colon")
ggsave(paste0("Figures/", today, "_antiTNF_colon.png"), plot = c_antiTNF)


df_remitter <- df_upa %>%
  group_by(BarCode) %>%
  summarise(remitter = if_else(any(pSES.CD < 5), "Y", "N"))

i_upa <- filter(df_upa, ContainerName == "Ileum") %>%
  ggplot(aes(remission, AU, col = remission,
             group = remission)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  facet_wrap(~`Target Name`, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(caption = "ileum")
print(i_upa)
c_upa <- filter(df_upa, ContainerName == "Colon") %>%
  ggplot(aes(remission, AU, col = remission,
             group = remission)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  facet_wrap(~`Target Name`, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(caption = "colon")
print(c_upa)

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

dodge <- position_dodge(width=0.9)
preplot %>%
  filter(Study != "C", !grepl("non-remiters", remission)) %>%
  ggplot(aes(remission, meanAU, col = Study, group = Study)) +
  geom_point() +
  expand_limits(y = 0) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2) +
  geom_line() +
  geom_hline(data = filter(preplot, Study == "C"), 
                aes(yintercept = ymin), linetype = "dotted") +
  geom_hline(data = filter(preplot, Study == "C"), 
                aes(yintercept = ymax), linetype = "dotted") +
  facet_wrap(Location ~ Target, scales = "free", ncol = 4) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  ylab("AU (mean\u00B1SEM)")

ile <- preplot %>%
  filter(Location == "ileum") %>%
  ggplot(aes(condition, AU, col = Study)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge()) + # https://stackoverflow.com/a/24019668/2886003
  facet_wrap( ~ Target, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  labs(title = "ileum") +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", tip.length = 0.05) +
  scale_color_manual(values = cols)

col + ile

# Table for reporting ####
bd3 <- bd %>%
  mutate(Sample_id = tolower(Sample_id)) %>%
  filter(str_detect(Sample_id, "-w")) %>%
  mutate(Sample_id = gsub(" reseq|rep", "", Sample_id),
         Pacient_id = as.character(as.numeric(gsub("-w.*", "", Sample_id))),
         Time = str_extract(Sample_id, "w[0-9]*"),
         Sample_id = paste0(Pacient_id, "-", Time)) %>%
  filter(Time %in% c("w0", "w14")) %>%
  as.data.frame() %>%
  select(Location, `Diagnostic age`, `Gender`, `DOB`, `CRP (mg/dL)`, `CDAI CD`,
         `Perianal disease`, `biopsied_segment`, `Sample_id`, `Treatment`,
         `Diagnostic age`, sample_date, CDEIS_partial)

df2 <- base::merge(bd3, reals2, all.x = FALSE, all.y = TRUE, by = "Sample_id")

# TODO Finish the tables
df2 %>%
  distinct(Pacient_id, .keep_all = TRUE) %>%
  NULL


