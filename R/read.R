library("xlsx")
library("readxl")
library("dplyr")
library("stringr")
library("ggplot2")

bd <- read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = c("n.a.", ""))
reals <- read_xls("data/test 2_20190502_100858_Results_Export.xls", range = "A16:H2284")

# Clean the original data
reals2 <- reals %>%
  filter(`Target Name` != "BETA ACTINA") %>%
  group_by(`Sample Name`, `Target Name`) %>%
  summarise(`ΔCт` = unique(`ΔCт Mean`)) %>%
  ungroup() %>%
  mutate(`Sample Name` = tolower(`Sample Name`),
         Pacient_id = gsub("w.*", "", `Sample Name`),
         Time = str_extract(`Sample Name`, "w[0-9]*"),
         Location = if_else(grepl("i$", `Sample Name`), "ileum", "colon"),
         Sample_id = paste0(Pacient_id, "-", Time),
         AU = 2^(-`ΔCт`)*1000) %>%
  select(-`ΔCт`)


# Clean the "database"
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

# Merge them together
df <- merge(bd2, reals2, all.x = FALSE, all.y = TRUE)

# Some errors
filter(df, Location != sample_location)

df$remission[df$Time == "w0"] <- "w0"
df$remission[df$Time == "w14"] <- ifelse(df$CD_endoscopic_remission[df$Time == "w14"] == "yes",
                                         "w14 remiters", "w14 non-remiters")

filter(df, sample_location == "ileum") %>%
  ggplot(aes(remission, AU, col = remission,
             group = remission)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  facet_wrap(~`Target Name`, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(caption = "ileum")

filter(df, sample_location == "colon") %>%
  ggplot(aes(remission, AU, col = remission,
             group = remission)) +
  geom_boxplot() +
  geom_jitter(width = 0.25) +
  facet_wrap(~`Target Name`, scales = "free", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  labs(caption = "colon")
