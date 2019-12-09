# Table for reporting ####
library("readxl")
library("dplyr")
library("stringr")
library("lubridate")
library("xlsx")

# Read tables
reals <- read_xlsx("processed/AU_markers.xlsx")
bd <- read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = c("n.a.", ""))


reals %>%
  select(-Target, -AU) %>%
  distinct() %>%
  as.data.frame() %>%
  write.xlsx(file = "processed/Samples.xlsx", row.names = FALSE)


# Extract the ids of the samples of the reals for the antiTNF
names_patients <- reals %>%
  filter(Study == "TNF") %>%
  pull(Patient) %>%
  unique() %>%
  sort()

access <- read_excel("data/20190702_consulta_BCN_bbddd_4.xlsx", na = c("N/A", "")) %>%
  mutate(Patient_id = as.character(as.numeric(`Patients.ID patient`)),
         Sample = paste0(as.numeric(Patient_id), "-w", as.numeric(Week))) %>%
  filter(Patient_id %in% names_patients) %>%
  arrange(as.numeric(Patient_id))
n <- sapply(access, function(x){sum(is.na(x))})
access <- access[, n != nrow(access)]

access <- access %>%
  filter(Sample %in% reals$Sample) %>%
  distinct()

# There are some dates of inclusions not know
# We take the date of the first visit as the inclusion date (likely it was as such)
access$`Inclusion date`[access$Patient_id == "169"] <- as.POSIXct("2017-01-10", tz = "UTC")
access$`Inclusion date`[access$Patient_id == "172"] <- as.POSIXct("2017-02-16", tz = "UTC")

access <- mutate(access,
                 Disease_duration = interval(ymd(`Birth Date`), ymd(`Inclusion date`))/years(1))


final <- merge(reals, access) %>%
  select(-AU, -Target) %>%
  distinct() %>%
  mutate(group = if_else(Time == "w14", remission, "w0")) %>%
  mutate(Loc = case_when(
    Location == "ileum" ~ "ileum",
    TRUE ~ "colon")) %>%
  distinct(Sample, .keep_all = TRUE)

final$`CD phenotype` <- factor(final$`CD phenotype`,
                               levels = c("Inflammatory", "Stricturing", "Penetrating",
                                          "Stricturing/Penetrating"))

write.csv(final, "processed/metadata.csv", row.names = FALSE)

final %>%
  filter(group == "w0") %>%
  count(Loc)

final %>%
  filter(group == "w0") %>%
  group_by(Gender, group, Loc) %>%
  count() %>%
  ungroup() %>%
  select(-group) %>%
  arrange(Loc)

final %>%
  filter(group == "w0") %>%
  group_by(`Age at diagnosis`, group, Loc) %>%
  count() %>%
  ungroup() %>%
  select(-group) %>%
  arrange(Loc, desc(`Age at diagnosis`))

final %>%
  filter(group == "w0") %>%
  group_by(Loc, group) %>%
  summarise(mean = mean(`Disease_duration`, na.rm = TRUE),
            median = median(`Disease_duration`, na.rm = TRUE),
            min = min(`Disease_duration`, na.rm = TRUE),
            max = max(`Disease_duration`, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-group) %>%
  mutate_if(is.numeric, round)

final %>%
  filter(group == "w0") %>%
  group_by(`CD location`, group, Loc) %>%
  count() %>%
  ungroup() %>%
  select(-group) %>%
  arrange(Loc, desc(`CD location`))

final %>%
  filter(group == "w0") %>%
  group_by(`CD phenotype`, group, Loc) %>%
  count() %>%
  ungroup() %>%
  select(-group) %>%
  arrange(Loc)

final %>%
  filter(group == "w0") %>%
  group_by(`CD perianal disease`, group, Loc) %>%
  count() %>%
  ungroup() %>%
  select(-group) %>%
  arrange(Loc)

final %>%
  filter(group == "w0") %>%
  group_by(Loc, group) %>%
  summarise(mean = mean(`SES-CD (global)`, na.rm = TRUE),
            median = median(`SES-CD (global)`, na.rm = TRUE),
            min = min(`SES-CD (global)`, na.rm = TRUE),
            max = max(`SES-CD (global)`, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-group)

final %>%
  filter(group == "w0") %>%
  group_by(Loc, group) %>%
  summarise(mean = mean(`CDAI CD`, na.rm = TRUE),
            median = median(`CDAI CD`, na.rm = TRUE),
            min = min(`CDAI CD`, na.rm = TRUE),
            max = max(`CDAI CD`, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-group)

final %>%
  filter(group == "w0") %>%
  group_by(Loc, group) %>%
  summarise(mean = mean(`CDEIS CD`, na.rm = TRUE),
            median = median(`CDEIS CD`, na.rm = TRUE),
            min = min(`CDEIS CD`, na.rm = TRUE),
            max = max(`CDEIS CD`, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-group)

final %>%
  filter(group == "w0") %>%
  group_by(Loc, group) %>%
  summarise(mean = mean(`CRP`, na.rm = TRUE),
          median = median(`CRP`, na.rm = TRUE),
          min = min(`CRP`, na.rm = TRUE),
          max = max(`CRP`, na.rm = TRUE)) %>%
  ungroup() %>%
  select(-group)

# Treatment drugs
drugs <- read_xlsx("data/20190704_Treatment.xlsx")
f <- drugs %>%
  filter(as.character(as.numeric(`ID Patient`)) %in%  reals$Patient,
         grepl("00$|14$", Visit)) %>%
  mutate(ID_patient = as.character(as.numeric(`ID Patient`))) %>%
  arrange(ID_patient, Drug) %>%
  filter(!is.na(Drug)) %>%
  select(Drug, Visit, ID_patient) %>%
  distinct() %>%
  group_by(Visit) %>%
  summarise(Treatment = paste(Drug, collapse = ", "),
            ID_patient = unique(ID_patient)) %>%
  mutate(Time = str_extract(Visit, "-w(.*$)"),
         Time = paste0("w", as.numeric(str_extract(Time, "[0-9]+")))) %>%
  rename(Patient = ID_patient)



loc_tnf <- reals %>%
  filter(Study == "TNF") %>%
  select(-AU, -Target) %>%
  distinct()

f2 <- merge(f, loc_tnf) %>%
  mutate(group = if_else(Time == "w0", "w0", remission)) %>%
  group_by(group, General_location, Treatment) %>%
  count() %>%
  arrange(General_location, group)
filter(f2, group == "w0") %>%
  write_xlsx(path = "processed/treatment_TNF.xlsx")


# Ask azu how to put this table inside the other? Classify just AZA+IFX/corticoides and others ?

segments <- read_xlsx("data/20190702_consulta_BCN_bbddd_4.xlsx")
seg2 <- segments %>%
  mutate(Patient = as.character(as.numeric(`Patients.ID patient`))) %>%
  rename(Location = Segmento) %>%
  mutate(Location = case_when(
    Location == "ascendente" ~ "ascending",
    Location == "descendente" ~ "descending",
    Location == "íleon" ~ "ileum",
    grepl("rect", Location) ~ "rectum",
    Location == "sigma" ~ "sigmoid",
    Location == "transverso" ~ "transverse"
    ),
    Sample = paste0(Patient, "-w", as.numeric(Week)))

reals2 <- select(reals, -AU, -Target) %>% distinct()
seg3 <- merge(seg2, reals2)
clean <- !sapply(seg3, function(x){all(is.na(x))})
seg3 <- seg3[, clean] %>%
  mutate(Location = if_else(grepl("sigm*", Location), "sigm", Location))

d <- sapply(unique(seg3$Location), grep, colnames(seg3), value = TRUE)
num <- sapply(seg3, is.numeric)
seg3[, num][seg3[, num] > 1000] <- NA
out <- seg3 %>%
  filter(Week %in% c("000", "014")) %>%
  mutate(group = if_else(Time == "w0", "w0", remission)) %>%
  group_by(group, General_location, Location) %>%
  select(Sample, unlist(d, use.names = FALSE)) %>%
  distinct(Sample, Location, .keep_all = TRUE) %>%
  ungroup()

# Check that each sample has just one time
stopifnot(!any(rowSums(table(seg3$Sample, seg3$Week)) != 1))

scores <- function(x){
  y <- x[grep(x["Location"], names(x))]
  if (length(y) == 2) {
  } else {
    y <- c(NA_integer_, y)
  }
  names(y) <- c("CDEIS", "SES-CD")
  y
}

m <- apply(out, 1, scores)
m <- apply(t(m), 2, as.numeric)
sem <- function(x){sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))}
out1 <- aggregate(m, out[, 1:2], mean, na.rm = TRUE)
out2 <- aggregate(m, out[, 1:2], sd, na.rm = TRUE)
out3 <- aggregate(m, out[, 1:2], sem)
colnames(out1)[3:4] <- paste(colnames(out1)[3:4], "mean")
colnames(out2)[3:4] <- paste(colnames(out2)[3:4], "sd")
colnames(out3)[3:4] <- paste(colnames(out3)[3:4], "sem")

# Global mean and sd
apply(m[out$group == "w0", ], 2, mean, na.rm = TRUE)
apply(m[out$group == "w0", ], 2, sd, na.rm = TRUE)


out <- merge(out1, out2)
out <- merge(out, out3)



writexl::write_xlsx(out, "processed/score_segments_TNF.xlsx")

# Table UPA ####

bd_upa <- read_xlsx("processed/M13-740_abbvie_database_210619.xlsx", na = c("n.a.", "NA"))
upa2 <- right_join(bd_upa, reals[reals$Study == "UPA", ],
            by = c("BarCode" = "Sample", "SubjectID" = "Patient",
                   "cd" = "General_location")) %>%
  mutate(group = if_else(Time == "w0", "w0", remission.y)) %>%
  select(remission.x, remission.y, BarCode, pSES.CD, cd, group) %>%
  distinct()

upa2$pSES.CD <- as.numeric(upa2$pSES.CD)
# More NA than expected.
meanSESCD <- aggregate(pSES.CD~group+cd, data = upa2, mean, na.rm = TRUE)
semSESCD <- aggregate(pSES.CD~group+cd, data = upa2, sem)
sdSESCD <- aggregate(pSES.CD~group+cd, data = upa2, sd, na.rm = TRUE)
colnames(meanSESCD)[3] <- paste(colnames(meanSESCD)[3], "mean")
colnames(semSESCD)[3] <- paste(colnames(semSESCD)[3], "sem")
colnames(sdSESCD)[3] <- paste(colnames(sdSESCD)[3], "sd")
m <- merge(meanSESCD, semSESCD)
m <- merge(m, sdSESCD)
writexl::write_xlsx(m, "processed/score_segments_UPA.xlsx")
