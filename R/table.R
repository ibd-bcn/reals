# Table for reporting ####
library("readxl")
library("dplyr")
library("stringr")
library("lubridate")

# Read tables
reals <- read.csv("processed/reals.csv", check.names = FALSE)
bd <- read_xls("data/bd_BCN_tnf_biopsies_110119.xls", na = c("n.a.", ""))

# Extract the ids of the samples of the reals for the antiTNF
names_patients <- reals %>%
  filter(grepl("[0-9]W",`Sample Name`, ignore.case = TRUE)) %>%
  distinct(`Sample Name`) %>%
  pull(`Sample Name`) %>%
  tolower() %>%
  str_split("w") %>%
  sapply(getElement, name = 1) %>%
  unique()

access <- read_excel("data/20190522_consulta_BCN_bbddd_3.xlsx", na = c("N/A", "")) %>%
  mutate(Patient_id = as.character(as.numeric(`id Patient`))) %>%
  filter(Patient_id %in% names_patients) %>%
  arrange(as.numeric(Patient_id))
# There are some dates of inclusions not know
# We take the date of the first visit as the inclusion date (likely it was as such)
access$`Inclusion date`[access$Patient_id == "169"] <- as.POSIXct("2017-01-10", tz = "UTC")
access$`Inclusion date`[access$Patient_id == "172"] <- as.POSIXct("2017-02-16", tz = "UTC")

access <- mutate(access,
                 Disease_duration = interval(ymd(`Birth Date`), ymd(`Inclusion date`))/years(1))


group_by(access, `Age at diagnosis`) %>% count()
mutate(access,
       Loc = case_when(
         Segmento == "íleon" ~ "ileum",
         !is.na(Segmento) ~ "colon")) %>%
  group_by(Loc) %>%
  count()

dff %>%
  filter(Study == "TNF") %>%
  distinct(Location, Time, Patient) %>%
  group_by(Location) %>%
  count()

group_by(access, Gender) %>% count()
group_by(access, `CD perianal disease`) %>% count()
group_by(access, `CD location`) %>% count()
group_by(access, tolower(`CD phenotype`)) %>% count()
summarise(access,
          mean = mean(`CDAI CD`, na.rm = TRUE),
          min = min(`CDAI CD`, na.rm = TRUE),
          max = max(`CDAI CD`, na.rm = TRUE))
summarise(access,
          mean = mean(`CRP`, na.rm = TRUE),
          min = min(`CRP`, na.rm = TRUE),
          max = max(`CRP`, na.rm = TRUE))
summarise(access,
          mean = mean(`SES-CD (global)`, na.rm = TRUE),
          min = min(`SES-CD (global)`, na.rm = TRUE),
          max = max(`SES-CD (global)`, na.rm = TRUE))
summarise(access,
          mean = mean(`Disease_duration`, na.rm = TRUE),
          min = min(`Disease_duration`, na.rm = TRUE),
          max = max(`Disease_duration`, na.rm = TRUE))
summarise(access,
          mean = mean(`CDEIS CD`, na.rm = TRUE),
          min = min(`CDEIS CD`, na.rm = TRUE),
          max = max(`CDEIS CD`, na.rm = TRUE))
summarise(access,
          mean = mean(`Disease_duration`, na.rm = TRUE),
          min = min(`Disease_duration`, na.rm = TRUE),
          max = max(`Disease_duration`, na.rm = TRUE))
group_by(access, tolower(`CD phenotype`), `CD perianal disease`) %>% count()
