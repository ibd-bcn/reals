## Real times Azu (UPA) review 2021 ##

library("xlsx")
library("gdata", warn.conflicts = FALSE)
library("ggplot2", warn.conflicts = FALSE)
library("ggpubr")
library("tidyr", quietly = TRUE, warn.conflicts = FALSE)
library("ggthemes", warn.conflicts = FALSE)
library("Rmisc", quietly = TRUE, warn.conflicts = FALSE)
library("dplyr", quietly = TRUE, warn.conflicts = FALSE)
library("gridExtra", quietly = TRUE, warn.conflicts = FALSE)
library("grid")


# Dades
new_files <- list.files(pattern = "*.xls", path = "data/reals2021/", full.names = TRUE)
dades0 <- lapply(new_files, readxl::read_xls, skip = 7, n_max = 192)
for (i in seq_along(dades0)) {
  dades0[[i]] <- cbind(Experiment = basename(new_files[i]), dades0[[i]])
}

merger <- function(x, y) {
    merge(x, y, all = TRUE, sort = FALSE)
}
dades <- Reduce(merger, dades0)
dades <- dades %>%
  filter(nchar(Well) <= 3, !is.na(Well), # Filter the ending of the tables
         !is.na(`Sample Name`))  # filter actual samples and do not use blanks

clinica <- readxl::read_xlsx("processed/M13-740_abbvie_database_210619.xlsx")
RT <- readxl::read_xlsx("data/RT UPA_healthy.xlsx")
clinica <- clinica[, c(
  "Week", "pSES.CD", "BarCode", "ContainerName",
  "SubjectID", "ulcers", "Biopsy_Location"
)]

# colnames(dades) <- c(
#     "Experiment", "Well", "Sample_name", "Target.Name",
#     "CT", "CT_Mean", "DcT_Mean", "DcT_SE", "Threshold",
#     "Auto.Baseline", "Baseline.Start", "Baseline.End"
# )

colnames(dades)[colnames(dades) == "Cт"] <- "CT"
colnames(dades)[colnames(dades) == "Sample Name"] <- "Sample_name"
colnames(dades)[colnames(dades) == "Target Name"] <- "Target.Name"
dades$CT[grep(pattern = "Und", dades$CT)] <- NA
dades$CT <- as.numeric(dades$CT)
dades$CTN <- as.numeric(dades$CT)

# Samples ####
# Samples missing both genes in the same well
bad_well <- dades %>%
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
bad_genes <- dades %>%
  filter(is.na(CT)) %>%
  select(Experiment, Well, Sample_name, Target.Name) %>%
  group_by(Sample_name, Target.Name) %>%
  summarize(n = n(), wells = list(Well), Exp = unique(Experiment)) %>%
  filter(n != 1, !is.na(Sample_name)) %>%
  select(-n)
bad_genes
stopifnot(nrow(bad_genes) == 8)

# Fill the missing values of ADCYAP with 40 as it is a lowly expressed gene
dades$CT[dades$Target.Name == "ADCYAP1" &
           is.na(dades$CT) &
           dades$Sample_name %in% tail(bad_genes$Sample_name, 6)] <- 40

# Missing sample
# Merge miriam and barcode

# Look for the beta actina of the same well as the genes ####
nested_well <- dades %>%
  filter(Target.Name != "BETA ACTINA") %>%
  group_by(Experiment, Target.Name, Sample_name) %>%
  nest(data = Well)

nested_well$data <- lapply(nested_well$data, unlist)

# Look for the wells how many have values
xy <- vector("list", length = nrow(nested_well))
for (i in seq_len(nrow(nested_well))) {
  x <- nested_well[i, ]
  wells_value <- dades %>%
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
  x[x$`Target Name` != "BETA ACTINA", "nas", drop = TRUE]})
table(gene)
####

dades$BarCode <- as.factor(dades$BarCode)
taula <- NULL

for (i in 1:length(levels(dades$BarCode))) {
  barcode <- levels(dades$BarCode)[i]
  if (barcode %in% RT$Sample.ID) {
    message("Present both in Azucena and Miram files")
  }
  if (!(barcode %in% RT$Sample.ID)) {
    message("NO RT")
  }
  if (barcode %in% clinica$BarCode) {
    message("Present both in Azucena and Lluís files")
  }
  if (!(barcode %in% clinica$BarCode)) {
    message("NO CLINICA for ", barcode)
  }

  if (barcode == "AC7980403") {
    message(i)
  }
  data <- dades[dades$BarCode == barcode, ]
  if (nrow(data) == 6) {
    data <- data[!is.na(data$CTN), ]
    bact <- mean(data$CTN[data$Target.Name == "BETA ACTINA"])
    targ <- mean(data$CTN[data$Target.Name == "KLRC2"])

    df <- data.frame(
      "sample" = barcode,
      "week" = unique(clinica$Week[clinica$BarCode == barcode]),
      "pSES.CD" = unique(clinica$pSES.CD[clinica$BarCode == barcode]),
      "ContainerName" = unique(clinica$ContainerName[clinica$BarCode == barcode]),
      "SubjectID" = unique(clinica$SubjectID[clinica$BarCode == barcode]),
      "ulcers" = unique(clinica$ulcers[clinica$BarCode == barcode]),
      "Biopsy_Location" = unique(clinica$Biopsy_Location[clinica$BarCode == barcode]),
      "CT_mean_B_ACTINA" = bact,
      "CT_mean_KLRC2" = targ,
      "DcT" = targ - bact,
      "AU" = 2^(bact - targ) * 1000
    )
    taula <- rbind(taula, df)
  }
  if (nrow(data) > 6) {
    message(barcode, " PRESENT IN TWO EXPERIMENTS")
    for (k in levels(as.factor(data$Experiment))) {
      datai <- data[data$Experiment == k, ]
      datai <- datai[!is.na(datai$CTN), ]
      bact <- mean(datai$CTN[datai$Target.Name == "BETA ACTINA"])
      targ <- mean(datai$CTN[datai$Target.Name == "KLRC2"])

      df <- data.frame(
        "sample" = paste(barcode, k),
        "week" = clinica$Week[clinica$BarCode == barcode],
        "pSES.CD" = clinica$pSES.CD[clinica$BarCode == barcode],
        "ContainerName" = clinica$ContainerName[clinica$BarCode == barcode],
        "SubjectID" = clinica$SubjectID[clinica$BarCode == barcode],
        "ulcers" = clinica$ulcers[clinica$BarCode == barcode],
        "Biopsy_Location" = clinica$Biopsy_Location[clinica$BarCode == barcode],
        "CT_mean_B_ACTINA" = bact,
        "CT_mean_KLRC2" = targ,
        "DcT" = targ - bact,
        "AU" = 2^(bact - targ) * 1000
      )
      taula <- rbind(taula, df)
    }
  }
}

## hi ha una mostra present a dos taules, suposo que per mirar tema batch.
## Ho parlo amb l'Azu, es tracta d'un error.

##
taula[grep("eds", taula$sample), ]
#     sample                                  week pSES.CD ContainerName
# 94  Y45570803 210619 UPA HEALTHY PLAT 1.eds   14       4  Biopsy Ileum
# 95  Y45570803 210619 UPA HEALTHY PLAT 3.eds   14       4  Biopsy Ileum

# SubjectID    ulcers Biopsy_Location CT_mean_B_ACTINA CT_mean_KLRC2
# 94     12501    yes           ileum         20.53750        31.294
# 95     12501    yes           ileum         20.72767        27.831

# DcT        AU
# 94 10.756499 0.5780575
# 95  7.103333 7.2724987


## alunes mostres tenen "NaN", mirem-les.

mostres <- as.character(taula$sample[is.nan(taula$DcT)])
dades_mostres <- dades[dades$BarCode %in% mostres, ]

## Aquestes 6 mostres ja estan marcades en groc a l'excel que
## m'ha passat l'Azu, han de tenir el CT mes alt de tot el seu grup.

taula$CT_mean_KLRC2[is.nan(taula$CT_mean_KLRC2)] <- max(taula$CT_mean_KLRC2[!is.nan(taula$CT_mean_KLRC2)])
taula$DcT <- taula$CT_mean_KLRC2 - taula$CT_mean_B_ACTINA
taula$AU <- 2^(-taula$DcT) * 1000
taula$week[taula$week == "52"] <- "14"
#  Arxiu original a 'O:/IBD_LAB/Personals/Ana/Azu/Dades_Mostres.xlsx'
write.xlsx(
  x = taula, file = "processed/Dades_mostres_reals21.xlsx",
  sheetName = "1", col.names = T, row.names = F
)

##### plots ######

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
