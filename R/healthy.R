##########################
## Real times Azu (UPA) ##
##########################

library(xlsx)
library(gdata)
library(ggplot2)
library(ggpubr)
library(plyr)
library(ggthemes)
library(Rmisc, quietly = TRUE)
library(gridExtra)
library(grid)


### 25/06/2019####
# Mostres "healthy", que finalment resulta que no totes ho son.

# Dades
dades <- readxl::read_xls("data/Healthy biopsies_20190625_105706_Results_Export.xls")
clinica <- readxl::read_xlsx("processed/M13-740_abbvie_database_210619.xlsx")
RT <- readxl::read_xlsx("data/RT UPA_healthy.xlsx")
clinica <- clinica[, c(
  "Week", "pSES.CD", "BarCode", "ContainerName",
  "SubjectID", "ulcers", "Biopsy_Location"
)]

colnames(dades) <- c(
  "Experiment", "Well", "Sample_name", "Target.Name",
  "CT", "CT_Mean", "DcT_Mean", "DcT_SE", "Threshold",
  "Auto.Baseline", "Baseline.Start", "Baseline.End"
)

dades[grep(pattern = "Und", dades$CT), "CT"] <- "NA"
dades$CTN <- as.numeric(dades$CT)

aa <- strsplit(dades$Sample_name, split = "H ") # splits miriam code and barcode
df <- data.frame(matrix(unlist(aa), nrow = dim(dades)[1], byrow = T)) #
dades$BarCode <- df$X2
inside <- clinica[clinica$BarCode %in% df$X2, ]
write_xlsx(inside, "processed/RT_healthy_metadata.xlsx")
rm(df)

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
    message(barcode " PRESENT IN TWO EXPERIMENTS")
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
  x = taula, file = "processed/Dades_mostres_healthy.xlsx",
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
pdf("processed/Reals_pSESCD_0.pdf")
plotcolon
plotileum

dev.off()
