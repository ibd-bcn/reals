##########################
## Real times Azu (UPA) ##
##########################

library(xlsx)
library(gdata)
library(ggplot2)
library(ggpubr)
library(plyr)
library(ggthemes)
library(Rmisc, quietly =T)
library(gridExtra)
library(grid)


###25/06/2019####
# Mostres "healthy", que finalment resulta que no totes ho son.

# Dades
dades <- readxl::read_xls('data/Healthy biopsies_20190625_105706_Results_Export.xls')
clinica <- readxl::read_xlsx('processed/M13-740_abbvie_database_210619.xlsx')
RT <- readxl::read_xlsx('data/RT UPA_healthy.xlsx')
clinica <- clinica[,c('Week', 'pSES.CD', 'BarCode', 'ContainerName',
                      'SubjectID', 'ulcers', 'Biopsy_Location')]

colnames(dades)  <- c('Experiment', 'Well', 'Sample_name', 'Target.Name',
                      'CT', 'CT_Mean', 'DcT_Mean','DcT_SE', "Threshold",
                      "Auto.Baseline",   "Baseline.Start",  "Baseline.End")

dades[grep(pattern = 'Und', dades$CT),'CT'] <- 'NA'
dades$CTN <- as.numeric(dades$CT)

aa<-strsplit(dades$Sample_name, split='H ') # splits miriam code and barcode
df<-data.frame(matrix(unlist(aa), nrow=dim(dades)[1], byrow=T)) #
dades$BarCode <- df$X2
rm(df)

dades$BarCode <- as.factor(dades$BarCode)
taula <- NULL

for(i in 1:length(levels(dades$BarCode))){
  if(levels(dades$BarCode)[i] %in% RT$Sample.ID){print("Present both in Azucena and Miram files")  }
  if(!(levels(dades$BarCode)[i] %in% RT$Sample.ID)){print('RED ALERT! NO RT')}
  if(levels(dades$BarCode)[i] %in% clinica$BarCode){    print("Present both in Azucena and Llu?s files")  }
  if(!(levels(dades$BarCode)[i] %in% clinica$BarCode)){print('RED ALERT! - NO CLINICA')}

  if(levels(dades$BarCode)[i] == 'AC7980403'){print(i)}
  data <- dades[dades$BarCode == levels(dades$BarCode)[i],]
  if(nrow(data) == 6){
    data <- data[!is.na(data$CTN),]
    bact <- mean(data$CTN[data$Target.Name == 'BETA ACTINA'])
    targ <- mean(data$CTN[data$Target.Name == 'KLRC2'])

    df <- data.frame('sample' = levels(dades$BarCode)[i],
                     'week' = unique(clinica$Week[clinica$BarCode == levels(dades$BarCode)[i]]),
                     'pSES.CD' = unique(clinica$pSES.CD[clinica$BarCode == levels(dades$BarCode)[i]]),
                     'ContainerName' = unique(clinica$ContainerName[clinica$BarCode == levels(dades$BarCode)[i]]),
                     'SubjectID'= unique(clinica$SubjectID[clinica$BarCode == levels(dades$BarCode)[i]]),
                     'ulcers'= unique(clinica$ulcers[clinica$BarCode == levels(dades$BarCode)[i]]),
                     'Biopsy_Location'= unique(clinica$Biopsy_Location[clinica$BarCode == levels(dades$BarCode)[i]]),
                     'CT_mean_B_ACTINA' = bact,
                     'CT_mean_KLRC2' = targ,
                     'DcT' = targ-bact,
                     'AU' = 2^(bact-targ)*1000)
    taula <- rbind(taula,df)
  }
  if(nrow(data)>6){
    print('RED ALERT - PRESENT IN TWO EXPERIMENTS')
    for(k in levels(as.factor(data$Experiment))){
      datai <- data[data$Experiment == k,]
      datai <- datai[!is.na(datai$CTN),]
      bact <- mean(datai$CTN[datai$Target.Name == 'BETA ACTINA'])
      targ <- mean(datai$CTN[datai$Target.Name == 'KLRC2'])

      df <- data.frame('sample' = paste(levels(dades$BarCode)[i], k),
                       'week' = clinica$Week[clinica$BarCode == levels(dades$BarCode)[i]],
                       'pSES.CD' = clinica$pSES.CD[clinica$BarCode == levels(dades$BarCode)[i]],
                       'ContainerName' = clinica$ContainerName[clinica$BarCode == levels(dades$BarCode)[i]],
                       'SubjectID'= clinica$SubjectID[clinica$BarCode == levels(dades$BarCode)[i]],
                       'ulcers'= clinica$ulcers[clinica$BarCode == levels(dades$BarCode)[i]],
                       'Biopsy_Location'= clinica$Biopsy_Location[clinica$BarCode == levels(dades$BarCode)[i]],
                       'CT_mean_B_ACTINA' = bact,
                       'CT_mean_KLRC2' = targ,
                       'DcT' = targ-bact,
                       'AU' = 2^(bact-targ)*1000)
      taula <- rbind(taula,df)
    }
  }

}

## hi ha una mostra present a dos taules, suposo que per mirar tema batch.
## Ho parlo amb l'Azu, es tracta d'un error.

##
taula[grep('eds', taula$sample),]
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
dades_mostres <- dades[dades$BarCode %in% mostres,]

## Aquestes 6 mostres ja estan marcades en groc a l'excel que
## m'ha passat l'Azu, han de tenir el CT mes alt de tot el seu grup.

taula$CT_mean_KLRC2[is.nan(taula$CT_mean_KLRC2)] <-  max(taula$CT_mean_KLRC2[!is.nan(taula$CT_mean_KLRC2)])
taula$DcT <- taula$CT_mean_KLRC2 - taula$CT_mean_B_ACTINA
taula$AU <- 2^(-taula$DcT)*1000
taula$week[taula$week == '52'] <- '14'

write.xlsx(x = taula, file = 'O:/IBD_LAB/Personals/Ana/Azu/Dades_Mostres.xlsx',
           sheetName = '1', col.names = T, row.names = F)

##### plots ######

## colon
taula_colon <- taula[taula$ContainerName == 'Biopsy Colon',]
taula_colon <- taula_colon[taula_colon$pSES.CD == '0',]

plot <- ggplot(data = taula_colon, aes(x = week, y=AU,fill = Biopsy_Location))+
  geom_boxplot(outlier.colour = 'white',
               fill = 'white')+
  geom_point(shape = 21,aes(size=pSES.CD),  position = position_jitterdodge(jitter.width = 0.1))+
  labs(title = 'KLRC2 - Colon Samples')+
  theme_classic()

plotcolon <- plot + stat_compare_means(method = 'wilcox.test',
                                       comparisons = list(c('0','14')),
                                       label = 'p.format',
                                       tip.length = 0.01)
plotcolon

## ileum
taula_ileum <- taula[taula$ContainerName == 'Biopsy Ileum',]
taula_ileum <- taula_ileum[taula_ileum$pSES.CD == '0',]

plot <- ggplot(data = taula_ileum, aes(x = week, y=AU,fill = pSES.CD))+
  geom_boxplot(outlier.colour = 'white',
               fill = 'white')+
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.98))+
  scale_fill_manual(values = c('#ffebee', '#ff8a80','#d32f2f', '#c62828','#b71c1c', '#ffcdd2', '#ef9a9a','#e57373', '#ef5350','#ff1744','#f44336','#e53935','#d50000','gray'),
                    labels = levels(taula$pSES.CD),
                    name = 'Partial SES-CD')+
  labs(title = 'KLRC2 - Ileum Samples')+
  theme_classic()

plotileum <- plot + stat_compare_means(method = 'wilcox.test',
                                       comparisons = list(c('0','14')),
                                       label = 'p.format',
                                       tip.length = 0.01)
plotileum


# pdf #
pdf('processed/Reals_pSESCD_0.pdf')
grid.arrange(j)
plotcolon
plotileum

dev.off()
