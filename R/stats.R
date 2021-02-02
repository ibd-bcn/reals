##########################
## Real times Azu (UPA) ##
##########################

library(gdata)
library("xlsx")
library(ggplot2)
library(ggpubr)
library(plyr)
library(ggthemes)
library(Rmisc, quietly =T)
library(gridExtra)
library(grid)

### 26/06/2019
# l'Azu em demana que faci l'analisi de les reals avui que en lluis no hi és (ell les ha fet
# la resta de dies). De fet ell m'ha deixat uns arxius a la O:/ que ara copiaré per comoditat al
# meu ordinador.
# la idea es comparar W0 vs W14R i W0 vs W14NR a UPA i TNF per separat i
# W0 UPA vs W0 TNF vs W0 Ctrl.

### 27/06/2019
# L'Azu ha tret tres valors. Repeteixo analisis.

### 27/06/2019####

# Dades
dades <- read.xlsx('processed/AU_markers.xlsx', sheetIndex = 1)

# Remove outliers
orig <- nrow(dades)
dades <- subset(dades, !(Target == "OSM" & Study == "TNF" &
                          Time == "w0" & AU > 35 & General_location == "colon"))
dades <- subset(dades, !(Target == "IL17A" & Study == "UPA" &
                          Time == "w0" & AU > 15 & General_location == "colon"))
dades <- subset(dades, !(Target == "OSM" & Study == "UPA" &
                          Time == "w0" & AU > 70 & General_location == "ileum"))
dades <- subset(dades, !(Target == "APQ8" & Time == "w0" &
                          AU > 130 & General_location == "ileum"))
dades <- subset(dades, !(Target == "RTNLB" & Time == "w0" & Study == "TNF" &
                          AU > 50 & General_location == "ileum"))
# In total it should be 6
stopifnot(orig-nrow(dades) == 6)

dades <- droplevels(dades)

clinica <- read.xlsx('processed/M13-740_abbvie_database_210619.xlsx', sheetIndex = 1)

clinica <- clinica[,c('Week', 'pSES.CD', 'BarCode', 'ContainerName',
                      'SubjectID', 'ulcers', 'Biopsy_Location')]

dades$myvar <- paste(dades$Time, dades$remission, sep='_')
dades$mytime <- mapvalues(dades$myvar, from = levels(as.factor(dades$myvar)),
                          to = c('Ctrl', 'w0', 'w0', 'w14_no', 'w14_yes'))
dades$myvar <- mapvalues(dades$myvar, from = levels(as.factor(dades$myvar)),
                         to = c('Ctrl', "w0_no", "w0_yes", 'w14_no', 'w14_yes'))
dades$Time <- mapvalues(dades$Time, from = 'C', to= 'Ctrl')
dades$Study <- mapvalues(dades$Study, from = 'C', to= 'Ctrl')
dades$summary <- paste(dades$Study, dades$General_location, dades$Target, dades$mytime, sep = '_')

dades$Target <- as.factor(dades$Target)


### dades eliminades ####

dades[dades$Target == 'OSM' & dades$Time == 'w0' &
        dades$General_location == 'ileum' & dades$AU > 71,'AU'] <- NA
dades[dades$Target == 'OSM' & dades$Time == 'w0' &
        dades$General_location == 'colon' & dades$Study == 'TNF' &
        dades$AU >38,'AU'] <- NA
dades[dades$Target == 'IL17A' & dades$Time == 'w0' &
        dades$General_location == 'colon' & dades$Study == 'UPA' &
        dades$AU > 16,'AU'] <- NA

dades[is.na(dades$AU),]
#         Sample Patient Time remission Target AU   Location Study General_location  myvar mytime
# 57   AB7379303   13407   w0        no    OSM NA      ileum   UPA            ileum  w0_no     w0
# 740  Y45576402   14101   w0        no  IL17A NA descending   UPA            colon  w0_no     w0
# 2520     21-w0      21   w0       yes    OSM NA     rectum   TNF            colon w0_yes     w0
#                 summary
# 57     UPA_ileum_OSM_w0
# 740  UPA_colon_IL17A_w0
# 2520   TNF_colon_OSM_w0


###### estadistiques #####


##### upadacitinib ####
statistics_upa_no <- NULL
statistics_upa_yes <- NULL

## ileum w0 vs w14 ####
for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'ileum',]
  mydata_upa <- mydata[mydata$Study == 'UPA',]

  a <- wilcox.test(mydata_upa$AU[mydata_upa$mytime == 'w0'],
                   mydata_upa$AU[mydata_upa$mytime == 'w14_yes'])
  df <- data.frame('Study' = 'UPA',
                   'Target' = i,
                   'General Location' = 'Ileum',
                   'Comparativa' = 'w0 vs w14 resmission',
                   'p.value' = a$p.value)

  statistics_upa_no <- rbind(statistics_upa_no, df)

  b <- wilcox.test(mydata_upa$AU[mydata_upa$mytime == 'w0'],
                   mydata_upa$AU[mydata_upa$mytime == 'w14_no'])
  df <- data.frame('Study' = 'UPA',
                   'Target' = i,
                   'Comparativa' = 'w0 vs w14 no resmission',
                   'General Location' = 'Ileum',
                   'p.value' = b$p.value)

  statistics_upa_yes <- rbind(statistics_upa_yes, df)


}

statistics_upa_no$fdr <- p.adjust(statistics_upa_no$p.value, method = 'fdr')
statistics_upa_yes$fdr <- p.adjust(statistics_upa_yes$p.value, method = 'fdr')
statistics_upa <- rbind(statistics_upa_no, statistics_upa_yes)

# colon w0 vs w14 ####
statistics_upa_no <- NULL
statistics_upa_yes <- NULL

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'colon',]
  mydata_upa <- mydata[mydata$Study == 'UPA',]

  a <- wilcox.test(mydata_upa$AU[mydata_upa$mytime == 'w0'],
                   mydata_upa$AU[mydata_upa$mytime == 'w14_yes'])
  df <- data.frame('Study' = 'UPA',
                   'Target' = i,
                   'General Location' = 'Colon',
                   'Comparativa' = 'w0 vs w14 resmission',
                   'p.value' = a$p.value)

  statistics_upa_no <- rbind(statistics_upa_no, df)

  b <- wilcox.test(mydata_upa$AU[mydata_upa$mytime == 'w0'],
                   mydata_upa$AU[mydata_upa$mytime == 'w14_no'])
  df <- data.frame('Study' = 'UPA',
                   'Target' = i,
                   'Comparativa' = 'w0 vs w14 no resmission',
                   'General Location' = 'Colon',
                   'p.value' = b$p.value)

  statistics_upa_yes <- rbind(statistics_upa_yes, df)
}

statistics_upa_no$fdr <- p.adjust(statistics_upa_no$p.value, method = 'fdr')
statistics_upa_yes$fdr <- p.adjust(statistics_upa_yes$p.value, method = 'fdr')

statistics_upa <- rbind(statistics_upa, statistics_upa_no)
statistics_upa <- rbind(statistics_upa, statistics_upa_yes)

write.xlsx(statistics_upa, file='processed/Statistics.xlsx', col.names = T, row.names = F,
           showNA = T, sheetName = 'Upadacitinib')


## w0: R vs NR ####

statistics_upa_no <- NULL
statistics_upa_yes <- NULL

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'ileum',]
  mydata_upa <- mydata[mydata$Study == 'UPA',]

  a <- wilcox.test(mydata_upa$AU[mydata_upa$myvar == 'w0_yes'],
                   mydata_upa$AU[mydata_upa$myvar == 'w0_no'])
  df <- data.frame('Study' = 'UPA',
                   'Target' = i,
                   'General Location' = 'Ileum',
                   'Comparativa' = 'w0: R vs NR',
                   'p.value' = a$p.value)

  statistics_upa_no <- rbind(statistics_upa_no, df)

  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'colon',]

  b <- wilcox.test(mydata_upa$AU[mydata_upa$myvar == 'w0_yes'],
                   mydata_upa$AU[mydata_upa$myvar == 'w0_no'])
  df <- data.frame('Study' = 'UPA',
                   'Target' = i,
                   'Comparativa' = 'w0: R vs NR',
                   'General Location' = 'Colon',
                   'p.value' = b$p.value)

  statistics_upa_yes <- rbind(statistics_upa_yes, df)

}

statistics_upa_no$fdr <- p.adjust(statistics_upa_no$p.value, method = 'fdr')
statistics_upa_yes$fdr <- p.adjust(statistics_upa_yes$p.value, method = 'fdr')
statistics_upa <- rbind(statistics_upa_no, statistics_upa_yes)

write.xlsx(statistics_upa, file='processed/Statistics.xlsx', col.names = T, row.names = F,
           showNA = T,  append = TRUE, sheetName = 'Upadacitinib w0- R vs NR')

##### TNF ####
statistics_tnf_no <- NULL
statistics_tnf_yes <- NULL

# ileum w0 vs w14 R ####
for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'ileum',]
  mydata_tnf <- mydata[mydata$Study == 'TNF',]

  a <- wilcox.test(mydata_tnf$AU[mydata_tnf$mytime == 'w0'],
                   mydata_tnf$AU[mydata_tnf$mytime == 'w14_yes'])
  df <- data.frame('Study' = 'TNF',
                   'Target' = i,
                   'General Location' = 'Ileum',
                   'Comparativa' = 'w0 vs w14 resmission',
                   'p.value' = a$p.value)

  statistics_tnf_yes <- rbind(statistics_tnf_yes, df)

  b <- wilcox.test(mydata_tnf$AU[mydata_tnf$mytime == 'w0'],
                   mydata_tnf$AU[mydata_tnf$mytime == 'w14_no'])
  df <- data.frame('Study' = 'TNF',
                   'Target' = i,
                   'Comparativa' = 'w0 vs w14 no resmission',
                   'General Location' = 'Ileum',
                   'p.value' = b$p.value)

  statistics_tnf_no <- rbind(statistics_tnf_no, df)


}

statistics_tnf_no$fdr <- p.adjust(statistics_tnf_no$p.value, method = 'fdr')
statistics_tnf_yes$fdr <- p.adjust(statistics_tnf_yes$p.value, method = 'fdr')
statistics_tnf <- rbind(statistics_tnf_no, statistics_tnf_yes)

# colon w0 vs W14 ####
statistics_tnf_no <- NULL
statistics_tnf_yes <- NULL

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'colon',]
  mydata_tnf <- mydata[mydata$Study == 'TNF',]

  a <- wilcox.test(mydata_tnf$AU[mydata_tnf$mytime == 'w0'],
                   mydata_tnf$AU[mydata_tnf$mytime == 'w14_yes'])
  df <- data.frame('Study' = 'TNF',
                   'Target' = i,
                   'General Location' = 'Colon',
                   'Comparativa' = 'w0 vs w14 resmission',
                   'p.value' = a$p.value)

  statistics_tnf_yes <- rbind(statistics_tnf_yes, df)

  b <- wilcox.test(mydata_tnf$AU[mydata_tnf$mytime == 'w0'],
                   mydata_tnf$AU[mydata_tnf$mytime == 'w14_no'])
  df <- data.frame('Study' = 'TNF',
                   'Target' = i,
                   'Comparativa' = 'w0 vs w14 no resmission',
                   'General Location' = 'Colon',
                   'p.value' = b$p.value)

  statistics_tnf_no <- rbind(statistics_tnf_no, df)
}

statistics_tnf_no$fdr <- p.adjust(statistics_tnf_no$p.value, method = 'fdr')
statistics_tnf_yes$fdr <- p.adjust(statistics_tnf_yes$p.value, method = 'fdr')

statistics_tnf <- rbind(statistics_tnf, statistics_tnf_no)
statistics_tnf <- rbind(statistics_tnf, statistics_tnf_yes)

write.xlsx(statistics_tnf, file='processed/Statistics.xlsx', col.names = T, row.names = F,
           showNA = T, append = T, sheetName = 'anti-TNF')

rm(statistics_tnf, statistics_tnf_no, statistics_tnf_yes, a, b, df)

## w0: R vs NR ####
statistics_tnf_no <- NULL
statistics_tnf_yes <- NULL

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'ileum',]
  mydata_tnf <- mydata[mydata$Study == 'TNF',]

  a <- wilcox.test(mydata_tnf$AU[mydata_tnf$myvar == 'w0_yes'],
                   mydata_tnf$AU[mydata_tnf$myvar == 'w0_no'])
  df <- data.frame('Study' = 'TNF',
                   'Target' = i,
                   'General Location' = 'Ileum',
                   'Comparativa' = 'w0: R vs NR',
                   'p.value' = a$p.value)

  statistics_tnf_no <- rbind(statistics_tnf_no, df)

  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'colon',]

  b <- wilcox.test(mydata_tnf$AU[mydata_tnf$myvar == 'w0_yes'],
                   mydata_tnf$AU[mydata_tnf$myvar == 'w0_no'])
  df <- data.frame('Study' = 'TNF',
                   'Target' = i,
                   'Comparativa' = 'w0: R vs NR',
                   'General Location' = 'Colon',
                   'p.value' = b$p.value)

  statistics_tnf_yes <- rbind(statistics_tnf_yes, df)

}

statistics_tnf_no$fdr <- p.adjust(statistics_tnf_no$p.value, method = 'fdr')
statistics_tnf_yes$fdr <- p.adjust(statistics_tnf_yes$p.value, method = 'fdr')
statistics_tnf <- rbind(statistics_tnf_no, statistics_tnf_yes)

write.xlsx(statistics_tnf, file='processed/Statistics.xlsx', col.names = T, row.names = F,
           showNA = T, append = TRUE, sheetName = 'anti-TNF w0- R vs NR')

##### weeks 0 #####

statistics_0_TNF <- NULL
statistics_0_UPA <- NULL
statistics_TNF_UPA <- NULL

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'ileum',]
  mydata <- mydata[mydata$Time == 'w0' | mydata$Time == 'Ctrl' ,]

  a <- wilcox.test(mydata$AU[mydata$Study == 'UPA'],
                   mydata$AU[mydata$Study == 'Ctrl'])
  df <- data.frame('Target' = i,
                   'General Location' = 'Ileum',
                   'Comparativa' = 'Ctrl vs UPA',
                   'p.value' = a$p.value)

  statistics_0_UPA <- rbind(statistics_0_UPA, df)

  b <- wilcox.test(mydata$AU[mydata$Study == 'TNF'],
                   mydata$AU[mydata$Study == 'Ctrl'])
  df <- data.frame('Target' = i,
                   'Comparativa' = 'Ctrl vs TNF',
                   'General Location' = 'Ileum',
                   'p.value' = b$p.value)

  statistics_0_TNF <- rbind(statistics_0_TNF, df)

  d <- wilcox.test(mydata$AU[mydata$Study == 'TNF'],
                   mydata$AU[mydata$Study == 'UPA'])
  df <- data.frame('Target' = i,
                   'Comparativa' = 'UPA vs TNF',
                   'General Location' = 'Ileum',
                   'p.value' = d$p.value)

  statistics_TNF_UPA <- rbind(statistics_TNF_UPA, df)


}

statistics_TNF_UPA$fdr <- p.adjust(statistics_TNF_UPA$p.value, method = 'fdr')
statistics_0_TNF$fdr <- p.adjust(statistics_0_TNF$p.value, method = 'fdr')
statistics_0_UPA$fdr <- p.adjust(statistics_0_UPA$p.value, method = 'fdr')

statistics_togh <- rbind(statistics_TNF_UPA, statistics_0_TNF, statistics_0_UPA)

statistics_0_TNF <- NULL
statistics_0_UPA <- NULL
statistics_TNF_UPA <- NULL

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'colon',]
  mydata <- mydata[mydata$Time == 'w0' | mydata$Time == 'Ctrl' ,]

  a <- wilcox.test(mydata$AU[mydata$Study == 'UPA'],
                   mydata$AU[mydata$Study == 'Ctrl'])
  df <- data.frame('Target' = i,
                   'General Location' = 'Colon',
                   'Comparativa' = 'Ctrl vs UPA',
                   'p.value' = a$p.value)

  statistics_0_UPA <- rbind(statistics_0_UPA, df)

  b <- wilcox.test(mydata$AU[mydata$Study == 'TNF'],
                   mydata$AU[mydata$Study == 'Ctrl'])
  df <- data.frame('Target' = i,
                   'Comparativa' = 'Ctrl vs TNF',
                   'General Location' = 'Colon',
                   'p.value' = b$p.value)

  statistics_0_TNF <- rbind(statistics_0_TNF, df)

  d <- wilcox.test(mydata$AU[mydata$Study == 'TNF'],
                   mydata$AU[mydata$Study == 'UPA'])
  df <- data.frame('Target' = i,
                   'Comparativa' = 'UPA vs TNF',
                   'General Location' = 'Colon',
                   'p.value' = d$p.value)

  statistics_TNF_UPA <- rbind(statistics_TNF_UPA, df)


}

statistics_TNF_UPA$fdr <- p.adjust(statistics_TNF_UPA$p.value, method = 'fdr')
statistics_0_TNF$fdr <- p.adjust(statistics_0_TNF$p.value, method = 'fdr')
statistics_0_UPA$fdr <- p.adjust(statistics_0_UPA$p.value, method = 'fdr')

statistics_togh <- rbind(statistics_togh, statistics_TNF_UPA, statistics_0_TNF,
                         statistics_0_UPA)

write.xlsx(statistics_togh, file='processed/Statistics.xlsx', col.names = T, row.names = F,
           showNA = T, sheetName = 'Weeks 0', append = T)

rm(statistics_togh, statistics_0_UPA, statistics_0_TNF, statistics_TNF_UPA, a, b, df, d,
   mydata, mydata_tnf, mydata_upa, i)


#####
##### plots #####
#####

pdf('Figures/Reals_Ileum.pdf')

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'ileum',]
  mydata_upa <- mydata[mydata$Study == 'UPA',]
  mydata_tnf <- mydata[mydata$Study == 'TNF',]
  mydata_0 <- mydata[mydata$Time == 'w0' | mydata$Time == 'Ctrl',]

  b <- ggplot(mydata_upa, aes(x=mytime, y=AU)) +
    ylab(i) +
    xlab('')+
    labs(title = paste('Ileum - UPA samples', i, sep = ' - '), y = 'AU') +
    geom_boxplot(position=position_dodge(),
                 outlier.colour = 'white',
                 fill = c(alpha('#336b87', 0.5), rep(alpha('#C0C0C0',0.4),1), rep(alpha('#708090',0.4),1)))+
    geom_point(shape = 21, size=2,position = position_jitterdodge(jitter.width = 0.11),aes(group = remission, fill = remission))+
    scale_fill_manual(values=c('#de7a22', '#6ab187'),
                      labels = c('No', 'Yes'),
                      name = 'Remission')+
    theme_classic() +
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size = 11), legend.position='none',
          axis.text.y = element_text(size = 11))


  labb <- quantile(mydata_upa[!is.na(mydata_upa[,'AU']),'AU'],0.99)
  labb2 <- labb + 0.1*max(mydata_upa[!is.na(mydata_upa[,'AU']),'AU'])

  upa <- b + stat_compare_means(method = "wilcox.test",
                                comparisons = list(c('w0', 'w14_no'),
                                                   c('w0', 'w14_yes')),
                                label="p.format",
                                label.y = c(labb, labb2),
                                tip.length = 0.01 )

  b <- ggplot(mydata_tnf, aes(x=mytime, y=AU)) +
    ylab(i) +
    xlab('')+
    labs(title = paste('Ileum - TNF samples', i, sep = ' - '), y = 'AU') +
    geom_boxplot(position=position_dodge(),
                 outlier.colour = 'white',
                 fill = c(alpha('#336b87', 0.5), rep(alpha('#C0C0C0',0.4),1), rep(alpha('#708090',0.4),1)))+
    geom_point(shape = 21, size=2,position = position_jitterdodge(jitter.width = 0.11),aes(group = remission, fill = remission))+
    scale_fill_manual(values=c('#de7a22', '#6ab187'),
                      labels = c('No', 'Yes'),
                      name = 'Remission')+
    theme_classic() +
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size = 11), legend.position='none',
          axis.text.y = element_text(size = 11))


  labb <- quantile(mydata_tnf[!is.na(mydata_tnf[,'AU']),'AU'],0.99)
  labb2 <- labb + 0.1*max(mydata_tnf[!is.na(mydata_tnf[,'AU']),'AU'])

  tnf <- b + stat_compare_means(method = "wilcox.test",
                                comparisons = list(c('w0', 'w14_no'),
                                                   c('w0', 'w14_yes')),
                                label="p.format",
                                label.y = c(labb, labb2),
                                tip.length = 0.01 )


  b <- ggplot(mydata_0, aes(x=Study, y=AU)) +
    ylab(i) +
    xlab('')+
    labs(title = paste('Ileum - W0 samples', i, sep = ' - '), y = 'AU') +
    geom_boxplot(position=position_dodge(),
                 outlier.colour = 'white',
                 fill = c(alpha('#336b87', 0.5), rep(alpha('#C0C0C0',0.4),1), rep(alpha('#708090',0.4),1)))+
    geom_point(shape = 21, size=2,position = position_jitterdodge(jitter.width = 0.11),aes(group = remission, fill = remission))+
    scale_fill_manual(values=c('#de7a22', '#6ab187'),
                      labels = c('No', 'Yes'),
                      name = 'Remission')+
    theme_classic() +
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size = 11), legend.position='none',
          axis.text.y = element_text(size = 11))


  labb <- quantile(mydata_0[!is.na(mydata_0[,'AU']),'AU'],0.99)
  labb2 <- labb + 0.1*max(mydata_0[!is.na(mydata_0[,'AU']),'AU'])

  zeros <- b + stat_compare_means(method = "wilcox.test",
                                  comparisons = list(c('Ctrl', 'TNF'),
                                                     c('Ctrl', 'UPA')),
                                  label="p.format",
                                  label.y = c(labb, labb2),
                                  tip.length = 0.01 )

  b <- ggplot(mydata_upa, aes(x=mytime, y=AU)) +
    ylab(i) +
    xlab('')+
    labs(title = paste('Ileum - UPA samples', i, sep = ' - '), y = 'AU') +
    geom_boxplot(position=position_dodge(),
                 outlier.colour = 'white',
                 fill = c(alpha('#336b87', 0.5), rep(alpha('#C0C0C0',0.4),1), rep(alpha('#708090',0.4),1)))+
    geom_point(shape = 21, size=2,position = position_jitterdodge(jitter.width = 0.11),aes(group = remission, fill = remission))+
    scale_fill_manual(values=c('#de7a22', '#6ab187'),
                      labels = c('No', 'Yes'),
                      name = 'Remission')+
    theme_classic() +
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size = 11), legend.position='left',
          axis.text.y = element_text(size = 11))

  legend <- cowplot::get_legend(b)



  grid.arrange(upa, tnf, zeros, legend,
               nrow = 2, ncol=2)

}
dev.off()


pdf('Figures/Reals_Colon.pdf')

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'colon',]
  mydata_upa <- mydata[mydata$Study == 'UPA',]
  mydata_tnf <- mydata[mydata$Study == 'TNF',]
  mydata_0 <- mydata[mydata$Time == 'w0' | mydata$Time == 'Ctrl',]

  b <- ggplot(mydata_upa, aes(x=mytime, y=AU)) +
    ylab(i) +
    xlab('')+
    labs(title = paste('Colon - UPA samples', i, sep = ' - '), y = 'AU') +
    geom_boxplot(position=position_dodge(),
                 outlier.colour = 'white',
                 fill = c(alpha('#336b87', 0.5), rep(alpha('#C0C0C0',0.4),1), rep(alpha('#708090',0.4),1)))+
    geom_point(shape = 21, size=2,position = position_jitterdodge(jitter.width = 0.11),aes(group = remission, fill = remission))+
    scale_fill_manual(values=c('#de7a22', '#6ab187'),
                      labels = c('No', 'Yes'),
                      name = 'Remission')+
    theme_classic() +
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size = 11), legend.position='none',
          axis.text.y = element_text(size = 11))


  labb <- quantile(mydata_upa[!is.na(mydata_upa[,'AU']),'AU'],0.99)
  labb2 <- labb + 0.1*max(mydata_upa[!is.na(mydata_upa[,'AU']),'AU'])

  upa <- b + stat_compare_means(method = "wilcox.test",
                                comparisons = list(c('w0', 'w14_no'),
                                                   c('w0', 'w14_yes')),
                                label="p.format",
                                label.y = c(labb, labb2),
                                tip.length = 0.01 )

  b <- ggplot(mydata_tnf, aes(x=mytime, y=AU)) +
    ylab(i) +
    xlab('')+
    labs(title = paste('Colon - TNF samples', i, sep = ' - '), y = 'AU') +
    geom_boxplot(position=position_dodge(),
                 outlier.colour = 'white',
                 fill = c(alpha('#336b87', 0.5), rep(alpha('#C0C0C0',0.4),1), rep(alpha('#708090',0.4),1)))+
    geom_point(shape = 21, size=2,position = position_jitterdodge(jitter.width = 0.11),aes(group = remission, fill = remission))+
    scale_fill_manual(values=c('#de7a22', '#6ab187'),
                      labels = c('No', 'Yes'),
                      name = 'Remission')+
    theme_classic() +
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size = 11), legend.position='none',
          axis.text.y = element_text(size = 11))


  labb <- quantile(mydata_tnf[!is.na(mydata_tnf[,'AU']),'AU'],0.99)
  labb2 <- labb + 0.1*max(mydata_tnf[!is.na(mydata_tnf[,'AU']),'AU'])

  tnf <- b + stat_compare_means(method = "wilcox.test",
                                comparisons = list(c('w0', 'w14_no'),
                                                   c('w0', 'w14_yes')),
                                label="p.format",
                                label.y = c(labb, labb2),
                                tip.length = 0.01 )


  b <- ggplot(mydata_0, aes(x=Study, y=AU)) +
    ylab(i) +
    xlab('')+
    labs(title = paste('Colon - W0 samples', i, sep = ' - '), y = 'AU') +
    geom_boxplot(position=position_dodge(),
                 outlier.colour = 'white',
                 fill = c(alpha('#336b87', 0.5), rep(alpha('#C0C0C0',0.4),1), rep(alpha('#708090',0.4),1)))+
    geom_point(shape = 21, size=2,position = position_jitterdodge(jitter.width = 0.11),aes(group = remission, fill = remission))+
    scale_fill_manual(values=c('#de7a22', '#6ab187'),
                      labels = c('No', 'Yes'),
                      name = 'Remission')+
    theme_classic() +
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size = 11), legend.position='none',
          axis.text.y = element_text(size = 11))


  labb <- quantile(mydata_0[!is.na(mydata_0[,'AU']),'AU'],0.99)
  labb2 <- labb + 0.1*max(mydata_0[!is.na(mydata_0[,'AU']),'AU'])

  zeros <- b + stat_compare_means(method = "wilcox.test",
                                  comparisons = list(c('Ctrl', 'TNF'),
                                                     c('Ctrl', 'UPA')),
                                  label="p.format",
                                  label.y = c(labb, labb2),
                                  tip.length = 0.01 )

  b <- ggplot(mydata_upa, aes(x=mytime, y=AU)) +
    ylab(i) +
    xlab('')+
    labs(title = paste('Colon - UPA samples', i, sep = ' - '), y = 'AU') +
    geom_boxplot(position=position_dodge(),
                 outlier.colour = 'white',
                 fill = c(alpha('#336b87', 0.5), rep(alpha('#C0C0C0',0.4),1), rep(alpha('#708090',0.4),1)))+
    geom_point(shape = 21, size=2,position = position_jitterdodge(jitter.width = 0.11),aes(group = remission, fill = remission))+
    scale_fill_manual(values=c('#de7a22', '#6ab187'),
                      labels = c('No', 'Yes'),
                      name = 'Remission')+
    theme_classic() +
    theme(axis.text.x  = element_text(angle=0, vjust=0.5, size = 11), legend.position='left',
          axis.text.y = element_text(size = 11))

  legend <- cowplot::get_legend(b)



  grid.arrange(upa, tnf, zeros, legend,
               nrow = 2, ncol = 2)

}
dev.off()





## Compare against controls ####
l <- vector("list", length =  length(unique(dff$Target))*2*2*2)
i <- 1
for (gene in unique(dff$Target)) {
  for (site in unique(dff$Location)) {
    for (study in c("UPA", "TNF")) {
      for (rem in c("yes", "no")) {
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
        i <- i + 1
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

preplot <- dades %>%
  group_by(remission, Target, Location, Study) %>%
  summarise(meanAU = mean(AU), sem = sd(AU, na.rm = TRUE)/sqrt(n())) %>%
  mutate(ymax = meanAU + sem, ymin = meanAU - sem)


bars_controls <- dades %>%
  filter(Study == "Ctrl") %>%
  group_by(Target, General_location) %>%
  summarise(meanAU = mean(AU), sem = sd(AU, na.rm = TRUE)/sqrt(n())) %>%
  mutate(ymax = meanAU + sem, ymin = meanAU - sem) %>%
  ungroup() %>%
  as.data.frame()

write.xlsx(bars_controls, "processed/controls_maxmin.xlsx", row.names = FALSE)

dw <- merge(between_studies, d) %>%
  filter(fdr < 0.05)

ws <- within_studies %>%
  mutate(remission = "w14 remiters")

db <- preplot %>%
  filter(Study != "C", !grepl("non-remiters", remission)) %>%
  merge(ws) %>%
  unique() %>%
  filter(fdr < 0.05)



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
