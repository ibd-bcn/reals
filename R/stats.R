##########################
## Real times Azu (UPA) ##
##########################

library(gdata)
library(ggplot2)
library(ggpubr)
library(plyr)
library(ggthemes)
library(Rmisc, quietly =T)
library(gridExtra)
library(grid)

###26/06/2019
# l'Azu em demana que faci l'analisi de les reals avui que en lluis no hi és (ell les ha fet
# la resta de dies). De fet ell m'ha deixat uns arxius a la O:/ que ara copiaré per comoditat al
# meu ordinador.
# la idea es comparar W0 vs W14R i W0 vs W14NR a UPA i TNF per separat i
# W0 UPA vs W0 TNF vs W0 Ctrl.

### 27/06/2019
# L'Azu ha tret tres valors. Repeteixo analisis.

###27/06/2019####

# Dades
dades <- readxl::read_xlsx('processed/AU_markers.xlsx')
clinica <- readxl::read_xlsx('processed//M13-740_abbvie_database_210619.xlsx')

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

dades[dades$Target == 'OSM' & dades$Time == 'w0' & dades$General_location == 'ileum' & dades$AU > 71,'AU'] <- NA

dades[dades$Target == 'OSM' & dades$Time == 'w0' & dades$General_location == 'colon' & dades$Study == 'TNF' & dades$AU >38,'AU'] <- NA

dades[dades$Target == 'IL17A' & dades$Time == 'w0' & dades$General_location == 'colon' & dades$Study == 'UPA' & dades$AU > 16,'AU'] <- NA

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

write.xlsx(statistics_upa, file='Statistics.xlsx', col.names = T, row.names = F,
           showNA = T, sheetName = 'Upadacitinib')

rm(statistics_upa, statistics_upa_no, statistics_upa_yes, a, b, df)

##### TNF ####
statistics_tnf_no <- NULL
statistics_tnf_yes <- NULL

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

write.xlsx(statistics_tnf, file='Statistics.xlsx', col.names = T, row.names = F,
           showNA = T, append = T, sheetName = 'anti-TNF')

rm(statistics_tnf, statistics_tnf_no, statistics_tnf_yes, a, b, df)

##### weeks 0 #####

statistics_0_TNF <- NULL
statistics_0_UPA <- NULL
statistics_TNF_UPA <- NULL

for(i in levels(dades$Target)){
  mydata <- dades[dades$Target == i,]
  mydata <- mydata[mydata$General_location == 'ileum',]
  mydata <- mydata[mydata$Time == 'w0' |mydata$Time == 'Ctrl' ,]

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
  mydata <- mydata[mydata$Time == 'w0' |mydata$Time == 'Ctrl' ,]

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

write.xlsx(statistics_togh, file='Statistics.xlsx', col.names = T, row.names = F,
           showNA = T, sheetName = 'Weeks 0', append = T)

rm(statistics_togh, statistics_0_UPA, statistics_0_TNF, statistics_TNF_UPA, a, b, df, d,
   mydata, mydata_tnf, mydata_upa, i)


#####
##### plots #####
#####

pdf('Documents/azucena_salas/Reals UPA/Reals_Ileum.pdf')

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


pdf('Documents/azucena_salas/Reals UPA/Reals_Colon.pdf')

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
               nrow = 2, ncol=2)

}
dev.off()
