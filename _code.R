###################
###################
## PREPROCESSING ##
###################
###################

#################
# PRELIMINARIES #
#################

# SET WORKING DIRECTORY
wd <- "SELECT_WORKING_DIRECTORY"
setwd(wd)

# SET SUBDIRECTORY NAMES (you must have these subfolders)
analyses <- paste0(wd,'/analyses/')
imed <- paste0(wd,'/imed/')
rawData <- paste0(wd,'/rawData')
qa <- paste0(wd,'/qa')

# SOURCE FUNCTIONS AND LIBRARIES 
source("_formatPlates.R")
source("_functionsAcrossPlates.R")
source("_fAnalysis.R")

#######################
# PLATE PREPROCESSING #
#######################

load(file=paste0(imed,'dat.RData'))
if(!exists('dat')){

  # Raw plate reader files (all assumed to be in the rawData folder)
  plateReaderExports <- c(
    "201005_BY_GPH1-GFP_IRA2_Blocksandv1-8_Plate1A_BYtoRM_EXPORT.xlsx",
    "201006_BY_GPH1-GFP_IRA2_Blocksandv1-8_Plate2A_BYtoRM_EXPORT.xlsx",
    "201009_BY_GPH1-GFP_IRA2_Blocksandv1-8_Plate1B_BYtoRM_EXPORT.xlsx",
    "201010_BY_GPH1-GFP_IRA2_Blocksandv1-8_Plate2B_BYtoRM_EXPORT.xlsx",
    "201016_BY_GPH1-GFP_IRA2_Blocksandv1-8_Plate1E_BYtoRM_EXPORT.xlsx",
    "201017_BY_GPH1-GFP_IRA2_Blocksandv1-8_Plate2E_BYtoRM_EXPORT.xlsx",
    "200907_BY_GPH1-GFP_IRA2_allVariants_Plate1B6_BYtoRM_EXPORT.xlsx",
    "200910_BY_GPH1-GFP_IRA2_allVariants_Plate1B7_BYtoRM_EXPORT.xlsx",
    "200912_BY_GPH1-GFP_IRA2_allVariants_Plate2B6_BYtoRM_EXPORT.xlsx",
    "200914_BY_GPH1-GFP_IRA2_allVariants_Plate3B6_BYtoRM_EXPORT.xlsx",
    "200915_BY_GPH1-GFP_IRA2_allVariants_Plate1B8_BYtoRM_EXPORT.xlsx",
    "200916_BY_GPH1-GFP_IRA2_allVariants_Plate2B7_BYtoRM_EXPORT.xlsx",
    "200917_BY_GPH1-GFP_IRA2_allVariants_Plate3B7_BYtoRM_EXPORT.xlsx",
    "200918_BY_GPH1-GFP_IRA2_allVariants_Plate1B9_BYtoRM_EXPORT.xlsx",
    "200919_BY_GPH1-GFP_IRA2_allVariants_Plate2B8_BYtoRM_EXPORT.xlsx",
    "200920_BY_GPH1-GFP_IRA2_allVariants_Plate3B8_BYtoRM_EXPORT.xlsx",
    "200923_BY_GPH1-GFP_IRA2_allVariants_Plate2B9_BYtoRM_EXPORT.xlsx",
    "200924_BY_GPH1-GFP_IRA2_allVariants_Plate3B9_BYtoRM_EXPORT.xlsx",
    "200925_BY_GPH1-GFP_IRA2_allVariants_Plate1B10_BYtoRM_EXPORT.xlsx",
    "200926_BY_GPH1-GFP_IRA2_allVariants_Plate2B10_BYtoRM_EXPORT.xlsx",
    "200927_BY_GPH1-GFP_IRA2_allVariants_Plate3B10_BYtoRM_EXPORT.xlsx",
    "210124_BY_GPH1-GFP_IRA2_ALLBlocks_Random_2_EXPORT.xlsx",
    "210125_BY_GPH1-GFP_IRA2_ALLBlocks_Random_3_EXPORT.xlsx",
    "210126_BY_GPH1-GFP_IRA2_ALLBlocks_Random_4_EXPORT.xlsx",
    "210127_BY_GPH1-GFP_IRA2_ALLBlocks_Random_5_EXPORT.xlsx",
    "210130_BY_GPH1-GFP_IRA2_ALLBlocks_Random_6_EXPORT.xlsx",
    "210311_BYandRM_GPH1-GFP_Glucose_EXPORT.xlsx",
    "210312_BYandRM_GPH1-GFP_Glucose_EXPORT.xlsx",
    "210129_BY_GPH1-GFP_IRA2_Block1syn_Random_1_EXPORT.xlsx",
    "210131_BY_GPH1-GFP_IRA2_Block1syn_Random_2_EXPORT.xlsx"
  )
  
  # descriptions must also be in rawData and match order of export files above 
  plateReaderDescriptionFiles <- c(
    "201005_Blocksandv1-8_Plate1A_PlateData.xlsx",
    "201006_Blocksandv1-8_Plate2A_PlateData.xlsx",
    "201009_Blocksandv1-8_Plate1B_PlateData.xlsx",
    "201010_Blocksandv1-8_Plate2B_PlateData.xlsx",
    "201016_Blocksandv1-8_Plate1E_PlateData.xlsx",
    "201017_Blocksandv1-8_Plate2E_PlateData.xlsx",
    "200907_Plate1B6_IRA2_Allvariants_random_PlateData.xlsx",
    "200910_Plate1B7_IRA2_Allvariants_random_PlateData.xlsx",
    "200912_Plate2B6_IRA2_Allvariants_random_PlateData.xlsx",
    "200914_Plate3B6_IRA2_Allvariants_random_PlateData.xlsx",
    "200915_Plate1B8_IRA2_Allvariants_random_PlateData.xlsx",
    "200916_Plate2B7_IRA2_Allvariants_random_PlateData.xlsx",
    "200917_Plate3B7_IRA2_Allvariants_random_PlateData.xlsx",
    "200918_Plate1B9_IRA2_Allvariants_random_PlateData.xlsx",
    "200919_Plate2B8_IRA2_Allvariants_random_PlateData.xlsx",
    "200920_Plate3B8_IRA2_Allvariants_random_PlateData.xlsx",
    "200923_Plate2B9_IRA2_Allvariants_random_PlateData.xlsx",
    "200924_Plate3B9_IRA2_Allvariants_random_PlateData.xlsx",
    "200925_Plate1B10_IRA2_Allvariants_random_PlateData.xlsx",
    "200926_Plate2B10_IRA2_Allvariants_random_PlateData.xlsx",
    "200927_Plate3B10_IRA2_Allvariants_random_PlateData.xlsx",
    "210124_Allblock_old_Random_PlateData.xlsx",
    "210125_Allblock_old_Random_PlateData.xlsx",
    "210126_Allblock_old_Random_PlateData.xlsx",
    "210127_Allblock_old_Random_PlateData.xlsx",
    "210130_Allblock_old_Random_PlateData_modif.xlsx",
    "210311_BYandRM_GPH1-GFP_PlateData.xlsx",
    "210312_BYandRM_GPH1-GFP_PlateData.xlsx",
    "210129_IRA2_Block1_submap2_PlateData_modif.xlsx",
    "210131_IRA2_Block1_submap2_PlateData_modif.xlsx"
  )
  
  # bind together each export with its corresponding descriptor 
  inputData <- cbind(plateReaderExports, plateReaderDescriptionFiles)
  
  # get formatted plates 
  fPlates <- lapply(1:nrow(inputData), function(i){
    print(inputData[i, 1]) # print input data
    knitPlates(
      dataFile = inputData[i, 1],
      descriptionFile = inputData[i, 2],
      outTextFile=paste(strsplit(inputData[i, 1], "_")[[1]][1], "output.txt", sep="_"),
      outPlotFile=paste(strsplit(inputData[i, 1], "_")[[1]][1], "output.pdf", sep="_"),
      userOD=0.25, outputFolder=qa, inputFolder=rawData
    )
  })
  names(fPlates) <- sapply(plateReaderDescriptionFiles, function(x){strsplit(x, "_")[[1]][1]})
  
  # put them all into one table for stats
  dat <- c()
  for (plate in names(fPlates)){
    dat <- rbind(dat, cbind(fPlates[[plate]], plate))
  }
  
  # remove rows with blanks and rows with no data (i.e. empty wells)
  dat <- dat[is.na(dat$Blank),]
  colnames(dat)[colnames(dat) == "catDesc"] <- "clone"
  dat <- dat[dat$clone != "____",]
  
  # log-transform GFP readings, get simple Genotype, and compress genos
  dat <- dat %>% mutate_at(c("Ratios_Midlog", "Ratios_UserDef", 
                             "Ratios_Sat", "Ratios_TMid", 
                             "growthRate", "capacity"), log2) %>%
    mutate(Genotype = sapply(dat$clone,function(x)
      str_extract(string = x,pattern = "(?<=\\().*(?=\\))"))) %>%
    rename(Background = Desc_1) %>% 
    mutate(Background = factor(Background)) %>%
    mutate(Genotype = ifelse(test=str_detect(Genotype,"BY[1-4]"),
                         "BY", Genotype)) %>% 
    mutate(Genotype = ifelse(test=str_detect(Genotype,"RM[1-4]"),
                         "RM", Genotype)) %>%
    mutate(Genotype = ifelse(test=str_detect(Genotype,"BLK0.5"),
                             "halfB1", Genotype)) %>%
    drop_na(any_of('Ratios_TMid')) 
  
  # save the dat table 
  save(dat,file=paste0(imed,'dat.RData')) 
}

###################
###################
## BY VS RM ONLY ##
###################
###################

# select the relevant platEs
oDat <- filter(dat,plate %in% c('210311','210312')) 
  
# prepare summaries 
oSummary <- sapply(oModels,function(x)summary(x)$coefficients[2,]) %>% 
  t() %>% as.data.frame() %>% 
  mutate(`Pr(>F)`= (lapply(oModels,anova) %>% 
                      sapply(function(x)x$`Pr(>F)`))) %>%
  mutate(Genotype=factor(levels(oDat$Background))) %>%
  rename(Effect=Estimate) %>% rename(SE=`Std. Error`) %>%
  mutate(pStar = cut(`Pr(>F)`, 
                     breaks=c(-Inf,0.001,0.01,0.05,Inf),
                     labels=c("***","**","*","")))

write.csv(oSummary,paste0(analyses,file='BYRM_only_ANOVAs.csv'))

templot <- oneBars(oSummary,6)
ggsave(templot,
       filename = paste0(analyses,'_onlyBR.pdf'),width=3,height=2)

######################
######################
## BLOCK ANALYSIS ####
######################
######################


####################################
# SELECT PLATES AND PREPARE LEVELS #
####################################

# select relevant plates
bDat <- filter(dat,plate %in% c('210124','210125','210126','210127','210130')) 

# decide on the levels for you block genotypes
bLevels <-  permutations(n=2, r=4, v=c('B','R'), repeats.allowed=T) %>% 
  as.data.frame() %>%
  arrange(desc(V1),desc(V2),desc(V3),desc(V4)) %>%
  apply(1,paste,collapse = "") %>% 
  as.data.frame() %>%
  `colnames<-`(c('geno')) %>%
  mutate(counts = str_count(geno,"R")) %>% 
  arrange(counts)

# breakout genotype information and select genotype order
bDat <- bDat$Genotype %>% (function(x)sapply(x,function(y)strsplit(y,""))) %>% 
  (function(x)do.call(cbind,x)) %>% t %>% 
  `colnames<-`(c('b1','b2','b3','b4')) %>% 
  `rownames<-`(c()) %>% cbind(bDat) %>%
  mutate(Genotype = factor(Genotype,levels=bLevels$geno))

########################
# ONE WAY ANOVA METHOD #
########################

bAnovas <- oneMaker(bDat,bLevels$geno,'BBBB')
write.csv(bAnovas,paste0(analyses,file='blockAnovas.csv'))

templot <- oneBars(bAnovas,6)
ggsave(templot,filename=paste0(analyses,'bAnovaBar.pdf'),
       width=3,height=3.5)

############################
# PROPER EPISTASIS TESTING #
############################

bModel <- lmer(data=bDat, 
               Ratios_TMid ~ (1|plate) + (1|clone) + 
                 (b1 + b2 + b3 + b4)^4)


epistasisResults <- read.delim(file=paste0(analyses,'epistasisResults_table.csv'),
                   sep=',',row.names = 1)

if(!exists('epistasisResults')){
  
  SE <- summary(bModel)$coef[, 2, drop = FALSE]
  
  epistasisResults <- cbind(fixef(bModel),SE) %>% 
    as.data.frame() %>%
    (function(x)cbind(factor(rownames(SE),
                             levels=rev(as.factor(rownames(SE)))),x)) %>%
    `colnames<-`(c('Genotype','Effect','SE')) %>%
    (function(x)x[-1,]) %>% `rownames<-`(c()) %>%
    mutate(sumSq=anova(bModel,type='I')$'Sum Sq') %>%
    mutate(Fval=anova(bModel,type='I')$'F value') %>%
    mutate(pVal=anova(bModel,type='I')$`Pr(>F)`) %>%
    mutate(Epistasis = c('No', 'Yes')[1+str_detect(Genotype, ':')])
  
  bBoot <- confint.merMod(bModel,method=c("boot")) %>% as.data.frame()
  bBoot <- bBoot[rownames(bBoot) %in% m$Genotype,]
  colnames(bBoot) <- c('percent2.5','percent97.5')
  epistasisResults  <- cbind(epistasisResults,bBoot)
  
  epistasisResults %>% write.csv(file=paste0(analyses,'epistasisResults_table.csv'))
  
  rm(bBoot,SE)
  
}

p <- lapply(c('No','Yes'),function(x)
  ggplot(filter(epistasisResults ,Epistasis == x),
         aes(x=Genotype, y=Effect)) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_bar(stat="identity", width=0.5,fill='gray',fill=NA) +
    geom_hline(yintercept = 0, color = "black") +
    geom_errorbar(aes(ymin=Effect-SE,ymax=Effect+SE,width=0.3),color='black') +
    coord_flip()
)

mapply(function(p,n)
  ggsave(p,filename = n,width=6.5,height=4),p,
  c(paste0(analyses,'nonEpistasis_plot.pdf'),
    paste0(analyses,'Epistasis_plot.pdf')))


# test if an R in block 4 actually reduces expression
lmer(data=bDat[(bDat$Genotype %in% c('RRRB','RRRR')),], 
     Ratios_TMid ~ (1|plate) + (1|clone) + Genotype) %>%
  anova()

##################################
# BOXPLOTS FOR EPISTASIS ZOOM IN #
##################################

tbox <- bDat[c('b1','b2')] %>%
  mutate(pc=residuals(lmer(data=bDat, 
                           Ratios_TMid ~ (1|plate) + (1|clone))))

fbox <- sapply(names(table(tbox$b1)),function(x)
  sapply(names(table(tbox$b1)),function(y)
    median(tbox[(tbox$b1 == x & tbox$b2 == y),]$pc))) %>% 
  t() %>% as.data.frame() %>%
  mutate(b1 = names(table(tbox$b1))) %>% cbind(
    matrix(rep(names(table(tbox$b1)),
               length(table(tbox$b1))),
           nrow=length(table(tbox$b1)),
           ncol=length(table(tbox$b1))) %>% t %>%
      as.data.frame() %>% `names<-`(c(1:length(table(tbox$b1))))
    
  )

# could not get ggplot to place line segments correctly on x-plane
# sadly they and their color were hard-coded
templot <- bDat[c('b1','b2')] %>%
  mutate(pc=residuals(lmer(data=bDat, 
                         Ratios_TMid ~ (1|plate) + (1|clone))))  %>% 
  ggplot(aes(x=b2,y=pc,fill=b1,color=b1)) + 
  theme_bw() + 
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(size = 0.4, position = position_jitterdodge()) + 
  scale_fill_manual(values=c('#0000FF78','#FF000078'),name = "Block 1") +
  scale_color_manual(values=c('#0000FFFF','#FF0000FF'),name = "Block 1") +
  ylab(expression('Log'[2]*'(Gph1-GFP/OD) at inflection point')) + 
  xlab('Block 2') + 
  geom_segment(aes(x = 0.82, y = fbox$B[1],
                   xend = 1.82, yend = fbox$R[1]),color='blue') +
  geom_segment(aes(x = 1.18, y = fbox$B[2],
                   xend = 2.18, yend = fbox$R[2]),color='red') 

ggsave(paste0(analyses,'b1b2boxplot.pdf'),plot = templot,device = NULL,
       path = NULL,scale = 1,width = 3,height = 4,
       units = c("in"), useDingbats=F)

rm(tbox,fbox,templot)

#############################
#############################
## SINGLE VARIANT ANALYSIS ##
#############################
#############################
####################################
# SELECT PLATES AND PREPARE LEVELS #
####################################

vDat <- filter(dat,
               plate %in% c('201005','201006','201009','201010',
                            '201016','201017','200907','200910',
                            '200912','200914','200915','200916',
                            '200917','200918','200919','200920',
                            '200923','200924','200925','200926',
                            '200927'))

# get your factor order right
vLevels <-  c('BY','RM','RBBB','BRBB','BBRB','BBBR','halfB1',varEz(1:26))
deBlocks <- c('RBBB','BRBB','BBRB','BBBR','halfB1')

# correct factor order
vDat$Genotype <- factor(vDat$Genotype,levels=as.factor(vLevels))

##################
# ONE-WAY ANOVAS #
##################

# For all variants
vAnovas <- oneMaker(filter(vDat,Genotype %nin% deBlocks),
                    vLevels[vLevels %nin% deBlocks],'BY')

write.csv(vAnovas,paste0(analyses,file='variantAnovas.csv'))

templot <- oneBars(vAnovas,6)

ggsave(templot,filename = paste0(analyses,'vAnovaBar.pdf'),
       width=6.5,height=3.5)

######################
# VARIANT BOOTSTRAPS #
######################

# PARAMETERS
straps <- 10000

# get plot limits for a figure match 
lims <- summary(c(unlist(b1straps),unlist(hbstraps)))

# BLOCK 1 VARIANT TEST
load(paste0(imed,'b1straps.RData'))
if(!exists('b1straps')){
  b1straps<- varStrap(straps,vDat,'BY','RBBB',varEz(1:8)) 
  save(b1straps,file=paste0(imed,'b1straps.RData'))
}

strapSignif(b1straps)
templot <- bdPlot(b1straps,'Full Allele','Variants',
                  "purple","orange",xLo = lims['Min.'],xHi=lims['Max.'])
ggsave(paste0(analyses,'b1straps.pdf'),plot = templot,device = cairo_pdf,
       path = NULL,scale = 1,width = 3.5,height = 2.5,units = c("in"))

# determine the overlap
sum(b1straps$svSums >= quantile(b1straps$diffBYRM,0.025))/length(b1straps$svSums)

# HALF GENE VARIANT TEST
load(paste0(imed,'hbstraps.RData'))
if(!exists('hbstraps')){
  hbstraps <- varStrap(straps,vDat,'BY','halfB1',varEz(1:5))
  save(hbstraps,file=paste0(imed,'hbstraps.RData'))
}

strapSignif(hbstraps)
templot <- bdPlot(hbstraps,'Full Allele','Variants',
                  "purple","orange",xLo = lims['Min.'],xHi=lims['Max.'])
ggsave(paste0(analyses,'hbstraps.pdf'),plot = templot,device = cairo_pdf,
       path = NULL,scale = 1,width = 3.5,height = 2.5,units = c("in")) 

# determine the overlap
sum(hbstraps$svSums >= quantile(hbstraps$diffBYRM,0.025))/length(hbstraps$svSums)

# WHOLE GENE VARIANT TEST
load(paste0(imed,'wgstraps.RData'))
if(!exists('wgstraps')){
  wgstraps<- varStrap(straps,vDat,'BY','RM',varEz(1:26))
  save(wgstraps,file=paste0(imed,'wgstraps.RData'))
}

strapSignif(wgstraps)
templot <- bdPlot(wgstraps,'Full Allele','Variants',"purple","orange")
ggsave(paste0(analyses,'wgstraps.pdf'),plot = templot,device = cairo_pdf,
       path = NULL,scale = 1,width = 6.5,height = 5.5,units = c("in")) 

# determine the overlap
sum(wgstraps$svSums >= quantile(wgstraps$diffBYRM,0.025))/length(wgstraps$svSums)
sum(wgstraps$diffBYRM <= quantile(wgstraps$svSums,0.975))/length(wgstraps$diffBYRM)
# how would you determine central quantile overlap and get this to apply on other scales
# for example at other confidence intervals 


###############################################
# VARIANT BOOTSTRAPS USING ONLY SELECT PLATES #
###############################################

lims <- summary(c(unlist(alt_b1straps),unlist(alt_hbstraps)))

# BLOCK 1 VARIANT TEST
load(paste0(imed,'alt_b1straps.RData'))
if(!exists('alt_b1straps')){
  alt_b1straps<- varStrap(straps,zDat,'BY','RBBB',varEz(1:8)) 
  save(alt_b1straps,file=paste0(imed,'alt_b1straps.RData'))
}

strapSignif(alt_b1straps)
templot <- bdPlot(alt_b1straps,'Full Allele','Variants',
                  "purple","orange",xLo = lims['Min.'],xHi=lims['Max.'])
ggsave(paste0(analyses,'alt_b1straps.pdf'),plot = templot,device = cairo_pdf,
       path = NULL,scale = 1,width = 3.5,height = 2.5,units = c("in"))


# HALF BLOCK 1 VARIANT TEST
load(paste0(imed,'alt_hbstraps.RData'))
if(!exists('alt_hbstraps')){
  alt_hbstraps<- varStrap(straps,zDat,'BY','halfB1',varEz(1:5)) 
  save(alt_hbstraps,file=paste0(imed,'alt_hbstraps.RData'))
}

strapSignif(alt_hbstraps)
templot <- bdPlot(alt_hbstraps,'Full Allele','Variants',
                  "purple","orange",xLo = lims['Min.'],xHi=lims['Max.'])
ggsave(paste0(analyses,'alt_hbstraps.pdf'),plot = templot,device = cairo_pdf,
       path = NULL,scale = 1,width = 3.5,height = 2.5,units = c("in"))

#######################
#######################
## BLOCK 1 ZOOM IN ####
#######################
#######################

####################################
# SELECT PLATES AND PREPARE LEVELS #
####################################

zDat <- filter(dat,
               plate %in% c('201005','201006','201009',
                            '201010','201016','201017')) %>%
  filter(Genotype %nin% c('BRBB','BBRB','BBBR'))

# get your factor order right
zLevels <-  c('BY','RM','RBBB','halfB1',varEz(1:8))

# correct factor order
zDat$Genotype <- factor(zDat$Genotype,levels=as.factor(zLevels))

##################
# ONE-WAY ANOVAS #
##################

# For all variants
zAnovas <- oneMaker(vDat,zLevels,'BY')
write.csv(vAnovas,paste0(analyses,file='zoomAnovas.csv'))

templot <- oneBars(zAnovas,6)
ggsave(templot,filename = paste0(analyses,'zAnovaBar.pdf'),
       width=3,height=5)

##############################################
# ONE-WAY ANOVAS USING SELECT VARIANT PLATES #
##############################################

# For all variants
alt_zAnovas <- oneMaker(zDat,zLevels,'BY')
write.csv(alt_zAnovas,paste0(analyses,file='alt_zoomAnovas.csv'))

templot <- oneBars(alt_zAnovas,6)
ggsave(templot,filename = paste0(analyses,'alt_zAnovaBar.pdf'),
       width=3,height=5)

#############################
#############################
## NON-SYNONYMOUS ANALYSES ##
#############################
#############################

####################################
# SELECT PLATES AND FORMAT FACTORS #
####################################

nDat <- filter(dat,plate %in% c('210129','210131'))  %>% 
  mutate(Nonsynonymous = case_when(
    Genotype == "BY" ~ 'BY',
    Genotype == "RM" ~ 'RM',
    Genotype == "synRMnonsynBY" ~ 'BY',
    Genotype == "synBYnonsynRM" ~ 'RM')) %>%
  mutate(Synonymous = case_when(
    Genotype == "BY" ~ 'BY',
    Genotype == "RM" ~ 'RM',
    Genotype == "synRMnonsynBY" ~ 'RM',
    Genotype == "synBYnonsynRM" ~ 'BY'
  ))

#############################
# BUILD MIXED-EFFECTS MODEL #
#############################

nModel <- nInterxn <- lmer(data=nDat,
                 Ratios_TMid ~ (1|plate) + (1|clone) + 
                   (Nonsynonymous + Synonymous)^2) 

anova(nInterxn,type='I') %>% 
  write.csv(file=paste0(analyses,'nInterxn_Anova_typeI.csv'))

anova(nInterxn,type='III') %>% 
  write.csv(file=paste0(analyses,'nInterxn_Anova_typeIII.csv'))


