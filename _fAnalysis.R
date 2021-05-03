#############
# FUNCTIONS #
#############

library(gtools)
library(splitstackshape) # for stratified samplings
library(stringr)
library(tidyverse)
library(merTools)
library(lmerTest)

bdPlot <- function(straps,lab1,lab2,col1,col2,lw=1,xLo=NULL,xHi=NULL){
  
  tidy <- straps %>% 
    `names<-`(c(lab1, lab2)) %>%
    gather(key='Distribution',value='Value') 
  
  
  p <- tidy %>% 
    ggplot(aes(x=Value,fill=Distribution)) + 
    geom_density(alpha=0.6) + theme_bw() + 
    scale_fill_manual(values=c(col1, col2)) +
    theme(legend.position=c(0.2,0.8),
          legend.title = element_blank(),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black")) +
    geom_vline(
      xintercept = unlist(
        quantile(
          unlist(
            filter(tidy,Distribution==lab1)['Value']),
          0.975)),
      color=col1,lwd=lw) +
    geom_vline(
      xintercept = unlist(
        quantile(
          unlist(
            filter(tidy,Distribution==lab1)['Value']),
          0.025)),
      color=col1,lwd=lw) +
    geom_vline(
      xintercept = unlist(
        quantile(
          unlist(
            filter(tidy,Distribution==lab2)['Value']),
          0.975)),
      color=col2,lwd=lw) +
    geom_vline(
      xintercept = unlist(
        quantile(
          unlist(
            filter(tidy,Distribution==lab2)['Value']),
          0.025)),
      color=col2,lwd=lw) +
    xlab('Allele Effect')
  
  if(length(xLo) > 0 & length(xHi) > 0){
    p <- p + xlim(xLo,xHi)
  }
  
  return(p)
  
}

# CREATES A BARPLOT FROM AN LMER OBJECT
lmerBar <- function(m){
  
  q <- cbind(fixef(m),
             summary(m)$coef[,2, drop = FALSE],
             coef(summary(m))[,"Pr(>|t|)"]) %>% 
    as.data.frame() %>%
    (function(x)cbind(factor(names(fixef(m)),
                             levels=rev(as.factor(
                               names(
                                 fixef(m)[-1])))),x)) %>%
    `colnames<-`(c('Genotype','Effect','SE','prT')) %>%
    (function(x)x[-1,]) %>% `rownames<-`(c()) %>%
    mutate(Epistasis = c('No', 'Yes')[1+str_detect(Genotype, ':')]) %>%
    mutate(pStar=cut(prT, breaks=c(-Inf,0.001,0.01,0.05,Inf), 
                     labels=c("***","**","*","")))
  
  p <- lapply(c('No','Yes'),function(x)
    ggplot(filter(q,Epistasis == x),aes(x=Genotype, y=Effect)) +
      theme_bw() +
      theme(axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_bar(stat="identity", width=0.9,fill='grey50',fill=NA) +
      geom_hline(yintercept = 0, color = "black") +
      geom_errorbar(aes(ymin=Effect-SE,ymax=Effect+SE,width=0.3)) +
      geom_text(aes(x=Genotype, y=(Effect+((SE+0.02)*sign(Effect))), label = pStar)) +
      scale_x_continuous(n.breaks = 5) +
      coord_flip()
  ) %>% `names<-`(c('Non-Epistatic','Epistatic'))
  
}

# plots one-way ANOVAs summaries based on genotype
oneBars <- function(m,b,f='grey50'){
  ggplot(data=m,aes(x=Genotype, y=Effect)) +
    theme_bw() +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_bar(stat="identity", width=0.9,fill=f,fill=NA) +
    geom_hline(yintercept = 0, color = "black") +
    geom_errorbar(aes(ymin=Effect-SE,ymax=Effect+SE,width=0.3)) +
    geom_text(aes(x=Genotype, 
                  y=(Effect+((SE+0.02)*sign(Effect))), label = pStar)) +
    scale_y_continuous(n.breaks = b) +
    coord_flip()
}

# summarizes one-way ANOVAs based on genotype
oneMaker <- function(dat,lvls,baseline){
  
  # generate each one-way anova
  models <-sapply(as.character(lvls[-1]),function(x)
    lmer(data=filter(dat,Genotype %in% c(baseline,x)), 
         Ratios_TMid ~ (1|plate) + (1|clone) + Genotype, REML=FALSE))
  
  # add additional information and plot 
  sapply(models,function(x)anova(x)) %>% t() %>% 
    apply(c(1,2),function(x)as.numeric(unlist(x))) %>% as.data.frame() %>%
    mutate(Effect = sapply(models,function(x)fixef(x)[-1])) %>%
    mutate(Genotype = factor(lvls[-1],
                             levels=rev(lvls[-1]))) %>%
    mutate(SE = sapply(models,function(x)
      summary(x)$coef[-1,'Std. Error', drop = FALSE])) %>%
    mutate(pStar = ifelse(`Pr(>F)`<0.05,`Pr(>F)`,NA)) %>%
    return()
}

# NegatE the %in% command
`%nin%` <- Negate(`%in%`)

# FIGURE OUT IF THE BOOTSTRAPS ARE SIGNIFICANTLY DIFFERENT
# FALSE INDICATES THE RANGES OVERLAP AND YOU ARE NOT SIGNIFICANT
# NOTE: TO CHECK IF TWO WELL FORMED RANGES OVERLAP
# A WELL FORMED RANGE IS ONE WHERE X_1 <= X_2
# THIS IS SUFFICIENT, AND IS ADAPTED BELOW: x1 <= y2 && y1 <= x2
strapSignif <- function(straps){
  apply(straps,2,function(x)quantile(x,c(0.025,0.975))) %>%
    (function(x)!(x[1,1] <=x[2,2] && x[1,2] <= x[2,1]))
}

# AN EASY AWY TO GET YOUR SNV NAMES
varEz <- function(x){
  return(paste0('BYv',x,'RM'))
}

# FUNCTION FOR PERFORMING EPISTASIS BOOTSTRAPS 
varStrap <- function(straps, vDat, bsl, alt, SNVs){
  
  # full target list
  targets <- c(bsl,alt,SNVs)
  
  # initialize your diffBYRM and svSums vectors
  df <- as.data.frame(matrix(NA,nrow=straps,ncol=2)) %>%
    `colnames<-`(c('diffBYRM','svSums'))
  
  # PREPARE A LIGHTER VERSION OF THE DATA FOR PERMUTATION
  lite <- vDat[which(vDat$Genotype %in% targets),
               c('Ratios_TMid','Genotype','plate','clone')]
  
  # MAIN BOOTSTRAPPING LOOP
  for(j in 1:straps){
    
    strat <- stratified(
      lite,"Genotype",
      replace=T,
      size=c((table(lite$Genotype)[table(lite$Genotype)!=0]))) %>% 
      as.data.frame()
    
    strapped <- fixef(
      lmer(data = strat,
           Ratios_TMid ~ (1|plate) + (1|clone) + Genotype, REML=FALSE))
    
    # get the names correct for proper assignment of results
    names(strapped) <- (str_remove(names(strapped),'Genotype'))
    
    # bootstrapped difference between BY and RM
    df[j,'diffBYRM'] <- strapped[alt]
    
    # bootstrapped sum of SNVs
    df[j,'svSums'] <- sum(strapped[SNVs])
    
    print(j)
    
  }
  
  return(df)
  
}
