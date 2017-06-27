library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)

#function to add NA for missing col & row
napad<-function(m) {m<-as.matrix(m); m[,1]<-NA; m<-rbind(m,rep(NA,ncol(m)))}
#function to make matrix symmetrical & diag == 1
resym<-function(m,diagval=1) {m[is.na(m)]<-0; m <- m+ t(m);  diag(m)<-diagval;return(m)}


#specify a function for the linear model for site predicting each roi-roi pair
getresid_sitevisit <- function(site,value) lm(value~site,data=data.frame(site,value))$residuals


# function to remove residuals
# to modify change:
#  1. add new function getresid_* copy of getresid_sitevisit with new regress.
#  2. add variable to `corrlong` by modifying the `collapse` subfunction
#  3. change resids= in mutate that creates withresid
mkresidCorr <- function(filenames,sites) {
  #read in all filenames
  allcor_list=lapply(filenames,read.table)
  #specify site
  
  #specify function that you will use later
  # needs allcor_list to be created
  collapse<-function(i){
    d<-allcor_list[[i]]
    site<-sites[[i]]
    dm<-as.matrix(d)
    dm[lower.tri(d,diag=T)]<-NA

    dlong<-melt(dm)
    dlong<-dlong[!is.na(dlong$value),]
    dlong$site<-site
    dlong$i<-i
    return(dlong)
  }
  
  #collapse all matrices into one long one
  #need to have sites & allcor_list defined already
  #put correlation data in long form (value), specifies the roi corr by row (Var1) and column (Var2)
  corrlong<-lapply(1:length(sites),collapse) %>% rbind.fill
  
  
  #model each roi-roi pair, get residuals
  # == this takes a bit ==
  withresid<-corrlong %>% 
      group_by(Var1,Var2)%>%
      mutate(resids=getresid_sitevisit(site,value))
  
  # reshape per subject (i)
  d.reshape <- withresid %>%group_by(i) %>% do(dcast(.,Var1~Var2,value.var='resids')) 

  # nest that per subject
  d.nest <- d.reshape %>% nest 

  # finally make the nested matrix look like how we expect corr mats too look
  res <- d.nest %>% mutate(data=lapply(data, function(x) x %>% napad %>% resym ))
}

#### CHANGE CODE BELOW

#subj = read.table('/raid6/Maria/Neuroimaging/PNC/scripts/Aim01/txt/PNC_subj_info.txt', head = T)

filenames <- list("600009963128_Parcels_MNI_111_warped_rs_pearson_original.txt","600405811873_Parcels_MNI_111_warped_rs_pearson_original.txt","601004917174_Parcels_MNI_111_warped_rs_pearson_original.txt")

sites<-c(1,2,1)

# USAGE:
res <- mkresidCorr(filenames,sites)
# test
res$data[[3]][1:4,1:4]
res$data[[1]][1:4,1:4]

# TODO: save
if(!dir.exists('out')) dir.create('out')


