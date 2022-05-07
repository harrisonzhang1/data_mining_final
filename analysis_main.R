rm(list=ls())
library(data.table)
library(dplyr)
library(RODBC)
library(BSDA)
library(ggplot2)
library(pheatmap)
odbc <- odbcConnect("XXX", "XXX", .rs.api.askForPassword("password"))
dir.data="dbmart_rollup_pasc"
dir.repo="mining_final_project/"
setwd(dir.data)
pcr.flag <- readRDS("pcr_flag.RDS")
colnames(pcr.flag)=c("patient_num","encounter_num",
                     "date","concept_code","concept",
                     "value")
load("phe_str.Rdata")
pcr.pos=pcr.flag %>%
  filter(concept_code=="pcr_positive") %>%
  arrange(patient_num,date)%>%
  group_by(patient_num)%>%
  slice(1)%>%
  ungroup()
all.files=list.files(dir.data); all.files=all.files[grepl(".csv",all.files)]
mylist=lapply(all.files,FUN = function(ll){
  junk=fread(ll,colClasses = "character")
  return(filter(junk,patient_num %in% pcr.pos$patient_num))
})
data=rbindlist(mylist)
rm(mylist)
data.pasc=left_join(data,
                    select(pcr.pos,
                           patient_num,
                           date,
                           concept_code))
data.pasc = data.pasc %>%
  mutate(inf_date=date) %>%
  select(-date) %>%
  mutate(days_since_admission=as.Date(start_date)-as.Date(inf_date))
data.pasc=data.pasc[!duplicated(data.pasc),]
data.pasc=data.pasc %>%
  filter(days_since_admission>=-(365*1))
save(data.pasc,
     file="P:/Harrison/data_mining/data_pasc.Rdata")
dat.icd.f=data.pasc %>% arrange(patient_num, days_since_admission)
dat.icd.f=dat.icd.f %>% group_by(patient_num, concept_cd) %>% slice(1) %>% ungroup()
save(dat.icd.f,
     file="P:/Harrison/data_mining/dat_icd_f.Rdata")
icd.keep=dat.icd.f %>%
  filter(days_since_admission>0,
         days_since_admission<=90,
         grepl("PheCode",concept_cd))
icd.keep$ID=paste0(icd.keep$patient_num,"_",icd.keep$concept_cd)
id.keep=unique(icd.keep$ID)
data.pasc$ID=paste0(data.pasc$patient_num,"_",data.pasc$concept_cd)
data.pasc=data.pasc %>%
  filter(ID %in% id.keep)
phe.str=phe.str%>%
  mutate(concept_cd=paste0("PheCode:",phenotype)) %>%
  select(-phenotype)
data.pasc=left_join(data.pasc,
                    phe.str)
# data.pasc=data.pasc%>%
#   filter(!is.na(group))
junk=filter(data.pasc,!is.na(description))
### Exploration of the Data
length(unique(data.pasc$patient_num))
dim(data.pasc)
summary(as.Date(data.pasc$start_date))
util=data.pasc%>%
  group_by(patient_num)%>%
  summarize(n=n())
summary(util$n)
summary(as.numeric(data.pasc$days_since_admission))
visit=data.pasc%>%
  group_by(patient_num)%>%
  summarize(n=length(unique(as.Date(start_date))))
summary(visit$n)
group=data.pasc%>%
  group_by(group)%>%
  summarize(p=n()/nrow(data.pasc))
### Feature Engineering
hosp = readRDS("hospitalizations.RDS")
pcr.pos=pcr.pos %>%
  filter(!is.na(date))
pcr.pos=pcr.pos[!duplicated(pcr.pos),]
pat=unique(data.pasc$patient_num)
hosp=left_join(select(pcr.pos,patient_num,date),
               hosp)
hosp$hosp_flag=ifelse(abs(as.Date(hosp$date)-as.Date(hosp$start_date))<=7,1,0)
hosp$hosp_flag[is.na(hosp$hosp_flag)]=0
sum(hosp$hosp_flag)

### Create TF - IDF matrix for kmeans clustering
### feataures in columns, 60k patients in rows

icd.tmp=data.pasc %>%
  filter(days_since_admission>30) %>%
  select(patient_num,concept_cd)
icd.tmp$patient_num=as.character(icd.tmp$patient_num)
icd.tmp$concept_cd=as.character(icd.tmp$concept_cd)
icd.tmp=table(icd.tmp)
icd.tmp=data.frame(icd.tmp)
icd.tmp=arrange(icd.tmp,concept_cd,patient_num)
num.patients=length(unique(icd.tmp$patient_num))
num.patients
sum(icd.tmp$patient_num[1:num.patients]==icd.tmp$patient_num[(num.patients+1):(num.patients*2)])
phecodes.conf=as.character(unique(icd.tmp$concept_cd))
res.out.90=as.list(phecodes.conf)
for(pp in 1:length(phecodes.conf)){
  res.out.90[[pp]]=c(icd.tmp$Freq[((pp-1)*num.patients+1):(pp*num.patients)])
}
tf=do.call("cbind",res.out.90)
colnames(tf)=phecodes.conf
rownames(tf)=icd.tmp$patient_num[1:num.patients]
dim(tf)
sum(c(tf==0))/(nrow(tf)*ncol(tf))
junk=tf>0
id.keep=which(colMeans(junk)>=0.01)
tf=tf[,id.keep]
tf=log(1+tf)
id.rm=phe.str$concept_cd[phe.str$group %in%
                           c("genitourinary",
                             "pregnancy complications",
                             "neoplasms",
                             "injuries & poisonings",
                             "congenital anomalies")]
id.rm=c("PheCode:1000",id.rm)
tf=tf[,which(!colnames(tf)%in%id.rm)]
nt=sapply(colnames(tf),FUN=function(ll){
  return(length(unique(data.pasc$patient_num[data.pasc$concept_cd==ll])))
})
idf=log(nrow(tf)/nt)
junk=as.list(1:nrow(tf))
for(tt in 1:nrow(tf)){junk[[tt]]=idf}
idf=do.call("rbind",junk)
tfidf=tf*idf
mtfidf=apply(tfidf,MARGIN=2,max)
hist(mtfidf)
cut=max(mtfidf[grepl("PheCode:10",names(mtfidf))])
mtfidf=data.frame(mtfidf=mtfidf)
ggplot(mtfidf, aes(x=mtfidf)) + 
  geom_histogram( colour="black", fill="white") + 
  ggtitle("Distribution of Largest Diagnosis Code TF-IDF Value across Patients") +
  xlab("Maximum TF-IDF Value") +
  ylab("Number of Diagnosis Codes") + 
  geom_vline(xintercept = cut,
             color = "red", size=1.5)
mtfidf=apply(tfidf,MARGIN=2,max)
tfidf.m=tfidf[,names(mtfidf)[which(mtfidf>cut)]]
pat.rm=rowSums(tfidf.m)
tfidf.m=tfidf.m[which(pat.rm!=0),]
hosp=filter(hosp,patient_num %in% rownames(tfidf.m))
hosp$patient_num=as.character(hosp$patient_num)
hosp.flag=sapply(rownames(tfidf.m),FUN=function(ll){
  if(sum(hosp$hosp_flag[hosp$patient_num==ll])>0){
    return(1)
  }else(return(0))
})
tfidf.m=cbind(tfidf.m,hosp.flag)
k_max <- 20
outs <- matrix(NA, ncol=2, nrow=k_max)
for(k_guess in seq_len(k_max)){
  set.seed(281983981)
  km_out <- kmeans(tfidf.m, centers=k_guess)
  outs[k_guess, 1] <- km_out$betweenss
  outs[k_guess, 2] <- km_out$tot.withinss
}
plot(outs[,1]/outs[,2],
     type="b",
     xlab="Number of Clusters",
     ylab="Ratio: Between Cluster SS and W/in Cluster SS",
     main="Clustering Metric To Choose K Number of Clusters")

km_out <- kmeans(tfidf.m, centers=5)
clus=data.frame("patient_num"=rownames(tfidf.m),
                "cluster"=km_out$cluster)
table(factor(clus$cluster))

plot.res=NULL
size=50
for(cc in 1:5){
m=apply(tfidf[as.character(clus$patient_num[clus$cluster==cc]),],MARGIN=2,max)
phe=tail(sort(m), size)
phe=names(phe)
grp=phe.str$group[phe.str$concept_cd %in% phe]
plot.res=rbind.data.frame(plot.res,
                          cbind.data.frame("Cluster"=cc,
                                           "Body_System"=grp,
                                           "count"=1))

}
plot.res=plot.res%>%
  group_by(Cluster,Body_System)%>%
  summarize(n=sum(count))%>%
  mutate("p"=n/size)%>%
  arrange(Cluster,Body_System)
grp=as.character(unique(plot.res$Body_System))
mm.plot=matrix(NA,nrow=length(grp),ncol=5)
for(cc in 1:5){
  tmp=sapply(grp,FUN=function(ll){
    tt=filter(plot.res,Cluster==cc,Body_System==ll)
    if(nrow(tt)!=0){return(tt$p)
      }else(return(0))
  })
  mm.plot[,cc]=tmp
}
colnames(mm.plot)=paste0("Cluster_",1:5)
rownames(mm.plot)=grp
pheatmap(mm.plot,
         cluster_rows = F,
         cluster_cols = F,
         main = "Body Systems of 50 Highest TF-IDF Diagnosis Code Features (Proportion of 50 Codes)")
clus$hosp_flag=tfidf.m[,"hosp.flag"]

clus%>%
  group_by(cluster)%>%
  summarise("p.hosp"=sum(hosp_flag)/n())


### Boostrapping
mm.plot.final=matrix(0,nrow=17,ncol=5)
for(bb in 1:50){
  print(bb)
tfidf.s=tfidf.m[sample(nrow(tfidf.m),replace = T),]
km_out <- kmeans(tfidf.s, centers=5)
clus=data.frame("patient_num"=rownames(tfidf.s),
                "cluster"=km_out$cluster)
table(factor(clus$cluster))
plot.res=NULL
size=50
for(cc in 1:5){
  m=apply(tfidf[as.character(clus$patient_num[clus$cluster==cc]),],MARGIN=2,max)
  phe=tail(sort(m), size)
  phe=names(phe)
  grp=phe.str$group[phe.str$concept_cd %in% phe]
  plot.res=rbind.data.frame(plot.res,
                            cbind.data.frame("Cluster"=cc,
                                             "Body_System"=grp,
                                             "count"=1))
  
}
plot.res=plot.res%>%
  group_by(Cluster,Body_System)%>%
  summarize(n=sum(count))%>%
  mutate("p"=n/size)%>%
  arrange(Cluster,Body_System)
grp=as.character(unique(phe.str$group))
mm.plot=matrix(NA,nrow=length(grp),ncol=5)
for(cc in 1:5){
  tmp=sapply(grp,FUN=function(ll){
    tt=filter(plot.res,Cluster==cc,Body_System==ll)
    if(nrow(tt)!=0){return(tt$p)
    }else(return(0))
  })
  mm.plot[,cc]=tmp
}
colnames(mm.plot)=paste0("Cluster_",1:5)
rownames(mm.plot)=grp
index=c("circulatory system",
        "dermatologic",
        "digestive",
        "endocrine/metabolic",
        "hematopoietic",
        "mental disorders",
        "musculoskeletal",
        "neurological",
        "respiratory",
        "sense organs",
        "symptoms",
        "infectious diseases")
pheatmap(mm.plot[index,],
         cluster_rows = F,
         cluster_cols = F,
         main = "Bootstrapped - Body Systems of Diagnoses with 50 Highest TF-IDF (Proportion of 50 Codes)")


}





























