rm(list=ls())
setwd('/Users/haleoriu/Documents/Takakia绘图和文章修改/outputs/Figure7_environment_association_20261013')
.libPaths(c(file.path(getwd(),'R_library'),.libPaths()))
library(vegan)
library(ggplot2)
library(ggrepel)
library(cowplot)
set.seed(20261013)
root <- '/Users/haleoriu/Documents/Takakia绘图和文章修改'
write_tsv <- function(x,file) write.table(x,file.path('tables',paste0(file,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
soil_path <- file.path(root,'dbRDA_soil_physicochemical_2026-07-30/soil_physicochemical_input.tsv')
soil <- read.delim(soil_path,check.names=FALSE)
rownames(soil) <- sub('^PS','',soil$Sample)
soil$Elevation <- soil$elev
soil_vars <- setdiff(names(soil),c('Sample','elev','Elevation'))
env <- soil[,c(soil_vars,'Elevation')]
envz <- as.data.frame(scale(env))
stopifnot(nrow(env)==30,!anyDuplicated(rownames(env)),all(is.finite(as.matrix(envz))))
vifcalc <- function(x) sapply(names(x),function(v) {if(ncol(x)==1) return(1);1/(1-summary(lm(reformulate(setdiff(names(x),v),response=v),data=x))$r.squared)})
screen <- function(name,vars,protected,threshold) {
 cur<-vars; steps<-list(); hist<-list(); k<-0
 repeat {
  vv<-vifcalc(envz[,cur,drop=FALSE]); candidates<-setdiff(cur,protected)
  action<-if(max(vv)<threshold)'stop' else if(length(candidates)==0)'infeasible' else 'remove'
  drop<-if(action=='remove') names(which.max(vv[candidates])) else NA_character_
  hist[[k+1]]<-data.frame(Preset=name,Step=k,Variable=names(vv),VIF=unname(vv),Protected=names(vv)%in%protected,Action=action,Removed=drop)
  steps[[k+1]]<-data.frame(Preset=name,Step=k,Predictors=paste(cur,collapse=','),Max_VIF=max(vv),Action=action,Removed=drop)
  if(action!='remove') break
  cur<-setdiff(cur,drop);k<-k+1
 }
 write_tsv(do.call(rbind,hist),paste0('VIF_',name,'_history'))
 write_tsv(do.call(rbind,steps),paste0('VIF_',name,'_steps'))
 write_tsv(data.frame(Preset=name,Variable=names(vv),VIF=unname(vv)),paste0('VIF_',name,'_final'))
 list(vars=cur,threshold=threshold,feasible=action=='stop',protected=protected,maxvif=max(vv))
}
presets<-list()
presets$oldVIF10<-screen('oldVIF10',names(env), 'Elevation',10)
infeasible<-screen('protected_SOC_TN_Elev_VIF10',names(env),c('SOC','TN','Elevation'),10)
presets$protectedVIF20<-screen('protectedVIF20',names(env),c('SOC','TN','Elevation'),20)
presets$soilVIF10<-screen('soilVIF10',soil_vars,c('SOC','TN'),10)
soil20<-screen('soilVIF20',soil_vars,c('SOC','TN'),20)
if(!identical(soil20$vars,presets$soilVIF10$vars))presets$soilVIF20<-soil20
stopifnot(!infeasible$feasible, all(vapply(presets,function(x)x$feasible,logical(1))))
protected_vif<-vifcalc(envz[,c('SOC','TN','Elevation')])
write_tsv(data.frame(Variable=names(protected_vif),VIF=unname(protected_vif)),'VIF_protected_only')
correlations<-do.call(rbind,lapply(names(env),function(a)do.call(rbind,lapply(names(env),function(b)data.frame(Variable1=a,Variable2=b,Pearson=cor(env[[a]],env[[b]]),Spearman=cor(env[[a]],env[[b]],method='spearman'))))))
write_tsv(correlations,'environment_pair_correlations')
write_tsv(data.frame(Preset=names(presets),Predictors=vapply(presets,function(x)paste(x$vars,collapse=','),character(1)),Threshold=vapply(presets,function(x)x$threshold,numeric(1)),MaxVIF=vapply(presets,function(x)x$maxvif,numeric(1))),'predictor_presets')
input_paths<-c(soil=soil_path,bacteria='/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data/asv_bacteria.txt',fungi='/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data/asv_fungi.txt',KO='/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/tpm_renormalized.txt')
prepare <- function(path,prefix,asv) {
 x<-read.delim(path,row.names=1,check.names=FALSE,comment.char='');ss<-paste0(prefix,rownames(env))
 stopifnot(all(ss%in%names(x)));m<-t(as.matrix(x[,ss]));storage.mode(m)<-'numeric';rownames(m)<-rownames(env)
 m<-m[,colSums(m)>0,drop=FALSE];stopifnot(all(is.finite(m)),all(m>=0),all(rowSums(m)>0))
 if(asv)m<-decostand(m,'total');m
}
profiles<-list('R bacteria'=prepare(input_paths['bacteria'],'S',TRUE),'E bacteria'=prepare(input_paths['bacteria'],'P',TRUE),'R fungi'=prepare(input_paths['fungi'],'S',TRUE),'E fungi'=prepare(input_paths['fungi'],'P',TRUE),'R KO'=prepare(input_paths['KO'],'S',FALSE),'E KO'=prepare(input_paths['KO'],'P',FALSE))
distances<-lapply(profiles,vegdist,method='bray')
saveRDS(list(env=env,envz=envz,presets=presets,profiles=profiles,distances=distances,soil_vars=soil_vars), 'models/prepared_data.rds')
write_tsv(data.frame(Input=names(input_paths),Path=unname(input_paths),MD5=unname(tools::md5sum(input_paths))),'input_manifest_overall')
results<-list();statslist<-list();termlist<-list();axislist<-list();siteslist<-list();arrowslist<-list()
labels<-c(AvK='K_ex',AvCa='Ca_ex',AvMg='Mg_ex',AvFe='Fe_ex',Moisture='Soil moisture')
pretty_var<-function(x)ifelse(x%in%names(labels),labels[x],x)
plot_one <- function(m,st,profile) {
 ss<-as.data.frame(scores(m,display='sites',choices=1:2,scaling=2));names(ss)<-c('x','y')
 aa<-as.data.frame(scores(m,display='bp',choices=1:2,scaling=2));names(aa)<-c('x','y');aa$Variable<-rownames(aa)
 multiplier<-.67*max(abs(as.matrix(ss)))/max(sqrt(aa$x^2+aa$y^2));aa[,c('x','y')]<-aa[,c('x','y')]*multiplier;aa$Label<-pretty_var(aa$Variable)
 ss$PairID<-rownames(ss);ss$Profile<-profile;ss$Preset<-st$Preset;aa$Profile<-profile;aa$Preset<-st$Preset
 siteslist[[paste(st$Preset,profile)]]<<-ss;arrowslist[[paste(st$Preset,profile)]]<<-aa
 ggplot(ss,aes(x,y))+geom_hline(yintercept=0,color='grey85',linewidth=.25)+geom_vline(xintercept=0,color='grey85',linewidth=.25)+geom_segment(data=aa,aes(x=0,y=0,xend=x,yend=y),color='grey45',linewidth=.45,arrow=grid::arrow(length=grid::unit(1.8,'mm')))+geom_point(shape=21,fill=if(substr(profile,1,1)=='R')'#82B3D7' else '#CF6C72',color='grey25',size=2.3,stroke=.3)+geom_text_repel(data=aa,aes(label=Label),size=3,family='Helvetica',seed=20261013,max.overlaps=Inf,min.segment.length=0)+labs(title=profile,subtitle=sprintf('adj. R² = %.3f; p = %.3f',st$Adjusted_R2,st$p),x=sprintf('db-RDA1 (%.1f%%)',st$Axis1pct),y=sprintf('db-RDA2 (%.1f%%)',st$Axis2pct))+theme_bw(base_size=11,base_family='Helvetica')+theme(panel.grid=element_blank(),plot.title=element_text(face='bold',hjust=.5),plot.subtitle=element_text(size=9,hjust=.5))
}
for(preset in names(presets)) {
 plots<-list();dataenv<-envz[,presets[[preset]]$vars,drop=FALSE]
 for(profile in names(profiles)) {
  cat('dbRDA',preset,profile,'\n');set.seed(20261013)
  d<-distances[[profile]];m<-capscale(d~.,data=dataenv,add='lingoes')
  overall<-anova.cca(m,permutations=999);axis<-anova.cca(m,by='axis',permutations=999);terms<-as.data.frame(anova.cca(m,by='margin',permutations=999));terms$Term<-rownames(terms);terms<-terms[terms$Term%in%names(dataenv),]
  terms<-data.frame(Preset=preset,Profile=profile,Term=terms$Term,Df=terms$Df,SumOfSqs=if('SumOfSqs'%in%names(terms))terms$SumOfSqs else terms$Variance,F=terms$F,p=terms$`Pr(>F)`)
  terms$q_BH_within_model<-p.adjust(terms$p,'BH');termlist[[paste(preset,profile)]]<-terms
  axs<-as.data.frame(axis);axs$Axis<-rownames(axs);axs$Preset<-preset;axs$Profile<-profile;axislist[[paste(preset,profile)]]<-axs
  er<-RsquareAdj(m);pct<-100*eigenvals(m,model='constrained')[1:2]/m$tot.chi
  st<-data.frame(Preset=preset,Profile=profile,N=30,N_features=ncol(profiles[[profile]]),Predictors=paste(names(dataenv),collapse=','),Raw_R2=er$r.squared,Adjusted_R2=er$adj.r.squared,F=overall$F[1],p=overall$`Pr(>F)`[1],Axis1pct=pct[1],Axis2pct=pct[2],MaxVIF=max(vif.cca(m)))
  stopifnot(st$MaxVIF<presets[[preset]]$threshold)
  statslist[[paste(preset,profile)]]<-st;results[[paste(preset,profile)]]<-m;plots[[profile]]<-plot_one(m,st,profile)
 }
 fig<-plot_grid(plotlist=plots,ncol=2,labels=letters[1:6],label_size=15)
 ggsave(file.path('figures',paste0('Figure7_dbRDA_',preset,'.pdf')),fig,width=9,height=11.8,device="pdf",family="Helvetica",useDingbats=FALSE)
 ggsave(file.path('figures',paste0('Figure7_dbRDA_',preset,'.png')),fig,width=9,height=11.8,dpi=300)
}
write_tsv(do.call(rbind,statslist),'dbRDA_statistics');write_tsv(do.call(rbind,termlist),'dbRDA_marginal_terms');write_tsv(do.call(rbind,axislist),'dbRDA_axes');write_tsv(do.call(rbind,siteslist),'dbRDA_sites');write_tsv(do.call(rbind,arrowslist),'dbRDA_arrows')
saveRDS(results,'models/dbRDA_models.rds')
sink('sessionInfo_01.txt');print(sessionInfo());sink()
cat('STAGE 01 COMPLETE\n')
