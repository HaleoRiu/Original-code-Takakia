rm(list=ls())
setwd(file.path(Sys.getenv("TAKAKIA_PROJECT_DIR", unset = "~/Documents/Takakia"), "outputs/Figure7_environment_association_20261013/evidence_priority_no_AN_v1"))
.libPaths(c(file.path(Sys.getenv("TAKAKIA_PROJECT_DIR", unset = "~/Documents/Takakia"), "outputs/Figure7_environment_association_20261013/R_library"),.libPaths()))
library(vegan)
library(ggplot2)
library(ggrepel)
library(cowplot)
set.seed(20261013)
root <- Sys.getenv("TAKAKIA_PROJECT_DIR", unset = "~/Documents/Takakia")
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
# Evidence ranking uses existing response-based results; exploratory selection, not independent validation.
mantel_source<-file.path(root,'controlled_altitude_2026-07-30/partial_Mantel_results_all.tsv')
s4_source<-file.path(root,'enriched_pathway_soil_partial_spearman_2026-07-30/Spearman_summed_pathway_soil_BH_results.tsv')
mantel<-read.delim(mantel_source);s4<-read.delim(s4_source)
rank_table<-do.call(rbind,lapply(soil_vars,function(v) {
 a<-mantel[mantel$env==v,];r<-s4[s4$Soil_variable==v&s4$Direction=='Rhizomatousphere enriched',];e<-s4[s4$Soil_variable==v&s4$Direction=='Endophytic enriched',]
 stopifnot(nrow(a)==6,nrow(r)==10,nrow(e)==9)
 data.frame(Variable=v,Mantel_significant_profiles=sum(a$p_BH<.05),S4_R_significant=sum(r$p_BH_within_pathway<.05),S4_E_significant=sum(e$p_BH_within_pathway<.05),S4_equal_habitat_significant_fraction=(mean(r$p_BH_within_pathway<.05)+mean(e$p_BH_within_pathway<.05))/2,Mantel_mean_abs_r=mean(abs(a$r)),S4_equal_habitat_mean_abs_rho=(mean(abs(r$rho))+mean(abs(e$rho)))/2,Input_order=match(v,soil_vars),Excluded_by_user=v=='AN')
}))
rank_table<-rank_table[order(-rank_table$Mantel_significant_profiles,-rank_table$S4_equal_habitat_significant_fraction,-rank_table$Mantel_mean_abs_r,-rank_table$S4_equal_habitat_mean_abs_rho,rank_table$Input_order),];rank_table$Evidence_rank<-seq_len(nrow(rank_table))
write_tsv(rank_table,'environment_evidence_priority_rank')
write_tsv(data.frame(Source=c('partial_Mantel_saved_results','S4_saved_results'),Path=c(mantel_source,s4_source),MD5=unname(tools::md5sum(c(mantel_source,s4_source)))),'input_manifest_evidence_rank')
rankvalue<-setNames(rank_table$Evidence_rank,rank_table$Variable);rankvalue['Elevation']<-0
all_pairs<-function(vars) {
 ij<-t(combn(vars,2));data.frame(Variable1=ij[,1],Variable2=ij[,2],Pearson=sapply(seq_len(nrow(ij)),function(i)cor(env[[ij[i,1]]],env[[ij[i,2]]])),Spearman=sapply(seq_len(nrow(ij)),function(i)cor(env[[ij[i,1]]],env[[ij[i,2]]],method='spearman')))
}
screen <- function(name,vars,protected,threshold) {
 cur<-vars;steps<-list();hist<-list();candidates_history<-list();k<-0
 repeat {
  vv<-vifcalc(envz[,cur,drop=FALSE]);candidates<-setdiff(cur,protected);drop<-NA_character_;pair1<-pair2<-NA_character_;strength<-NA_real_
  cp<-if(length(cur)>1)all_pairs(cur)else data.frame();action<-'stop';reason<-'all_pair_strength_lt_0.7_and_joint_VIF_lt10'
  if(nrow(cp)>0) {
   cp$Strength<-pmax(abs(cp$Pearson),abs(cp$Spearman));cp<-cp[order(-cp$Strength,match(cp$Variable1,names(env)),match(cp$Variable2,names(env))),]
   if(cp$Strength[1]>=.7) {
    pair1<-cp$Variable1[1];pair2<-cp$Variable2[1];strength<-cp$Strength[1];ends<-c(pair1,pair2);eligible<-setdiff(ends,protected)
    if(length(eligible)==0){action<-'infeasible';reason<-'protected_pair_strong'}else{drop<-eligible[which.max(rankvalue[eligible])];action<-'remove';reason<-'strong_pair_remove_weaker_evidence_endpoint'}
   }
  }
  if(action=='stop'&&max(vv)>=threshold) {
   after<-sapply(candidates,function(v)max(vifcalc(envz[,setdiff(cur,v),drop=FALSE])))
   improve<-max(vv)-after;eligible<-candidates[improve>1e-8]
   if(length(eligible)==0){action<-'infeasible';reason<-'no_VIF_improving_candidate'}else{drop<-eligible[which.max(rankvalue[eligible])];action<-'remove';reason<-'joint_VIF_remove_weakest_improving_candidate'}
   candidates_history[[length(candidates_history)+1]]<-data.frame(Preset=name,Step=k,Candidate=candidates,Evidence_rank=unname(rankvalue[candidates]),MaxVIF_before=max(vv),MaxVIF_after=unname(after),Improvement=unname(improve),Eligible_improvement_gt1e_8=improve>1e-8,Selected=candidates%in%drop)
  }
  hist[[k+1]]<-data.frame(Preset=name,Step=k,Variable=names(vv),VIF=unname(vv),Evidence_rank=unname(rankvalue[names(vv)]),Protected=names(vv)%in%protected,Action=action,Reason=reason,Removed=drop)
  steps[[k+1]]<-data.frame(Preset=name,Step=k,Predictors=paste(cur,collapse=','),Max_VIF=max(vv),Action=action,Reason=reason,Conflict_variable1=pair1,Conflict_variable2=pair2,Conflict_strength=strength,Removed=drop)
  if(action!='remove')break
  cur<-setdiff(cur,drop);k<-k+1
 }
 write_tsv(do.call(rbind,hist),paste0('VIF_',name,'_history'));write_tsv(do.call(rbind,steps),paste0('VIF_',name,'_steps'));write_tsv(data.frame(Preset=name,Variable=names(vv),VIF=unname(vv)),paste0('VIF_',name,'_final'))
 if(length(candidates_history))write_tsv(do.call(rbind,candidates_history),paste0('VIF_',name,'_candidate_deletion_diagnostics'))
 list(vars=cur,threshold=threshold,feasible=action=='stop',protected=protected,maxvif=max(vv))
}
presets<-list()
presets$Elevation_evidencePriorityVIF10<-screen('Elevation_evidencePriorityVIF10',setdiff(names(env),'AN'),'Elevation',10)
presets$soil_evidencePriorityVIF10<-screen('soil_evidencePriorityVIF10',setdiff(soil_vars,'AN'),character(),10)
stopifnot(all(vapply(presets,function(x)x$feasible,logical(1))))
correlations<-do.call(rbind,lapply(names(env),function(a)do.call(rbind,lapply(names(env),function(b)data.frame(Variable1=a,Variable2=b,Pearson=cor(env[[a]],env[[b]]),Spearman=cor(env[[a]],env[[b]],method='spearman'))))))
write_tsv(correlations,'environment_pair_correlations')
# Diagnose all 13 starting variables, not only protected SOC/TN.
unique_pairs<-correlations[match(correlations$Variable1,names(env))<match(correlations$Variable2,names(env)),]
unique_pairs$Abs_Pearson<-abs(unique_pairs$Pearson);unique_pairs$Abs_Spearman<-abs(unique_pairs$Spearman)
unique_pairs$Pearson_flag<-ifelse(unique_pairs$Abs_Pearson>=.8,'high_abs_r_ge_0.8',ifelse(unique_pairs$Abs_Pearson>=.7,'strong_abs_r_ge_0.7','below_0.7'))
unique_pairs$Spearman_flag<-ifelse(unique_pairs$Abs_Spearman>=.8,'high_abs_rho_ge_0.8',ifelse(unique_pairs$Abs_Spearman>=.7,'strong_abs_rho_ge_0.7','below_0.7'))
unique_pairs<-unique_pairs[order(-unique_pairs$Abs_Pearson),];stopifnot(nrow(unique_pairs)==78)
write_tsv(unique_pairs,'environment_unique_pairs_sorted')
ep<-unique_pairs[unique_pairs$Variable1=='Elevation'|unique_pairs$Variable2=='Elevation',];write_tsv(ep,'Elevation_all12_soil_pair_correlations')
finalpairs<-do.call(rbind,lapply(names(presets),function(p) {x<-unique_pairs[unique_pairs$Variable1%in%presets[[p]]$vars&unique_pairs$Variable2%in%presets[[p]]$vars,];x$Preset<-p;x}))
write_tsv(finalpairs,'preset_final_pair_correlations')
write_tsv(do.call(rbind,lapply(names(presets),function(p){x<-finalpairs[finalpairs$Preset==p,];data.frame(Preset=p,Max_abs_Pearson=max(x$Abs_Pearson),Max_abs_Spearman=max(x$Abs_Spearman),N_strong_Pearson_ge07=sum(x$Abs_Pearson>=.7),N_high_Pearson_ge08=sum(x$Abs_Pearson>=.8))})),'preset_final_correlation_diagnostics')
pdf('figures/Environmental_correlations_all13.pdf',width=10,height=9,family='Helvetica',useDingbats=FALSE)
for(method in c('Pearson','Spearman')) {
 x<-correlations;x$Variable1<-factor(x$Variable1,levels=names(env));x$Variable2<-factor(x$Variable2,levels=rev(names(env)));x$Value<-x[[method]]
 p<-ggplot(x,aes(Variable1,Variable2,fill=Value))+geom_tile(color='white',linewidth=.2)+geom_text(aes(label=sprintf('%.2f',Value)),size=2.6)+scale_fill_gradient2(low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-1,1))+coord_equal()+labs(title=paste('All soil variables and elevation:',method,'correlations'),x=NULL,y=NULL,caption='Models require max(|Pearson r|, |Spearman rho|) < 0.7 and joint VIF < 10.\nAN excluded; conflicts resolved by the fixed evidence ranking; elevation retained in the elevation model.')+theme_minimal(base_size=10,base_family='Helvetica')+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=45,hjust=1),plot.caption=element_text(hjust=0,size=8))
 print(p)
}
dev.off()

write_tsv(data.frame(Preset=names(presets),Predictors=vapply(presets,function(x)paste(x$vars,collapse=','),character(1)),Threshold=vapply(presets,function(x)x$threshold,numeric(1)),MaxVIF=vapply(presets,function(x)x$maxvif,numeric(1))),'predictor_presets')
if('screen_only'%in%commandArgs(trailingOnly=TRUE))quit(save='no')
input_paths<-c(soil=soil_path,bacteria='/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data/asv_bacteria.txt',fungi='/Volumes/Haleo_DK01/WD_Takakia_F/amplicon-data/asv_fungi.txt',KO='/Volumes/Haleo_DK01/WD_Takakia_F/metagenome-data/tpm_renormalized.txt')
old_prepared<-readRDS('../models/prepared_data.rds')
stopifnot(identical(env,old_prepared$env),identical(envz,old_prepared$envz))
profiles<-old_prepared$profiles;distances<-old_prepared$distances
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
