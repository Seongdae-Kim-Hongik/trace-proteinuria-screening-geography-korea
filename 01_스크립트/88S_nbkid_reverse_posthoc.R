# 88S_nbkid_reverse_posthoc.R - Post hoc analyses (after results were seen) of the association opposite to the prespecified direction (public assistance recipient rate).
#  P0 principal model; P1 plus sensitivity-only variables; P2 alone with fixed adjustments; P3 plus log population density; P4 between/within
#  decomposition; P6 without physicians; P7 2013-2024 only; P8 expected counts conditional on the kidney-disease-suspected judgement; P9 without uptake.
#  Usage: TAG=NBKID ST_VER=std3 ST_GRAPH=<graph> Rscript 88S_nbkid_reverse_posthoc.R. Output: results_std3/S46_nbkid_posthoc.csv
suppressPackageStartupMessages({library(INLA);library(sf);library(jsonlite)})
TAG<-Sys.getenv("TAG","NBKID");VER<-Sys.getenv("ST_VER","std3");inla.setOption(num.threads="1:1")
OUTD<-file.path("..","04_분석",TAG);RES<-file.path(OUTD,paste0("results_",VER))
SHP<-Sys.getenv("ST_SHP",unset="data/boundary/final.shp")
S<-fromJSON(file.path(OUTD,sprintf("config_%s_%s.json",TAG,VER)),simplifyVector=TRUE)
DD<-read.csv(file.path(OUTD,sprintf("panel_%s_%s.csv",TAG,VER)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);g<-inla.read.graph(Sys.getenv("ST_GRAPH"));stopifnot(g$n==nrow(G))
DD$area<-match(DD$fbd_region,G$region);DD$time<-DD$year-min(DD$year)+1;DD$yiid<-DD$time
FORCED<-S$forced;HYP<-as.character(S$sets_after_vif$S1_narrative_full)
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
z<-function(x)(x-mean(x))/sd(x)
fit<-function(vars){f<-as.formula(paste("cases ~ 1 +",paste(sprintf("`%s`",vars),collapse="+"),"+ f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc) + f(time,model='rw1',scale.model=TRUE,hyper=pcp) + f(yiid,model='iid',hyper=pcp)"))
  inla(f,family="nbinomial",data=DD,E=DD$E,control.fixed=list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001),control.compute=list(waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}
out<-NULL
add<-function(an,m,vs){fx<-m$summary.fixed;for(v in vs){if(!v%in%rownames(fx))next;r<-fx[v,];out<<-rbind(out,data.frame(analysis=an,cov=v,IRR=round(exp(r$mean),4),lo=round(exp(r$`0.025quant`),4),hi=round(exp(r$`0.975quant`),4),cred=(exp(r$`0.025quant`)>1|exp(r$`0.975quant`)<1),WAIC=round(m$waic$waic,1)))};cat(an,"done\n")}
z<-function(x)(x-mean(x))/sd(x)
F1<-"basic_livelihood_rate";REV<-c(F1)
DD$popdens_log<-z(log(DD$population/(DD$area_m2/1e6)))
for(v in REV){m<-ave(DD[[v]],DD$fbd_region);DD[[paste0(v,"__B")]]<-(m-mean(m))/sd(m);w<-DD[[v]]-m;DD[[paste0(v,"__W")]]<-w/sd(w)}
MAIN<-unique(c(HYP,FORCED));SENS<-intersect(c("tertiary_educated_share","older_single_household_share","moderate_pa_rate","public_sports_facilities_per_100k"),names(DD))
add("P0_main",fit(MAIN),REV)
add("P1_plus_sensitivity_vars",fit(c(MAIN,SENS)),c(REV,SENS))
add("P2a_livelihood_alone",fit(c(F1,FORCED)),F1)
add("P3_plus_popdens",fit(c(MAIN,"popdens_log")),c(REV,"popdens_log"))
BW<-c(setdiff(MAIN,REV),paste0(REV,"__B"),paste0(REV,"__W"));add("P4_between_within",fit(BW),c(paste0(REV,"__B"),paste0(REV,"__W")))
add("P6_drop_physicians",fit(setdiff(MAIN,"physicians_per_1000")),REV)
dx<-DD;dx$cases[dx$year%in%2011:2012]<-NA;DD0<-DD;DD<-dx;add("P7_excl_2011_2012",fit(MAIN),REV);DD<-DD0
K1<-read.csv(file.path("..","04_분석","KIDS","panel_KIDS.csv"),stringsAsFactors=FALSE);K1$fbd_region<-enc2utf8(K1$fbd_region);k<-paste(DD$year,DD$fbd_region)
ks<-K1$cases[match(k,paste(K1$year,K1$fbd_region))];cat("P8 kids missing",sum(is.na(ks)),"\n")
h<-ks/DD$offset;h[is.na(h)]<-ave(h,DD$year,FUN=function(x)median(x,na.rm=TRUE))[is.na(h)];Ec<-DD$E*(1-h);Ec<-Ec*ave(DD$cases,DD$year,FUN=sum)/ave(Ec,DD$year,FUN=sum)
cat("P8 kids share range",round(range(h),3),"\n");DD0<-DD;DD$E<-Ec;add("P8_hierarchy_conditional",fit(MAIN),REV);DD<-DD0
add("P9_drop_screening_uptake",fit(setdiff(MAIN,"screening_uptake")),REV)
write.csv(out,file.path(RES,"S46_nbkid_posthoc.csv"),row.names=FALSE);print(out)
