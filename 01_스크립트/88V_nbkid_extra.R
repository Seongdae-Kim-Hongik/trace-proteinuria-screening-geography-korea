# 88V_nbkid_extra.R - Additional sensitivity analyses of the principal model (10 hypothesised characteristics + 5 fixed adjustments).
#  (1) one sensitivity-only variable added at a time: annual mean temperature, screening-measured suspected hypertension share,
#  survey obesity rate, period-standardised screening obesity-management share; (2) without Sejong; without districts with boundary
#  changes (Sejong, Cheongju, Gunwi); (3) without 2020-2021; from 2015; (4) BYM2 mixing-parameter priors PC(0.5, 2/3), Beta(1,1), Jeffreys Beta(0.5,0.5).
#  Usage: ST_VER=std3 ST_SHP=<boundary file> Rscript 88V_nbkid_extra.R NBKID (requires the 81 outputs, incl. results_std3/nb_NBKID.graph).
#  Output: results_std3/S60_extra.csv (rate ratios, 95% CrI, WAIC and mixing parameter for each analysis; resumable)
suppressPackageStartupMessages({library(INLA);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);TAG<-args[1];ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
VER<-Sys.getenv("ST_VER",unset="std3");SUF<-paste0("_",VER);OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF))
inla.setOption(num.threads="1:1");inla.setOption(inla.timeout=900)
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
SHP<-Sys.getenv("ST_SHP",unset="data/boundary/final.shp")
reg<-enc2utf8(sf::st_read(SHP,quiet=TRUE)$region);DD$fbd_region<-enc2utf8(DD$fbd_region)
DD$area<-match(DD$fbd_region,reg);stopifnot(!anyNA(DD$area))
g<-inla.read.graph(file.path(RES,sprintf("nb_%s.graph",TAG)))
V<-unique(c(S$sets_after_vif$S1_narrative_full,S$forced));EXTRA<-c("annual_mean_temperature_v","screening_suspected_hypertension_share_v","obesity_rate_sr","ext__scr_obesity_pstd")
stopifnot(all(c(V,EXTRA)%in%names(DD)))
zs<-function(x)(x-mean(x))/sd(x)
pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001);bt<-function(v)sprintf("`%s`",v)
phiprior<-list(pc=list(prior="pc",param=c(0.5,2/3)),beta11=list(prior="logitbeta",param=c(1,1)),jeffreys=list(prior="logitbeta",param=c(0.5,0.5)))
fit<-function(d,vars,phi="pc"){d<-d;for(v in vars)d[[v]]<-zs(d[[v]]);Y0<-min(d$year);d$time<-d$year-Y0+1;d$yiid<-d$time
  hy<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=phiprior[[phi]])
  f<-as.formula(paste("cases ~",paste(c("1",bt(vars),"f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=hy)","f(time,model='rw1',scale.model=TRUE,hyper=pcp)","f(yiid,model='iid',hyper=pcp)"),collapse="+")))
  inla(f,family="nbinomial",data=d,E=d$E,control.fixed=CF,control.compute=list(waic=TRUE),control.inla=list(control.vb=list(emergency=50)))}
OUT<-file.path(RES,"S60_extra.csv");R<-if(file.exists(OUT))read.csv(OUT,stringsAsFactors=FALSE) else NULL
add<-function(lab,m,vars,n){if(is.null(m))return();s<-m$summary.fixed[vars,];ph<-if("Phi for area"%in%rownames(m$summary.hyperpar))m$summary.hyperpar["Phi for area","mean"] else NA
  R<<-rbind(R,data.frame(analysis=lab,cov=vars,IRR=round(exp(s$mean),4),lo=round(exp(s$`0.025quant`),4),hi=round(exp(s$`0.975quant`),4),cred=(s$`0.025quant`>0|s$`0.975quant`<0),WAIC=round(m$waic$waic,1),phi=round(ph,3),n_obs=n));write.csv(R,OUT,row.names=FALSE);cat(lab,"done\n")}
done<-function(lab)!is.null(R)&&lab%in%R$analysis
run<-function(lab,d,vars,phi="pc"){if(done(lab))return();m<-tryCatch(fit(d,vars,phi),error=function(e){message(lab," FAILED: ",conditionMessage(e));NULL});add(lab,m,vars,nrow(d))}
for(x in EXTRA)run(paste0("plus_",sub("ext__","",x)),DD,c(V,x))
run("drop_sejong",DD[!grepl("세종",DD$fbd_region),],V)
run("drop_boundary_change",DD[!grepl("세종|청주|군위",DD$fbd_region),],V)
run("drop_2020_2021",DD[!DD$year%in%c(2020,2021),],V)
run("from_2015",DD[DD$year>=2015,],V)
for(p in names(phiprior))run(paste0("phi_prior_",p),DD,V,p)
print(R[R$cov%in%c("basic_livelihood_rate","walking_practice_rate")|grepl("^plus_",R$analysis)&R$cov%in%EXTRA,c("analysis","cov","IRR","lo","hi","WAIC")])
