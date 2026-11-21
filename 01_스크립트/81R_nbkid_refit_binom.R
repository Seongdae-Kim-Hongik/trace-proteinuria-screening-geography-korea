# 81R_nbkid_refit_binom.R - Reproducibility and likelihood comparison. (A) refits the principal model twice and records WAIC and coefficient differences;
#  (B) compares five likelihoods (Poisson, negative binomial, binomial, beta-binomial, binomial + observation iid) with the same covariates.
#  Usage: ST_VER=std3 ST_GRAPH=<graph> Rscript 81R_nbkid_refit_binom.R NBKID
#  Output: results_std3/S1b_refit_twice.csv, S1c_likelihood5.csv, S1c_likelihood5_coef.csv, S0_repro_meta.csv
suppressPackageStartupMessages({library(INLA);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE);TAG<-args[1];ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
VER<-Sys.getenv("ST_VER",unset="std3");SUF<-paste0("_",VER);OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF))
inla.setOption(num.threads="1:1");inla.setOption(inla.timeout=900)
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE)
SHP<-Sys.getenv("ST_SHP",unset="data/boundary/final.shp")
reg<-enc2utf8(sf::st_read(SHP,quiet=TRUE)$region);DD$fbd_region<-enc2utf8(DD$fbd_region)
DD$area<-match(DD$fbd_region,reg);stopifnot(!anyNA(DD$area));Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time;DD$obs<-seq_len(nrow(DD))
g<-inla.read.graph(file.path(RES,sprintf("nb_%s.graph",TAG)))
pc<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)),phi=list(prior="pc",param=c(0.5,2/3)));pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
V<-unique(c(S$sets_after_vif$S1_narrative_full,S$forced));bt<-function(v)sprintf("`%s`",v)
RE<-"f(area,model='bym2',graph=g,scale.model=TRUE,constr=TRUE,hyper=pc)+f(time,model='rw1',scale.model=TRUE,hyper=pcp)+f(yiid,model='iid',hyper=pcp)"
stopifnot(all(DD$cases<=DD$offset),all(DD$E>0))
DD$pE<-DD$E/DD$offset;stopifnot(all(DD$pE>0&DD$pE<1));DD$off<-qlogis(DD$pE)
fit<-function(fam,olre=FALSE){
  rhs<-paste(c("1",bt(V),RE,if(olre)"f(obs,model='iid',hyper=pcp)"),collapse="+")
  if(fam%in%c("poisson","nbinomial")){f<-as.formula(paste("cases ~",rhs));m<-inla(f,family=fam,data=DD,E=DD$E,control.fixed=CF,control.compute=list(waic=TRUE,dic=TRUE),control.inla=list(control.vb=list(emergency=50)))
    fv<-m$summary.fitted.values$mean;stopifnot(abs(sum(fv*DD$E)/sum(DD$cases)-1)<0.02)
  } else {f<-as.formula(paste("cases ~",rhs,"+ offset(off)"));m<-inla(f,family=fam,data=DD,Ntrials=DD$offset,control.fixed=CF,control.compute=list(waic=TRUE,dic=TRUE),control.inla=list(control.vb=list(emergency=50)))
    fv<-m$summary.fitted.values$mean;stopifnot(all(fv>0&fv<1))}
  m}
co<-function(m,lab){r<-m$summary.fixed[V,];data.frame(model=lab,cov=V,est=round(exp(r$mean),4),lo=round(exp(r$`0.025quant`),4),hi=round(exp(r$`0.975quant`),4),cred=(r$`0.025quant`>0|r$`0.975quant`<0))}
F0<-file.path(RES,"S0_repro_meta.csv")
write.csv(data.frame(INLA_version=as.character(packageVersion("INLA")),R_version=R.version.string,threads="1:1",seed_multiverse=20261019,seed_moran_resid=1,seed_moran_year=2,graph=Sys.getenv("ST_GRAPH"),n_obs=nrow(DD),date=format(Sys.time())),F0,row.names=FALSE)
F1<-file.path(RES,"S1b_refit_twice.csv")
if(!file.exists(F1)){a<-fit("nbinomial");b<-fit("nbinomial");ca<-co(a,"fit1");cb<-co(b,"fit2")
  write.csv(data.frame(WAIC1=round(a$waic$waic,2),WAIC2=round(b$waic$waic,2),dWAIC=round(b$waic$waic-a$waic$waic,3),max_abs_dlogcoef=signif(max(abs(log(ca$est)-log(cb$est))),3),
    phi1=round(a$summary.hyperpar["Phi for area","mean"],4),phi2=round(b$summary.hyperpar["Phi for area","mean"],4)),F1,row.names=FALSE);cat("refit done\n")}
F2<-file.path(RES,"S1c_likelihood5.csv")
if(!file.exists(F2)){L<-NULL;C<-NULL
  for(x in list(c("poisson","F"),c("nbinomial","F"),c("binomial","F"),c("betabinomial","F"),c("binomial","T"))){lab<-if(x[2]=="T")"binomial_olre" else x[1]
    m<-tryCatch(fit(x[1],x[2]=="T"),error=function(e){message(lab," failed: ",conditionMessage(e));NULL});if(is.null(m))next
    L<-rbind(L,data.frame(family=lab,effect=if(x[1]%in%c("poisson","nbinomial"))"prevalence ratio" else "odds ratio",WAIC=round(m$waic$waic,1),DIC=round(m$dic$dic,1)));C<-rbind(C,co(m,lab));cat(lab,"done\n")
    write.csv(L,F2,row.names=FALSE);write.csv(C,file.path(RES,"S1c_likelihood5_coef.csv"),row.names=FALSE)}}
print(read.csv(F1));print(read.csv(F2))
