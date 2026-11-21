# 88k_struct_kh_anem.R - Knorr-Held space-time interaction models (Types II-IV) for a 14-year panel with one district unobserved allowed.
#  Usage: TAG=<TAG> ST_VER=<version> Rscript 88k_struct_kh_anem.R kh. Uses kh_constraints_anem.R.
#  Output: results_<version>/S43_part_KH_typeII.csv, S43_part_KH_typeIII.csv, S43_part_KH_typeIV.csv
#                   maps: results_std/S20_obs_fit_resid.csv, figs_std/S20_obs_fit_resid_pooled.png · _windows.png
suppressPackageStartupMessages({library(INLA);library(sf);library(spdep);library(jsonlite);library(Matrix)})
args<-commandArgs(trailingOnly=TRUE);GROUP<-args[1];stopifnot(!is.na(GROUP))
TAG<-Sys.getenv("TAG",unset="NBKID");FOCAL<-Sys.getenv("FOCAL",unset="ext__walking_rate");VER<-Sys.getenv("ST_VER",unset="std3");SUF<-paste0("_",VER)
ROOT<-normalizePath(file.path(dirname(sub("--file=","",grep("--file=",commandArgs(),value=TRUE))),".."))
source(file.path(ROOT,"01_스크립트","kh_constraints_anem.R"))
SHP<-Sys.getenv("ST_SHP",unset="data/boundary/final.shp")
OUTD<-file.path(ROOT,"04_분석",TAG);RES<-file.path(OUTD,paste0("results",SUF));FIG<-file.path(OUTD,paste0("figs",SUF))
inla.setOption(num.threads="1:1")
LOG<-function(...)cat(sprintf("[88 %s %s] ",GROUP,format(Sys.time(),"%H:%M:%S")),sprintf(...),"\n",sep="")
S<-fromJSON(file.path(OUTD,sprintf("config_%s%s.json",TAG,SUF)),simplifyVector=TRUE)
INFO<-fromJSON(file.path(OUTD,sprintf("panel_INFO_%s%s.json",TAG,SUF)))
META<-read.csv(file.path(OUTD,sprintf("panel_META_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);LAB<-setNames(META$label_en,META$variable)
DD<-read.csv(file.path(OUTD,sprintf("panel_%s%s.csv",TAG,SUF)),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);DD$fbd_region<-enc2utf8(DD$fbd_region)
HB<-readLines(file.path(OUTD,sprintf("HM_best_%s%s.txt",TAG,SUF)));kv<-function(k)sub(paste0("^",k,"="),"",grep(paste0("^",k,"="),HB,value=TRUE)[1])
MAINFAM<-kv("main_family");BEST<-kv("best_set");FORCED<-S$forced;SETS<-lapply(S$sets_after_vif,as.character);HYPV<-SETS[[BEST]];MAINV<-unique(c(HYPV,FORCED))
if(file.exists(file.path(RES,"S20_panel_extra.csv"))){EX<-read.csv(file.path(RES,"S20_panel_extra.csv"),fileEncoding="UTF-8-BOM",stringsAsFactors=FALSE);EX$fbd_region<-enc2utf8(EX$fbd_region)
DD<-merge(DD,EX,by=c("year","fbd_region"),all.x=TRUE,sort=FALSE);stopifnot(nrow(DD)==INFO$n_rows)}
WIN<-INFO$windows

G<-st_read(SHP,quiet=TRUE);G$region<-enc2utf8(G$region);n<-nrow(G)
nb<-poly2nb(G,queen=TRUE);isl<-which(card(nb)==0);cen<-st_coordinates(st_centroid(st_geometry(G)))
if(length(isl)){kn<-knearneigh(cen,k=1)$nn;for(i in isl){j<-kn[i,1];nb[[i]]<-sort(unique(c(nb[[i]][nb[[i]]>0],j)));nb[[j]]<-sort(unique(c(nb[[j]][nb[[j]]>0],i)))}}
if(nzchar(Sys.getenv("ST_GRAPH"))){sg<-inla.read.graph(Sys.getenv("ST_GRAPH"));stopifnot(sg$n==nrow(G));for(i in seq_len(sg$n))nb[[i]]<-as.integer(sort(sg$nbs[[i]]));attr(nb,"sym")<-TRUE}
g<-inla.read.graph(file.path(RES,sprintf("nb_%s.graph",TAG)));stopifnot(g$n==n,all(g$nnbs==card(nb)))
DD<-DD[order(DD$year,match(DD$fbd_region,G$region)),];DD$area<-match(DD$fbd_region,G$region);Y0<-min(DD$year);DD$time<-DD$year-Y0+1;DD$yiid<-DD$time;DD$id<-seq_len(nrow(DD))
pcU<-function(u=0.5,a=0.01,pu=0.5,pa=2/3)list(prec=list(prior="pc.prec",param=c(u,a)),phi=list(prior="pc",param=c(pu,pa)))
pc<-pcU();pcp<-list(prec=list(prior="pc.prec",param=c(0.5,0.01)));if(nzchar(Sys.getenv("KH_INIT")))pcp$prec$initial<-as.numeric(Sys.getenv("KH_INIT"))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=0.001)
SPf<-function(gr="g",hyper="pc",extra="")sprintf("f(area,model='bym2',graph=%s,scale.model=TRUE,constr=TRUE,hyper=%s%s)",gr,hyper,extra)
TI<-"f(time,model='rw1',scale.model=TRUE,hyper=pcp)";IT<-"f(yiid,model='iid',hyper=pcp)"
bt<-function(v)sprintf("`%s`",v)
fit<-function(vars,data=DD,sp=SPf(),extra=character(0),ecol="E",robust=FALSE,env=parent.frame()){
  rhs<-c("1",bt(vars),extra,sp,TI,IT);rhs<-rhs[nzchar(rhs)];form<-as.formula(paste("cases ~",paste(rhs,collapse="+")),env=env)
  ci<-if(robust)list(control.vb=list(emergency=50),int.strategy="eb",strategy="gaussian") else list(control.vb=list(emergency=50))
  data$EOFF<-data[[ecol]]
  inla(form,family=MAINFAM,data=data,E=EOFF,control.fixed=CF,control.compute=list(dic=TRUE,waic=TRUE,cpo=FALSE,config=FALSE),control.inla=ci,verbose=FALSE)}
safe<-function(...){a<-list(...);tryCatch(fit(...),error=function(e){message("  fit failed: ",conditionMessage(e)," -> robust");tryCatch(fit(...,robust=TRUE),error=function(e2){message("  robust failed: ",conditionMessage(e2));NULL})})}
irr<-function(m,v){r<-m$summary.fixed[v,,drop=FALSE];data.frame(cov=v,IRR=round(exp(r$mean),3),lo=round(exp(r$`0.025quant`),3),hi=round(exp(r$`0.975quant`),3),cred=(r$`0.025quant`>0|r$`0.975quant`<0),stringsAsFactors=FALSE)}
hp<-function(m,pat){h<-m$summary.hyperpar;k<-grep(pat,rownames(h));if(length(k))round(h[k[1],"mean"],3) else NA}
save_part<-function(an,grp,m,vars,data,t0,note="",cmp="",added=character(0)){
  if(is.null(m)){write.csv(data.frame(analysis=an,group=grp,cov=NA,IRR=NA,lo=NA,hi=NA,cred=NA,WAIC=NA,note=paste("FIT FAILED",note)),file.path(RES,sprintf("S43_part_%s.csv",an)),row.names=FALSE);LOG("%s FAILED",an);return(invisible())}
  r<-do.call(rbind,lapply(vars,function(v)irr(m,v)));r$role<-ifelse(r$cov%in%HYPV,"hypothesis",ifelse(r$cov%in%FORCED,"forced","added"))
  out<-data.frame(analysis=an,group=grp,r,WAIC=round(m$waic$waic,1),DIC=round(m$dic$dic,1),n_obs=sum(!is.na(data$cases)),n_regions=length(unique(data$area[!is.na(data$cases)])),
    years=sprintf("%d-%d",min(data$year[!is.na(data$cases)]),max(data$year[!is.na(data$cases)])),phi=hp(m,"^Phi for area"),comparable_group=cmp,
    minutes=round(as.numeric(difftime(Sys.time(),t0,units="mins")),1),note=note,stringsAsFactors=FALSE)
  write.csv(out,file.path(RES,sprintf("S43_part_%s.csv",an)),row.names=FALSE,fileEncoding="UTF-8")
  LOG("%-22s WAIC %9.1f  cred hyp %2d/%d  walking %.3f (%.3f-%.3f)  · %.1f min",an,out$WAIC[1],sum(r$cred[r$role=="hypothesis"]),length(HYPV),r$IRR[r$cov==FOCAL],r$lo[r$cov==FOCAL],r$hi[r$cov==FOCAL],out$minutes[1])}
done<-function(an)file.exists(file.path(RES,sprintf("S43_part_%s.csv",an)))
run<-function(an,grp,expr_fun,vars=MAINV,data=DD,note="",cmp=""){if(done(an)){LOG("%s exists, skip",an);return(invisible())};t0<-Sys.time();m<-expr_fun();save_part(an,grp,m,vars,data,t0,note,cmp);invisible(m)}

if(GROUP=="nb"){
  m<-run("main_queen","1_neighbour",function()safe(MAINV),note="주모형 재적합(S3 와 같은 사양)",cmp="neighbour")
  if(!is.null(m)){mu<-m$summary.fitted.values$mean
    fx<-m$summary.fixed$mean;names(fx)<-gsub("`","",rownames(m$summary.fixed));X<-as.matrix(DD[,MAINV]);lpf<-fx["(Intercept)"]+X%*%fx[MAINV]
    write.csv(data.frame(year=DD$year,fbd_region=DD$fbd_region,window=DD$window,cases=DD$cases,E=DD$E,fitted_RR=mu,fixed_lp=as.numeric(lpf)),file.path(RES,"S20_fitted_main.csv"),row.names=FALSE,fileEncoding="UTF-8")}
  for(k in 3:5){an<-sprintf("knn%d",k);if(done(an))next
    nbk<-make.sym.nb(knn2nb(knearneigh(cen,k=k)));gf<-file.path(RES,sprintf("S20_nb_knn%d.graph",k));nb2INLA(gf,nbk);gk<-inla.read.graph(gf)
    run(an,"1_neighbour",function()safe(MAINV,sp=SPf("gk")),note=sprintf("KNN k=%d(중심점, 대칭화), 평균 이웃 %.2f, 부분그래프 %d",k,mean(card(nbk)),n.comp.nb(nbk)$nc),cmp="neighbour")}
}

if(GROUP=="struct"){
  run("icar_besag","2_structure",function()safe(MAINV,sp="f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=pcp)"),note="ICAR(besag)만",cmp="structure")
  run("iid_only","2_structure",function()safe(MAINV,sp="f(area,model='iid',hyper=pcp)"),note="공간 iid 만",cmp="structure")
  run("no_spatial","2_structure",function()safe(MAINV,sp=""),note="공간효과 없음(RW1+연도 iid 만)",cmp="structure")
  pcA<-pcU(1,0.01,0.5,0.5);pcB<-pcU(2,0.05,0.5,0.5);pcC<-pcU(0.5,0.01,0.5,0.5)
  run("pc_phi_flat","2_structure",function()safe(MAINV,sp=SPf(hyper="pcC")),note="BYM2 PC: prec(0.5,0.01) 유지, phi P(phi<0.5)=0.5(중립)",cmp="structure")
  run("pc_relaxA","2_structure",function()safe(MAINV,sp=SPf(hyper="pcA")),note="BYM2 PC 완화 A: prec P(sd>1)=0.01, phi P(phi<0.5)=0.5",cmp="structure")
  run("pc_relaxB","2_structure",function()safe(MAINV,sp=SPf(hyper="pcB")),note="BYM2 PC 완화 B: prec P(sd>2)=0.05, phi P(phi<0.5)=0.5",cmp="structure")
}

## ───────── 3. Knorr-Held II·III·IV ─────────
if(GROUP=="kh"){
  khfit<-function(type,d) kh_fit_repaired(type,d,MAINV)
  for(ty in c("II","III","IV")){an<-sprintf("KH_type%s",ty);if(done(an))next
    inla.setOption(inla.timeout=if(ty=="IV")5400 else 0)
    t0<-Sys.time();m<-khfit(ty,DD);note<-sprintf("Knorr-Held Type %s(BYM2+RW1+iid 위에 공간×시간 상호작용, 합 0 제약)",ty);d<-DD
    save_part(an,"3_spacetime",m,MAINV,d,t0,note,if(identical(d,DD))"spacetime" else "spacetime_2016")}
  inla.setOption(inla.timeout=0)
}


if(GROUP=="conf"){
  if(!done("RSR")){Xb<-apply(as.matrix(DD[,MAINV]),2,function(x)tapply(x,DD$area,mean))[as.character(1:n),];Xb<-scale(Xb,scale=FALSE)
    RSRC<<-list(A=cbind(t(Xb),matrix(0,ncol(Xb),n)),e=rep(0,ncol(Xb)))
    run("RSR","4_spatial_confounding",function()safe(MAINV,sp=SPf(extra=",extraconstr=RSRC"),env=.GlobalEnv),note=sprintf("RSR: BYM2 총효과 ⟂ 시군구 평균 공변량 %d개(제약)",ncol(Xb)),cmp="main_data")}
  miss<-function(an,cols){m<-setdiff(cols,names(DD));if(length(m))LOG("%s skipped: column(s) not in the panel: %s",an,paste(m,collapse=", "));length(m)>0}
  if(!miss("negctl_road","road_pavement_rate"))run("negctl_road","5_negcontrol",function()safe(c(MAINV,"road_pavement_rate")),vars=c(MAINV,"road_pavement_rate"),note="도로포장률(음성대조)을 통제로 추가",cmp="main_data")
  if(!miss("urban_quintile, urban_rw2",c("urban_q5","urban_rw2_bin"))){
  DD$urban_q5f<-factor(DD$urban_q5);for(q in 2:5)DD[[sprintf("urbQ%d",q)]]<-as.integer(DD$urban_q5==q)
  V0<-setdiff(MAINV,"urban_population_share")
  run("urban_quintile","5_negcontrol",function()safe(c(V0,sprintf("urbQ%d",2:5))),vars=c(V0,sprintf("urbQ%d",2:5)),note="도시지역 인구 비율 선형항 → 5분위 범주(기준 Q1)",cmp="main_data",data=DD)
  run("urban_rw2","5_negcontrol",function()safe(V0,extra="f(urban_rw2_bin,model='rw2',scale.model=TRUE,constr=TRUE,hyper=pcp)",data=DD),vars=V0,note="도시지역 인구 비율 선형항 → 20분위 rw2",cmp="main_data",data=DD)
  if(!miss("urban_q5_plus_road","road_pavement_rate"))run("urban_q5_plus_road","5_negcontrol",function()safe(c(V0,sprintf("urbQ%d",2:5),"road_pavement_rate")),vars=c(V0,sprintf("urbQ%d",2:5),"road_pavement_rate"),note="도시성 5분위 + 도로포장률",cmp="main_data",data=DD)
  }
  if(!miss("access_outside","outside_care_share"))run("access_outside","6_care_use",function()safe(c(MAINV,"outside_care_share")),vars=c(MAINV,"outside_care_share"),note="관외진료 비율(거주지 기준 전체 질환 관외/(관내+관외) 진료실인원, 1 SD) 추가",cmp="main_data")
  if(!miss("access_clinics","clinics_per_100k"))run("access_clinics","6_care_use",function()safe(c(MAINV,"clinics_per_100k")),vars=c(MAINV,"clinics_per_100k"),note="의원 수/10만(1 SD) 추가",cmp="main_data")
  if(!miss("access_both",c("outside_care_share","clinics_per_100k")))run("access_both","6_care_use",function()safe(c(MAINV,"outside_care_share","clinics_per_100k")),vars=c(MAINV,"outside_care_share","clinics_per_100k"),note="관외진료 비율 + 의원 밀도 추가",cmp="main_data")
}

if(GROUP=="maps"){suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(patchwork)})
  F<-read.csv(file.path(RES,"S20_fitted_main.csv"),fileEncoding="UTF-8",stringsAsFactors=FALSE);F$fbd_region<-enc2utf8(F$fbd_region);F$mu<-F$fitted_RR*F$E
  per<-c(pooled="2011-2024",setNames(sapply(names(WIN),function(k)sprintf("%d-%d",WIN[[k]][1],WIN[[k]][2])),names(WIN)))
  agg<-do.call(rbind,lapply(names(per),function(p){x<-if(p=="pooled")F else F[F$window==p,]
    x$pf<-exp(x$fixed_lp);x$pf<-x$pf*sum(x$cases)/sum(x$E*x$pf);x$Epf<-x$E*x$pf
    a<-x%>%group_by(fbd_region)%>%summarise(cases=sum(cases),Esum=sum(E),mu=sum(mu),Epf=sum(Epf),.groups="drop")%>%rename(E=Esum)
    data.frame(period=p,years=per[[p]],fbd_region=a$fbd_region,cases=a$cases,E=round(a$E,1),SIR_obs=round(a$cases/a$E,4),RR_fit=round(a$mu/a$E,4),RR_cov=round(a$Epf/a$E,4),
      resid_fit=round(a$cases/a$mu,4),resid_cov=round(a$cases/a$Epf,4),stringsAsFactors=FALSE)}))
  write.csv(agg,file.path(RES,"S20_obs_fit_resid.csv"),row.names=FALSE,fileEncoding="UTF-8")
  sm<-agg%>%group_by(period)%>%summarise(cor_obs_fit=round(cor(log(SIR_obs),log(RR_fit)),3),cor_obs_cov=round(cor(log(SIR_obs),log(RR_cov)),3),resid_fit_min=min(resid_fit),resid_fit_max=max(resid_fit),resid_cov_p10=round(quantile(resid_cov,.1),3),resid_cov_p90=round(quantile(resid_cov,.9),3),.groups="drop")
  write.csv(sm,file.path(RES,"S20_obs_fit_resid_summary.csv"),row.names=FALSE,fileEncoding="UTF-8");print(as.data.frame(sm))
  RRB<-c(-Inf,0.67,0.8,0.9,1,1.1,1.25,1.5,Inf);RRLAB<-c("<0.67","0.67-0.8","0.8-0.9","0.9-1","1-1.1","1.1-1.25","1.25-1.5",">=1.5")
  PAL8<-setNames(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","#FDDBC7","#F4A582","#D6604D","#B2182B"),RRLAB)
  RB<-c(-Inf,0.95,0.97,0.99,1,1.01,1.03,1.05,Inf);RLAB<-c("<0.95","0.95-0.97","0.97-0.99","0.99-1","1-1.01","1.01-1.03","1.03-1.05",">=1.05");PALR<-setNames(unname(PAL8),RLAB)
  G$sido<-substr(as.character(G$sigungu_cd),1,2);SIDO<-suppressWarnings(G%>%group_by(sido)%>%summarise(.groups="drop"))
  thm<-theme_void(base_size=10)+theme(strip.text=element_text(face="bold",size=9,hjust=0),legend.title=element_text(face="bold",size=8),legend.text=element_text(size=7),plot.margin=margin(4,4,4,4))
  lv<-function(SP,labs,col){d<-SP[seq_along(labs),];d[[col]]<-factor(labs,levels=labs);geom_sf(data=d,aes(fill=.data[[col]]),alpha=0,color=NA,show.legend=TRUE,inherit.aes=FALSE)}
  mk<-function(pp){SP<-merge(G[,c("region","geometry")],agg[agg$period%in%pp,],by.x="region",by.y="fbd_region")
    L<-data.frame(var=c("SIR_obs","RR_fit","RR_cov","resid_fit"),lab=c("Observed SIR","Fitted RR (full model)","Covariate-predicted RR (fixed effects)","Residual ratio, observed / fitted"),stringsAsFactors=FALSE)
    plots<-list()
    for(p in pp)for(i in 1:4){v<-L$var[i];S2<-SP[SP$period==p,];res<-v=="resid_fit"
      S2$cl<-factor(cut(S2[[v]],if(res)RB else RRB,labels=if(res)RLAB else RRLAB),levels=if(res)RLAB else RRLAB);S2$panel<-sprintf("%s, %s",L$lab[i],per[[p]])
      plots[[length(plots)+1]]<-ggplot(S2)+geom_sf(aes(fill=cl),color="grey70",linewidth=0.05)+lv(S2,if(res)RLAB else RRLAB,"cl")+geom_sf(data=SIDO,fill=NA,color="grey35",linewidth=0.15,inherit.aes=FALSE)+
        facet_wrap(~panel)+scale_fill_manual(values=if(res)PALR else PAL8,breaks=rev(if(res)RLAB else RRLAB),limits=if(res)RLAB else RRLAB,drop=FALSE,name=if(res)"Observed /\nfitted" else "SIR or RR")+
        guides(fill=guide_legend(nrow=2,override.aes=list(alpha=1,color="grey40")))+thm+theme(legend.position="bottom")}
    wrap_plots(plots,ncol=4)}
  ggsave(file.path(FIG,"S20_obs_fit_resid_pooled.png"),mk("pooled"),width=15,height=5.4,dpi=300,bg="white")
  ggsave(file.path(FIG,"S20_obs_fit_resid_windows.png"),mk(names(WIN)),width=15,height=15,dpi=300,bg="white")
  LOG("maps saved")
}

## ───────── combine ─────────
if(GROUP=="review"){
  F50<-file.path(RES,"S50_resid_moran_by_spec.csv");F50p<-file.path(RES,"S50_parsimonious.csv");F50r<-file.path(RES,"S50_region_slopes.csv")
  OBSR<-G$region%in%unique(DD$fbd_region[!is.na(DD$cases)])
  moran_obs<-function(x,regions){k<-!is.na(x);if(sum(k)<10)return(list(statistic=NA,p.value=NA));kk<-G$region%in%regions[k];nbk<-subset(nb,kk);lwk<-nb2listw(nbk,style="W",zero.policy=TRUE);moran.mc(x[k][match(G$region[kk],regions[k])],lwk,999,zero.policy=TRUE,alternative="two.sided")}
  rmoran<-function(m){mu<-m$summary.fitted.values$mean*DD$E;sz<-if(MAINFAM=="nbinomial")m$summary.hyperpar[grep("size",rownames(m$summary.hyperpar))[1],"mean"] else Inf
    pr<-(DD$cases-mu)/sqrt(mu+mu^2/sz);pra<-tapply(pr,DD$area,mean);set.seed(1);mm<-moran_obs(as.numeric(pra[as.character(which(OBSR))]),G$region[OBSR]);c(I=round(unname(mm$statistic),4),p=round(mm$p.value,4))}
  S7s<-read.csv(file.path(RES,"S7_multiverse_summary.csv"),stringsAsFactors=FALSE);ROBV<-S7s$cov[S7s$robust%in%c(TRUE,"TRUE","True")]
  if(nzchar(Sys.getenv("ROBV_OVERRIDE")))ROBV<-strsplit(Sys.getenv("ROBV_OVERRIDE"),",")[[1]]
  LOG("robust: %s",paste(ROBV,collapse=", "))
  specs<-list(principal=list(v=MAINV,sp=SPf()),forced_only=list(v=FORCED,sp=SPf()),parsimonious=list(v=unique(c(ROBV,FORCED)),sp=SPf()),
              no_spatial=list(v=MAINV,sp=""),iid_district=list(v=MAINV,sp="f(area,model='iid',hyper=pcp)"),icar_district=list(v=MAINV,sp="f(area,model='besag',graph=g,scale.model=TRUE,constr=TRUE,hyper=pcp)"),
              typeI_interaction=list(v=MAINV,sp=SPf(),ex="f(id,model='iid',hyper=pcp)"))
  if(!file.exists(F50)){rows<-list();for(nm in names(specs)){s<-specs[[nm]];t0<-Sys.time();m<-safe(s$v,sp=s$sp,extra=if(is.null(s$ex))character(0) else s$ex)
      if(is.null(m)){rows[[nm]]<-data.frame(spec=nm,n_cov=length(s$v),WAIC=NA,resid_moran_I=NA,p_perm=NA,alternative="two.sided",expected_I=round(-1/(sum(OBSR)-1),4));next}
      rm<-rmoran(m);rows[[nm]]<-data.frame(spec=nm,n_cov=length(s$v),WAIC=round(m$waic$waic,1),resid_moran_I=rm["I"],p_perm=rm["p"],alternative="two.sided",expected_I=round(-1/(sum(OBSR)-1),4))
      if(nm=="parsimonious"){r<-do.call(rbind,lapply(s$v,function(v)irr(m,v)));r$role<-ifelse(r$cov%in%FORCED,"forced","robust");r$WAIC<-round(m$waic$waic,1);write.csv(r,F50p,row.names=FALSE,fileEncoding="UTF-8")}
      LOG("%-18s WAIC %9.1f  resid Moran %.3f (p %.3f)  %.1f min",nm,rows[[nm]]$WAIC,rm["I"],rm["p"],as.numeric(difftime(Sys.time(),t0,units="mins")))}
    write.csv(do.call(rbind,rows),F50,row.names=FALSE,fileEncoding="UTF-8")}
  if(!file.exists(F50r)){
    RG<-c("서울시"="Capital","인천시"="Capital","경기도"="Capital","대전시"="Chungcheong","세종시"="Chungcheong","충청북도"="Chungcheong","충청남도"="Chungcheong","광주시"="Honam","전라북도"="Honam","전라남도"="Honam",
          "부산시"="Yeongnam","대구시"="Yeongnam","울산시"="Yeongnam","경상북도"="Yeongnam","경상남도"="Yeongnam","강원도"="Gangwon-Jeju","제주도"="Gangwon-Jeju")
    DD$rg<-unname(RG[enc2utf8(DD$sido_std)]);if(anyNA(DD$rg))stop("권역 미매칭 시도: ",paste(unique(DD$sido_std[is.na(DD$rg)]),collapse=","))
    m0<-safe(MAINV);w0<-m0$waic$waic;out<-list()
    for(v in ROBV){D2<-DD;nv<-character(0);for(r in unique(DD$rg)){cn<-paste0(v,"__",gsub("-","_",r));D2[[cn]]<-D2[[v]]*(D2$rg==r);nv<-c(nv,cn)}
      m<-safe(c(setdiff(MAINV,v),nv),data=D2);if(is.null(m))next
      r<-do.call(rbind,lapply(nv,function(x)irr(m,x)));r$covariate<-v;r$region<-sub(paste0("^",v,"__"),"",r$cov);r$n_district_years<-sapply(r$region,function(x)sum(gsub("-","_",DD$rg)==x))
      r$dWAIC_vs_common<-round(m$waic$waic-w0,1);r$common_IRR<-irr(m0,v)$IRR;out[[v]]<-r;LOG("region slopes %s dWAIC %.1f",v,r$dWAIC_vs_common[1])}
    write.csv(do.call(rbind,out),F50r,row.names=FALSE,fileEncoding="UTF-8")}
  F50a<-file.path(RES,"S50_access_sensitivity.csv")
  if(!file.exists(F50a)&&"outside_care_share"%in%names(DD)){m<-safe(c(MAINV,"outside_care_share"))
    if(!is.null(m)){r<-do.call(rbind,lapply(c(MAINV,"outside_care_share"),function(v)irr(m,v)));r$role<-ifelse(r$cov%in%HYPV,"hypothesis",ifelse(r$cov%in%FORCED,"forced","added"));r$WAIC<-round(m$waic$waic,1)
      m0a<-safe(MAINV);r$main_IRR<-sapply(r$cov,function(v)if(v%in%MAINV)irr(m0a,v)$IRR[1] else NA);write.csv(r,F50a,row.names=FALSE,fileEncoding="UTF-8");LOG("access sensitivity WAIC %.1f",r$WAIC[1])}}
}

if(GROUP=="combine"){
  fs<-list.files(RES,pattern="^S43_part_.*\\.csv$",full.names=TRUE);P<-do.call(rbind,lapply(fs,function(f)read.csv(f,fileEncoding="UTF-8",stringsAsFactors=FALSE)))
  ad<-function(f,an,grp,note,cmp){x<-read.csv(file.path(RES,f),stringsAsFactors=FALSE);data.frame(analysis=an,group=grp,cov=x$cov,IRR=x$IRR,lo=x$lo,hi=x$hi,cred=x$cred,
    role=ifelse(x$cov%in%HYPV,"hypothesis",ifelse(x$cov%in%FORCED,"forced","added")),WAIC=x$WAIC,DIC=NA,n_obs=INFO$n_rows,n_regions=INFO$n_sgg,years="2011-2024",phi=NA,comparable_group=cmp,minutes=NA,note=note,stringsAsFactors=FALSE)}
  P<-rbind(P,ad("S9_typeI_IRR.csv","KH_typeI_S9","3_spacetime","Knorr-Held Type I(기존 S9, 81)","spacetime"),ad("S10_offset_sensitivity_IRR.csv","offset_sensitivity_S10","8_offset",sprintf("S10 (81) offset %s",INFO$offset_sensitivity),"E_alt"))
  P$cmp2<-ifelse(P$comparable_group%in%c("neighbour","structure","spacetime","main_data"),"main_data",P$comparable_group)
  w<-unique(P[,c("analysis","comparable_group","cmp2","WAIC")]);mainW<-w$WAIC[w$analysis=="main_queen"]
  P$dWAIC_vs_main<-ifelse(P$cmp2=="main_data",round(P$WAIC-mainW,1),NA)
  P$min_WAIC_neighbour<-P$comparable_group=="neighbour"&P$WAIC==min(w$WAIC[w$comparable_group=="neighbour"],na.rm=TRUE)
  sw<-w[w$comparable_group=="structure"|w$analysis=="main_queen",];P$min_WAIC_structure<-P$analysis==sw$analysis[which.min(sw$WAIC)]
  M<-read.csv(file.path(RES,"S3_main_IRR.csv"),stringsAsFactors=FALSE);mm<-setNames(log(M$IRR),M$cov);mc<-setNames(M$cred,M$cov)
  P$same_sign_as_main<-ifelse(P$cov%in%M$cov,sign(log(P$IRR))==sign(mm[P$cov])|P$IRR==1,NA);P$main_cred<-mc[P$cov]
  ROB<-read.csv(file.path(RES,"S7_multiverse_summary.csv"),stringsAsFactors=FALSE);ROBV<-ROB$cov[ROB$robust];P$robust8<-P$cov%in%ROBV
  P<-P[order(P$group,P$analysis,match(P$role,c("hypothesis","forced","added")),match(P$cov,c(HYPV,FORCED))),]
  P$cmp2<-NULL;write.csv(P,file.path(RES,"S43_sensitivity_IRR.csv"),row.names=FALSE,fileEncoding="UTF-8")
  sm<-do.call(rbind,lapply(split(P,P$analysis),function(x){h<-x[x$role=="hypothesis",];r<-h[h$robust8,];wk<-h[h$cov==FOCAL,]
    data.frame(analysis=x$analysis[1],group=x$group[1],WAIC=x$WAIC[1],dWAIC_vs_main=x$dWAIC_vs_main[1],comparable_group=x$comparable_group[1],n_obs=x$n_obs[1],years=x$years[1],phi=x$phi[1],
      n_cred_hyp=sum(h$cred,na.rm=TRUE),robust8_cred_same_sign=sum(r$cred&r$same_sign_as_main,na.rm=TRUE),robust8_lost=paste(r$cov[!(r$cred&r$same_sign_as_main)],collapse=";"),
      walking=if(nrow(wk))sprintf("%.3f (%.3f-%.3f)%s",wk$IRR,wk$lo,wk$hi,if(isTRUE(wk$cred))"*" else "") else NA,
      sign_flips_vs_main=paste(h$cov[!is.na(h$same_sign_as_main)&!h$same_sign_as_main&h$cred],collapse=";"),
      min_WAIC_neighbour=any(x$min_WAIC_neighbour),min_WAIC_structure=any(x$min_WAIC_structure),note=x$note[1],stringsAsFactors=FALSE)}))
  sm<-sm[order(sm$group,sm$analysis),];write.csv(sm,file.path(RES,"S43_sensitivity_fit.csv"),row.names=FALSE,fileEncoding="UTF-8");print(sm[,c("analysis","WAIC","dWAIC_vs_main","n_cred_hyp","robust8_cred_same_sign","walking","robust8_lost")])
}
LOG("done")
