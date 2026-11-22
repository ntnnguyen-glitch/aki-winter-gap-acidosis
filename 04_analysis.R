suppressMessages({library(data.table);library(survival);library(splines);library(sandwich);library(WeightIt);library(MatchIt);library(mice);library(EValue)})
set.seed(2026); M_IMP<-40
D<-as.data.table(readRDS("D2.rds"))
prep<-function(d){d<-copy(d)
 d[,hd28:=as.integer(!is.na(hdeath_h)&hdeath_h<=672)]; d[,t28:=pmax(pmin(fifelse(is.na(hdeath_h),disch_h,hdeath_h),672),0.5)/24]
 d[,d28:=as.integer(!is.na(death_h)&death_h<=672)]; d[,d90:=as.integer(!is.na(death_h)&death_h<=2160)]; d[,t90:=pmax(pmin(fifelse(is.na(death_h),2160,death_h),2160),0.5)/24]; d[,d365:=as.integer(!is.na(death_h)&death_h<=8760)]; d[,t365:=pmax(pmin(fifelse(is.na(death_h),8760,death_h),8760),0.5)/24]
 d[,c_death:=as.integer(!is.na(hdeath_h)&hdeath_h<=720)]; d[,c_rrt:=as.integer(!is.na(rrt_h)&rrt_h<=720)]; d[,c_prd:=fifelse(is.na(cr_last),NA_integer_,as.integer(cr_last>=2*base))]
 d[,make30:=as.integer(c_death==1|c_rrt==1|(!is.na(c_prd)&c_prd==1))]; d[,int48:=as.integer(!is.na(intub_h)&intub_h<=48)]
 d[,gpos:=pmax(gap,0)/5]; d[,gneg:=pmax(-gap,0)/5]; d[,dbf:=factor(db)]; d[,hosp:=as.character(hosp)]; d}
P0<-prep(D[supp=="None"])
M2v<-c("age","male","chf","ckd","copd","dm","liver","cancer","sepsis")
M3v<-c(M2v,"sofa_nr","lact","hco3","alb","bun","k","cl","glu","hb","wbc","plt","map","hr","vaso","stage_idx")
impv<-c("lact","alb","bun","k","cl","glu","hb","wbc","plt","map","hr","gcs","pf")
mkimp<-function(d,m){x<-d[,c("hd28","make30","int48","gap","dbf",M3v,"gcs","pf"),with=FALSE]; x<-as.data.frame(x); meth<-make.method(x); meth[!names(x)%in%impv]<-""
 im<-mice(x,m=m,method=meth,maxit=10,printFlag=FALSE,seed=2026); lapply(1:m,function(i){z<-copy(d); cc<-complete(im,i); for(v in impv) set(z,j=v,value=cc[[v]]); z})}
IMP<-mkimp(P0,M_IMP); cat("imputed\n")
rub<-function(est,var){m<-length(est); q<-mean(est); U<-mean(var); Bv<-if(m>1) var(est) else 0; Tt<-U+(1+1/m)*Bv; r<-(1+1/m)*Bv/U; df<-if(Bv>0) (m-1)*(1+1/r)^2 else 1e6; se<-sqrt(Tt); p<-2*pt(-abs(q/se),df); c(est=exp(q),lo=exp(q-qt(.975,df)*se),hi=exp(q+qt(.975,df)*se),p=p)}
multi<-function(d) length(unique(d$hosp))>1
VC<-function(f,d,clus=TRUE){ if(clus&&multi(d)) vcovCL(f,cluster=d$hosp) else if(inherits(f,"coxph")) vcov(f) else vcovHC(f,type="HC0")}
fo<-function(y,v,pool,x="grp") as.formula(paste(y,"~",paste(c(x,v,if(pool) "dbf"),collapse="+")))
gl<-function(dl,y,v,x="grp",pat="^grp"){pool<-length(unique(dl[[1]]$db))>1; res<-lapply(dl,function(d){f<-glm(fo(y,v,pool,x),binomial,d); V<-VC(f,d); i<-grep(pat,names(coef(f))); list(b=coef(f)[i],v=diag(V)[i])})
 nm<-names(res[[1]]$b); t(sapply(seq_along(nm),function(j) rub(sapply(res,function(r) r$b[j]),sapply(res,function(r) r$v[j])))) |> `rownames<-`(nm)}
cx<-function(dl,tt,ev,v,x="grp",pat="^grp"){pool<-length(unique(dl[[1]]$db))>1; res<-lapply(dl,function(d){f<-as.formula(paste0("Surv(",tt,",",ev,")~",paste(c(x,v,if(pool) "strata(dbf)"),collapse="+")))
 fit<-if(multi(d)) coxph(f,data=d,cluster=hosp) else coxph(f,data=d); i<-grep(pat,names(coef(fit))); list(b=coef(fit)[i],v=diag(vcov(fit))[i])})
 nm<-names(res[[1]]$b); t(sapply(seq_along(nm),function(j) rub(sapply(res,function(r) r$b[j]),sapply(res,function(r) r$v[j])))) |> `rownames<-`(nm)}
fgr<-function(dl,tcol,v,hor){pool<-length(unique(dl[[1]]$db))>1; res<-lapply(dl,function(d){d<-copy(d); te<-fifelse(is.na(d[[tcol]]),Inf,d[[tcol]]); td<-fifelse(is.na(d$hdeath_h),Inf,d$hdeath_h); ft<-pmin(te,td,d$disch_h,hor); ft<-pmax(ft,0.01)
 st<-fifelse(te<=ft&te<=hor,1L,fifelse(td<=ft&td<=hor,2L,0L)); d[,`:=`(ft=ft,st=factor(st,0:2,c("cens","ev","death")))]; fg<-finegray(Surv(ft,st)~.,data=d[,c("ft","st","grp",v,"dbf","hosp"),with=FALSE],etype="ev")
 f<-as.formula(paste("Surv(fgstart,fgstop,fgstatus)~",paste(c("grp",v,if(pool) "strata(dbf)"),collapse="+"))); fit<-if(multi(d)) coxph(f,data=fg,weights=fgwt,cluster=hosp) else coxph(f,data=fg,weights=fgwt); i<-grep("^grp",names(coef(fit))); list(b=coef(fit)[i],v=diag(vcov(fit))[i],n=sum(st=="ev"))})
 nm<-names(res[[1]]$b); out<-t(sapply(seq_along(nm),function(j) rub(sapply(res,function(r) r$b[j]),sapply(res,function(r) r$v[j])))); rownames(out)<-nm; attr(out,"nev")<-res[[1]]$n; out}
prepfg<-function(d){d<-copy(d); dh<-fifelse(is.na(d$hdeath_h),Inf,d$hdeath_h); ft<-pmax(pmin(dh,d$disch_h,672),0.5); st<-fifelse(dh<=ft&dh<=672,1L,fifelse(d$disch_h<=ft&d$disch_h<672,2L,0L)); d[,`:=`(ftd=ft/24,std=factor(st,0:2,c("cens","death","disch")))]; d}
fgfit<-function(d,x,v,pool,w=NULL){d<-prepfg(d); if(!is.null(w)) d[,ww:=w] else d[,ww:=1]; d[,rid:=.I]; cols<-unique(c("ftd","std",x,v,"dbf","hosp","ww","rid")); fg<-finegray(Surv(ftd,std)~.,data=d[,cols,with=FALSE],etype="death"); f<-as.formula(paste("Surv(fgstart,fgstop,fgstatus)~",paste(c(x,v,if(pool) "strata(dbf)"),collapse="+"))); fg$wt<-fg$fgwt*fg$ww; if(length(unique(d$hosp))>1) coxph(f,data=fg,weights=wt,cluster=hosp) else coxph(f,data=fg,weights=wt,cluster=rid)}
fgd<-function(dl,v,x="grp",pat="^grp"){pool<-length(unique(dl[[1]]$db))>1; res<-lapply(dl,function(d){fit<-fgfit(d,x,v,pool); i<-grep(pat,names(coef(fit))); list(b=coef(fit)[i],v=diag(vcov(fit))[i])})
 nm<-names(res[[1]]$b); t(sapply(seq_along(nm),function(j) rub(sapply(res,function(r) r$b[j]),sapply(res,function(r) r$v[j])))) |> `rownames<-`(nm)}
psm_iptw<-function(dl){pool<-length(unique(dl[[1]]$db))>1; cv<-c(M3v,if(pool) "dbf"); fps<-as.formula(paste("grp~",paste(cv,collapse="+")))
 R<-lapply(dl,function(d){d<-copy(d); w<-weightit(fps,d,method="glm",estimand="ATE"); d[,w:=pmin(w$weights,quantile(w$weights,.99))]; wo<-weightit(fps,d,method="glm",estimand="ATO"); d[,wo:=wo$weights]
  f1<-glm(hd28~grp,quasibinomial,d,weights=w); V1<-VC(f1,d); f2<-glm(hd28~grp,quasibinomial,d,weights=wo); V2<-VC(f2,d)
  c1<-fgfit(d,"grp",NULL,FALSE,d$w)
  d2<-d[grp!="Resp. alkalosis"]; d2[,tr:=as.integer(grp=="Resp. acidosis")]; mt<-matchit(as.formula(paste("tr~",paste(cv,collapse="+"))),d2,method="nearest",caliper=0.2); md<-as.data.table(match.data(mt))
  f3<-glm(hd28~tr,binomial,md,weights=weights); V3<-vcovCL(f3,cluster=md$subclass)
  list(ip=list(b=coef(f1)[2:3],v=diag(V1)[2:3]),ow=list(b=coef(f2)[2:3],v=diag(V2)[2:3]),iph=list(b=coef(c1)[1:2],v=diag(vcov(c1))[1:2]),ps=list(b=coef(f3)[2],v=V3[2,2]),np=nrow(md)/2)})
 pk<-function(nm,j) rub(sapply(R,function(r) r[[nm]]$b[j]),sapply(R,function(r) r[[nm]]$v[j]))
 list(iptw=rbind(acid=pk("ip",1),alk=pk("ip",2)),ow=rbind(acid=pk("ow",1),alk=pk("ow",2)),iptwhr=rbind(acid=pk("iph",1),alk=pk("iph",2)),psm=c(pk("ps",1),pairs=round(mean(sapply(R,function(r) r$np)))))}
rcs<-function(dl){pool<-length(unique(dl[[1]]$db))>1; Bm<-ns(dl[[1]]$gap,df=4); gg<-seq(-15,25,0.25); C<-predict(Bm,gg)-matrix(predict(Bm,0),length(gg),4,byrow=TRUE)
 R<-lapply(dl,function(d){d<-copy(d); X<-predict(Bm,d$gap); for(j in 1:4) set(d,j=paste0("b",j),value=X[,j]); f<-glm(fo("hd28",M3v,pool,paste0("b",1:4)),binomial,d); f0<-glm(fo("hd28",M3v,pool,"1"),binomial,d); fl<-glm(fo("hd28",M3v,pool,"gap"),binomial,d)
  i<-paste0("b",1:4); V<-VC(f,d)[i,i]; list(f=f,f0=f0,fl=fl,lp=as.vector(C%*%coef(f)[i]),vv=rowSums((C%*%V)*C))})
 po<-D1(as.mira(lapply(R,`[[`,"f")),as.mira(lapply(R,`[[`,"f0"))); pn<-D1(as.mira(lapply(R,`[[`,"f")),as.mira(lapply(R,`[[`,"fl")))
 cur<-t(sapply(seq_along(gg),function(k) rub(sapply(R,function(r) r$lp[k]),sapply(R,function(r) r$vv[k])))); list(p=c(po$result[1,"P(>F)"],pn$result[1,"P(>F)"]),curve=data.table(gap=gg,or=cur[,1],lo=cur[,2],hi=cur[,3]))}
run<-function(dl){o<-list(); d1<-dl[[1]]; o$n<-nrow(d1); o$t1<-d1[,.(n=.N,ev=sum(hd28)),keyby=grp]
 o$crude<-gl(dl,"hd28",NULL); o$m2<-gl(dl,"hd28",M2v); o$m3<-gl(dl,"hd28",M3v); o$hr<-fgd(dl,M3v)
 o$pw<-gl(dl,"hd28",M3v,c("gpos","gneg"),"^gpos|^gneg"); o$pwhr<-fgd(dl,M3v,c("gpos","gneg"),"^gpos|^gneg")
 o$km<-{q<-prepfg(d1); cmprsk::cuminc(q$ftd,as.integer(q$std)-1L,q$grp,cencode=0)$Tests[1,"pv"]}; o$rcs<-rcs(dl); o$ps<-psm_iptw(dl)
 o$ev<-evalues.OR(o$m3[1,1],o$m3[1,2],o$m3[1,3],rare=FALSE)
 o$make<-gl(dl,"make30",M3v); o$c_death<-gl(dl,"c_death",M3v); o$c_rrt<-gl(dl,"c_rrt",M3v); o$c_prd<-gl(lapply(dl,function(d) d[!is.na(c_prd)]),"c_prd",M3v)
 o$int<-fgr(dl,"intub_h",M3v,48); o$rrt7<-fgr(dl,"rrt_h",M3v,168)
 o$makeparts<-d1[,.(make=sum(make30),death=sum(c_death),rrt=sum(c_rrt),prd=sum(c_prd,na.rm=TRUE),prd_unknown=sum(is.na(c_prd)&c_death==0&c_rrt==0),cr_last_h_med=median(cr_last_h,na.rm=TRUE)/24),by=grp]
 if(all(d1$db=="MIMIC-IV")){o$d28<-cx(dl,"t28m","d28",M3v); o$d90<-cx(dl,"t90","d90",M3v); o$d365<-cx(dl,"t365","d365",M3v)}
 o}
for(i in seq_along(IMP)) IMP[[i]][,t28m:=pmax(pmin(fifelse(is.na(death_h),672,death_h),672),0.5)/24]
sets<-list("MIMIC-IV"=lapply(IMP,function(d) d[db=="MIMIC-IV"]),"eICU-CRD"=lapply(IMP,function(d) d[db=="eICU-CRD"]),"Pooled"=IMP)
R<-list(); for(k in names(sets)){R[[k]]<-run(sets[[k]]); cat(k,"\n"); print(round(R[[k]]$m3,3))}
# ---- sensitivity (pooled and per db) ----
S<-list(); for(k in names(sets)){dl<-sets[[k]]; d0<-P0[db%in%unique(dl[[1]]$db)]; s<-list()
 cc<-d0[complete.cases(d0[,c(M3v),with=FALSE])]; s[["Complete-case analysis"]]<-c(gl(list(cc),"hd28",M3v)[1,],n=nrow(cc))
 md<-copy(d0); for(v in impv) md[,(v):=fifelse(is.na(get(v)),median(get(v),na.rm=TRUE),get(v)),by=db]; s[["Median imputation"]]<-c(gl(list(md),"hd28",M3v)[1,],n=nrow(md))
 s[["Model 3 + Glasgow Coma Scale"]]<-c(gl(dl,"hd28",c(M3v,"gcs"))[1,],n=nrow(dl[[1]]))
 s[["Model 3 + Glasgow Coma Scale + PaO2/FiO2"]]<-c(gl(dl,"hd28",c(M3v,"gcs","pf"))[1,],n=nrow(dl[[1]]))
 sub<-function(nm,f,v=M3v){dd<-lapply(dl,function(d) d[eval(f,d)]); s[[nm]]<<-c(tryCatch(gl(dd,"hd28",v)[1,],error=function(e) rep(NA,4)),n=nrow(dd[[1]]))}
 sub("AKI established at or before the index gas",quote(aki_pre==1)); sub("Glasgow Coma Scale >8",quote(!is.na(gcs)&gcs>8)); sub("pH <7.35",quote(ph<7.35)); sub("Blood gas bicarbonate only",quote(src=="ABG"))
 sub("Deaths within 24 h excluded",quote(is.na(hdeath_h)|hdeath_h>=24)); sub("Chronic kidney disease excluded",quote(ckd==0),setdiff(M3v,"ckd")); sub("Diabetic ketoacidosis excluded",quote(dka==0))
 d5<-lapply(dl,function(d){d<-copy(d); d[,grp:=factor(fifelse(gap< -5,"Resp. alkalosis",fifelse(gap>5,"Resp. acidosis","Adequate")),levels=levels(d$grp))]; d}); s[["Adequate band -5 to +5 mmHg"]]<-c(gl(d5,"hd28",M3v)[1,],n=nrow(dl[[1]]))
 # adding NIV / all supported patients (median-imputed, cluster SE)
 a<-prep(D[db%in%unique(dl[[1]]$db)]); for(v in impv) a[,(v):=fifelse(is.na(get(v)),median(get(v),na.rm=TRUE),get(v)),by=db]
 s[["NIV patients added (NIV as covariate)"]]<-c(gl(list(a[supp!="IMV"]),"hd28",c(M3v,"niv"))[1,],n=nrow(a[supp!="IMV"]))
 s[["All patients incl. invasive ventilation (support as covariate)"]]<-c(gl(list(a),"hd28",c(M3v,"imv","niv"))[1,],n=nrow(a))
 S[[k]]<-s; cat("sens",k,"\n")}
# ---- subgroups (pooled) ----
sgdef<-list("Age >=65 years"=quote(age>=65),"Male sex"=quote(male==1),"Sepsis"=quote(sepsis==1),"KDIGO stage 3"=quote(stage==3),"Chronic kidney disease"=quote(ckd==1),"Lactate >2 mmol/L"=quote(lact>2),"Vasopressor at index"=quote(vaso==1),"COPD"=quote(copd==1),"eICU-CRD (vs MIMIC-IV)"=quote(db=="eICU-CRD"))
rmv<-list("Age >=65 years"="age","Male sex"="male","Sepsis"="sepsis","KDIGO stage 3"="stage_idx","Chronic kidney disease"="ckd","Lactate >2 mmol/L"="lact","Vasopressor at index"="vaso","COPD"="copd","eICU-CRD (vs MIMIC-IV)"=character(0)); SG<-list()
for(s in names(sgdef)){dl<-IMP; f<-sgdef[[s]]; vr<-setdiff(M3v,rmv[[s]]); y<-lapply(dl,function(d) d[eval(f,d)]); nn<-lapply(dl,function(d) d[!eval(f,d)])
 r1<-gl(y,"hd28",vr); r0<-gl(nn,"hd28",vr)
 fa<-lapply(dl,function(d){d<-copy(d); d[,z:=as.integer(eval(f,d))]; glm(fo("hd28",c(vr,"z"),TRUE),binomial,d)}); fb<-lapply(dl,function(d){d<-copy(d); d[,z:=as.integer(eval(f,d))]; glm(fo("hd28",c(vr,"z","grp:z"),TRUE),binomial,d)})
 pint<-D1(as.mira(fb),as.mira(fa))$result[1,"P(>F)"]; SG[[s]]<-list(yes=r1,no=r0,nyes=nrow(y[[1]]),nno=nrow(nn[[1]]),pint=pint)}
# balance (imputation 1, pooled)
d<-IMP[[1]]; cv<-c(M3v,"dbf"); w<-weightit(as.formula(paste("grp~",paste(cv,collapse="+"))),d,method="glm",estimand="ATE"); d[,w:=pmin(w$weights,quantile(w$weights,.99))]; d[,eicu:=as.integer(db=="eICU-CRD")]
smd<-function(x,g,wt,a,b){i<-g==a;j<-g==b; abs(weighted.mean(x[i],wt[i])-weighted.mean(x[j],wt[j]))/sqrt((var(x[i])+var(x[j]))/2)}
BT<-rbindlist(lapply(c(M3v,"eicu"),function(v){x<-d[[v]];g<-d$grp;data.table(var=v,acid_un=smd(x,g,rep(1,nrow(d)),"Resp. acidosis","Adequate"),acid_w=smd(x,g,d$w,"Resp. acidosis","Adequate"),alk_un=smd(x,g,rep(1,nrow(d)),"Resp. alkalosis","Adequate"),alk_w=smd(x,g,d$w,"Resp. alkalosis","Adequate"))}))
miss<-P0[,lapply(.SD,function(x) round(100*mean(is.na(x)),1)),by=db,.SDcols=impv]
recl<-NULL; if(file.exists("D_v1.rds")){v1<-as.data.table(readRDS("D_v1.rds"))[imv==0]; recl<-data.table(db=c("MIMIC-IV","eICU-CRD"),old=c(v1[db=="MIMIC-IV",.N],v1[db=="eICU-CRD",.N]),new=c(P0[db=="MIMIC-IV",.N],P0[db=="eICU-CRD",.N]),kept=c(sum(v1[db=="MIMIC-IV"]$id%in%P0[db=="MIMIC-IV"]$id),sum(v1[db=="eICU-CRD"]$id%in%P0[db=="eICU-CRD"]$id)))}
saveRDS(list(R=R,S=S,SG=SG,BT=BT,miss=miss,recl=recl,fl=readRDS("fl2.rds"),m=M_IMP),"R2.rds"); saveRDS(P0,"P0.rds"); cat("done\n")
