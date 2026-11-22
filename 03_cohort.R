suppressMessages({library(DBI);library(data.table)})
con<-dbConnect(duckdb::duckdb(shared_home=FALSE),"w.duckdb",read_only=TRUE); Q<-function(s) as.data.table(dbGetQuery(con,s))
H<-function(a,b) as.numeric(difftime(a,b,units="hours")); fl<-list()
near<-function(L,key,ids,ref,lo,hi,a,z,fun="near"){x<-L[v>=lo&v<=hi]; x[,rh:=h-ref[match(get(key),ids)]]; x<-x[rh>=a&rh<=z]
 if(fun=="near"){x[,d:=abs(rh)]; setorderv(x,c(key,"d")); x<-x[,.SD[1],by=key]} else x<-x[,.(v=if(fun=="max") max(v) else min(v)),by=key]; x$v[match(ids,x[[key]])]}
akipair<-function(cr,key,upper){ # cr: key,h,v ; returns keys with KDIGO creatinine AKI at a value with h in [0,upper] (upper may be vector per key via column 'up')
 a<-cr[,.(k=get(key),hj=h,vj=v)]; b<-cr[,.(k=get(key),hi=h,vi=v)]; a<-a[hj>=-48]
 x<-merge(a,b,by="k",allow.cartesian=TRUE)[hi<hj & hi>=-168]
 x<-x[(vj/vi>=1.5 & hj-hi<=168)|(vj-vi>=0.3 & hj-hi<=48)]; x}
# ================= MIMIC =================
m<-Q("select * from st"); fl$m0<-nrow(m)
dx<-Q("select hadm_id,icd_code from dx"); has<-function(p) as.integer(m$hadm_id %in% dx[grepl(p,icd_code)]$hadm_id)
m[,esrd:=has("^(5856|N186|Z992|V4511|V560|V561|V568|Z4901|Z4931|Z4932)")]
L<-Q("select stay_id,specimen_id,itemid,charttime,v,value from lab"); L<-merge(L,m[,.(stay_id,intime)],by="stay_id"); L[,h:=H(charttime,intime)]
cr<-L[itemid==50912&v>0.1&v<25,.(stay_id,h,v)]
pr<-akipair(cr,"stay_id"); ak48<-pr[hj>=0&hj<=48,unique(k)]
b<-cr[h>=-168&h<=6,.(base=min(v)),by=stay_id]; pk<-cr[h>=0&h<=48,.(peak=max(v)),by=stay_id]
m<-merge(merge(m,b,by="stay_id",all.x=TRUE),pk,by="stay_id",all.x=TRUE)
pe<-Q("select * from pe"); pe<-pe[itemid!=225436]; pe<-merge(pe,m[,.(stay_id,intime)],by="stay_id"); pe[,`:=`(hs=H(starttime,intime),he=H(endtime,intime))]
rr_<-pe[itemid!=225792]; r48<-rr_[hs<=48,unique(stay_id)]
m[,aki:=as.integer(stay_id%in%ak48)]
m[,stage:=fifelse(aki==0,0L,fifelse(peak/base>=3|peak>=4|stay_id%in%r48,3L,fifelse(peak/base>=2,2L,1L)))]
fl$m_noesrd<-m[esrd==0,.N]; fl$m_aki<-m[esrd==0&aki==1,.N]
sp<-L[itemid==52033,.(specimen_id,st=value)]
g<-dcast(L[itemid %in% c(50818,50820,50803,50821)&h>=-6&h<=24],stay_id+specimen_id+h~itemid,value.var="v",fun.aggregate=function(x) x[1])
setnames(g,c("50818","50820","50803","50821"),c("paco2","ph","hco3a","pao2g"),skip_absent=TRUE)
g<-g[!is.na(paco2)&!is.na(ph)]; g[,st:=sp$st[match(specimen_id,sp$specimen_id)]]; g<-g[st=="ART."]
se<-L[itemid==50882&v>2&v<60,.(stay_id,hh=h,hs=v)]
gs<-merge(g[,.(stay_id,specimen_id,h)],se,by="stay_id",allow.cartesian=TRUE)[abs(hh-h)<=2][,d:=abs(hh-h)][order(d)][,.SD[1],by=specimen_id]
g[,hco3s:=gs$hs[match(specimen_id,gs$specimen_id)]]; g[,ok:=!is.na(hco3a)&hco3a>2&hco3a<60]; g[,hco3:=fifelse(ok,hco3a,hco3s)]; g[,src:=fifelse(ok,"ABG","serum")]; g[,hdiff:=fifelse(ok,0,gs$hh[match(specimen_id,gs$specimen_id)]-h)]
g<-g[!is.na(hco3)&paco2>=10&paco2<=120&ph>=6.6&ph<=7.8]
fl$m_abg<-m[esrd==0&aki==1&stay_id%in%g$stay_id,.N]
g<-g[hco3<22&ph<7.40][order(h)][,.SD[1],by=stay_id]
m<-merge(m,g[,.(stay_id,hidx=h,paco2,ph,hco3,src,hdiff,pao2g)],by="stay_id"); m<-m[esrd==0&aki==1]; fl$m_ma<-nrow(m)
m<-m[!stay_id%in%merge(rr_,m[,.(stay_id,hidx)],by="stay_id")[hs<hidx]$stay_id]; fl$m_norrt<-nrow(m)
# ventilatory support at index
p2<-Q("select * from pe2"); p2<-merge(p2,m[,.(stay_id,intime,hidx)],by="stay_id"); p2[,`:=`(hs=H(starttime,intime),he=H(endtime,intime))]
c2<-Q("select * from ce2"); c2<-merge(c2,m[,.(stay_id,intime,hidx)],by="stay_id"); c2[,rh:=H(charttime,intime)-hidx]; c2w<-c2[rh>=-2&rh<=2]
imv_ids<-unique(c(p2[itemid==225792&hs<=hidx+2&(is.na(he)|he>=hidx-1)]$stay_id, c2w[itemid==223900&grepl("ETT|Trach",value,ignore.case=TRUE)]$stay_id, c2w[itemid==226732&grepl("Endotracheal|Trach",value,ignore.case=TRUE)]$stay_id, c2w[itemid %in% c(220339,224684,224695,224696,224686)]$stay_id, c2w[itemid==223849&!grepl("CPAP|Bipap|NIV|Non",value,ignore.case=TRUE)]$stay_id))
niv_ids<-unique(c(p2[itemid==225794&hs<=hidx+2&(is.na(he)|he>=hidx-1)]$stay_id, c2w[itemid==226732&grepl("Bipap|CPAP",value,ignore.case=TRUE)]$stay_id))
m[,imv:=as.integer(stay_id%in%imv_ids)]; m[,niv:=as.integer(stay_id%in%niv_ids&imv==0)]
iv<-p2[itemid==225792&hs>hidx+2,.(t=min(hs)),by=stay_id]; m[,intub_h:=iv$t[match(stay_id,iv$stay_id)]-hidx]
rf<-merge(rr_,m[,.(stay_id,hidx)],by="stay_id")[hs>=hidx,.(t=min(hs)),by=stay_id]; m[,rrt_h:=rf$t[match(stay_id,rf$stay_id)]-hidx]
dt<-as.POSIXct(m$deathtime,tz="UTC"); m[,hdeath_h:=H(dt,intime)-hidx]; m[hdeath_h<0,hdeath_h:=0]
dd<-as.POSIXct(ifelse(is.na(m$dod),NA,paste(as.character(m$dod),"12:00:00")),tz="UTC"); dall<-dt; dall[is.na(dall)]<-dd[is.na(dall)]
m[,death_h:=H(dall,intime)-hidx]; m[death_h<0,death_h:=0]; m[,disch_h:=H(dischtime,intime)-hidx]
# AKI established at index
pre<-pr[merge(pr[,.(k)],m[,.(k=stay_id,hidx)],by="k")[,.I[1],by=k]$V1,][0]
prm<-merge(pr,m[,.(k=stay_id,hidx)],by="k")[hj<=hidx]; m[,aki_pre:=as.integer(stay_id%in%prm$k)]
bp<-merge(cr,m[,.(stay_id,hidx)],by="stay_id")[h>=-168&h<=hidx,.(bmin=min(v),bmax=max(v)),by=stay_id]
m[,`:=`(bpre=bp$bmin[match(stay_id,bp$stay_id)],cpre=bp$bmax[match(stay_id,bp$stay_id)])]
m[,stage_idx:=fifelse(aki_pre==0,0L,fifelse(cpre/bpre>=3|cpre>=4,3L,fifelse(cpre/bpre>=2,2L,1L)))]
ids<-m$stay_id; ref<-m$hidx; L2<-L[stay_id%in%ids,.(stay_id,itemid,h,v)]
gv<-function(it,lo,hi,a=-24,z=2,fun="near") near(L2[itemid%in%it],"stay_id",ids,ref,lo,hi,a,z,fun)
m[,`:=`(lact=gv(50813,.1,30),alb=gv(50862,1,6),bun=gv(51006,1,250),creat=gv(50912,.1,25),k=gv(50971,1.5,9),cl=gv(50902,60,140),glu=gv(50931,20,1500),hb=gv(51222,3,20),wbc=gv(51301,.1,200),plt=gv(51265,1,1500),ag=gv(50868,1,60),bili=gv(50885,.1,60,-24,2,"max"),pltmin=gv(51265,1,1500,-24,2,"min"),crmax=gv(50912,.1,25,-24,2,"max"),fio2l=gv(50816,21,100,-2,2))]
crl<-merge(cr,m[,.(stay_id,hidx,disch_h)],by="stay_id")[h>hidx&h<=hidx+720&(h-hidx)<=disch_h+24][order(h)][,.SD[.N],by=stay_id]; m[,`:=`(cr_last=crl$v[match(stay_id,crl$stay_id)],cr_last_h=crl$h[match(stay_id,crl$stay_id)]-m$hidx)]
ce<-Q("select * from ce"); ce<-merge(ce[stay_id%in%ids],m[,.(stay_id,intime)],by="stay_id"); ce[,h:=H(charttime,intime)]; ce[itemid==223761,`:=`(v=(v-32)/1.8,itemid=223762)]
cv<-function(it,lo,hi,a=-6,z=2,fun="near") near(ce[itemid%in%it],"stay_id",ids,ref,lo,hi,a,z,fun)
m[,`:=`(hr=cv(220045,20,250,-6,6),map=cv(c(220052,220181,225312),20,200,-6,6),mapmin=cv(c(220052,220181,225312),20,200,-24,2,"min"),fio2c=cv(223835,21,100,-2,2))]
gc<-ce[itemid%in%c(220739,223900,223901),.(s=sum(v),n=.N),by=.(stay_id,h)][n==3]; gc[,v:=s]; m[,gcs:=near(gc,"stay_id",ids,ref,3,15,-6,6)]
inp<-Q("select * from inp"); inp<-merge(inp[stay_id%in%ids],m[,.(stay_id,intime,hidx)],by="stay_id"); inp[,`:=`(hs=H(starttime,intime),he=H(endtime,intime))]
m[,vaso:=as.integer(stay_id%in%inp[hs<=hidx+1&he>=hidx-1]$stay_id)]
m[,fio2:=fifelse(!is.na(fio2l),fio2l,fifelse(!is.na(fio2c),fio2c,NA_real_))]; m[,pf:=pao2g/(fio2/100)]
sc<-function(x,br) fifelse(is.na(x),0L,as.integer(findInterval(x,br)))
m[,sofa_nr:=fifelse(is.na(pltmin),0L,4L-sc(pltmin,c(20,50,100,150)))+sc(bili,c(1.2,2,6,12))+fifelse(vaso==1,3L,fifelse(!is.na(mapmin)&mapmin<70,1L,0L))+sc(crmax,c(1.2,2,3.5,5))]
cx<-function(p) as.integer(m$hadm_id%in%dx[grepl(p,icd_code)]$hadm_id)
m[,`:=`(chf=cx("^(428|I50)"),ckd=cx("^(585[1-5]|5859|N18[1-5]|N189)"),copd=cx("^(49[0-6]|J4[0-4]|J47)"),dm=cx("^(250|E1[0-4])"),liver=cx("^(571|K70|K7[2-6])"),cancer=cx("^(1[4-9][0-9]|20[0-8]|C[0-9])"),sepsis=cx("^(99591|99592|78552|A40|A41|R652)"),pna=cx("^(48[0-6]|J1[2-8])"),dka=cx("^(2501|E1[0-4]1)"))]
m[,`:=`(male=as.integer(gender=="M"),db="MIMIC-IV",hosp="BIDMC",id=stay_id)]
# ================= eICU =================
e<-Q("select cast(patientunitstayid as bigint) s,uniquepid,gender,age,hospitalid,try_cast(hospitaldischargeoffset as int) hdo,hospitaldischargestatus hds,try_cast(unitvisitnumber as int) uvn,hospitaldischargeyear yr from ept")
e<-e[uvn==1][order(yr,s)][,.SD[1],by=uniquepid]; e[,age:=as.numeric(ifelse(age=="> 89","90",age))]; e<-e[!is.na(age)&age>=18]; fl$e0<-nrow(e)
EL<-Q("select * from elab"); EL<-EL[s%in%e$s&!is.na(v)]; EL[,h:=o/60]
edx<-Q("select s,d from edx"); eph<-Q("select s,d from eph"); fe<-function(p) as.integer(e$s%in%c(edx[grepl(p,d)]$s,eph[grepl(p,d)]$s))
aps<-Q("select cast(patientunitstayid as bigint) s,try_cast(dialysis as int) dial from eaps")
e[,esrd:=as.integer(fe("esrd|end stage renal|end-stage renal|hemodialysis|peritoneal dialysis|renal failure - dialysis|chronic kidney disease.*stage 5")==1|s%in%aps[dial==1]$s)]
ecr<-EL[n=="creatinine"&v>0.1&v<25,.(s,h,v)]; epr<-akipair(ecr,"s"); eak<-epr[hj>=0&hj<=48,unique(k)]
b<-ecr[h>=-168&h<=6,.(base=min(v)),by=s]; pk<-ecr[h>=0&h<=48,.(peak=max(v)),by=s]; e<-merge(merge(e,b,by="s",all.x=TRUE),pk,by="s",all.x=TRUE)
tr<-Q("select * from etr"); tr<-tr[s%in%e$s]; tr[,h:=o/60]; err<-tr[grepl("dialysis|crrt|cvvh|ultrafiltration",d)]
rc<-Q("select * from erc"); rc<-rc[s%in%e$s]; vst<-rbind(rc[,.(s,h=vs/60,he=ve/60)],tr[grepl("mechanical ventilation|intubation|endotracheal",d),.(s,h,he=NA_real_)])
e[,aki:=as.integer(s%in%eak)]; e[,stage:=fifelse(aki==0,0L,fifelse(peak/base>=3|peak>=4|s%in%err[h<=48]$s,3L,fifelse(peak/base>=2,2L,1L)))]
fl$e_noesrd<-e[esrd==0,.N]; fl$e_aki<-e[esrd==0&aki==1,.N]
pc<-EL[n=="paCO2"&h>=-6&h<=24&v>=10&v<=120]
prx<-function(nm,lo,hi,w){x<-merge(pc[,.(s,o)],EL[n==nm&v>=lo&v<=hi,.(s,o2=o,v)],by="s",allow.cartesian=TRUE)[abs(o2-o)<=w][,d:=abs(o2-o)][order(d)][,.SD[1],by=.(s,o)]; list(v=x$v[match(paste(pc$s,pc$o),paste(x$s,x$o))],dt=(x$o2-x$o)[match(paste(pc$s,pc$o),paste(x$s,x$o))]/60)}
a1<-prx("pH",6.6,7.8,15); a2<-prx("HCO3",2,60,15); a3<-prx("bicarbonate",2,60,120); a4<-prx("paO2",20,600,15); a5<-prx("FiO2",21,100,60)
pc[,`:=`(ph=a1$v,hco3a=a2$v,hco3s=a3$v,pao2g=a4$v,fio2l=a5$v)]; pc[,hco3:=fifelse(!is.na(hco3a),hco3a,hco3s)]; pc[,src:=fifelse(!is.na(hco3a),"ABG","serum")]; pc[,hdiff:=fifelse(!is.na(hco3a),a2$dt,a3$dt)]
pc<-pc[!is.na(ph)&!is.na(hco3)]; fl$e_abg<-e[esrd==0&aki==1&s%in%pc$s,.N]
pc<-pc[hco3<22&ph<7.40][order(h)][,.SD[1],by=s]
e<-merge(e,pc[,.(s,hidx=h,paco2=v,ph,hco3,src,hdiff,pao2g,fio2l)],by="s"); e<-e[esrd==0&aki==1]; fl$e_ma<-nrow(e)
e<-e[!s%in%merge(err,e[,.(s,hidx)],by="s")[h<hidx]$s]; fl$e_norrt<-nrow(e)
vi<-merge(vst,e[,.(s,hidx)],by="s")
r2<-Q("select * from erc2"); r2<-merge(r2,e[,.(s,hidx)],by="s"); r2[,rh:=o/60-hidx]; r2w<-r2[rh>=-2&rh<=2]
n2<-Q("select * from enc2"); n2<-merge(n2,e[,.(s,hidx)],by="s"); n2[,rh:=o/60-hidx]; n2w<-n2[rh>=-2&rh<=2]
ventlab<-c("PEEP","Vent Rate","Tidal Volume (set)","Peak Insp. Pressure","Plateau Pressure","Pressure Support","Mean Airway Pressure","Pressure Control","Adult Con Setting Set Vt","Set Vt (Servo,LTV)","Inspiratory Pressure, Set","Endotracheal Tube Placement","Secured at-ETT","Tube Size","Exhaled TV (machine)","RT Vent On/Off","Set Vt (Drager)","Adult Con Setting Set RR")
nivpat<-"bipap|bi-pap|cpap|niv|non-invasive|noninvasive|nippv|mask ventilation"; imvpat<-"ett|endotrach|trach|ventilator|vent "
e_imv<-unique(c(vi[h<=hidx+2&(is.na(he)|he>=hidx-1)]$s, r2w[l%in%ventlab]$s, r2w[l%in%c("O2 Device","Oxygen Delivery Method")&grepl(imvpat,tolower(v))]$s, n2w[n=="O2 Admin Device"&grepl(imvpat,tolower(v))]$s, n2w[n=="Verbal"&grepl("intubat|ett|trach",tolower(v))]$s))
e_niv<-unique(c(r2w[l%in%c("O2 Device","Oxygen Delivery Method","CPAP","PEEP/CPAP")&(l%in%c("CPAP","PEEP/CPAP")|grepl(nivpat,tolower(v)))]$s, n2w[n=="O2 Admin Device"&grepl(nivpat,tolower(v))]$s))
e[,imv:=as.integer(s%in%e_imv)]; e[,niv:=as.integer(s%in%e_niv&imv==0)]
iv<-vi[h>hidx+2,.(t=min(h)),by=s]; e[,intub_h:=iv$t[match(s,iv$s)]-hidx]
rf<-merge(err,e[,.(s,hidx)],by="s")[h>=hidx,.(t=min(h)),by=s]; e[,rrt_h:=rf$t[match(s,rf$s)]-hidx]
e[,disch_h:=hdo/60-hidx]; e[,hdeath_h:=fifelse(hds=="Expired",pmax(disch_h,0),NA_real_)]; e[,death_h:=NA_real_]
eprm<-merge(epr,e[,.(k=s,hidx)],by="k")[hj<=hidx]; e[,aki_pre:=as.integer(s%in%eprm$k)]
bp<-merge(ecr,e[,.(s,hidx)],by="s")[h>=-168&h<=hidx,.(bmin=min(v),bmax=max(v)),by=s]; e[,`:=`(bpre=bp$bmin[match(s,bp$s)],cpre=bp$bmax[match(s,bp$s)])]
e[,stage_idx:=fifelse(aki_pre==0,0L,fifelse(cpre/bpre>=3|cpre>=4,3L,fifelse(cpre/bpre>=2,2L,1L)))]
ids<-e$s; ref<-e$hidx; E2<-EL[s%in%ids,.(s,n,h,v)]
ge<-function(nm,lo,hi,a=-24,z=2,fun="near") near(E2[n%in%nm],"s",ids,ref,lo,hi,a,z,fun)
e[,`:=`(lact=ge("lactate",.1,30),alb=ge("albumin",1,6),bun=ge("BUN",1,250),creat=ge("creatinine",.1,25),k=ge("potassium",1.5,9),cl=ge("chloride",60,140),glu=ge("glucose",20,1500),hb=ge("Hgb",3,20),wbc=ge("WBC x 1000",.1,200),plt=ge("platelets x 1000",1,1500),ag=ge("anion gap",1,60),bili=ge("total bilirubin",.1,60,-24,2,"max"),pltmin=ge("platelets x 1000",1,1500,-24,2,"min"),crmax=ge("creatinine",.1,25,-24,2,"max"))]
crl<-merge(ecr,e[,.(s,hidx,disch_h)],by="s")[h>hidx&h<=hidx+720&(h-hidx)<=disch_h+24][order(h)][,.SD[.N],by=s]; e[,`:=`(cr_last=crl$v[match(s,crl$s)],cr_last_h=crl$h[match(s,crl$s)]-e$hidx)]
vp<-Q("select * from evp"); vp<-vp[s%in%ids]; vp[,h:=o/60]; va<-Q("select * from evap"); va<-va[s%in%ids]; va[,h:=o/60]
vhr<-vp[!is.na(hr),.(s,h,v=hr)]; vmap<-rbind(vp[!is.na(map),.(s,h,v=map)],va[!is.na(map),.(s,h,v=map)])
e[,`:=`(hr=near(vhr,"s",ids,ref,20,250,-6,6),map=near(vmap,"s",ids,ref,20,200,-6,6),mapmin=near(vmap,"s",ids,ref,20,200,-24,2,"min"))]
gn<-n2[n=="GCS Total",.(s,h=o/60,v=suppressWarnings(as.numeric(v)))][!is.na(v)]; e[,gcs:=near(gn,"s",ids,ref,3,15,-6,6)]
inf<-Q("select * from einf"); inf<-merge(inf[s%in%ids],e[,.(s,hidx)],by="s"); e[,vaso:=as.integer(s%in%inf[o/60>=hidx-6&o/60<=hidx+1]$s)]
e[,fio2:=fifelse(!is.na(fio2l),fio2l,NA_real_)]; e[,pf:=pao2g/(fio2/100)]
e[,sofa_nr:=fifelse(is.na(pltmin),0L,4L-sc(pltmin,c(20,50,100,150)))+sc(bili,c(1.2,2,6,12))+fifelse(vaso==1,3L,fifelse(!is.na(mapmin)&mapmin<70,1L,0L))+sc(crmax,c(1.2,2,3.5,5))]
e[,`:=`(chf=fe("congestive heart failure|heart failure"),ckd=fe("chronic kidney|chronic renal|renal insufficiency"),copd=fe("copd|chronic obstructive|emphysema"),dm=fe("diabetes"),liver=fe("cirrhosis|chronic liver|hepatic failure"),cancer=fe("cancer|malignan|carcinoma|lymphoma|leukemia|metasta"),sepsis=fe("sepsis|septic"),pna=fe("pneumonia"),dka=fe("ketoacidosis"))]
e[,`:=`(male=as.integer(gender=="Male"),db="eICU-CRD",hosp=as.character(hospitalid),id=s)]
cols<-c("id","db","hosp","age","male","base","peak","stage","stage_idx","aki_pre","hidx","paco2","ph","hco3","src","hdiff","imv","niv","intub_h","rrt_h","hdeath_h","death_h","disch_h","cr_last","cr_last_h","lact","alb","bun","creat","k","cl","glu","hb","wbc","plt","ag","hr","map","gcs","pf","vaso","sofa_nr","chf","ckd","copd","dm","liver","cancer","sepsis","pna","dka")
D<-rbind(m[,..cols],e[,..cols]); D[,gap:=paco2-(1.5*hco3+8)]
D[,grp:=factor(fifelse(gap< -2,"Resp. alkalosis",fifelse(gap>2,"Resp. acidosis","Adequate")),levels=c("Adequate","Resp. acidosis","Resp. alkalosis"))]
D[,supp:=fifelse(imv==1,"IMV",fifelse(niv==1,"NIV","None"))]
fl$m_imv<-m[imv==1,.N]; fl$m_niv<-m[niv==1,.N]; fl$e_imv<-e[imv==1,.N]; fl$e_niv<-e[niv==1,.N]
saveRDS(D,"D2.rds"); saveRDS(fl,"fl2.rds"); print(unlist(fl))
print(D[supp=="None",.(n=.N,gcs3=mean(gcs==3,na.rm=T),gap=median(gap),acid=mean(grp=="Resp. acidosis"),hd28=mean(!is.na(hdeath_h)&hdeath_h<=672),akipre=mean(aki_pre)),by=db])
print(D[supp=="None",.(n=.N,hd28=mean(!is.na(hdeath_h)&hdeath_h<=672)),by=.(db,grp)][order(db,grp)])
