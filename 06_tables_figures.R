suppressMessages({library(data.table);library(survival);library(ggplot2);library(survminer);library(patchwork);library(jsonlite)})
X<-readRDS("R2.rds"); P<-as.data.table(readRDS("P0.rds")); R<-X$R; fl<-X$fl
P[,grp:=factor(grp,levels=c("Adequate","Resp. acidosis","Resp. alkalosis"))]
f1<-function(x,d=1) formatC(x,format="f",digits=d)
mi<-function(x,d=1){q<-quantile(x,c(.5,.25,.75),na.rm=TRUE); sprintf("%s (%s, %s)",f1(q[1],d),f1(q[2],d),f1(q[3],d))}
np<-function(x){x<-x[!is.na(x)]; sprintf("%d (%s)",sum(x),f1(100*mean(x)))}
fp<-function(p) ifelse(p<0.001,"<0.001",sprintf("%.3f",p))
mkT1<-function(P){G<-levels(P$grp); pk<-function(v) kruskal.test(P[[v]]~P$grp)$p.value; pc<-function(v) chisq.test(table(P[[v]],P$grp))$p.value
 row<-function(lab,f,v,p) data.table(Characteristic=lab,Overall=f(P[[v]]),A=f(P[grp==G[1]][[v]]),B=f(P[grp==G[2]][[v]]),C=f(P[grp==G[3]][[v]]),P=fp(p))
 cont<-list(age=c("Age, years",1),ph=c("pH",2),paco2=c("PaCO2, mmHg",1),hco3=c("Bicarbonate, mmol/L",1),gap=c("Winter gap, mmHg",1),ag=c("Anion gap, mmol/L",1),lact=c("Lactate, mmol/L",1),creat=c("Creatinine, mg/dL",1),bun=c("Blood urea nitrogen, mg/dL",0),k=c("Potassium, mmol/L",1),cl=c("Chloride, mmol/L",0),glu=c("Glucose, mg/dL",0),hb=c("Haemoglobin, g/dL",1),wbc=c("White cell count, 10^9/L",1),plt=c("Platelets, 10^9/L",0),alb=c("Albumin, g/dL",1),hr=c("Heart rate, beats/min",0),map=c("Mean arterial pressure, mmHg",0),gcs=c("Glasgow Coma Scale",0),pf=c("PaO2/FiO2, mmHg",0),sofa_nr=c("SOFA, non-respiratory non-neurological",0))
 bin<-list(male="Male sex",chf="Heart failure",ckd="Chronic kidney disease",copd="COPD",dm="Diabetes",liver="Chronic liver disease",cancer="Malignancy",sepsis="Sepsis",pna="Pneumonia",vaso="Vasopressor infusion at index gas",aki_pre="AKI present at or before index gas")
 rows<-list(); if(length(unique(P$db))>1) rows[[1]]<-row("Source database, eICU-CRD",np,"eicu",pc("eicu"))
 for(v in c("age","male","chf","ckd","copd","dm","liver","cancer","sepsis","pna","ph","paco2","hco3","gap","ag","lact","creat","bun","k","cl","glu","hb","wbc","plt","alb","hr","map","gcs","pf","sofa_nr","vaso","aki_pre")) rows[[length(rows)+1]]<-if(!is.null(cont[[v]])) row(cont[[v]][1],function(x) mi(x,as.numeric(cont[[v]][2])),v,pk(v)) else row(bin[[v]],np,v,pc(v))
 for(s in 1:3){P[,tmp:=as.integer(stage==s)]; rows[[length(rows)+1]]<-row(paste0("KDIGO AKI stage ",s),np,"tmp",chisq.test(table(P$stage,P$grp))$p.value)}
 for(v in list(c("hd28","In-hospital death within 28 days"),c("make30","MAKE30"),c("int48","Intubation within 48 h"),c("rrt7","KRT within 7 days"))) rows[[length(rows)+1]]<-row(v[2],np,v[1],pc(v[1]))
 T<-rbindlist(rows); nn<-table(P$grp); setnames(T,c("Characteristic",sprintf("Overall (n = %s)",format(nrow(P),big.mark=",")),sprintf("Adequate compensation (n = %s)",format(nn[[1]],big.mark=",")),sprintf("Concomitant respiratory acidosis (n = %s)",format(nn[[2]],big.mark=",")),sprintf("Concomitant respiratory alkalosis (n = %s)",format(nn[[3]],big.mark=",")),"P value")); T}
P[,eicu:=as.integer(db=="eICU-CRD")]; P[,rrt7:=as.integer(!is.na(rrt_h)&rrt_h<=168)]
T1<-mkT1(P); saveRDS(T1,"T1.rds"); S2<-list(M=mkT1(P[db=="MIMIC-IV"]),E=mkT1(P[db=="eICU-CRD"])); saveRDS(S2,"S2.rds")
cnt<-P[,.(n=.N,hd=sum(hd28),mk=sum(make30),int=sum(int48),rrt=sum(rrt7)),keyby=.(db,grp)]; cntp<-P[,.(n=.N,hd=sum(hd28),mk=sum(make30),int=sum(int48),rrt=sum(rrt7)),keyby=grp]
EX<-list(hosp_e=uniqueN(P[db=="eICU-CRD"]$hosp),tidx=mi(P$hidx),abg=P[,mean(src=="ABG"),by=db],cnt=cnt,cntp=cntp,disch28_e=P[db=="eICU-CRD"&hd28==0&disch_h<672,.N],dead_after_disch_m=P[db=="MIMIC-IV"&hd28==0&d28==1,.N],agep=P[,.(pre=sum(aki_pre)),by=db])
saveRDS(EX,"EX.rds"); print(EX)
# ---------- Figures ----------
th<-theme_classic(base_size=9,base_family="Arial")+theme(plot.background=element_rect(fill="white",colour=NA))
pal<-c("#000000","#D55E00","#0072B2"); labs3<-c("Adequate compensation (-2 to +2 mmHg)","Concomitant respiratory acidosis (> +2 mmHg)","Concomitant respiratory alkalosis (< -2 mmHg)")
nm<-R[["MIMIC-IV"]]$n; ne<-R[["eICU-CRD"]]$n
mk<-function(m,e) sprintf("MIMIC-IV n = %s | eICU-CRD n = %s",format(m,big.mark=","),format(e,big.mark=","))
bx<-data.table(lab=c(paste0("First ICU stay, age \u226518 years\n",mk(fl$m0,fl$e0)),paste0("Excluded: end-stage kidney disease\n",mk(fl$m0-fl$m_noesrd,fl$e0-fl$e_noesrd)),paste0("KDIGO AKI (creatinine criteria) within 48 h\n",mk(fl$m_aki,fl$e_aki)),paste0("Arterial blood gas with PaCO2, pH and HCO3\n(6 h before to 24 h after ICU admission)\n",mk(fl$m_abg,fl$e_abg)),paste0("Metabolic acidosis (HCO3 <22 mmol/L, pH <7.40)\n",mk(fl$m_ma,fl$e_ma)),
 paste0("Excluded: kidney replacement therapy before index gas\n",mk(fl$m_ma-fl$m_norrt,fl$e_ma-fl$e_norrt)),paste0("Excluded: ventilatory support at index gas\nInvasive: ",mk(fl$m_imv,fl$e_imv),"\nNoninvasive: ",mk(fl$m_niv,fl$e_niv)),paste0("Analytic cohort, unassisted breathing at index gas\n",mk(nm,ne),"\nPooled n = ",format(nm+ne,big.mark=","))),
 x=c(0,1.3,0,0,0,1.3,1.3,0),y=c(10,8.85,7.7,6.4,5.1,3.95,2.65,1.3),w=c(1.2,1.2,1.2,1.2,1.2,1.25,1.25,1.2),h=c(.42,.42,.42,.55,.42,.42,.6,.55))
g<-0.13  # gap between arrowhead and box border
v<-data.table(y0=c(10-.42,7.7-.42,6.4-.55,5.1-.42),y1=c(7.7+.42,6.4+.55,5.1+.42,1.3+.55)+g)
hz<-data.table(y=c(8.85,3.95,2.65),xe=c(1.3-0.6,1.3-0.625,1.3-0.625)-0.04)
pf<-ggplot()+geom_rect(data=bx,aes(xmin=x-w/2,xmax=x+w/2,ymin=y-h,ymax=y+h),fill="white",colour="black",linewidth=0.4)+geom_text(data=bx,aes(x=x,y=y,label=lab),family="Arial",size=2.5,lineheight=0.95)+
 geom_segment(data=v,aes(x=0,xend=0,y=y0,yend=y1),arrow=arrow(length=unit(1.5,"mm"),type="closed"),linewidth=0.4)+
 geom_segment(data=hz,aes(x=0,xend=xe,y=y,yend=y),arrow=arrow(length=unit(1.5,"mm"),type="closed"),linewidth=0.4)+
 coord_cartesian(xlim=c(-0.65,1.97),ylim=c(0.6,10.5),expand=FALSE)+theme_void()+theme(plot.background=element_rect(fill="white",colour=NA))
ggsave("output/Figure1.tiff",pf,width=170,height=155,units="mm",dpi=300,compression="lzw",bg="white"); ggsave("output/Figure1.png",pf,width=170,height=155,units="mm",dpi=150,bg="white")
cifp<-function(d,ttl,yl){d<-copy(d); dh<-fifelse(is.na(d$hdeath_h),Inf,d$hdeath_h); ft<-pmax(pmin(dh,d$disch_h,672),0.5); st<-fifelse(dh<=ft&dh<=672,1L,fifelse(d$disch_h<=ft&d$disch_h<672,2L,0L)); d[,`:=`(ftd=ft/24,st=st)]
 ci<-cmprsk::cuminc(d$ftd,d$st,d$grp,cencode=0); p<-ci$Tests[1,"pv"]; pv<-if(p<0.001) "Gray's test P < 0.001" else sprintf("Gray's test P = %.3f",p)
 cd<-rbindlist(lapply(seq_along(levels(d$grp)),function(i){e<-ci[[paste(levels(d$grp)[i],1)]]; data.table(t=e$time,y=e$est,g=labs3[i])})); cd[,g:=factor(g,levels=labs3)]
 a<-ggplot(cd,aes(t,y,colour=g))+geom_step(linewidth=0.8)+scale_colour_manual(values=pal)+scale_x_continuous(limits=c(0,28),breaks=seq(0,28,7))+scale_y_continuous(limits=c(0,0.5))+annotate("text",x=0.5,y=0.47,label=pv,hjust=0,family="Arial",size=2.6)+labs(title=ttl,x="Days",y=yl)+th+theme(legend.position="none",plot.title=element_text(face="bold",size=9))
 rk<-rbindlist(lapply(seq_along(levels(d$grp)),function(i) data.table(t=seq(0,28,7),g=labs3[i],n=sapply(seq(0,28,7),function(tt) sum(d$grp==levels(d$grp)[i]&d$ftd>=tt)))))
 rk[,g:=factor(g,levels=rev(labs3))]
 b<-ggplot(rk,aes(t,g,label=n,colour=g))+geom_text(family="Arial",size=2.2)+scale_colour_manual(values=rev(pal))+scale_x_continuous(limits=c(-1,29),breaks=seq(0,28,7))+labs(title="Number at risk (in hospital)",x=NULL,y=NULL)+theme_void(base_family="Arial")+theme(legend.position="none",plot.title=element_text(size=7),plot.background=element_rect(fill="white",colour=NA))
 a/b+plot_layout(heights=c(3,1))}
k1<-cifp(P[db=="MIMIC-IV"],"A  MIMIC-IV","Cumulative incidence of in-hospital death"); k2<-cifp(P[db=="eICU-CRD"],"B  eICU-CRD",""); k3<-cifp(P,"C  Pooled cohort","")
dd<-data.table(x=1:3,y=1,g=factor(labs3,levels=labs3)); lg<-cowplot::get_legend(ggplot(dd,aes(x,y,colour=g))+geom_line(linewidth=0.8)+scale_colour_manual(values=pal,name="Winter gap:",labels=c("Adequate (-2 to +2 mmHg)","Respiratory acidosis (> +2 mmHg)","Respiratory alkalosis (< -2 mmHg)"))+th+theme(legend.position="bottom",legend.text=element_text(size=7)))
pk2<-(k1|k2|k3)/wrap_elements(lg)+plot_layout(heights=c(12,1))
ggsave("output/Figure2.tiff",pk2,width=170,height=105,units="mm",dpi=300,compression="lzw",bg="white"); ggsave("output/Figure2.png",pk2,width=170,height=105,units="mm",dpi=150,bg="white")
cv<-R$Pooled$rcs$curve[gap>=-8]; rp<-R$Pooled$rcs$p
pr3<-ggplot(cv,aes(gap,or))+geom_ribbon(aes(ymin=lo,ymax=hi),fill="grey82")+geom_line(linewidth=0.7)+geom_hline(yintercept=1,linetype=2,colour="grey40")+geom_vline(xintercept=c(-2,2),linetype=3,colour="grey40")+scale_y_log10(breaks=c(0.5,1,1.5,2,3,4))+
 annotate("text",x=max(cv$gap),y=0.52,hjust=1,vjust=0,family="Arial",size=3,label=sprintf("P for overall %s\nP for non-linearity %s",ifelse(rp[1]<0.001,"< 0.001",sprintf("= %.3f",rp[1])),ifelse(rp[2]<0.001,"< 0.001",sprintf("= %.3f",rp[2]))))+labs(x="Winter gap, measured minus expected PaCO2 (mmHg)",y="Adjusted OR, in-hospital death within 28 days")+th
ggsave("output/Figure3.tiff",pr3,width=120,height=90,units="mm",dpi=300,compression="lzw",bg="white"); ggsave("output/Figure3.png",pr3,width=120,height=90,units="mm",dpi=200,bg="white")
SG<-X$SG; lv<-list("Age >=65 years"=c("Age, years","\u226565","<65"),"Male sex"=c("Sex","Male","Female"),"Sepsis"=c("Sepsis","Yes","No"),"KDIGO stage 3"=c("KDIGO AKI stage","3","1-2"),"Chronic kidney disease"=c("Chronic kidney disease","Yes","No"),"Lactate >2 mmol/L"=c("Lactate, mmol/L",">2","\u22642"),"Vasopressor at index"=c("Vasopressor use","Yes","No"),"COPD"=c("COPD","Yes","No"),"eICU-CRD (vs MIMIC-IV)"=c("Database","eICU-CRD","MIMIC-IV"))
fd<-rbindlist(lapply(names(lv),function(s){x<-SG[[s]]; rbind(data.table(lab=lv[[s]][1],or=NA,lo=NA,hi=NA,n=NA,pint=x$pint,hdr=TRUE),data.table(lab=paste0("   ",lv[[s]][2:3]),or=c(x$yes[1,1],x$no[1,1]),lo=c(x$yes[1,2],x$no[1,2]),hi=c(x$yes[1,3],x$no[1,3]),n=c(x$nyes,x$nno),pint=NA,hdr=FALSE))}))
fd[,yy:=.N:1]; fd[,txt:=fifelse(hdr,"",sprintf("%.2f (%.2f-%.2f)",or,lo,hi))]; fd[,pt:=fifelse(hdr,ifelse(pint<0.001,"<0.001",sprintf("%.3f",pint)),"")]; fd[,nt:=fifelse(hdr,"",format(n,big.mark=","))]
p4<-ggplot(fd,aes(y=yy))+geom_point(aes(x=or),shape=15,size=1.8,na.rm=TRUE)+geom_errorbar(aes(xmin=lo,xmax=hi),width=0.25,orientation="y",na.rm=TRUE,linewidth=0.4)+geom_segment(aes(x=1,xend=1,y=0.5,yend=max(fd$yy)+0.5),linetype=2,colour="grey40",linewidth=0.3)+
 scale_x_log10(limits=c(0.08,20),breaks=c(0.5,1,2,4))+geom_text(aes(x=0.085,label=lab,fontface=ifelse(hdr,"bold","plain")),hjust=0,family="Arial",size=2.7)+geom_text(aes(x=0.36,label=nt),hjust=1,family="Arial",size=2.7)+geom_text(aes(x=4.6,label=txt),hjust=0,family="Arial",size=2.7)+geom_text(aes(x=19.5,label=pt),hjust=1,family="Arial",size=2.7)+
 annotate("text",x=c(0.085,0.36,4.6,19.5),y=max(fd$yy)+1.5,label=c("Subgroup","n","OR (95% CI)","P for\ninteraction"),hjust=c(0,1,0,1),fontface="bold",family="Arial",size=2.7)+coord_cartesian(ylim=c(0.5,max(fd$yy)+2.3))+labs(x="Adjusted OR, respiratory acidosis vs adequate compensation",y=NULL)+th+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
ggsave("output/Figure4.tiff",p4,width=170,height=140,units="mm",dpi=300,compression="lzw",bg="white"); ggsave("output/Figure4.png",p4,width=170,height=140,units="mm",dpi=200,bg="white")
cat("done\n")
