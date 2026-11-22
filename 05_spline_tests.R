# patch: RCS overall/non-linear and gap-category LR tests, pooled with D2 (Li et al. 1991)
src<-readLines("R/04_analysis.R"); stop_at<-grep("^run<-function",src)-1
eval(parse(text=src[1:stop_at]))
for(i in seq_along(IMP)) IMP[[i]][,t28m:=pmax(pmin(fifelse(is.na(death_h),672,death_h),672),0.5)/24]
sets<-list("MIMIC-IV"=lapply(IMP,function(d) d[db=="MIMIC-IV"]),"eICU-CRD"=lapply(IMP,function(d) d[db=="eICU-CRD"]),"Pooled"=IMP)
D2p<-function(dv,k){m<-length(dv); r<-(1+1/m)*var(sqrt(dv)); Dst<-(mean(dv)/k-(m+1)/(m-1)*r)/(1+r); df2<-k^(-3/m)*(m-1)*(1+1/r)^2; pf(max(Dst,0),k,df2,lower.tail=FALSE)}
rcs2<-function(dl){pool<-length(unique(dl[[1]]$db))>1; Bm<-ns(dl[[1]]$gap,df=4); gg<-seq(-15,25,0.25); C<-predict(Bm,gg)-matrix(predict(Bm,0),length(gg),4,byrow=TRUE)
 R<-lapply(dl,function(d){d<-copy(d); X<-predict(Bm,d$gap); for(j in 1:4) set(d,j=paste0("b",j),value=X[,j])
  f<-glm(fo("hd28",M3v,pool,paste0("b",1:4)),binomial,d); f0<-glm(fo("hd28",M3v,pool,"1"),binomial,d); fl<-glm(fo("hd28",M3v,pool,"gap"),binomial,d); fg<-glm(fo("hd28",M3v,pool),binomial,d)
  i<-paste0("b",1:4); V<-VC(f,d)[i,i]; list(ov=deviance(f0)-deviance(f),nl=deviance(fl)-deviance(f),inc=deviance(f0)-deviance(fg),lp=as.vector(C%*%coef(f)[i]),vv=rowSums((C%*%V)*C))})
 cur<-t(sapply(seq_along(gg),function(k) rub(sapply(R,function(r) r$lp[k]),sapply(R,function(r) r$vv[k]))))
 list(p=c(overall=D2p(sapply(R,`[[`,"ov"),4),nonlin=D2p(sapply(R,`[[`,"nl"),3)),inc=D2p(sapply(R,`[[`,"inc"),2),curve=data.table(gap=gg,or=cur[,1],lo=cur[,2],hi=cur[,3]))}
X<-readRDS("R2.rds"); for(k in names(sets)){r<-rcs2(sets[[k]]); X$R[[k]]$rcs<-r[c("p","curve")]; X$R[[k]]$inc<-r$inc; cat(k,signif(r$p,3),"inc",signif(r$inc,3),"\n")}
saveRDS(X,"R2.rds"); cat("patched\n")
