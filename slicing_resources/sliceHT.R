source("/mnt/wigstore4/user/krasnitz/rfunctions/peaks.valleys.AK1.ssc")
source("/mnt/wigstore4/user/krasnitz/rfunctions/set.tiers.new.R")
source("/mnt/wigstore4/user/krasnitz/TCGA/OV/Agi244A/hmsEventCalls/newMatchedMasking")
source("/home/becerra/code/slicing/resources/src/helperFunctions.R")
chromrange<-1:22
mastcols<-
	c("start","end","chrom","median","error","gstart","gend","amps","dels")

samples <- load_samples(classes = c("T", "F", "M"), sampleList = "resources/meta/sampleList.csv")
dir <- "resources/prev_run_8_2_2018_3/"

segtables <- lapply(samples, function(sample){
		tryCatch({
			table <- read.table(paste0(dir, sample, "_segtable.tsv"),header=T,as.is=T, stringsAsFactors = FALSE)
			return(table)
		}, error = function(e) NULL)
	})
raw_htprep <- do.call(rbind, segtables)

#
# Convert appopriate columns to numeric
#
num_data <- data.frame(data.matrix(raw_htprep))
numeric_columns <- sapply(num_data,function(x){mean(as.numeric(is.na(x)))<0.5})
htprep <- data.frame(num_data[,numeric_columns], raw_htprep[,!numeric_columns])
htprep$chrom <- raw_htprep$chrom

htprep[htprep[,"chrom"]=="X","chrom"]<-"23"
htprep[htprep[,"chrom"]=="Y","chrom"]<-"24"
htprep[,"chrom"]<-as.numeric(htprep[,"chrom"])
tumoramps<-htprep[,"marginalprob"]<0.001&htprep[,"mediandev"]>0
tumordels<-htprep[,"marginalprob"]<0.001&htprep[,"mediandev"]<0
htprep<-cbind(htprep,tumoramps,tumordels)
unames<-unique(htprep[,"ID"])
jointaslices<-matrix(ncol=9,nrow=0)
dimnames(jointaslices)[[2]]<-
c("start","end","tier","parent","chrom","is.parent","exp.index","gstart","gend")
jointdslices<-matrix(ncol=9,nrow=0)
dimnames(jointdslices)[[2]]<-
c("start","end","tier","parent","chrom","is.parent","exp.index","gstart","gend")

print(unames)
for(i in 1:length(unames)){
	print(unames[i])
	a<-htprep[htprep[,"ID"]==unames[i],
		c("start","end","chrom","mediandev","segerr","abs.pos.start","abs.pos.end",
		"tumoramps","tumordels")]
	a<-a[a[,"chrom"]%in%chromrange,,drop=F]
	dimnames(a)[[2]]<-mastcols
	z<-a[,"end"]-a[,"start"]+1
	a[,"end"]<-cumsum(z)
	a[,"start"]<-a[,"end"]-z+1
	chromstarts<-a[match(unique(a[,"chrom"]),a[,"chrom"]),"start"]
  chromends<-c((chromstarts[-1]-1),max(a[,"end"]))
	astpos<-cbind(a[,"median"],z,rep(1,nrow(a)),a)
	astpos[!astpos[,"amps"],1]<-(-5)
	astpos<-simplify.shorty(astpos,boundaries=chromends)
	astpos[,"end"]<-cumsum(astpos[,2])
	astpos[,"start"]<-c(1,(astpos[-nrow(astpos),"end"]+1))
  astpos[,"median"]<-astpos[,1]
  mastpos<-as.matrix(astpos[astpos[,"amps"],4:8,drop=F])
	astneg<-cbind(a[,"median"],z,rep(1,nrow(a)),a)
	astneg[!astneg[,"dels"],1]<-5
	astneg<-simplify.shorty(astneg,boundaries=chromends)
	astneg[,"end"]<-cumsum(astneg[,2])
	astneg[,"start"]<-c(1,(astneg[-nrow(astneg),"end"]+1))
  astneg[,"median"]<-astneg[,1]
  mastneg<-as.matrix(astneg[astneg[,"dels"],4:8,drop=F])
	aslices<-make.mounds(mastpos,i,inflate=2,mysign=1)
	if(nrow(aslices)>0){
		gstart<-a[match(aslices[,"start"],a[,"start"]),"gstart"]
		gend<-a[match(aslices[,"end"],a[,"end"]),"gend"]
		jointaslices<-rbind(jointaslices,cbind(aslices,gstart,gend))
	}
 	 dslices<-make.mounds(mastneg,i,inflate=2,mysign=-1)
	if(nrow(dslices)>0){
		gstart<-a[match(dslices[,"start"],a[,"start"]),"gstart"]
		gend<-a[match(dslices[,"end"],a[,"end"]),"gend"]
		jointdslices<-rbind(jointdslices,cbind(dslices,gstart,gend))
	}
}
tmp<-rbind(jointaslices,jointdslices)
tmp<-cbind(tmp,c(rep(1,nrow(jointaslices)),rep(0,nrow(jointdslices))))
tmp<-tmp[order(tmp[,"exp.index"]),]
tmp[,"parent"]<-rev(cummin(rev((1:nrow(tmp))+(1-tmp[,"is.parent"])*nrow(tmp))))
tmp<-cbind(tmp,match(tmp[,"exp.index"],unique(tmp[,"exp.index"])))
dimnames(tmp)[[2]][(ncol(tmp)-1):ncol(tmp)]<-c("lesion.type","profile.index")
jointslices<-tmp
rm(tmp)
tsl<-set.tiers(jointslices)
jointslices<-data.frame(I(unames[tsl[,"exp.index"]]),tsl)
dimnames(jointslices)[[2]][1]<-"profID"
write.table(jointslices,
	"/home/becerra/code/slicing/output/organoidSlices.txt",col.names=T,
	row.names=F,quote=F,sep="\t")
