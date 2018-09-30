# sliced tables are initially generated with tiers numbered upside down: the
# top tier is number 1. Re-number tiers from bottom up; also set parents
set.tiers<-
	function(slicetab,startcol="start",endcol="end",profcol="profile.index"){
	a<-max(slicetab[,endcol])+1
	b<-max(slicetab[,endcol]-slicetab[,startcol])+1
	ordhead<-2*slicetab[,profcol]*a*b+2*slicetab[,startcol]*b+
		b-slicetab[,endcol]+slicetab[,startcol]-1
	ordtail<-2*slicetab[,profcol]*a*b+2*slicetab[,endcol]*b+
		b+slicetab[,endcol]-slicetab[,startcol]
	z<-cbind(c(ordhead,ordtail),
			c(rep(1,nrow(slicetab)),rep(-1,nrow(slicetab))),
			rbind(slicetab,slicetab))
	dimnames(z)[[2]]<-c("ord","pol",dimnames(slicetab)[[2]])
	z<-z[order(z[,"ord"]),-match("ord",dimnames(z)[[2]]),drop=F]
	tier<-cumsum(z[,"pol"])[z[,"pol"]==1]
	z<-z[z[,"pol"]==1,-match("pol",dimnames(z)[[2]]),drop=F]
	parent<-rep(0,nrow(z))
	if(max(tier)>1)for(i in max(tier):(min(tier)+1))
		parent[tier==i]<-cummax((1:nrow(z))*(tier==(i-1)))[tier==i]
	if("tier"%in%dimnames(z)[[2]])z[,"tier"]<-tier
	else{
		dz<-dimnames(z)[[2]]
		z<-cbind(z,tier)
		dimnames(z)[[2]]<-c(dz,"tier")
	}
	if("parent"%in%dimnames(z)[[2]])z[,"parent"]<-parent
	else{
		dz<-dimnames(z)[[2]]
		z<-cbind(z,parent)
		dimnames(z)[[2]]<-c(dz,"parent")
	}
	z
}
		
