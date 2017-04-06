# Calculates technical variation for a bisulfite sequencing experiment

# --------------------------------
# Computing standard error in methylation estimate
# for average CpG coverage of a certain size
# Assumes there is no bias towards coverage of some CpGs versus others

p <- seq(0,1,0.05); # %M in steps of 5%
maxR <- 100
# rows are expected number of reads for a CpG, column is 
# true methylation.
se <- matrix(NA,nrow=maxR, ncol=length(p))
for (n in 1:maxR) { se[n,] <- sqrt((p*(1-p))/n)}
colnames(se) <- paste("p",seq(0,100,5))
rownames(se) <- paste("nr",1:nrow(se),sep="")

# plot standard error at M%. One series per expected number of 
# reads per CpG.
library(RColorBrewer)
pos <- c(10,30,50,100) #seq(5,30,5)
plot(0,0,xlim=c(0,20),ylim=c(0,0.2), 
     bty='n', type='n', xaxt='n', 
     xlab="M/(M+U)",ylab="SD(sample~Bi(p_i))", 
     main="SD in %M calls for varying CpG coverage\n(WGBS model)"); 

pal <- brewer.pal(n=7,name="Dark2"); pal <- pal[-1]
ctr <- 1; 
for (k in pos) { 
	lines(se[k,],col=pal[ctr],type="o",pch=16); 
	ctr <- ctr+1
}
legend("topleft", bty='n',cex=0.8, border=NA, 
	legend=paste(pos,"rd/CpG"),fill=pal[1:(ctr-1)])
axis(1,at=1:20, labels=seq(5,100,5))
abline(h=c(0.02,0.05,0.1),col='red',lty=3)


#### plot variation in %M for a coverage of interest.
###cvg_level <- 4
###plot(0,0, type='n',ylim=c(0,110),xlim=c(0,110),xlab="Reported %M", ylab="Bounds of measurement error",
###     main=sprintf("Technical variation for %iX coverage", cvg_level))
###xpos <- seq(5,100,5); 
###rect(xleft=xpos-0.5, xright=xpos+0.5, 
###     ybottom=xpos-(se[cvg_level,]*100), 
###     ytop=xpos + (100*se[cvg_level,]), 
###     col="orange",border=NA)

