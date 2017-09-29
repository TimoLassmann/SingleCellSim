
// plot fit
fit = read.csv("GIGP_fit.csv",header=F,row.names=1)
plot(fit[,1], log = "xy")
lines(fit[,2],col = c("red"))

// add sim data at same depth.

samples = read.table("GIGP_sampled_2214851.csv",header = F)
points(samples[,2],col = c("green"),pch= 20)




;; plot sim single cell data
library(ggplot2);
library(gplots);

nzmean <- function(x) {
if (all(x==0)) 0 else mean(x)
}
nzcount <- function(x) {
sum(x !=0);
}

for(i in 1:length(files)){

	mat = read.table(files[i],header = F)
	x = colSums(mat[,1:ncol(mat)])
	y= mat;
	for( j in 1:ncol(mat)){
		y[,j] = mat[,j] / x[j]*1000000;
	}

	mean = apply(y,1,nzmean)
	count <- apply(y,1, nzcount)
	me = as.data.frame(cbind(mean,count))
	c <- ggplot(me, aes(me$mean, me$count/96 * 100))
	c + geom_point() + scale_x_log10(lim = c(0.1,1000))+labs(x="Mean t.p.m",y="Fraction Detected (% of 96 )",title="Detection vs tpm")
	ggsave(names[i])
}
