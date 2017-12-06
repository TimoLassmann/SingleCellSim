
library(ggplot2);
library(gplots);


stats = read.csv("GIGP_parameters.csv", header = T, row.names = 1) 
stats = t(stats) 
cond =   ifelse(grepl("Normal", row.names(stats)),"Normal","Tumour")
stats = cbind(stats,cond) 
stats = as.data.frame(stats) 
stats$S = as.numeric(as.character(stats$S))
ggplot(as.data.frame(stats), aes(x=cond, y= S)) + geom_violin() +  ylab("Estimated unique number of transcipts")+  labs(title = "BLCA")
ggsave("violin_S.pdf", width = 20, height = 20, units = "cm" )
ggplot(as.data.frame(stats), aes(x=cond, y= S)) + geom_boxplot() + geom_jitter(width = 0.2) +  ylab("Estimated unique number of transcipts")+  labs(title = "BLCA") 
ggsave("bar_jitter_S.pdf", width = 20, height = 20, units = "cm")





theme_grey(base_size = 18) 

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
	c <- ggplot(me, aes(me$mean, me$count/ncol(mat) * 100))
  ylabel = paste("Fraction Detected (% of",ncol(mat),")");
	c + geom_point() + scale_x_log10(lim = c(0.1,1000))+labs(x="Mean t.p.m",y=ylabel,title="Detection vs tpm")
	ggsave(names[i])
}

mat = read.csv("extable_CPhII001_CA8WCANXX_ACATTA_L006_R1.fastq.bam_96_80000_CV0.750000.csv",header = T,row.names = 1) 


cat = gsub("Blankgene[[:digit:]]+","Blank",rownames(mat) )



y = mat;
mean = apply(y,1,nzmean)
count <- apply(y,1, nzcount)
count = count / ncol(mat) * 100
me = as.data.frame(cbind(mean,count,cat))

c <- ggplot(me, aes(me$mean, me$count/ncol(mat) * 100))
  ylabel = paste("Fraction Detected (% of",ncol(mat),")");
c =	c + geom_point() + scale_x_log10(lim = c(0.1,1000))+labs(x="Mean t.p.m",y=ylabel,title="Detection vs tpm")
c +facet_wrap(~cat)
c += geom_point(data = gene, stat="bin") 



ggplot(me) + geom_point(aes(x = as.numeric(as.character(mean)), y = as.numeric(as.character(count))  , colour = cat)) + scale_x_log10(lim = c(0.1,1000))+labs(x="Mean t.p.m",y=ylabel,title="Detection vs tpm") + facet_wrap(~cat)


'

// good from here... 

library(viridis)
library(ggplot2);
library(gplots);
library(tibble);
library(dplyr);

nzmean <- function(x) {
if (all(x==0)) 0 else mean(x)
}
nzcount <- function(x) {
sum(x !=0);
}


mat = read.csv("extable_CPhII001_CA8WCANXX_ACATTA_L006_R1.fastq.bam_96_80000_CV0.750000.csv",header = T,row.names = 1) 
y = mat;
mean = apply(y,1,nzmean)
count <- apply(y,1, nzcount)
count = count / ncol(mat) * 100

cat = gsub("Blankgene[[:digit:]]+","Blank",rownames(mat) )


scheme = viridis_pal(option = "D")(length(unique(cat))-1) 

scheme = c("#999999",scheme)

ylabel = paste("Fraction Detected (% of",ncol(mat),")");



tib = tibble(mean = mean, count = count , cat = cat) 

gene = tib %>%  filter(cat != "Blank")

c = ggplot(tib, aes (mean,count))
c = c+ geom_point(aes(colour = cat))
c = c + labs(x="Mean counts",y=ylabel,title="Detection vs counts")
c = c + scale_color_manual(values= scheme)
c = c + scale_x_log10(lim = c(0.1,1000))
c + geom_text(aes(label = cat),check_overlap = TRUE,data = gene,col =c("black")) + theme(legend.position = c(0.8, 0.2))

ggsave("Simulation80Kreadsprecell.pdf", width = 40, height = 40, units = "cm")

//reformat ex table..

mat = read.csv("genecounts.csv",header =T , row.names = 1) 
tmp <-with(mat, aggregate(x = mat[,1:111] , by=list(Symbol), FUN=median))

out = as.data.frame(apply(tmp[,2:ncol(tmp)],1, median, na.rm = TRUE))
rownames(out) = tmp[,1] 

colnames(out) = c("Mean")



write.csv(out,"mean.csv",quote =FALSE ) 



//ing gagag 

IFNA1,IFNA2,IFNA4,IFNA5,IFNA6,IFNA7,IFNA8,IFNA10,IFNA14,IFNA16,IFNA17,IFNA21,IFNB1,IFNL1,IFNL2,IFNL3,IFNW1

CMPK2,DDX58,DHX58,GBP1,GBP3,GBP5,HELZ2,HERC5,IFI35,IFI44,IFI44L,IFI6,IFIH1,IFIT1,IFIT2,IFIT3,IFITM1,IFITM3,IL33,IRF7,IRF8,ISG15,ISG20,LILRA5,MX1,MX2,NT5C3A,OAS1,OAS2,OAS3,OASL,OSM,PARP14,PARP9,PHF11,PLSCR1,RSAD2,RTP4,SAMD9L,STAT1,STAT2,TAP1,TDRD7,TLR2,TLR4,TLR5,TNFAIP3,TNFSF13B,TOR1B,TRIM22,TRIM5,UBE2L6,USP18,XAF1,ZBP1,ZC3HAV1


