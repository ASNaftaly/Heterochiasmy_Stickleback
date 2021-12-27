library(ggplot2)
library(cowplot)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#reading in data and creating dataframes to create figures

l1f_observed <- read.table("L1F_filteredpeaks_10bin_percents.txt",header=T,sep="\t")
l2f_observed <-read.table("L2F_filteredpeaks_10bin_percents.txt",header=T,sep="\t")
o1_observed <- read.table("O1_filteredpeaks_10bin_percents.txt",header=T,sep="\t")
o2_observed <-read.table("O2_filteredpeaks_10bin_percents.txt",header=T,sep="\t")

l1m_observed <- read.table("L1M_filteredpeaks_10bin_percents.txt",header=T,sep="\t")
l2m_observed <-read.table("L2M_filteredpeaks_10bin_percents.txt",header=T,sep="\t")
t1_observed <- read.table("T1_filteredpeaks_10bin_percents.txt",header=T,sep="\t")
t2_observed <-read.table("T2_filteredpeaks_10bin_percents.txt",header=T,sep="\t")

l1f_expected <- read.table("Random_peaks_10bins_percents_L1F.txt",header=F, sep="\t")
colnames(l1f_expected)[1] <- "Bin.Num"
l2f_expected <- read.table("Random_peaks_10bins_percents_L2F.txt",header=F, sep="\t")
colnames(l2f_expected)[1] <- "Bin.Num"
o1_expected <- read.table("Random_peaks_10bins_percents_O1.txt",header=F, sep="\t")
colnames(o1_expected)[1] <- "Bin.Num"
o2_expected <- read.table("Random_peaks_10bins_percents_O2.txt",header=F, sep="\t")
colnames(o2_expected)[1] <- "Bin.Num"

l1m_expected <- read.table("Random_peaks_10bins_percents_L1M.txt",header=F, sep="\t")
colnames(l1m_expected)[1] <- "Bin.Num"
l2m_expected <- read.table("Random_peaks_10bins_percents_L2M.txt",header=F, sep="\t")
colnames(l2m_expected)[1] <- "Bin.Num"
t1_expected <- read.table("Random_peaks_10bins_percents_T1.txt",header=F, sep="\t")
colnames(t1_expected)[1] <- "Bin.Num"
t2_expected <- read.table("Random_peaks_10bins_percents_T2.txt",header=F, sep="\t")
colnames(t2_expected)[1] <- "Bin.Num"

#combining samples
average_fl_observed_percents <- (l1f_observed$Peak.Proportion + l2f_observed$Peak.Proportion)/2
averaged_female_liver_observed <- data.frame(Bin_Num=l1f_observed$Bin.Num, Percent=average_fl_observed_percents)

average_ml_observed_percents <- (l1m_observed$Peak.Proportion + l2m_observed$Peak.Proportion)/2
averaged_male_liver_observed <- data.frame(Bin_Num=l1m_observed$Bin.Num, Percent=average_ml_observed_percents)

average_o_observed_percents <- (o1_observed$Peak.Proportion + o2_observed$Peak.Proportion)/2
averaged_ovary_observed <- data.frame(Bin_Num=o1_observed$Bin.Num, Percent=average_o_observed_percents)

average_t_observed_percents <- (t1_observed$Peak.Proportion + t2_observed$Peak.Proportion)/2
averaged_testis_observed <- data.frame(Bin_Num=t1_observed$Bin.Num, Percent=average_t_observed_percents)

combined_female_liver_expected_percents <- merge(l1f_expected, l2f_expected, by="Bin.Num")
combined_male_liver_expected_percents <- merge(l1m_expected, l2m_expected, by="Bin.Num")
combined_ovary_expected_percents <- merge(o1_expected, o2_expected, by="Bin.Num")
combined_testis_expected_percents <- merge(t1_expected, t2_expected, by="Bin.Num")

#calculating mean, minimum, and maximum percentage

l1f_expected_stats <- data.frame(Bin=l1f_expected[,1],Means=rowMeans(l1f_expected[,-1]))
l1f_expected_stats$Min <- apply(subset(l1f_expected,select=V2:V10001),1,min)
l1f_expected_stats$Max <- apply(subset(l1f_expected,select=V2:V10001),1,max)

l2f_expected_stats <- data.frame(Bin=l2f_expected[,1],Means=rowMeans(l2f_expected[,-1]))
l2f_expected_stats$Min <- apply(subset(l2f_expected,select=V2:V10001),1,min)
l2f_expected_stats$Max <- apply(subset(l2f_expected,select=V2:V10001),1,max)

fl_expected_stats <- data.frame(Bin=combined_female_liver_expected_percents[,1],Means=rowMeans(combined_female_liver_expected_percents[,-1]))
fl_expected_stats$Min <- apply(subset(combined_female_liver_expected_percents,select=V2.x:V10001.y),1,min)
fl_expected_stats$Max <- apply(subset(combined_female_liver_expected_percents,select=V2.x:V10001.y),1,max)

l1m_expected_stats <- data.frame(Bin=l1m_expected[,1],Means=rowMeans(l1m_expected[,-1]))
l1m_expected_stats$Min <- apply(subset(l1m_expected,select=V2:V10001),1,min)
l1m_expected_stats$Max <- apply(subset(l1m_expected,select=V2:V10001),1,max)

l2m_expected_stats <- data.frame(Bin=l2m_expected[,1],Means=rowMeans(l2m_expected[,-1]))
l2m_expected_stats$Min <- apply(subset(l2m_expected,select=V2:V10001),1,min)
l2m_expected_stats$Max <- apply(subset(l2m_expected,select=V2:V10001),1,max)

ml_expected_stats <- data.frame(Bin=combined_male_liver_expected_percents[,1],Means=rowMeans(combined_male_liver_expected_percents[,-1]))
ml_expected_stats$Min <- apply(subset(combined_male_liver_expected_percents,select=V2.x:V10001.y),1,min)
ml_expected_stats$Max <- apply(subset(combined_male_liver_expected_percents,select=V2.x:V10001.y),1,max)

o1_expected_stats <- data.frame(Bin=o1_expected[,1],Means=rowMeans(o1_expected[,-1]))
o1_expected_stats$Min <- apply(subset(o1_expected,select=V2:V10001),1,min)
o1_expected_stats$Max <- apply(subset(o1_expected,select=V2:V10001),1,max)

o2_expected_stats <- data.frame(Bin=o2_expected[,1],Means=rowMeans(o2_expected[,-1]))
o2_expected_stats$Min <- apply(subset(o2_expected,select=V2:V10001),1,min)
o2_expected_stats$Max <- apply(subset(o2_expected,select=V2:V10001),1,max)

o_expected_stats <- data.frame(Bin=combined_ovary_expected_percents[,1],Means=rowMeans(combined_ovary_expected_percents[,-1]))
o_expected_stats$Min <- apply(subset(combined_ovary_expected_percents,select=V2.x:V10001.y),1,min)
o_expected_stats$Max <- apply(subset(combined_ovary_expected_percents,select=V2.x:V10001.y),1,max)

t1_expected_stats <- data.frame(Bin=t1_expected[,1],Means=rowMeans(t1_expected[,-1]))
t1_expected_stats$Min <- apply(subset(t1_expected,select=V2:V10001),1,min)
t1_expected_stats$Max <- apply(subset(t1_expected,select=V2:V10001),1,max)

t2_expected_stats <- data.frame(Bin=t2_expected[,1],Means=rowMeans(t2_expected[,-1]))
t2_expected_stats$Min <- apply(subset(t2_expected,select=V2:V10001),1,min)
t2_expected_stats$Max <- apply(subset(t2_expected,select=V2:V10001),1,max)

t_expected_stats <- data.frame(Bin=combined_testis_expected_percents[,1],Means=rowMeans(combined_testis_expected_percents[,-1]))
t_expected_stats$Min <- apply(subset(combined_testis_expected_percents,select=V2.x:V10001.y),1,min)
t_expected_stats$Max <- apply(subset(combined_testis_expected_percents,select=V2.x:V10001.y),1,max)


#Main text figure 1
jpeg(file="Binned_peak_percents_maintext_Individual1_V3.jpeg",width=1800,height=1000,quality=100)

fl_plot <- ggplot(data = l1f_observed, aes(x=Bin.Num, y=Peak.Proportion)) + geom_bar(stat="identity",fill="tan4",alpha=0.75) + ylab("Percent of Accessible\n Regions in Bin") + xlab("") + ylim(0,20) + theme(axis.text=element_text(size=32),axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),size=36),plot.margin=unit(c(1.5,1,1,1),"cm")) + geom_hline(yintercept=10, color="black",size=2,linetype="dashed")

ml_plot <- ggplot(data = l1m_observed, aes(x=Bin.Num, y=Peak.Proportion)) + geom_bar(stat="identity",fill="tan4",alpha=0.75) + ylab("Percent of Accessible\n Regions in Bin") + xlab("") + ylim(0,20) + theme(axis.text=element_text(size=32),axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),size=36),plot.margin=unit(c(1.5,1,1,1),"cm")) + geom_hline(yintercept=10, color="black",size=2,linetype="dashed")

o_plot <- ggplot(data = o1_observed, aes(x=Bin.Num, y=Peak.Proportion)) + geom_bar(stat="identity",fill="purple",alpha=0.75) + ylab("Percent of Accessible\n Regions in Bin") + xlab("") + ylim(0,20) + theme(axis.text=element_text(size=32),axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),size=36),plot.margin=unit(c(1.5,1,1,1),"cm")) + geom_hline(yintercept=10, color="black",size=2,linetype="dashed")

t_plot <- ggplot(data = t1_observed, aes(x=Bin.Num, y=Peak.Proportion)) + geom_bar(stat="identity",fill="blue",alpha=0.75) + ylab("Percent of Accessible\n Regions in Bin") + xlab("") + ylim(0,20) + theme(axis.text=element_text(size=32),axis.title.y=element_text(margin=margin(t=0,r=20,b=0,l=0),size=36),plot.margin=unit(c(1.5,1,1,1),"cm")) + geom_hline(yintercept=10, color="black",size=2,linetype="dashed")

plot_grid(fl_plot, ml_plot,o_plot, t_plot,labels=c("A","B","C", "D"), ncol=2, nrow=2,label_size=30)

dev.off()
