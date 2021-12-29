l1f_lw_20bins <- read.table("L1F_LW_averecombrates_perpeak_20Bins.txt",header=T, sep='\t')
l2f_lw_20bins <- read.table("L2F_LW_averecombrates_perpeak_20Bins.txt",header=T, sep='\t')
l1m_lw_20bins <- read.table("L1M_LW_averecombrates_perpeak_20Bins.txt",header=T, sep='\t')
l2m_lw_20bins <- read.table("L2M_LW_averecombrates_perpeak_20Bins.txt",header=T, sep='\t')

o1_lw_20bins <- read.table("O1_LW_averecombrates_perpeak_20Bins.txt",header=T, sep='\t')
o2_lw_20bins <- read.table("O2_LW_averecombrates_perpeak_20Bins.txt",header=T, sep='\t')
t1_lw_20bins <- read.table("T1_LW_averecombrates_perpeak_20Bins.txt",header=T, sep='\t')
t2_lw_20bins <- read.table("T2_LW_averecombrates_perpeak_20Bins.txt",header=T, sep='\t')


l1f_last_bins_rates <- subset(l1f_lw_20bins, Bin.Num == "0.05" | Bin.Num=="1.0",select=c("Bin.Num","Ave.Peak.Recomb.Rate"))
l2f_last_bins_rates <- subset(l2f_lw_20bins, Bin.Num == "0.05" | Bin.Num=="1.0",select=c("Bin.Num","Ave.Peak.Recomb.Rate"))
l1m_last_bins_rates <- subset(l1m_lw_20bins, Bin.Num == "0.05" | Bin.Num=="1.0",select=c("Bin.Num","Ave.Peak.Recomb.Rate"))
l2m_last_bins_rates <- subset(l2m_lw_20bins, Bin.Num == "0.05" | Bin.Num=="1.0",select=c("Bin.Num","Ave.Peak.Recomb.Rate"))

o1_last_bins_rates <- subset(o1_lw_20bins, Bin.Num == "0.05" | Bin.Num=="1.0",select=c("Bin.Num","Ave.Peak.Recomb.Rate"))
o2_last_bins_rates <- subset(o2_lw_20bins, Bin.Num == "0.05" | Bin.Num=="1.0",select=c("Bin.Num","Ave.Peak.Recomb.Rate"))
t1_last_bins_rates <- subset(t1_lw_20bins, Bin.Num == "0.05" | Bin.Num=="1.0",select=c("Bin.Num","Ave.Peak.Recomb.Rate"))
t2_last_bins_rates <- subset(t2_lw_20bins, Bin.Num == "0.05" | Bin.Num=="1.0",select=c("Bin.Num","Ave.Peak.Recomb.Rate"))

l1f_last_bins_rates$Tissue <- rep("L1F",times=length(l1f_last_bins_rates$Bin.Num))
l2f_last_bins_rates$Tissue <- rep("L2F",times=length(l2f_last_bins_rates$Bin.Num))
l1m_last_bins_rates$Tissue <- rep("L1M",times=length(l1m_last_bins_rates$Bin.Num))
l2m_last_bins_rates$Tissue <- rep("L2M",times=length(l2m_last_bins_rates$Bin.Num))
o1_last_bins_rates$Tissue <- rep("O1",times=length(o1_last_bins_rates$Bin.Num))
o2_last_bins_rates$Tissue <- rep("O2",times=length(o2_last_bins_rates$Bin.Num))
t1_last_bins_rates$Tissue <- rep("T1",times=length(t1_last_bins_rates$Bin.Num))
t2_last_bins_rates$Tissue <- rep("T2",times=length(t2_last_bins_rates$Bin.Num))


combined_rates_bins_df <- data.frame(rates = c(l1f_last_bins_rates$Ave.Peak.Recomb.Rate, l2f_last_bins_rates$Ave.Peak.Recomb.Rate, l1m_last_bins_rates$Ave.Peak.Recomb.Rate, l2m_last_bins_rates$Ave.Peak.Recomb.Rate, o1_last_bins_rates$Ave.Peak.Recomb.Rate, o2_last_bins_rates$Ave.Peak.Recomb.Rate, t1_last_bins_rates$Ave.Peak.Recomb.Rate, t2_last_bins_rates$Ave.Peak.Recomb.Rate), tissue= c(l1f_last_bins_rates$Tissue, l2f_last_bins_rates$Tissue, l1m_last_bins_rates$Tissue, l2m_last_bins_rates$Tissue, o1_last_bins_rates$Tissue, o2_last_bins_rates$Tissue, t1_last_bins_rates$Tissue, t2_last_bins_rates$Tissue))

end_bins_test <- ggplot(data=combined_rates_bins_df, aes(x=tissue, y=rates)) + geom_boxplot(notch=TRUE) + ylim(0,0.001)
