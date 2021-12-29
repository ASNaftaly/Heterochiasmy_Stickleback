library(ggplot2)

#LW rates
#read files

l1f_lw_rates <- read.table("L1F_LW_autosomes_averecombrates_per_peak.txt",header=TRUE, sep="\t")
l2f_lw_rates <- read.table("L2F_LW_autosomes_averecombrates_per_peak.txt",header=TRUE, sep="\t")
l1m_lw_rates <- read.table("L1M_LW_autosomes_averecombrates_per_peak.txt",header=TRUE, sep="\t")
l2m_lw_rates <- read.table("L2M_LW_autosomes_averecombrates_per_peak.txt",header=TRUE, sep="\t")

o1_lw_rates <- read.table("O1_LW_autosomes_averecombrates_per_peak.txt",header=TRUE, sep="\t")
o2_lw_rates <- read.table("O2_LW_autosomes_averecombrates_per_peak.txt",header=TRUE, sep="\t")
t1_lw_rates <- read.table("T1_LW_autosomes_averecombrates_per_peak.txt",header=TRUE, sep="\t")
t2_lw_rates <- read.table("T2_LW_autosomes_averecombrates_per_peak.txt",header=TRUE, sep="\t")

l1f_lw_rates$Tissue <- rep("L1F",times=length(l1f_lw_rates$Chr))
l2f_lw_rates$Tissue <- rep("L2F",times=length(l2f_lw_rates$Chr))
l1m_lw_rates$Tissue <- rep("L1M",times=length(l1m_lw_rates$Chr))
l2m_lw_rates$Tissue <- rep("L2M",times=length(l2m_lw_rates$Chr))
o1_lw_rates$Tissue <- rep("O1",times=length(o1_lw_rates$Chr))
o2_lw_rates$Tissue <- rep("O2",times=length(o2_lw_rates$Chr))
t1_lw_rates$Tissue <- rep("T1",times=length(t1_lw_rates$Chr))
t2_lw_rates$Tissue <- rep("T2",times=length(t2_lw_rates$Chr))

combined_rates_df <- data.frame(rates = c(l1f_lw_rates$Ave.Recomb.Rate, l2f_lw_rates$Ave.Recomb.Rate, l1m_lw_rates$Ave.Recomb.Rate, l2m_lw_rates$Ave.Recomb.Rate, o1_lw_rates$Ave.Recomb.Rate, o2_lw_rates$Ave.Recomb.Rate, t1_lw_rates$Ave.Recomb.Rate, t1_lw_rates$Ave.Recomb.Rate, t2_lw_rates$Ave.Recomb.Rate), tissue= c(l1f_lw_rates$Tissue, l2f_lw_rates$Tissue, l1m_lw_rates$Tissue, l2m_lw_rates$Tissue, o1_lw_rates$Tissue, o2_lw_rates$Tissue, t1_lw_rates$Tissue, t1_lw_rates$Tissue, t2_lw_rates$Tissue))

plot_test <- ggplot(data=combined_rates_df, aes(x=tissue, y=rates)) + geom_boxplot(notch=TRUE) + ylim(0,0.001)
