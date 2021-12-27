#Autosomes
#read files
l1f_autosomes <- read.table("L1F_TPM_ATACcov_correlation.txt", header=F, sep="\t")
l2f_autosomes <- read.table("L2F_TPM_ATACcov_correlation.txt", header=F, sep="\t")

l1m_autosomes <- read.table("L1M_TPM_ATACcov_correlation.txt", header=F, sep="\t")
l2m_autosomes <- read.table("L2M_TPM_ATACcov_correlation.txt", header=F, sep="\t")

o1_autosomes <- read.table("O1_TPM_ATACcov_correlation.txt", header=F, sep="\t")
o2_autosomes <- read.table("O2_TPM_ATACcov_correlation.txt", header=F, sep="\t")

t1_autosomes <- read.table("T1_TPM_ATACcov_correlation.txt", header=F, sep="\t")
t2_autosomes <- read.table("T2_TPM_ATACcov_correlation.txt", header=F, sep="\t")

#correlations

cor(l1f_autosomes$V2, l1f_autosomes$V3, method="spearman")
cor.test(l1f_autosomes$V2, l1f_autosomes$V3, method="spearman")
#p =
#rho =

cor(l2f_autosomes$V2, l2f_autosomes$V3, method="spearman")
cor.test(l2f_autosomes$V2, l2f_autosomes$V3, method="spearman")
#p =
#rho =

cor(l1m_autosomes$V2, l1m_autosomes$V3, method="spearman")
cor.test(l1m_autosomes$V2, l1m_autosomes$V3, method="spearman")
#p =
#rho =

cor(l2m_autosomes$V2, l2m_autosomes$V3, method="spearman")
cor.test(l2m_autosomes$V2, l2m_autosomes$V3, method="spearman")
#p =
#rho =

cor(o1_autosomes$V2, o1_autosomes$V3, method="spearman")
cor.test(o1_autosomes$V2, o1_autosomes$V3, method="spearman")
#p =
#rho =

cor(o2_autosomes$V2, o2_autosomes$V3, method="spearman")
cor.test(o2_autosomes$V2, o2_autosomes$V3, method="spearman")
#p =
#rho =

cor(t1_autosomes$V2, t1_autosomes$V3, method="spearman")
cor.test(t1_autosomes$V2, t1_autosomes$V3, method="spearman")
#p =
#rho =

cor(t2_autosomes$V2, t2_autosomes$V3, method="spearman")
cor.test(t2_autosomes$V2, t2_autosomes$V3, method="spearman")
#p =
#rho =

jpeg(file="ATAC_readcov_vs_TPM_V1.jpeg",width=1400,height=2000,quality=100)

l1m_plot <- ggplot(data=l1m_autosomes, aes(x=V2, y=V3)) + geom_point(stat="identity") + xlab("Average ATAC-seq\nRead Coverage") + ylab("TPM") + ylim(0,5000) +theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), text=element_text(size=30),panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(color="black"), axis.title.y=element_text(size=35, margin=margin(t=1,r=20,b=1,l=1)),axis.title.x=element_text(size=35,margin=margin(t=20,r=1,b=1,l=1)), plot.margin=unit(c(2,0.5,0.5,0.5),"cm")) + xlim(0,100)

l1m_plot_annotated <- l1m_plot + annotate("text",label="r = -0.027", x=80, y=4750, size=10) + annotate("text",label="p < 0.001",x=80, y=4250, size=10)

l2m_plot <- ggplot(data=l2m_autosomes, aes(x=V2, y=V3)) + geom_point(stat="identity") + xlab("Average ATAC-seq\nRead Coverage") + ylab("TPM") + ylim(0,5000) + theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), text=element_text(size=30),panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(color="black"), axis.title.y=element_text(size=35, margin=margin(t=1,r=20,b=1,l=1)),axis.title.x=element_text(size=35,margin=margin(t=20,r=1,b=1,l=1)), plot.margin=unit(c(2,0.5,0.5,0.5),"cm"))  + xlim(0,175)

l2m_plot_annotated <- l2m_plot + annotate("text",label="r = -0.065",x=150, y=4750, size=10) + annotate("text",label="p < 0.001",x=150, y=4250, size=10)

l1f_plot <- ggplot(data=l1f_autosomes, aes(x=V2, y=V3)) + geom_point(stat="identity") + xlab("Average ATAC-seq\nRead Coverage") + ylab("TPM") + ylim(0,5000) + theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), text=element_text(size=30),panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(color="black"), axis.title.y=element_text(size=35, margin=margin(t=1,r=20,b=1,l=1)),axis.title.x=element_text(size=35,margin=margin(t=20,r=1,b=1,l=1)), plot.margin=unit(c(2,0.5,0.5,0.5),"cm"))  + xlim(0,200)

l1f_plot_annotated <- l1f_plot + annotate("text",label="r = 0.026",x=155, y=3750, size=10) + annotate("text",label="p < 0.001",x=155, y=3250, size=10)

l2f_plot <- ggplot(data=l2f_autosomes, aes(x=V2, y=V3)) + geom_point(stat="identity") + xlab("Average ATAC-seq\nRead Coverage") + ylab("TPM") + ylim(0,5000) + theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), text=element_text(size=30),panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(color="black"), axis.title.y=element_text(size=35,margin=margin(t=1,r=20,b=1,l=1)),axis.title.x=element_text(size=35,margin=margin(t=20,r=1,b=1,l=1)), plot.margin=unit(c(2,0.5,0.5,0.5),"cm"))  + xlim(0,200)

l2f_plot_annotated <- l2f_plot + annotate("text",label="r = -0.009",x=155, y=3750, size=10) + annotate("text",label="p = 0.143",x=155, y=3250, size=10)

o1_plot <- ggplot(data=o1_autosomes, aes(x=V2, y=V3)) + geom_point(stat="identity") + xlab("Average ATAC-seq\nRead Coverage") + ylab("TPM") + ylim(0,5000) +theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), text=element_text(size=30),panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(color="black"), axis.title.y=element_text(size=35, margin=margin(t=1,r=20,b=1,l=1)),axis.title.x=element_text(size=35,margin=margin(t=20,r=1,b=1,l=1)), plot.margin=unit(c(2,0.5,0.5,0.5),"cm")) + xlim(0,100)

o1_plot_annotated <- o1_plot + annotate("text",label="r = 0.034", x=80, y=4750, size=10) + annotate("text",label="p < 0.001",x=80, y=4250, size=10)

o2_plot <- ggplot(data=o2_autosomes, aes(x=V2, y=V3)) + geom_point(stat="identity") + xlab("Average ATAC-seq\nRead Coverage") + ylab("TPM") + ylim(0,5000) +theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), text=element_text(size=30),panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(color="black"), axis.title.y=element_text(size=35, margin=margin(t=1,r=20,b=1,l=1)),axis.title.x=element_text(size=35,margin=margin(t=20,r=1,b=1,l=1)), plot.margin=unit(c(2,0.5,0.5,0.5),"cm")) + xlim(0,100)

o2_plot_annotated <- o2_plot + annotate("text",label="r = -0.057", x=80, y=4750, size=10) + annotate("text",label="p < 0.001",x=80, y=4250, size=10)

t1_plot <- ggplot(data=t1_autosomes, aes(x=V2, y=V3)) + geom_point(stat="identity") + xlab("Average ATAC-seq\nRead Coverage") + ylab("TPM") + ylim(0,5000) +theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), text=element_text(size=30),panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(color="black"), axis.title.y=element_text(size=35, margin=margin(t=1,r=20,b=1,l=1)),axis.title.x=element_text(size=35,margin=margin(t=20,r=1,b=1,l=1)), plot.margin=unit(c(2,0.5,0.5,0.5),"cm")) + xlim(0,100)

t1_plot_annotated <- t1_plot + annotate("text",label="r = -0.0002", x=80, y=4750, size=10) + annotate("text",label="p = 0.978",x=80, y=4250, size=10)

t2_plot <- ggplot(data=t2_autosomes, aes(x=V2, y=V3)) + geom_point(stat="identity") + xlab("Average ATAC-seq\nRead Coverage") + ylab("TPM") + ylim(0,5000) +theme_bw() + theme(axis.text.x=element_text(size=40), axis.text.y=element_text(size=40), text=element_text(size=30),panel.border=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(color="black"), axis.title.y=element_text(size=35, margin=margin(t=1,r=20,b=1,l=1)),axis.title.x=element_text(size=35,margin=margin(t=20,r=1,b=1,l=1)), plot.margin=unit(c(2,0.5,0.5,0.5),"cm")) + xlim(0,100)

t2_plot_annotated <- t2_plot + annotate("text",label="r = -0.057", x=80, y=4750, size=10) + annotate("text",label="p < 0.001",x=80, y=4250, size=10)

row_1 <- plot_grid(l1m_plot_annotated, l2m_plot_annotated,ncol=2, nrow=1)
row_2 <- plot_grid(l1f_plot_annotated, l2f_plot_annotated, ncol=2, nrow=1)
row_3 <- plot_grid(o1_plot_annotated, o2_plot_annotated, ncol=2, nrow=1)
row_4 <- plot_grid(t1_plot_annotated, t2_plot_annotated, ncol=2, nrow=1)

plot_grid(row_1, row_2, row_3, row_4, ncol=1,labels=c("A","B", "C", "D"),label_size=40)

dev.off()
