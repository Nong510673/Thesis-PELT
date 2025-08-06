# library(devtools)


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# BiocManager::install("S4Vectors")

#install_github("poyuliu/KTU")

# 必要的包
library(Biostrings)
library(KTU)
library(tidyr) 
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pwalign)


# 设置工作目录
wd <- "E:/轉KTU" # 输入文件存储位置
setwd(wd)

# 导入 ASV 表，保留 #OTU ID 为数据框的一列
asv <- read.csv("feature-table_M_D.csv", header = TRUE, check.names = FALSE)

# 将 #OTU ID 修改为与 metadata 一致的 X.OTU.ID
colnames(asv)[1] <- "X.OTU.ID"

# 检查列名
colnames(asv)

# 导入 fasta 文件
sequences <- readDNAStringSet("dna-sequences_M_D.fasta")

# 导入 metadata 表
metadata <- read.csv("metadata.csv")
#metadata <- metadata[-1, ] # 去除无用的第一行

# 分离 Taxon 列为 Kingdom, Phylum 等分类列
metadata <- metadata %>%
  separate(Taxon, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), 
           sep = "; ", 
           fill = "right", 
           extra = "merge")

# 按 Phylum 筛选 metadata 子集
metadata_p__Chloroflexota <- metadata %>% filter(Phylum == "p__Chloroflexota")
metadata_p__Bacillota <- metadata %>% filter(Phylum == "p__Bacillota")
metadata_p__Pseudomonadota <- metadata %>% filter(Phylum == "p__Pseudomonadota")
metadata_p__Planctomycetota <- metadata %>% filter(Phylum == "p__Planctomycetota")
metadata_p__Bacteroidota <- metadata %>% filter(Phylum == "p__Bacteroidota")

# 剔除特定 Phylum
ids_to_remove <- c("p__Chloroflexota", "p__Bacillota", "p__Pseudomonadota", "p__Planctomycetota", "p__Bacteroidota")
metadata_p__resudual <- metadata %>% filter(!Phylum %in% ids_to_remove)

# 从 ASV 表中筛选子集
asv_p__Chloroflexota <- asv %>% filter(X.OTU.ID %in% metadata_p__Chloroflexota$X.OTU.ID)
asv_p__Bacillota <- asv %>% filter(X.OTU.ID %in% metadata_p__Bacillota$X.OTU.ID)
asv_p__Pseudomonadota <- asv %>% filter(X.OTU.ID %in% metadata_p__Pseudomonadota$X.OTU.ID)
asv_p__Planctomycetota <- asv %>% filter(X.OTU.ID %in% metadata_p__Planctomycetota$X.OTU.ID)
asv_p__Bacteroidota <- asv %>% filter(X.OTU.ID %in% metadata_p__Bacteroidota$X.OTU.ID)
asv_p__resudual <- asv %>% filter(X.OTU.ID %in% metadata_p__resudual$X.OTU.ID)

# 检查列的类型并转换为数值类型（忽略 X.OTU.ID 列）
convert_to_numeric <- function(df) {
  character_cols <- colnames(df)[sapply(df, is.character)]
  character_cols <- setdiff(character_cols, "X.OTU.ID")
  df[character_cols] <- lapply(df[character_cols], function(x) as.numeric(as.character(x)))
  return(df)
}

asv_p__Chloroflexota <- convert_to_numeric(asv_p__Chloroflexota)
asv_p__Bacillota <- convert_to_numeric(asv_p__Bacillota)
asv_p__Pseudomonadota <- convert_to_numeric(asv_p__Pseudomonadota)
asv_p__Planctomycetota <- convert_to_numeric(asv_p__Planctomycetota)
asv_p__Bacteroidota <- convert_to_numeric(asv_p__Bacteroidota)
asv_p__resudual <- convert_to_numeric(asv_p__resudual)

# 从 sequences 中筛选对应的序列
sequences_p__Chloroflexota <- sequences[names(sequences) %in% metadata_p__Chloroflexota$X.OTU.ID]
sequences_p__Bacillota <- sequences[names(sequences) %in% metadata_p__Bacillota$X.OTU.ID]
sequences_p__Pseudomonadota <- sequences[names(sequences) %in% metadata_p__Pseudomonadota$X.OTU.ID]
sequences_p__Planctomycetota <- sequences[names(sequences) %in% metadata_p__Planctomycetota$X.OTU.ID]
sequences_p__Bacteroidota <- sequences[names(sequences) %in% metadata_p__Bacteroidota$X.OTU.ID]
sequences_p__resudual <- sequences[names(sequences) %in% metadata_p__resudual$X.OTU.ID]

# 将筛选后的序列保存为新 fasta 文件
writeXStringSet(sequences_p__Chloroflexota, "dna-sequences-Chloroflexota.fasta")
writeXStringSet(sequences_p__Bacillota, "dna-sequences-Bacillota.fasta")
writeXStringSet(sequences_p__Pseudomonadota, "dna-sequences-Pseudomonadota.fasta")
writeXStringSet(sequences_p__Planctomycetota, "dna-sequences-Planctomycetota.fasta")
writeXStringSet(sequences_p__Bacteroidota, "dna-sequences-Bacteroidota.fasta")
writeXStringSet(sequences_p__resudual, "dna-sequences-resudual.fasta")

# 检查结果
print("FASTA 文件保存完成。")

#Run KTU

#Chloroflexi
kluster.Chloroflexota <- klustering(repseq = "dna-sequences-Chloroflexota.fasta",
                                    feature.table = asv_p__Chloroflexota,
                                    write.fasta = TRUE,cores = 10)
saveRDS(kluster.Chloroflexota,file="kluster-Chloroflexota.RDS") 
ktu_eval.Chloroflexota <- KTUsim.eval(klusterRDS = "kluster-Chloroflexota.RDS", ASVfasta = "dna-sequences-Chloroflexota.fasta")
write.csv(kluster.Chloroflexota$KTU.table,file="KTU.table.Chloroflexota.csv")
write.csv(kluster.Chloroflexota$ReqSeq ,file="rep-seq.Chloroflexota.csv")
write.csv(kluster.Chloroflexota$clusters ,file="cluster.Chloroflexota.csv")

#Firmicutes
kluster.Bacillota <- klustering(repseq = "dna-sequences-Bacillota.fasta",
                                feature.table = asv_p__Bacillota,
                                write.fasta = TRUE,cores = 10)
saveRDS(kluster.Bacillota,file="kluster-Bacillota.RDS") 
ktu_eval.Bacillota <- KTUsim.eval(klusterRDS = "kluster-Bacillota.RDS", ASVfasta = "dna-sequences-Bacillota.fasta")
write.csv(kluster.Bacillota$KTU.table,file="KTU.table.Bacillota.csv")
write.csv(kluster.Bacillota$ReqSeq ,file="rep-seq.Bacillota.csv")
write.csv(kluster.Bacillota$clusters ,file="cluster.Bacillota.csv")

#Proteobacteria
kluster.Pseudomonadota <- klustering(repseq = "dna-sequences-Pseudomonadota.fasta",
                                     feature.table = asv_p__Pseudomonadota,
                                     write.fasta = TRUE,cores = 10)
saveRDS(kluster.Pseudomonadota,file="kluster-Pseudomonadota.RDS") 
ktu_eval.Pseudomonadota <- KTUsim.eval(klusterRDS = "kluster-Pseudomonadota.RDS", ASVfasta = "dna-sequences-Pseudomonadota.fasta")
write.csv(kluster.Pseudomonadota$KTU.table,file="KTU.table.Pseudomonadota.csv")
write.csv(kluster.Pseudomonadota$ReqSeq ,file="rep-seq.Pseudomonadota.csv")
write.csv(kluster.Pseudomonadota$clusters ,file="cluster.Pseudomonadota.csv")

#Planctomycetota
kluster.Planctomycetota <- klustering(repseq = "dna-sequences-Planctomycetota.fasta",
                                      feature.table = asv_p__Planctomycetota,
                                      write.fasta = TRUE,cores = 10)
saveRDS(kluster.Planctomycetota,file="kluster-Planctomycetota.RDS") 
ktu_eval.Planctomycetota <- KTUsim.eval(klusterRDS = "kluster-Planctomycetota.RDS", ASVfasta = "dna-sequences-Planctomycetota.fasta")
write.csv(kluster.Planctomycetota$KTU.table,file="KTU.table.Planctomycetota.csv")
write.csv(kluster.Planctomycetota$ReqSeq ,file="rep-seq.Planctomycetota.csv")
write.csv(kluster.Planctomycetota$clusters ,file="cluster.Planctomycetota.csv")

#Bacteroidota
kluster.Bacteroidota <- klustering(repseq = "dna-sequences-Bacteroidota.fasta",
                                   feature.table = asv_p__Bacteroidota,
                                   write.fasta = TRUE,cores = 10)
saveRDS(kluster.Bacteroidota,file="kluster-Bacteroidota.RDS") 
ktu_eval.Bacteroidota <- KTUsim.eval(klusterRDS = "kluster-Bacteroidota.RDS", ASVfasta = "dna-sequences-Bacteroidota.fasta")
write.csv(kluster.Bacteroidota$KTU.table,file="KTU.table.Bacteroidota.csv")
write.csv(kluster.Bacteroidota$ReqSeq ,file="rep-seq.Bacteroidota.csv")
write.csv(kluster.Bacteroidota$clusters ,file="cluster.Bacteroidota.csv")

#resudual
kluster.resudual <- klustering(repseq = "dna-sequences-resudual.fasta",
                               feature.table = asv_p__resudual,
                               write.fasta = TRUE,cores = 10)
saveRDS(kluster.resudual,file="kluster-resudual.RDS") 
ktu_eval.resudual <- KTUsim.eval(klusterRDS = "kluster-resudual.RDS", ASVfasta = "dna-sequences-resudual.fasta")
write.csv(kluster.resudual$KTU.table,file="KTU.table.resudual.csv")
write.csv(kluster.resudual$ReqSeq ,file="rep-seq.resudual.csv")
write.csv(kluster.resudual$clusters ,file="cluster.resudual.csv")




#Merge the sub result of the KTU 
KTU.table.summary=rbind(kluster.Chloroflexota$KTU.table,kluster.Bacillota$KTU.table,kluster.Pseudomonadota$KTU.table,
                        kluster.Planctomycetota$KTU.table,kluster.Bacteroidota$KTU.table,kluster.resudual$KTU.table)

kmer.table.summary=cbind(kluster.Chloroflexota$kmer.table,kluster.Bacillota$kmer.table,kluster.Pseudomonadota$kmer.table,
                         kluster.Planctomycetota$kmer.table,kluster.Bacteroidota$kmer.table,kluster.resudual$kmer.table)

ktu_eval.summary=rbind(ktu_eval.Chloroflexota$eachmean,ktu_eval.Bacillota$eachmean,ktu_eval.Pseudomonadota$eachmean,
                       ktu_eval.Planctomycetota$eachmean,ktu_eval.Bacteroidota$eachmean,ktu_eval.resudual$eachmean)

write.csv(KTU.table.summary,file="KTU.table.summary.csv")
write.csv(kmer.table.summary,file="kmer.table.summary.csv")
write.csv(ktu_eval.summary,file="ktu_evaluate.summary.csv")


# Convert count KTU table to relative abundance table
KTU.table.summary_relative_abundance <- KTU.table.summary %>%
  mutate(across(everything(), ~ . / sum(.)))
#write.csv(KTU.table.summary_relative_abundance,file="KTU.table.summary_relative_abundance.csv")


#histogram
#how many ASVs are clustered in to a KTU
p.n.feature.histogram=
  ggplot(ktu_eval.summary, aes(x = n.feature)) +
  geom_histogram(binwidth = 1, fill = "navyblue", color = "black") +
  labs(x = "# of ASVs / KTU", y = "Frequency") +
  #scale_x_continuous(limits = c(0, 15),breaks=seq(0,15,1))+
  theme(axis.text.y = element_text(colour = "black", size = 18), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.text = element_text(size = 20, colour ="black"), 
        #legend.position = c(0.6,0.88),
        legend.position = "bottom",
        legend.direction = "horizontal", ##Vertical  #horizontal
        legend.box = "vertical", ##Vertical  #horizontal
        axis.title.y = element_text( size = 22), 
        axis.title.x = element_text( size = 20, colour = "black"),
        #axis.title.x = element_blank(),#hide label of x axis
        legend.title = element_text(size = 20, colour = "black"), 
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2)#框線
  )

#mean ASV similarity within a KTU
p.similarity.pct.histogram=
  ggplot(ktu_eval.summary, aes(x = similarity.pct)) +
  geom_histogram(binwidth = 0.5, fill = "navyblue", color = "black") +
  labs(x = "% of ASVs within KTU", y = "Frequency") +
  theme(axis.text.y = element_text(colour = "black", size = 18), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.text = element_text(size = 20, colour ="black"), 
        #legend.position = c(0.6,0.88),
        legend.position = "bottom",
        legend.direction = "horizontal", ##Vertical  #horizontal
        legend.box = "vertical", ##Vertical  #horizontal
        axis.title.y = element_text( size = 22), 
        axis.title.x = element_text( size = 20, colour = "black"),
        #axis.title.x = element_blank(),#hide label of x axis
        legend.title = element_text(size = 20, colour = "black"), 
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2)#框線
  )

#Divergence of ASVs within KTU 
p.divergence.histogram=
  ggplot(ktu_eval.summary, aes(x = divergence)) +
  geom_histogram(binwidth = 0.002, fill = "navyblue", color = "black") +
  labs(x = "Divergence of ASVs within KTU", y = "Frequency") +
  theme(axis.text.y = element_text(colour = "black", size = 18), #, 1face = "bold"
        axis.text.x = element_text(colour = "black", size = 18), 
        legend.text = element_text(size = 20, colour ="black"), 
        #legend.position = c(0.6,0.88),
        legend.position = "bottom",
        legend.direction = "horizontal", ##Vertical  #horizontal
        legend.box = "vertical", ##Vertical  #horizontal
        axis.title.y = element_text(size = 22), 
        axis.title.x = element_text(size = 20, colour = "black"),
        #axis.title.x = element_blank(),#hide label of x axis
        legend.title = element_text(size = 20, colour = "black"), 
        axis.ticks=element_line(colour="black",size=1,linetype=1), #Modify the size and color of the line label on the axis
        axis.ticks.length=unit(0.4,"lines"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 2)#框線
  )


P.merge.KTU.histogram=
  ggarrange(p.n.feature.histogram,p.similarity.pct.histogram,p.divergence.histogram,
            align = "hv",
            labels = c("A", "B", "C", "D"),
            ncol=1, nrow=3,
            #common.legend = TRUE, legend="bottom",
            font.label = list(size = 30))

plot(p.n.feature.histogram)
plot(p.similarity.pct.histogram)
plot(p.divergence.histogram)
plot(P.merge.KTU.histogram)

ggsave(plot = P.merge.KTU.histogram,"01_P.merge.KTU.histogram.png", width = 5.5, height = 11)



# #how many ASVs are clustered in to a KTU
# histNdistri(histdata = ktu_eval$eachmean$n.feature,las=1,breaks=30,
#             coll = colset.d.5[1],
#             xlab.text = "# of ASVs / KTU",
#             legend.posit = "topright",legend.text = paste0(round(mean(ktu_eval$eachmean$n.feature),2)," ASVs/KTU on average"))
# #mean ASV similarity within a KTU
# histNdistri(histdata = ktu_eval$eachmean$similarity.pct,las=1,breaks=30,
#             coll = colset.d.5[2],
#             xlab.text = "% of ASVs within KTU",
#             legend.posit = "topleft",legend.text = paste0("Mean similarity: ",round(ktu_eval$globalmean,2),"%"))
# #Divergence of ASVs within KTU  
# histNdistri(histdata = ktu_eval$eachmean$divergence[!is.na(ktu_eval$eachmean$divergence)],las=1,breaks=30,
#             coll = colset.d.5[3],
#             xlab.text = "Divergence of ASVs within KTU",
#             legend.posit = "topright",legend.text = paste0("Mean divergence: ",round(ktu_eval$globaldivergence,2)))
