library(tidyverse)
library(reshape2)
library(DESeq2)
library(ggpubr)
library(TTR)
library(magrittr)


# basic count data --------------------------------------------------------------

count_matrix <- read.table('Osy63_parents.matrix', header = T, row.names = 1)
snp_genes <- read.table('snp_genes.bed', header = F)
count_matrix <- count_matrix %>% rownames_to_column('ID') %>% filter(ID %in% snp_genes$V4) %>% column_to_rownames('ID')
ordered_tissues <- c('leaf', 'panicle', 'root')
rna_expression_matrix <- read.table('refOmh63_rna.gene.TPM.not_cross_norm.avg_reps.byLog.matrix', header = T, row.names = 1)


# leaf delta--------------------------------------------------------------
leaf_contrast_count_matrix <- count_matrix %>% select(contains('leaf'))
leaf_contrast_count_matrix <- leaf_contrast_count_matrix[which(rowSums(leaf_contrast_count_matrix)>=1),]
leaf_coldata <- data.frame(SampleID=names(leaf_contrast_count_matrix))
leaf_coldata$Condition <- str_split(leaf_coldata$SampleID, pattern = '_') %>% map_chr(~.x[length(.x)])
leaf_coldata$Condition <- factor(leaf_coldata$Condition, levels = c('MH63', 'ZS97'))
leaf_coldata$SeqType <- str_split(leaf_coldata$SampleID, pattern = '_') %>% map_chr(~.x[3])
leaf_coldata$SeqType <- factor(leaf_coldata$SeqType, levels = c('RNA','Ribo'))
leaf_coldata$Batch <- rep(1,8)

leaf_ddsMat <- DESeqDataSetFromMatrix(countData = leaf_contrast_count_matrix, colData = leaf_coldata, design =~ Condition + SeqType + Condition:SeqType)
leaf_ddsMat$SeqType <- relevel(leaf_ddsMat$SeqType,"RNA")
leaf_ddsMat <- DESeq(leaf_ddsMat)
# resultsNames(leaf_ddsMat)
# [1] "Intercept"                 "Condition_ZS97_vs_MH63"    "SeqType_Ribo_vs_RNA"      
# [4] "ConditionZS97.SeqTypeRibo"
leaf_res <- results(object = leaf_ddsMat, name = resultsNames(leaf_ddsMat)[4])
leaf_res$padj[is.na(leaf_res$padj)]  <- 1
leaf_res$log2FoldChange[is.na(leaf_res$log2FoldChange)] <- 0
write.table(rownames(leaf_res)[which(leaf_res$padj <= 0.05 & abs(leaf_res$log2FoldChange)>=1)],"phasing/leaf/gene_lists/DTEGs.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(leaf_res[which(leaf_res$padj <= 0.05 & abs(leaf_res$log2FoldChange)>=1),],"phasing/leaf/fold_changes/deltaTE.txt",quote=F,sep="\t",col.names = T,row.names = T)

## DESeq2 for Ribo-seq
ind = which(leaf_coldata$SeqType == "Ribo")
leaf_coldata_ribo = leaf_coldata[ind,]
leaf_contrast_ribo_count_matrix <- leaf_contrast_count_matrix %>% select(contains('Ribo'))
leaf_ddsMat_ribo <- DESeqDataSetFromMatrix(countData = leaf_contrast_ribo_count_matrix, colData = leaf_coldata_ribo, design = ~ Condition)
leaf_ddsMat_ribo <- DESeq(leaf_ddsMat_ribo)
leaf_res_ribo <- results(leaf_ddsMat_ribo, contrast = c('Condition', 'ZS97', 'MH63'))
leaf_res_ribo <- lfcShrink(dds = leaf_ddsMat_ribo, coef = resultsNames(leaf_ddsMat_ribo)[2], res = leaf_res_ribo)
leaf_res_ribo$padj[is.na(leaf_res_ribo$padj)] <- 1
leaf_res_ribo$log2FoldChange[is.na(leaf_res_ribo$log2FoldChange)] <- 0
write.table(leaf_res_ribo[which(leaf_res_ribo$padj <= 0.05 & abs(leaf_res_ribo$log2FoldChange)>=1),],"phasing/leaf/fold_changes/deltaRibo.txt",quote=F,sep="\t",col.names = T,row.names = T)

## DESeq2 for RNA-seq
ind = which(leaf_coldata$SeqType == "RNA")
leaf_coldata_rna = leaf_coldata[ind,]
leaf_contrast_rna_count_matrix <- leaf_contrast_count_matrix %>% select(contains('RNA'))
leaf_ddsMat_rna <- DESeqDataSetFromMatrix(countData = leaf_contrast_rna_count_matrix, colData = leaf_coldata_rna, design = ~ Condition)
leaf_ddsMat_rna <- DESeq(leaf_ddsMat_rna)
leaf_res_rna <- results(leaf_ddsMat_rna, contrast = c('Condition', 'ZS97', 'MH63'))
leaf_res_rna <- lfcShrink(dds = leaf_ddsMat_rna, coef = resultsNames(leaf_ddsMat_rna)[2], res = leaf_res_rna)
leaf_res_rna$padj[is.na(leaf_res_rna$padj)]  <- 1
leaf_res_rna$log2FoldChange[is.na(leaf_res_rna$log2FoldChange)] <- 0
length(which(leaf_res_rna$padj <= 0.05 & abs(leaf_res_rna$log2FoldChange)>=1))
write.table(rownames(leaf_res_rna)[which(leaf_res_rna$padj <= 0.05 & abs(leaf_res_rna$log2FoldChange)>=1)],"phasing/leaf/gene_lists/DTG.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(leaf_res_rna[which(leaf_res_rna$padj<=0.05 & abs(leaf_res_rna$log2FoldChange)>=1),],"phasing/leaf/fold_changes/deltaRNA.txt",quote=F,sep="\t",col.names = T,row.names = T)

# leaf gene RNA TE FCplot -------------------------------------------------
leaf_aste_res <- leaf_res[which(leaf_res$padj <= 0.05 & abs(leaf_res$log2FoldChange)>=1),]
leaf_aste_genes <- rownames(leaf_aste_res)
leaf_asribo_res <- leaf_res_ribo[which(leaf_res_ribo$padj <= 0.05 & abs(leaf_res_ribo$log2FoldChange)>=1),]
leaf_asribo_genes <- rownames(leaf_asribo_res)
leaf_ase_res <- leaf_res_rna[which(leaf_res_rna$padj <= 0.05 & abs(leaf_res_rna$log2FoldChange)>=1),]
leaf_ase_genes <- rownames(leaf_ase_res)

leaf_mRNAOnly <- rownames(leaf_res)[which((leaf_res$padj > 0.05 | abs(leaf_res$log2FoldChange)<1) & (leaf_res_rna$padj <= 0.05 & abs(leaf_res_rna$log2FoldChange)>=1))]
leaf_TEOnly <- rownames(leaf_res)[which((leaf_res$padj <= 0.05 & abs(leaf_res$log2FoldChange)>=1) & (leaf_res_rna$padj > 0.05 | abs(leaf_res_rna$log2FoldChange)<1))]
leaf_FCboth <- which((leaf_res$padj <= 0.05 & abs(leaf_res$log2FoldChange)>=1) & (leaf_res_rna$padj <= 0.05 & abs(leaf_res_rna$log2FoldChange)>=1))
leaf_compensatory <- rownames(leaf_res)[leaf_FCboth[which(leaf_res[leaf_FCboth,2]*leaf_res_rna[leaf_FCboth,2] <0 )]]
leaf_reinforcing <- rownames(leaf_res)[leaf_FCboth[which(leaf_res[leaf_FCboth,2]*leaf_res_rna[leaf_FCboth,2] >=0 )]]
leaf_FCothers <- rownames(leaf_res)[! rownames(leaf_res) %in% c(leaf_mRNAOnly,leaf_TEOnly,leaf_reinforcing, leaf_compensatory)]

leaf_fc <- data.frame(fold_te=leaf_res[,2], fold_rna=leaf_res_rna[,2]) %>% na.omit()
nrow(leaf_fc)
leaf_fc$class[rownames(leaf_fc) %in% leaf_mRNAOnly] <- 'mRNA_only'
leaf_fc$class[rownames(leaf_fc) %in% leaf_TEOnly] <- 'TE_only'
leaf_fc$class[rownames(leaf_fc) %in% leaf_compensatory] <- 'Compensatory'
leaf_fc$class[rownames(leaf_fc) %in% leaf_reinforcing] <- 'Reinforcing'
leaf_fc$class[rownames(leaf_fc) %in% leaf_FCothers] <- 'Others'
leaf_fc$class <- factor(leaf_fc$class, levels = c('mRNA_only','TE_only','Compensatory','Reinforcing','Others'))
leaf_fc$class %>% table()
p1 <- leaf_fc %>% ggplot(aes(x=fold_rna, y=fold_te, color=class)) + geom_point(size=1.5) + theme_bw(base_size = 15) + scale_color_manual(values = c('#4169E1', '#CC524F', '#DAA520', '#3CB371', 'grey')) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text = element_text(color='black'), plot.margin = unit(rep(2,4),'lines')) + labs(x=expression(paste(Log[2],'(',mRNA[ZS97],'/',mRNA[MH63],')')), y=expression(paste(Log[2],'(',TE[ZS97],'/',TE[MH63],')'))) + 
  geom_vline(xintercept = 0, size=0.8, linetype='dashed') + 
  geom_vline(xintercept = -1, size=0.8, linetype='dashed', color='#696969') +
  geom_vline(xintercept = 1, size=0.8, linetype='dashed', color='#696969') + 
  geom_hline(yintercept = 0, size=0.8, linetype='dashed') +
  geom_hline(yintercept = -1, size=0.8, linetype='dashed', color='#696969') + 
  geom_hline(yintercept = 1, size=0.8, linetype='dashed', color='#696969')
p2 <- leaf_fc %>% select(class) %>% table() %>% as.data.frame() %>% magrittr::set_colnames(c('class','Freq')) %>% filter(class!='Others') %>% ggplot(aes(x=class,y=Freq)) + geom_bar(aes(fill=class), stat = 'identity', width = 0.4) + theme_bw(base_size = 15) + coord_flip() + scale_x_discrete(limits=rev(c('mRNA_only','TE_only','Compensatory','Reinforcing'))) + scale_fill_manual(values = c('#4169E1', '#CC524F', '#DAA520', '#3CB371')) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text = element_text(color='black'), plot.margin = unit(rep(2,4),'lines')) + labs(x='', y='Number of gene') + scale_y_continuous(expand = c(0,0))
ggsave(ggarrange(p1, p2, ncol = 2, nrow = 1), filename = 'phasing/leaf/leaf_FCplot.pdf', width = 14, height = 7)

write.table(leaf_mRNAOnly,"phasing/leaf/gene_lists/leaf_mRNAOnly.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(leaf_TEOnly,"phasing/leaf/gene_lists/leaf_TEOnly.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(leaf_compensatory,"phasing/leaf/gene_lists/leaf_compensatory.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(leaf_reinforcing,"phasing/leaf/gene_lists/leaf_reinforcing.txt",quote=F,sep="\t",col.names = F,row.names = F)


## leaf total mRNA expression level ----------------------------------------
leaf_expression <- rna_expression_matrix %>% select(contains('Osy63_leaf')) %>% rownames_to_column('ID') %>% filter(ID %in% rownames(leaf_contrast_count_matrix)) %>% column_to_rownames('ID')
names(leaf_expression) <- 'mRNA'
leaf_expression$group[rownames(leaf_expression) %in% leaf_mRNAOnly] <- 'mRNAOnly'
leaf_expression$group[rownames(leaf_expression) %in% leaf_TEOnly] <- 'TEOnly'
leaf_expression$group[! rownames(leaf_expression) %in% c(leaf_mRNAOnly,leaf_TEOnly)] <- 'Others'
leaf_expression$group <- factor(leaf_expression$group, levels = c('Others','mRNAOnly','TEOnly'))
leaf_expression %>% mutate(m=log10(mRNA+1)) %>% ggboxplot(x='group',y='m', width = 0.5,color='grey', fill = 'group', palette = c('#808080',rgb(208,128,91,maxColorValue = 255), rgb(94,160,205,maxColorValue = 255))) + stat_compare_means(comparisons = list(c('mRNAOnly','Others'),c('mRNAOnly','TEOnly'),c('TEOnly','Others'))) + labs(x='',y=expression(paste(Log[10],'(Total mRNA TPM+1)'))) + theme_bw(base_size = 15) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text = element_text(color = 'black')) + stat_compare_means(label.y = 8) + ggsave(filename = 'phasing/leaf/leaf_divExpression.pdf', width = 7, height = 7)


# panicle delta -----------------------------------------------------------------

panicle_contrast_count_matrix <- count_matrix %>% select(contains('panicle'))
panicle_contrast_count_matrix <- panicle_contrast_count_matrix[which(rowSums(panicle_contrast_count_matrix)>=1),]
panicle_coldata <- data.frame(SampleID=names(panicle_contrast_count_matrix))
panicle_coldata$Condition <- str_split(panicle_coldata$SampleID, pattern = '_') %>% map_chr(~.x[length(.x)])
panicle_coldata$Condition <- factor(panicle_coldata$Condition, levels = c('MH63', 'ZS97'))
panicle_coldata$SeqType <- str_split(panicle_coldata$SampleID, pattern = '_') %>% map_chr(~.x[3])
panicle_coldata$SeqType <- factor(panicle_coldata$SeqType, levels = c('RNA','Ribo'))
panicle_coldata$Batch <- rep(1,8)
panicle_ddsMat <- DESeqDataSetFromMatrix(countData = panicle_contrast_count_matrix, colData = panicle_coldata, design =~ Condition + SeqType + Condition:SeqType)
panicle_ddsMat$SeqType <- relevel(panicle_ddsMat$SeqType,"RNA")
panicle_ddsMat <- DESeq(panicle_ddsMat)
# resultsNames(panicle_ddsMat)
# [1] "Intercept"                 "Condition_ZS97_vs_MH63"    "SeqType_Ribo_vs_RNA"      
# [4] "ConditionZS97.SeqTypeRibo"
panicle_res <- results(object = panicle_ddsMat, name = resultsNames(panicle_ddsMat)[4])
panicle_res$padj[is.na(panicle_res$padj)] <- 1
panicle_res$log2FoldChange[is.na(panicle_res$log2FoldChange)] <- 0
length(which(panicle_res$padj <= 0.05 & abs(panicle_res$log2FoldChange)>=1))
write.table(rownames(panicle_res)[which(panicle_res$padj <= 0.05 & abs(panicle_res$log2FoldChange)>=1)],"phasing/panicle/gene_lists/DTEGs.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(panicle_res[which(panicle_res$padj <= 0.05 & abs(panicle_res$log2FoldChange)>=1),],"phasing/panicle/fold_changes/deltaTE.txt",quote=F,sep="\t",col.names = T,row.names = T)
panicle_aste_genes <- rownames(panicle_res)[which(panicle_res$padj <= 0.05 & abs(panicle_res$log2FoldChange)>=1)]

## DESeq2 for Ribo-seq
ind = which(panicle_coldata$SeqType == "Ribo")
panicle_coldata_ribo = panicle_coldata[ind,]
panicle_contrast_ribo_count_matrix <- panicle_contrast_count_matrix %>% select(contains('Ribo'))
panicle_ddsMat_ribo <- DESeqDataSetFromMatrix(countData = panicle_contrast_ribo_count_matrix, colData = panicle_coldata_ribo, design = ~ Condition)
panicle_ddsMat_ribo <- DESeq(panicle_ddsMat_ribo)
panicle_res_ribo <- results(panicle_ddsMat_ribo, contrast = c('Condition', 'ZS97', 'MH63'))
panicle_res_ribo <- lfcShrink(dds = panicle_ddsMat_ribo, coef = resultsNames(panicle_ddsMat_ribo)[2], res = panicle_res_ribo)
panicle_res_ribo$padj[is.na(panicle_res_ribo$padj)] <- 1
panicle_res_ribo$log2FoldChange[is.na(panicle_res_ribo$log2FoldChange)] <- 0
length(which(panicle_res_ribo$padj <= 0.05 & abs(panicle_res_ribo$log2FoldChange)>=1))
write.table(panicle_res_ribo[which(panicle_res_ribo$padj <= 0.05 & abs(panicle_res_ribo$log2FoldChange)>=1),],"phasing/panicle/fold_changes/deltaRibo.txt",quote=F,sep="\t",col.names = T,row.names = T)
panicle_asribo_genes <- rownames(panicle_res_ribo)[which(panicle_res_ribo$padj <= 0.05 & abs(panicle_res_ribo$log2FoldChange)>=1)]

## DESeq2 for RNA-seq
ind = which(panicle_coldata$SeqType == "RNA")
panicle_coldata_rna = panicle_coldata[ind,]
panicle_contrast_rna_count_matrix <- panicle_contrast_count_matrix %>% select(contains('RNA'))
panicle_ddsMat_rna <- DESeqDataSetFromMatrix(countData = panicle_contrast_rna_count_matrix, colData = panicle_coldata_rna, design = ~ Condition)
panicle_ddsMat_rna <- DESeq(panicle_ddsMat_rna)
panicle_res_rna <- results(panicle_ddsMat_rna, contrast = c('Condition', 'ZS97', 'MH63'))
panicle_res_rna <- lfcShrink(dds = panicle_ddsMat_rna, coef = resultsNames(panicle_ddsMat_rna)[2], res = panicle_res_rna)
panicle_res_rna$padj[is.na(panicle_res_rna$padj)] <- 1
panicle_res_rna$log2FoldChange[is.na(panicle_res_rna$log2FoldChange)] <- 0
length(which(panicle_res_rna$padj <= 0.05 & abs(panicle_res_rna$log2FoldChange)>=1))
write.table(rownames(panicle_res_rna)[which(panicle_res_rna$padj <= 0.05 & abs(panicle_res_rna$log2FoldChange)>=1)],"phasing/panicle/gene_lists/DTG.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(panicle_res_rna[which(panicle_res_rna$padj<=0.05 & abs(panicle_res_rna$log2FoldChange)>=1),],"parent/panicle/fold_changes/deltaRNA.txt",quote=F,sep="\t",col.names = T,row.names = T)
panicle_ase_genes <- rownames(panicle_res_rna)[which(panicle_res_rna$padj <= 0.05 & abs(panicle_res_rna$log2FoldChange)>=1)]
#

## panicle RNA TE FCplot ---------------------------------------------------
panicle_aste_res <- panicle_res[which(panicle_res$padj <= 0.05 & abs(panicle_res$log2FoldChange)>=1),]
panicle_aste_genes <- rownames(panicle_aste_res)
panicle_asribo_res <- panicle_res_ribo[which(panicle_res_ribo$padj <= 0.05 & abs(panicle_res_ribo$log2FoldChange)>=1),]
panicle_asribo_genes <- rownames(panicle_asribo_res)
panicle_ase_res <- panicle_res_rna[which(panicle_res_rna$padj <= 0.05 & abs(panicle_res_rna$log2FoldChange)>=1),]
panicle_ase_genes <- rownames(panicle_ase_res)

panicle_mRNAOnly <- rownames(panicle_res)[which((panicle_res$padj > 0.05 | abs(panicle_res$log2FoldChange)<1) & (panicle_res_rna$padj <= 0.05 & abs(panicle_res_rna$log2FoldChange)>=1))]
panicle_TEOnly <- rownames(panicle_res)[which((panicle_res$padj <= 0.05 & abs(panicle_res$log2FoldChange)>=1) & (panicle_res_rna$padj > 0.05 | abs(panicle_res_rna$log2FoldChange)<1))]
panicle_FCboth <- which((panicle_res$padj <= 0.05 & abs(panicle_res$log2FoldChange)>=1) & (panicle_res_rna$padj <= 0.05 & abs(panicle_res_rna$log2FoldChange)>=1))
panicle_compensatory <- rownames(panicle_res)[panicle_FCboth[which(panicle_res[panicle_FCboth,2]*panicle_res_rna[panicle_FCboth,2] <0 )]]
panicle_reinforcing <- rownames(panicle_res)[panicle_FCboth[which(panicle_res[panicle_FCboth,2]*panicle_res_rna[panicle_FCboth,2] >=0 )]]
panicle_FCothers <- rownames(panicle_res)[! rownames(panicle_res) %in% c(panicle_mRNAOnly,panicle_TEOnly,panicle_reinforcing, panicle_compensatory)]

panicle_fc <- data.frame(fold_te=panicle_res[,2], fold_rna=panicle_res_rna[,2]) %>% na.omit()
nrow(panicle_fc)
panicle_fc$class[rownames(panicle_fc) %in% panicle_mRNAOnly] <- 'mRNA_only'
panicle_fc$class[rownames(panicle_fc) %in% panicle_TEOnly] <- 'TE_only'
panicle_fc$class[rownames(panicle_fc) %in% panicle_compensatory] <- 'Compensatory'
panicle_fc$class[rownames(panicle_fc) %in% panicle_reinforcing] <- 'Reinforcing'
panicle_fc$class[rownames(panicle_fc) %in% panicle_FCothers] <- 'Others'
panicle_fc$class <- factor(panicle_fc$class, levels = c('mRNA_only','TE_only','Compensatory','Reinforcing','Others'))
panicle_fc$class %>% table()
p1 <- panicle_fc %>% ggplot(aes(x=fold_rna, y=fold_te, color=class)) + geom_point(size=1.5) + theme_bw(base_size = 15) + scale_color_manual(values = c('#4169E1', '#CC524F', '#DAA520', '#3CB371', 'grey')) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text = element_text(color='black'), plot.margin = unit(rep(2,4),'lines')) + labs(x=expression(paste(Log[2],'(',mRNA[ZS97],'/',mRNA[MH63],')')), y=expression(paste(Log[2],'(',TE[ZS97],'/',TE[MH63],')'))) + 
  geom_vline(xintercept = 0, size=0.8, linetype='dashed') + 
  geom_vline(xintercept = -1, size=0.8, linetype='dashed', color='#696969') +
  geom_vline(xintercept = 1, size=0.8, linetype='dashed', color='#696969') + 
  geom_hline(yintercept = 0, size=0.8, linetype='dashed') +
  geom_hline(yintercept = -1, size=0.8, linetype='dashed', color='#696969') + 
  geom_hline(yintercept = 1, size=0.8, linetype='dashed', color='#696969')
p2 <- panicle_fc %>% select(class) %>% table() %>% as.data.frame() %>% magrittr::set_colnames(c('class','Freq')) %>% filter(class!='Others') %>% ggplot(aes(x=class,y=Freq)) + geom_bar(aes(fill=class), stat = 'identity', width = 0.4) + theme_bw(base_size = 15) + coord_flip() + scale_x_discrete(limits=rev(c('mRNA_only','TE_only','Compensatory','Reinforcing'))) + scale_fill_manual(values = c('#4169E1', '#CC524F', '#DAA520', '#3CB371')) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text = element_text(color='black'), plot.margin = unit(rep(2,4),'lines')) + labs(x='', y='Number of gene') + scale_y_continuous(expand = c(0,0))
ggsave(ggarrange(p1, p2, ncol = 2, nrow = 1), filename = 'phasing/panicle/panicle_FCplot.pdf', width = 14, height = 7)

write.table(panicle_mRNAOnly,"phasing/panicle/gene_lists/panicle_mRNAOnly.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(panicle_TEOnly,"phasing/panicle/gene_lists/panicle_TEOnly.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(panicle_compensatory,"phasing/panicle/gene_lists/panicle_compensatory.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(panicle_reinforcing,"phasing/panicle/gene_lists/panicle_reinforcing.txt",quote=F,sep="\t",col.names = F,row.names = F)


## panicle total mRNA expression level -------------------------------------
panicle_expression <- rna_expression_matrix %>% select(contains('Osy63_panicle')) %>% rownames_to_column('ID') %>% filter(ID %in% rownames(panicle_contrast_count_matrix)) %>% column_to_rownames('ID')
names(panicle_expression) <- 'mRNA'
panicle_expression$group[rownames(panicle_expression) %in% panicle_mRNAOnly] <- 'mRNAOnly'
panicle_expression$group[rownames(panicle_expression) %in% panicle_TEOnly] <- 'TEOnly'
panicle_expression$group[! rownames(panicle_expression) %in% c(panicle_mRNAOnly,panicle_TEOnly)] <- 'Others'
panicle_expression$group <- factor(panicle_expression$group, levels = c('Others','mRNAOnly','TEOnly'))
panicle_expression %>% mutate(m=log10(mRNA+1)) %>% ggboxplot(x='group',y='m', width = 0.5,color='grey', fill = 'group', palette = c('#808080',rgb(208,128,91,maxColorValue = 255), rgb(94,160,205,maxColorValue = 255))) + stat_compare_means(comparisons = list(c('mRNAOnly','Others'),c('mRNAOnly','TEOnly'),c('TEOnly','Others'))) + labs(x='',y=expression(paste(Log[10],'(Total mRNA TPM+1)'))) + theme_bw(base_size = 15) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text = element_text(color = 'black')) + stat_compare_means(label.y = 5.5) + ggsave(filename = 'phasing/panicle/panicle_divExpression.pdf', width = 7, height = 7)


# root delta --------------------------------------------------------------
root_contrast_count_matrix <- count_matrix %>% select(contains('root'))
root_contrast_count_matrix <- root_contrast_count_matrix[which(rowSums(root_contrast_count_matrix)>=1),]
root_coldata <- data.frame(SampleID=names(root_contrast_count_matrix))
root_coldata$Condition <- str_split(root_coldata$SampleID, pattern = '_') %>% map_chr(~.x[length(.x)])
root_coldata$Condition <- factor(root_coldata$Condition, levels = c('MH63', 'ZS97'))
root_coldata$SeqType <- str_split(root_coldata$SampleID, pattern = '_') %>% map_chr(~.x[3])
root_coldata$SeqType <- factor(root_coldata$SeqType, levels = c('RNA','Ribo'))
root_coldata$Batch <- rep(1,8)
root_ddsMat <- DESeqDataSetFromMatrix(countData = root_contrast_count_matrix, colData = root_coldata, design =~ Condition + SeqType + Condition:SeqType)
root_ddsMat$SeqType <- relevel(root_ddsMat$SeqType,"RNA")
root_ddsMat <- DESeq(root_ddsMat)
# resultsNames(root_ddsMat)
# [1] "Intercept"                 "Condition_ZS97_vs_MH63"    "SeqType_Ribo_vs_RNA"      
# [4] "ConditionZS97.SeqTypeRibo"
root_res <- results(object = root_ddsMat, name = resultsNames(root_ddsMat)[4])
root_res$padj[is.na(root_res$padj)]  <- 1
root_res$log2FoldChange[is.na(root_res$log2FoldChange)] <- 0
length(which(root_res$padj <= 0.05 & abs(root_res$log2FoldChange)>=1))
write.table(rownames(root_res)[which(root_res$padj <= 0.05 & abs(root_res$log2FoldChange)>=1)],"phasing/root/gene_lists/DTEGs.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(root_res[which(root_res$padj <= 0.05 & abs(root_res$log2FoldChange)>=1),],"phasing/root/fold_changes/deltaTE.txt",quote=F,sep="\t",col.names = T,row.names = T)
root_aste_genes <- rownames(root_res)[which(root_res$padj <= 0.05 & abs(root_res$log2FoldChange)>=1)]

## DESeq2 for Ribo-seq
ind = which(root_coldata$SeqType == "Ribo")
root_coldata_ribo = root_coldata[ind,]
root_contrast_ribo_count_matrix <- root_contrast_count_matrix %>% select(contains('Ribo'))
root_ddsMat_ribo <- DESeqDataSetFromMatrix(countData = root_contrast_ribo_count_matrix, colData = root_coldata_ribo, design = ~ Condition)
root_ddsMat_ribo <- DESeq(root_ddsMat_ribo)
root_res_ribo <- results(root_ddsMat_ribo, contrast = c('Condition', 'ZS97', 'MH63'))
root_res_ribo <- lfcShrink(dds = root_ddsMat_ribo, coef = resultsNames(root_ddsMat_ribo)[2], res = root_res_ribo)
root_res_ribo$padj[is.na(root_res_ribo$padj)]  <- 1
root_res_ribo$log2FoldChange[is.na(root_res_ribo$log2FoldChange)] <- 0
length(which(root_res_ribo$padj <= 0.05 & abs(root_res_ribo$log2FoldChange)>=1))
write.table(root_res_ribo[which(root_res_ribo$padj <= 0.05 & abs(root_res_ribo$log2FoldChange)>=1),],"phasing/root/fold_changes/deltaRibo.txt",quote=F,sep="\t",col.names = T,row.names = T)
root_asribo_genes <- rownames(root_res_ribo)[which(root_res_ribo$padj <= 0.05 & abs(root_res_ribo$log2FoldChange)>=1)]

## DESeq2 for RNA-seq
ind = which(root_coldata$SeqType == "RNA")
root_coldata_rna = root_coldata[ind,]
root_contrast_rna_count_matrix <- root_contrast_count_matrix %>% select(contains('RNA'))
root_ddsMat_rna <- DESeqDataSetFromMatrix(countData = root_contrast_rna_count_matrix, colData = root_coldata_rna, design = ~ Condition)
root_ddsMat_rna <- DESeq(root_ddsMat_rna)
root_res_rna <- results(root_ddsMat_rna, contrast = c('Condition', 'ZS97', 'MH63'))
root_res_rna <- lfcShrink(dds = root_ddsMat_rna, coef = resultsNames(root_ddsMat_rna)[2], res = root_res_rna)
root_res_rna$padj[is.na(root_res_rna$padj)] <- 1
root_res_rna$log2FoldChange[is.na(root_res_rna$log2FoldChange)] <- 0
length(which(root_res_rna$padj <= 0.05 & abs(root_res_rna$log2FoldChange)>=1))
write.table(rownames(root_res_rna)[which(root_res_rna$padj <= 0.05 & abs(root_res_rna$log2FoldChange)>=1)],"phasing/root/gene_lists/DTG.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(root_res_rna[which(root_res_rna$padj<=0.05 & abs(root_res_rna$log2FoldChange)>=1),],"phasing/root/fold_changes/deltaRNA.txt",quote=F,sep="\t",col.names = T,row.names = T)
root_ase_genes <- rownames(root_res_rna)[which(root_res_rna$padj <= 0.05 & abs(root_res_rna$log2FoldChange)>=1)]


## root RNA TE FCplot ------------------------------------------------------
root_aste_res <- root_res[which(root_res$padj <= 0.05 & abs(root_res$log2FoldChange)>=1),]
root_aste_genes <- rownames(root_aste_res)
root_asribo_res <- root_res_ribo[which(root_res_ribo$padj <= 0.05 & abs(root_res_ribo$log2FoldChange)>=1),]
root_asribo_genes <- rownames(root_asribo_res)
root_ase_res <- root_res_rna[which(root_res_rna$padj <= 0.05 & abs(root_res_rna$log2FoldChange)>=1),]
root_ase_genes <- rownames(root_ase_res)

root_mRNAOnly <- rownames(root_res)[which((root_res$padj > 0.05 | abs(root_res$log2FoldChange)<1) & (root_res_rna$padj <= 0.05 & abs(root_res_rna$log2FoldChange)>=1))]
root_TEOnly <- rownames(root_res)[which((root_res$padj <= 0.05 & abs(root_res$log2FoldChange)>=1) & (root_res_rna$padj > 0.05 | abs(root_res_rna$log2FoldChange)<1))]
root_FCboth <- which((root_res$padj <= 0.05 & abs(root_res$log2FoldChange)>=1) & (root_res_rna$padj <= 0.05 & abs(root_res_rna$log2FoldChange)>=1))
root_compensatory <- rownames(root_res)[root_FCboth[which(root_res[root_FCboth,2]*root_res_rna[root_FCboth,2] <0 )]]
root_reinforcing <- rownames(root_res)[root_FCboth[which(root_res[root_FCboth,2]*root_res_rna[root_FCboth,2] >=0 )]]
root_FCothers <- rownames(root_res)[! rownames(root_res) %in% c(root_mRNAOnly,root_TEOnly,root_reinforcing, root_compensatory)]

root_fc <- data.frame(fold_te=root_res[,2], fold_rna=root_res_rna[,2]) %>% na.omit()
nrow(root_fc)
root_fc$class[rownames(root_fc) %in% root_mRNAOnly] <- 'mRNA_only'
root_fc$class[rownames(root_fc) %in% root_TEOnly] <- 'TE_only'
root_fc$class[rownames(root_fc) %in% root_compensatory] <- 'Compensatory'
root_fc$class[rownames(root_fc) %in% root_reinforcing] <- 'Reinforcing'
root_fc$class[rownames(root_fc) %in% root_FCothers] <- 'Others'
root_fc$class <- factor(root_fc$class, levels = c('mRNA_only','TE_only','Compensatory','Reinforcing','Others'))
root_fc$class %>% table()
p1 <- root_fc %>% ggplot(aes(x=fold_rna, y=fold_te, color=class)) + geom_point(size=1.5) + theme_bw(base_size = 15) + scale_color_manual(values = c('#4169E1', '#CC524F', '#DAA520', '#3CB371', 'grey')) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text = element_text(color='black'), plot.margin = unit(rep(2,4),'lines')) + labs(x=expression(paste(Log[2],'(',mRNA[ZS97],'/',mRNA[MH63],')')), y=expression(paste(Log[2],'(',TE[ZS97],'/',TE[MH63],')'))) + 
  geom_vline(xintercept = 0, size=0.8, linetype='dashed') + 
  geom_vline(xintercept = -1, size=0.8, linetype='dashed', color='#696969') +
  geom_vline(xintercept = 1, size=0.8, linetype='dashed', color='#696969') + 
  geom_hline(yintercept = 0, size=0.8, linetype='dashed') +
  geom_hline(yintercept = -1, size=0.8, linetype='dashed', color='#696969') + 
  geom_hline(yintercept = 1, size=0.8, linetype='dashed', color='#696969')
p2 <- root_fc %>% select(class) %>% table() %>% as.data.frame() %>% magrittr::set_colnames(c('class','Freq')) %>% filter(class!='Others') %>% ggplot(aes(x=class,y=Freq)) + geom_bar(aes(fill=class), stat = 'identity', width = 0.4) + theme_bw(base_size = 15) + coord_flip() + scale_x_discrete(limits=rev(c('mRNA_only','TE_only','Compensatory','Reinforcing'))) + scale_fill_manual(values = c('#4169E1', '#CC524F', '#DAA520', '#3CB371')) + theme(panel.grid = element_blank(), legend.position = 'none', axis.text = element_text(color='black'), plot.margin = unit(rep(2,4),'lines')) + labs(x='', y='Number of gene') + scale_y_continuous(expand = c(0,0))
ggsave(ggarrange(p1, p2, ncol = 2, nrow = 1), filename = 'phasing/root/root_FCplot.pdf', width = 14, height = 7)

write.table(root_mRNAOnly,"phasing/root/gene_lists/root_mRNAOnly.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(root_TEOnly,"phasing/root/gene_lists/root_TEOnly.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(root_compensatory,"phasing/root/gene_lists/root_compensatory.txt",quote=F,sep="\t",col.names = F,row.names = F)
write.table(root_reinforcing,"phasing/root/gene_lists/root_reinforcing.txt",quote=F,sep="\t",col.names = F,row.names = F)

root_mRNAOnly_GO <- read.table('phasing/root/gene_lists/root_mRNAOnly.txt.GOseq.enriched.plot', header = T, sep = "\t")
root_mRNAOnly_GO$group <- 'mRNA_Only'
root_TEOnly_GO <- read.table('phasing/root/gene_lists/root_TEOnly.txt.GOseq.enriched.plot', header = T, sep = "\t")
root_TEOnly_GO$group <- 'TE_Only'
bind_rows(root_mRNAOnly_GO, root_TEOnly_GO) %>% mutate(group=factor(group, levels = c('TE_Only', 'mRNA_Only'))) %>% ggplot(aes(x=numDEInCat, y=term)) + geom_point(aes(size=-log10(over_represented_pvalue)),color=rgb(52,145,191, maxColorValue = 255)) + facet_grid(group~ontology, scales = 'free_y', space = 'free') + theme_bw(base_size = 12) + theme(legend.position = 'top', axis.text = element_text(color = 'black')) + labs(y='', x='Number of genes', size=expression(paste(-Log[10],'(p-value)'))) + ggsave(filename = 'phasing/root/gene_lists/GO.pdf', width = 14, height = 7)

# root total mRNA expression level ----------------------------------------
root_expression <- rna_expression_matrix %>% select(contains('Osy63_root')) %>% rownames_to_column('ID') %>% filter(ID %in% rownames(root_contrast_count_matrix)) %>% column_to_rownames('ID')
names(root_expression) <- 'mRNA'
root_expression$group[rownames(root_expression) %in% root_mRNAOnly] <- 'mRNAOnly'
root_expression$group[rownames(root_expression) %in% root_TEOnly] <- 'TEOnly'
root_expression$group[! rownames(root_expression) %in% c(root_mRNAOnly,root_TEOnly)] <- 'Others'
root_expression$group <- factor(root_expression$group, levels = c('Others','mRNAOnly','TEOnly'))
root_expression %>% mutate(m=log10(mRNA+1)) %>% ggboxplot(x='group',y='m', width = 0.5,color='grey', fill = 'group', palette = c('#808080',rgb(208,128,91,maxColorValue = 255), rgb(94,160,205,maxColorValue = 255))) + stat_compare_means(comparisons = list(c('mRNAOnly','Others'),c('mRNAOnly','TEOnly'),c('TEOnly','Others'))) + labs(x='',y=expression(paste(Log[10],'(Total mRNA TPM+1)'))) + theme_bw(base_size = 15) + theme(legend.position = 'none', panel.grid = element_blank(), axis.text = element_text(color = 'black')) + stat_compare_means(label.y = 6) + ggsave(filename = 'phasing/root/root_divExpression.pdf', width = 7, height = 7)


# results analysis --------------------------------------------------------

TEOnly_union <- Reduce(union, list(leaf_TEOnly,panicle_TEOnly,root_TEOnly))
Reduce(union, list(leaf_TEOnly,panicle_TEOnly,root_TEOnly)) %>% sort %>% write.table(file = 'phasing/TE_only_genes.list', quote = F, sep = "\t", row.names = F, col.names = F)