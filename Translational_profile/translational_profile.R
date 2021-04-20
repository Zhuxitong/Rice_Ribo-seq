library(tidyverse)
library(reshape2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(Cairo)


# load basic data ---------------------------------------------------------

ribo_expression_matrix <- read.table('refOmh63_ribo.gene.TPM.not_cross_norm.avg_reps.byLog.matrix', header = T, row.names = 1)
rna_expression_matrix <- read.table('refOmh63_rna.gene.TPM.not_cross_norm.avg_reps.byLog.matrix', header = T, row.names = 1)

ordered_samples <- c("Omh63_leaf", "Omh63_panicle", "Omh63_root", "Osy63_leaf", "Osy63_panicle", "Osy63_root", "Ozs97_leaf", "Ozs97_panicle", "Ozs97_root")
ordered_tissues <- c('leaf', 'panicle', 'root')
ordered_specie <- c('Omh63','Osy63','Ozs97')
all_genes <- rownames(ribo_expression_matrix)

merged_expression_matrix <- bind_cols(ribo_expression_matrix, rna_expression_matrix)

# expression var within tissue --------------------------------------------

variance_table <- list()
for(i in ordered_samples){
  specie <- str_split(i,pattern = '_')[[1]][1]
  tissue <- str_split(i,pattern = '_')[[1]][2]
  
  a <- merged_expression_matrix %>% select(contains(i) & contains('Ribo')) %>% filter(.[[1]]>=0.1 & .[[1]]<=10000) %>% rownames_to_column('gene_id')
  names(a)[2] <- 'value'
  a$type <- 'Ribo'
  
  b <- merged_expression_matrix %>% select(contains(i) & contains('RNA')) %>% filter(.[[1]]>=0.1 & .[[1]]<=10000) %>% rownames_to_column('gene_id')
  names(b)[2] <- 'value'
  b$type <- 'RNA'
  
  my_variance <- bind_rows(a, b) %>% group_by(type) %>% summarise(var=var(log10(value)))
  my_variance$specie <- specie
  my_variance$tissue <- tissue
  variance_table[[i]] <- my_variance
  
  ribo_var <- my_variance[1,'var'] %>% as.numeric()
  rna_var <- my_variance[2,'var'] %>% as.numeric()
  
  p <- bind_rows(a, b) %>% ggplot(aes(x=log10(value), fill=type)) + geom_density(alpha=0.7, color=NA) + theme_classic(base_size = 20) + theme(legend.title = element_blank(), legend.position = c(0.8,0.9), plot.title = element_text(hjust = 0.5), plot.margin = unit(rep(1.5,4),'lines')) + labs(x=expression(paste('Expression value(',log[10],'(TPM))' )), y='Density', title = i) + scale_x_continuous(breaks = -1:5, labels = -1:5, limits = c(-1,5)) + scale_y_continuous(breaks = seq(0,0.6,0.1), labels = seq(0,0.6,0.1), limits = c(0,0.6)) + scale_fill_manual(values = c('#CC524F','#04A7BE'), labels = c(paste0('Ribo var=',round(ribo_var, 4)), paste0('RNA var=', round(rna_var,4))))
  output <- paste0(i,'.expression.pdf')
  ggsave(filename = output, plot = p, width = 8, height = 8)
}

variance_table <- bind_rows(variance_table)
variance_table$specie <- factor(variance_table$specie, levels = ordered_specie)
variance_table %>% ggplot(aes(x=tissue, y=var, color=type, shape=tissue)) + geom_point(size=4) + theme_classic(base_size = 20) + facet_grid(. ~ specie, switch = 'x') + theme(strip.placement = 'outside', strip.background = element_blank(), panel.spacing = unit(2,'lines'), legend.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(rep(1.5,4), 'lines')) + labs(x="", y='Expression variance') + scale_color_manual(values = c('#CC524F','#04A7BE')) + scale_y_continuous(breaks = seq(0.5,0.9,0.1), labels = seq(0.5,0.9,0.1), limits = c(0.5, 0.9))


variance_table %>% reshape2::dcast(specie + tissue ~ type, value.var = 'var') %>% mutate(diff=Ribo-RNA, diff_pct=diff/RNA*100) %>% arrange(diff_pct)
variance_table %>% reshape2::dcast(specie + tissue ~ type, value.var = 'var') %>% mutate(diff=Ribo-RNA, diff_pct=diff/RNA*100) %>% arrange(tissue, diff_pct)
variance_table %>% reshape2::dcast(specie + tissue ~ type, value.var = 'var') %>% mutate(diff=Ribo-RNA, diff_pct=diff/RNA*100) %>% group_by(tissue) %>% summarise(ave=mean(diff_pct))



# expressed genes in each tissue ------------------------------------------

shared_list <- list()
for( i in ordered_samples){
  shared <- merged_expression_matrix %>% select(contains(i)) %>% filter(.[[1]] >= 0.1 & .[[2]] >= 0.1) %>% count()
  shared_data <- data.frame(sample=i, num=shared)
  sup <- str_split_fixed(i, pattern = '_', 2) %>% as.data.frame()
  names(sup) <- c('Specie','tissue')
  shared_data <- bind_cols(shared_data, sup)
  shared_list[[i]] <- shared_data
}

shared_list <- bind_rows(shared_list)
shared_list$total <- 39406
shared_list$pct <- shared_list$n/shared_list$total * 100

expressed_genes <- merged_expression_matrix %>% map_df(~{.x=.x[which(.x >= 0.1 & .x <= 10000)]; length(.x)}) %>% t() %>% as.data.frame() %>% rownames_to_column('Sample')
names(expressed_genes)[2] <- c('num')
addition_info <- str_split_fixed(expressed_genes$Sample, pattern = '_', n = 3) %>% as.data.frame()
names(addition_info) <- c('Specie', 'tissue', 'type')
expressed_genes <- bind_cols(expressed_genes, addition_info)
expressed_genes$total <- 39406
expressed_genes$Specie <- factor(expressed_genes$Specie, levels = ordered_specie)
expressed_genes %>% mutate(pct=num/total*100) %>% ggplot(aes(x=pct, y=Specie)) + geom_point(aes(color=type), size=4) + facet_grid(~tissue) + geom_point(data = shared_list, aes(x=pct, y=Specie), size=4,shape=9) + theme_bw(base_size = 20) + theme(legend.position = 'top', legend.title = element_blank(), strip.placement = 'outside', strip.background = element_blank(), panel.grid = element_blank()) + labs(x='% expressed genes', y="") + scale_color_manual(values = c('#CC524F','#04A7BE')) + scale_y_discrete(limits=rev(ordered_specie))


expressed_genes %>% mutate(pct=num/total*100) %>%  arrange(type, tissue)
expressed_genes %>% group_by(type, tissue) %>% summarise(n=mean(num)) %>% mutate(pct=n/39406*100)
shared_list %>% arrange(tissue, Specie)

a <- ribo_expression_matrix %>% melt()
a$TPM <- cut(a$value,breaks=c(0,0.1,1,10,100,1000,10000,Inf), right = F, dig.lab = 5)
lab_x <- as.factor(unique(a$variable) %>% str_split(pattern = '_') %>% map_chr(.f = function(x) paste(x[1], x[2], sep = "_")))
p1 <- a %>% group_by(variable, TPM) %>% count() %>% ggplot(aes(x=variable,y=n, fill=TPM)) + geom_bar(stat = 'identity', width = 0.5) + coord_flip() + scale_fill_viridis_d(alpha = 0.8) + theme_classic() + xlab("") + ylab('Number of genes') + scale_y_continuous(expand = c(0,0), limits = c(0, 45000)) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0, size = rel(1.2)), axis.text = element_text(color = 'black')) + labs(title = 'Ribo-seq') + scale_x_discrete(labels=rev(lab_x))

b <- rna_expression_matrix %>% melt()
b$TPM <- cut(b$value,breaks=c(0,0.1,1,10,100,1000,10000,Inf), right = F, dig.lab = 5)
p2 <- b %>% group_by(variable, TPM) %>% count() %>% ggplot(aes(x=variable,y=n, fill=TPM)) + geom_bar(stat = 'identity', width = 0.5) + coord_flip() + scale_fill_viridis_d(alpha = 0.8) + theme_classic() + xlab("") + ylab('Number of genes') + scale_y_continuous(expand = c(0,0), limits = c(0, 45000)) + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y =element_blank(), axis.text = element_text(color = 'black')) + labs(title = 'RNA-seq')

p <- ggarrange(p1,p2,common.legend = T, legend = 'right', widths = c(1.3,1))
ggsave(filename = 'expression_profile.pdf', plot = p, width = 8, height = 8, family='Times')

a_file <- a %>% group_by(variable, TPM) %>% count() 
b_file <- b %>% group_by(variable, TPM) %>% count()

write.table(x = a_file, file = 'ribo_expression_profile.txt', quote = F, row.names = F, col.names = F, sep = "\t")
write.table(x = b_file, file = 'rna_expression_profile.txt', quote = F, row.names = F, col.names = F, sep = "\t")



# calculate TEs -----------------------------------------------------------

TE_list <- list()
TE_ranges <- list()
for(i in ordered_samples){
  TE_dt <- merged_expression_matrix %>% select(contains(i)) %>% filter(.[[1]] >= 0.1 & .[[2]] >= 0.1) %>% rownames_to_column('ID') %>% mutate(TE=.[[2]]/.[[3]])
  TE_limits <- quantile(TE_dt$TE, c(0.025,0.975))
  #TE_limits <- quantile(TE_dt$TE, c(0.01,0.99))
  TE_dt <- TE_dt %>% filter(TE >= TE_limits[1] & TE <= TE_limits[2])
  names(TE_dt)[2:3] <- c('Ribo', 'RNA')
  TE_dt$sample <- i
  sup <- str_split_fixed(i, pattern = '_', 2) %>% as.data.frame()
  names(sup) <- c('Specie','tissue')
  TE_dt <- bind_cols(TE_dt, sup)
  TE_list[[i]] <- TE_dt
  
  TE_range <- data.frame(sample=i, TE_min=as.numeric(TE_limits[1]), TE_max=as.numeric(TE_limits[2]), TE_range=TE_limits[2]/TE_limits[1])
  TE_range <- bind_cols(TE_range,sup)
  TE_ranges[[i]] <- TE_range
}

TE_list <- bind_rows(TE_list)
p <- TE_list %>% ggplot(aes(x=log2(TE))) + stat_ecdf(aes(color=Specie, linetype=tissue),size=1.5) + scale_linetype_manual(values = c(2,3,6)) + theme_classic(base_size = 20) + labs(x=expression(paste(Log[2],'(TE)')), y='Cumulative frequency') + theme(legend.title = element_blank(), legend.position = c(0.1,0.8), plot.margin = unit(rep(1.5,4), 'lines'), legend.background = element_blank()) + scale_color_manual(values = c("#3498BF", '#F48959', "#52A563")) + scale_x_continuous(breaks = seq(-4,4,1), labels = seq(-4,4,1))
ggsave(filename = 'TE_cumulative_freq.pdf', plot = p, width = 8, height = 8, family='Times')
write.table(x = TE_list, file = 'TE_dis.txt', quote = F, sep = "\t", row.names = F, col.names = F)

TE_list %>% ggplot(aes(x=sample, y=log2(TE))) + geom_boxplot(aes(fill=sample))


# ks.test for TE ecdf
combn(ordered_samples,2, simplify = FALSE) %>% map(.f = function(x){
  p <- ks.test(jitter(TE_list[TE_list$sample==x[1], 'TE']), ecdf(TE_list[TE_list$sample==x[2], 'TE']))
  c(sample1=x[1],sample2=x[2],p_value=p$p.value)
}) %>% bind_rows() %>% as.data.frame()


# TE ranges
TE_ranges <- bind_rows(TE_ranges)
p <- TE_ranges %>% ggplot(aes(x=tissue, y=TE_range)) + geom_bar(aes(fill=tissue), color='white', stat = 'identity', width = 0.7) + facet_grid(~Specie, switch = 'x') + scale_fill_manual(values = c('#70BF73','#FFCE38','#D58F4E')) + theme_classic(base_size = 20) + theme(panel.grid = element_blank(), strip.background = element_blank(), strip.placement = 'outside', axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.background = element_blank(), legend.position = 'top', legend.title = element_blank(), plot.margin = unit(rep(1.5,4), 'lines'), panel.spacing = unit(2,'lines')) + scale_y_continuous(expand = c(0,0), limits = c(0,200)) + geom_text(aes(y=TE_range, label=round(TE_range,0)), size=5, vjust=-0.5) + labs(x="", y='TE range')
ggsave(filename = 'TE_range.pdf', plot = p, width = 8, height = 8, family='Times')

TE_ranges %>% mutate(TE_min_log=log2(TE_min), TE_max_log=log2(TE_max), TE_range_log=TE_max_log-TE_min_log) %>% ggplot(aes(x=tissue, y=TE_range_log)) + geom_bar(aes(fill=tissue), color='white', stat = 'identity', width = 0.7) + facet_grid(~Specie, switch = 'x') + scale_fill_manual(values = c('#70BF73','#FFCE38','#D58F4E')) + theme_classic(base_size = 20) + theme(panel.grid = element_blank(), strip.background = element_blank(), strip.placement = 'outside', axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.background = element_blank(), legend.position = 'top', legend.title = element_blank(), plot.margin = unit(rep(1.5,4), 'lines'), panel.spacing = unit(2,'lines')) + scale_y_continuous(expand = c(0,0), limits = c(0,8)) + geom_text(aes(y=TE_range_log, label=round(TE_range_log,2)), size=5, vjust=-0.5) + labs(x="", y='TE range')

# TE and RNA aboundance
p <- TE_list %>% ggplot(aes(x=log10(RNA), y=log2(TE))) + geom_hex() + facet_wrap(~sample) + scale_fill_viridis_c() + theme_bw(base_size = 20) + theme(panel.grid = element_blank(), plot.margin = unit(rep(1.5,4), 'lines')) + labs(x=expression(paste("RNA: ",log[10],"(TPM)")), y=expression(paste(Log[2],"(TE)")), fill="Gene\ncount")
ggsave(filename = 'TE_RNAabundance.pdf', plot = p, width = 8, height = 8, family='Times')

p <- TE_list %>% group_by(Specie, tissue) %>% summarise(c=cor(log10(RNA),log2(TE))) %>% ggplot(aes(x=tissue, y=c)) + geom_point(aes(color=tissue), size=4) + facet_grid(~Specie, switch = 'x') + theme_classic(base_size = 20) + scale_y_continuous(limits = c(-0.25,0.25)) + scale_color_manual(values = c('#70BF73','#FFCE38','#D58F4E')) + theme(panel.grid = element_blank(), strip.background = element_blank(), strip.placement = 'outside', axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.background = element_blank(), legend.position = 'top', legend.title = element_blank(), plot.margin = unit(rep(1.5,4), 'lines'), panel.spacing = unit(2,'lines')) + labs(x="", y=expression(paste(rho, ' between RNA aboundance and TE'))) + geom_hline(yintercept = 0, linetype=3, size=1)
ggsave(filename = 'TE_RNAabundance_cor.pdf', plot = p, width = 8, height = 8, family='Times')


# distance between expression levels --------------------------------------

expressed_genes <- list()
for(i in ordered_tissues){
  expressed_gene <- merged_expression_matrix %>% select(contains(paste(i,'RNA', sep = "_"))) %>% filter_all(all_vars(.>=0.1)) %>% rownames_to_column('ID') %>% select('ID')
  expressed_genes[[i]] <- expressed_gene
}
expressed_genes <- bind_rows(expressed_genes)
expressed_genes <- unique(expressed_genes$ID)
length(expressed_genes)
write.table(x = expressed_genes, file = 'expressed_genes-ForGo.txt', quote = F, sep = "\t", row.names = F, col.names = F)

#
expressed_gene_matrix <- merged_expression_matrix %>% rownames_to_column('ID') %>% filter(ID %in% expressed_genes) %>% column_to_rownames('ID')
tissue_comb <- read.table('tissue_comb.txt', header = F)
dis_matrix <- list()
for(i in 1:nrow(tissue_comb)){
    a <- tissue_comb[i,] %>% as.character()
    a <- paste0(a,'_RNA')
    a_dis <- expressed_gene_matrix %>% select(a) %>% mutate(x=log10(.[[1]]+1), y=log10(.[[2]]+1), z=(x-y)^2) %>% summarise(dis=sqrt(sum(z)/nrow(expressed_gene_matrix))) %>% as.numeric()
    
    b <- tissue_comb[i,] %>% as.character()
    b <- paste0(b,'_Ribo')
    b_dis <- expressed_gene_matrix %>% select(b) %>% mutate(x=log10(.[[1]]+1), y=log10(.[[2]]+1), z=(x-y)^2) %>% summarise(dis=sqrt(sum(z)/nrow(expressed_gene_matrix))) %>% as.numeric()
    dis <- data.frame(type=c('RNA','Ribo'), tissue1=c(a[1],b[1]), tissue2=c(a[2],b[2]), dis=c(a_dis,b_dis))
    dis_matrix[[i]] <- dis
}
dis_matrix <- bind_rows(dis_matrix)
dis_matrix %>% arrange(type) %>% mutate(t1=str_split(tissue1, pattern = '_', n = 3) %>% map_chr(., ~paste(.x[1],.x[2], sep = '_')), t2=str_split(tissue2, pattern = '_', n = 3) %>% map_chr(., ~paste(.x[1],.x[2], sep = '_'))) %>% select(t1, t2, type, dis) %>% dcast(formula = t1 + t2 ~ type, value.var = 'dis') %>% mutate(ti1=str_split(t1, pattern = '_', n = 2) %>% map_chr(., ~.x[2]), ti2=str_split(t2, pattern = '_', n = 2) %>% map_chr(., ~.x[2])) %>% filter(ti1==ti2) %>% arrange(ti1) %>% write_tsv('distance_intratissue.txt')
dis_matrix %>% arrange(type) %>% mutate(t1=str_split(tissue1, pattern = '_', n = 3) %>% map_chr(., ~paste(.x[1],.x[2], sep = '_')), t2=str_split(tissue2, pattern = '_', n = 3) %>% map_chr(., ~paste(.x[1],.x[2], sep = '_'))) %>% select(t1, t2, type, dis) %>% dcast(formula = t1 + t2 ~ type, value.var = 'dis') %>% mutate(ti1=str_split(t1, pattern = '_', n = 2) %>% map_chr(., ~.x[2]), ti2=str_split(t2, pattern = '_', n = 2) %>% map_chr(., ~.x[2])) %>% filter(ti1!=ti2) %>% arrange(ti1) %>% write_tsv('distance_intertissue.txt')

#
dis_matrix_plot <- read.table('tissue_distance.txt', header = T, row.names = 1)
dis_matrix_plot[dis_matrix_plot==0] <- NA
pdf(file = 'tissue_distance.pdf', width = 8, height = 8, family = 'Times')
col_fun <- colorRamp2(breaks = c(0,0.7), colors = c('white', '#4169E1'))
ht <- dis_matrix_plot %>% as.matrix %>% 
  Heatmap(cluster_rows = F, 
          cluster_columns = F,
          row_names_side = 'left',
          column_names_side = 'top',
          column_names_rot = 45,
          
          name = 'Distance',
          
          na_col = 'white',
          top_annotation = HeatmapAnnotation(Specie=rep(x = ordered_specie, each=3), Tissue=rep(ordered_tissues,3),col = list(Specie=c(Omh63="#3498BF", Osy63='#F48959', Ozs97="#52A563"), Tissue=c(leaf='#70BF73',panicle='#FFCE38',root='#D58F4E')), show_legend = F, simple_anno_size = unit(0.7,'cm'), annotation_label = c(' ',' ')),
          left_annotation = HeatmapAnnotation(which = 'row',Specie=rep(x = ordered_specie, each=3), Tissue=rep(ordered_tissues,3),col = list(Specie=c(Omh63="#3498BF", Osy63='#F48959', Ozs97="#52A563"), Tissue=c(leaf='#70BF73',panicle='#FFCE38',root='#D58F4E')), simple_anno_size = unit(0.7,'cm'), annotation_name_side = 'top'),
          
          col = col_fun,
          rect_gp = gpar(type = "none"),
          cell_fun = function(j, i, x, y, w, h, fill){
            grid.rect(x,y,w, h, gp = gpar(col = 'grey', fill=NA))
            if(j==i | is.na(dis_matrix_plot[i,j])){
              grid.text(label = '*', x = x, y = y)
            }
            else if(i > j){
              grid.circle(x=x,y=y,r = 0.05, gp = gpar(fill=col_fun(dis_matrix_plot[i,j]), col=NA))
              grid.text(sprintf("%.2f", dis_matrix_plot[i, j]), x, y, gp = gpar(fontsize = 10))
            }
            else {
              #grid.polygon(x=c(0:3)/10,y=, id.lengths = 1, gp = gpar(col=NA, fill=col_fun(dis_matrix_plot[i,j])))
              grid.rect(x,y,unit(0.08,'npc'),unit(0.08,'npc'), gp = gpar(col=NA, fill=col_fun(dis_matrix_plot[i,j])))
              grid.text(sprintf("%.2f", dis_matrix_plot[i, j]), x, y, gp = gpar(fontsize = 10))
              #grid.polygon(x=unit.c(x - 0.5*width, x + 0.5*width, x + 0.5*width),y=unit.c(y + 0.5*height, y + 0.5*height, y - 0.5*height), gp = gpar(fill = col_fun(dis_matrix_plot[i,j]), col = NA))
              #grid.polygon(x=unit.c(x+0.5*w,x+0.5*w,y+0.5*h),y=unit.c(y,y,x), gp = gpar(fill = col_fun(dis_matrix_plot[i,j]), col = NA))
            }
            #grid.text(sprintf("%.2f", dis_matrix_plot[i, j]), x, y, gp = gpar(fontsize = 10))
          },
          
          heatmap_width = unit(6.5,'inches'),
          heatmap_height = unit(6.5, 'inches')
)

lgd <- Legend(labels = c('Ribo', 'RNA'), title = 'Type', background = 'white', type = 'points', pch = c(16,15), size = unit(4,'mm'))
draw(ht, annotation_legend_list = lgd)
dev.off()