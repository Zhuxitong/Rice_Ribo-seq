library(tidyverse)
library(reshape2)
library(riboWaltz)
library(psych)


# Test --------------------------------------------------------------------

## annotation file
annotation_dt <- create_annotation(gtfpath = 'Omh63.maker.gtf')
## 
reads_list <- bamtolist(bamfolder = ".", annotation = annotation_dt)
filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom",length_filter_vector = 24:30)
psite_offset <- psite(filtered_list, plot = TRUE, cl = 100, log_file = 'test.log', log_file_dir = './')
write.table(x = psite_offset, file = 'psite_offset_temp.txt', quote = F, row.names = F, sep = "\t")

reads_psite_list <- psite_info(data = filtered_list, offset = psite_offset, site = c("psite",'asite','esite'), fastapath = 'Omh63.mRNA.fa', fasta_genome = F, refseq_sep = ' ')

write.table(x = reads_psite_list[[1]], file = 'test.txt', quote = F, row.names = F, sep = "\t")

example_frames_stratified <- frame_psite_length(reads_psite_list, region = "all", cl = 100)
example_frames_stratified$dt %>% head()

example_frames <- frame_psite(reads_psite_list, sample = names(reads_psite_list), region = "all")
example_frames$plot
example_frames$dt %>% head()

example_metaprofile <- metaprofile_psite(reads_psite_list, annotation_dt, sample = names(reads_psite_list), utr5l = 20, cdsl = 100, utr3l = 20, plot_title = "sample.transcript", transcripts = test_trans$transcript)
example_metaprofile$dt %>% distinct(reg)
example_metaprofile$plot_Ozs97_root_Ribo_B2_R1_refOmh63_rsem.transcript.sorted


example_psite_region <- region_psite(reads_psite_list, annotation_dt, sample = names(reads_psite_list))

print('Done')

# summary -----------------------------------------------------------------

summary_psite_info <- read.table('psite_offset.txt', header = T)
summary_psite_info %>% select(c("length","corrected_offset_from_5","sample")) %>% reshape2::dcast(formula = sample ~ length, value.var = 'corrected_offset_from_5')
summary_psite_info %>% select(c("length","offset_from_5","sample")) %>% reshape2::dcast(formula = sample ~ length, value.var = 'offset_from_5')
mode_fun <- function(x){
  return(as.numeric(names(table(x))[table(x) == max(table(x))])[1])
}
summary_psite_info %>% select(c("length","corrected_offset_from_5","sample")) %>% group_by(length) %>% summarise(most=mode_fun(corrected_offset_from_5)) %>% write.table(file = 'psite_offset_most.txt', row.names = F, quote = F, sep = "\t")


# frame -------------------------------------------------------------------

ribo_frame <- read.table('ribo_frame.txt', header = F, quote = "", sep = "\t")
rna_frame <- read.table('rna_frame.txt', header = F, quote = "", sep = "\t")
frame_analysis <- bind_rows(ribo_frame, rna_frame)
names(frame_analysis) <- c('region','frame','count','percentage','sample','type')
frame_analysis$region <- factor(frame_analysis$region, levels = c("5' UTR", "CDS", "3' UTR"))
frame_analysis$frame <- factor(frame_analysis$frame, levels = c(0:2))

frame_analysis %>% group_by(type, region, frame) %>% summarise(n=sum(count)) %>% mutate(freq=n/sum(n)*100)
frame_analysis %>% group_by(type, region, frame) %>% summarise(n=sum(count)) %>% mutate(freq=n/sum(n)*100) %>% ggplot(aes(x = frame, y = freq, fill=type)) + geom_bar(stat = "identity",position = 'dodge', width = 0.7) + facet_grid(~region) + theme_bw(base_size = 20) + labs(x = "Frame", y = "P-site signal (%)") + theme(legend.title = element_blank(), legend.position = c(0.9,0.9), legend.background = element_blank(), plot.margin = unit(rep(1.5,4), 'lines'), panel.grid = element_blank(), axis.title.y = element_text(vjust = 3)) + scale_fill_manual(values = c('#F07368','#1DB0E6')) + scale_y_continuous(limits = c(0,80))


# metaprofile -------------------------------------------------------------

normalizeSide <- function(x) {
  normalizeFactor <- mean(x[1:20])
  #normalizeFactor <- mean(x)
  return(x/normalizeFactor)
}

metaprofile_analysis <- list()
for(i in dir('./', pattern = '*metaprofile.txt$')){
  metaprofile <- read.table(i, header = F, sep = "\t", quote = "")
  metaprofile_nor <- metaprofile %>% group_by(V5,V1,V2) %>% summarise(n=sum(V3)) %>% arrange(V2) %>% map_at(vars(n), ~normalizeSide(.x)) %>% map_df(~.x)
  metaprofile_analysis[[i]] <- metaprofile_nor
}

metaprofile_analysis <- bind_rows(metaprofile_analysis)
names(metaprofile_analysis) <- c('type', 'distance', 'reg', 'nor')

metaprofile_analysis$reg <- factor(metaprofile_analysis$reg, levels = c('Distance from start (nt)','Distance from stop (nt)'))

metaprofile_analysis %>% ggplot(aes(x = distance, y = nor, color = type)) + geom_line(size = 1.2) + facet_grid( .~reg, scales = 'free', switch = "x") + theme_bw(base_size = 20) + theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), strip.placement = "outside", legend.position = c(1, 1), legend.justification = c(1, 1), legend.title = element_blank(), legend.background = element_blank(), plot.margin = unit(rep(1.5,4), 'lines')) + scale_color_manual(values = c('#F07368','#1DB0E6')) + ylab("Normalized read count")


# metaprofile for sample and tissues --------------------------------------

ribo_metaprofile <- read.table('ribo_metaprofile.txt', header = F, sep = "\t", quote = "")
names(ribo_metaprofile) <- c('distance','reg','count','sample','type')
ribo_metaprofile$reg <- factor(ribo_metaprofile$reg, levels = c('Distance from start (nt)','Distance from stop (nt)'))
ribo_metaprofile$sample_chr <- ribo_metaprofile$sample %>% str_split(pattern = '_') %>% map_chr(~paste0(.x[1],.x[2],.x[4],.x[5]))

ribo_metaprofile %>% filter(str_detect(sample,'leaf')) %>% ggplot(aes(x = distance, y = count)) + geom_line(size = 1.2) + facet_grid( sample_chr ~reg, scales = 'free', switch = 'x') + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), strip.placement = "outside", strip.text.y = element_text(size = 8), plot.margin = unit(rep(1.5,4), 'lines'))
ribo_metaprofile %>% filter(str_detect(sample,'root')) %>% ggplot(aes(x = distance, y = count)) + geom_line(size = 1.2) + facet_grid( sample_chr ~reg, scales = 'free', switch = 'x') + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), strip.placement = "outside",strip.text.y = element_text(size = 8), plot.margin = unit(rep(1.5,4), 'lines'))
ribo_metaprofile %>% filter(str_detect(sample,'panicle')) %>% ggplot(aes(x = distance, y = count)) + geom_line(size = 1.2) + facet_grid( sample_chr ~reg, scales = 'free', switch = 'x') + theme_bw(base_size = 12) + theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), strip.placement = "outside", strip.text.y = element_text(size = 8),plot.margin = unit(rep(1.5,4), 'lines'))


omh63_CDS_info <- read.table('Omh63.CDS.info', header = F)
omh63_CDS_info %>% ggplot(aes(x=V3)) + geom_histogram(binwidth = 100) + scale_x_continuous(breaks = seq(0,3000,500), labels = seq(0,3000,500), limits = c(0,3000)) + theme_classic(base_size = 20)
summary(omh63_CDS_info$V3)

test_trans <- annotation_dt %>% filter(l_utr5>0 & l_cds>=300 & l_utr3>0) %>% select(transcript)
write.table(test_trans$transcript, "utr_cds300.trans.list", quote = F, sep = "\t", row.names = F, col.names = F)



# codon -------------------------------------------------------------------

first_50_codon <- read.table('utr_cds300.trans.matrix', header = F, row.names = 1)
# first_50_codon_plot <- first_50_codon %>% map_dbl(mean) %>% as.data.frame()
first_50_codon_plot <- first_50_codon %>% map_dbl(geometric.mean) %>% as.data.frame()
names(first_50_codon_plot) <- 'usage'
first_50_codon_plot$x <- 1:50
first_50_codon_plot$type <- 'Distance from start(codon)'
ggplot(first_50_codon_plot, aes(x=x,y=usage)) + geom_point()


last_50_codon <- read.table('utr_cds300.trans.last50.matrix', header = F, row.names = 1)
# last_50_codon_plot <- last_50_codon %>% map_dbl(mean) %>% as.data.frame()
last_50_codon_plot <- last_50_codon %>% map_dbl(geometric.mean) %>% as.data.frame()
names(last_50_codon_plot) <- 'usage'
last_50_codon_plot$x <- -50:-1
last_50_codon_plot$type <- 'Distance from stop(codon)'
ggplot(last_50_codon_plot, aes(x=x,y=usage)) + geom_point()

codon_analysis <- bind_rows(first_50_codon_plot, last_50_codon_plot)
ggplot(codon_analysis, aes(x=x,y=usage)) + geom_point(color='#A67DBA',size=2) + facet_grid(~type, scales = 'free_x', switch = 'x') + theme_bw(base_size = 20) + theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.background = element_blank(), strip.placement = "outside", legend.position = c(1, 1), legend.justification = c(1, 1), legend.title = element_blank(), legend.background = element_blank(), plot.margin = unit(rep(1.5,4), 'lines'), axis.title.y = element_text(vjust=3)) + ylab('Relative codon usage') + ggsave('codon_usage.pdf', width = 7, height = 4)
