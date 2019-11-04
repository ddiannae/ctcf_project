library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
setwd("~/Workspace/ctcf_project/R/")
mycolors <- 

annot <- read_tsv("/labs/csbig/subtipos_mama_2018/Biomart_EnsemblG94_GRCh38_p12_karyotype.txt",
                    col_names = c("ensemblID", "chr", "start", "end", "karyoband"),
                    col_types = cols(col_character(), col_character(), col_integer(), 
                                     col_integer(), col_character()),
                    skip = 1)
genes <- read_tsv("/labs/csbig/subtipos_mama_2018/genes_in_exp_matrix.txt", col_names = "ensemblID")

annot <- annot %>% semi_join(genes)

ctcfs <- read_tsv("../data/non_gene_overlaping_ctcfbs.tsv", 
                  col_types = cols(col_character(), col_integer(), col_integer()))

gene_counts <- annot %>% group_by(chr) %>% tally()
gene_counts$chr <- factor(gene_counts$chr, levels = as.character(c(seq(1:22), "X")))

png(filename = "../figures/gene_count.png", width = 1200, height = 600)
ggplot(gene_counts, aes(x = chr, y = n, fill = chr)) +
  geom_bar(stat = "identity") +
  theme_few(base_size = 20) +
  theme(legend.position = "none") +
  ylab("Genes") +
  xlab("Chromosome") +
  scale_fill_manual(values = getPalette(23))
dev.off()

ctcfs_counts <- ctcfs %>% group_by(chr) %>% tally()
ctcfs_counts$chr <- factor(ctcfs_counts$chr, levels = as.character(c(seq(1:22), "X")))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))

png(filename = "../figures/ctcf_count.png", width = 1200, height = 600)
ggplot(ctcfs_counts, aes(x = chr, y = n, fill = chr)) +
  geom_bar(stat = "identity") +
  theme_few(base_size = 20) +
  theme(legend.position = "none") +
  ylab("CTCF binding sites") +
  xlab("Chromosome") +
  scale_fill_manual(values = getPalette(23))
dev.off()

chrs <- as.character(c(seq(1:22), "X"))

gb_ctfs <- parallel::mclapply(X = chrs, mc.cores = 5, FUN = function(c){
  gs_in_c <- annot %>% filter(chr == c)
  cs_in_c <- ctcfs %>% filter(chr == c) 
  cs_in_c <- add_row(cs_in_c, chr = c, start = 0, end = 0) %>% arrange(start)
  
  genes <- lapply(seq(1:(nrow(cs_in_c)-1)), function(i){
    s <- unlist(cs_in_c[i, "end"])
    e <- unlist(cs_in_c[i+1, "start"])
    
    ngenes <- gs_in_c %>% filter(start >= s & end <= e) %>% nrow()
    return(list(i = i, ngenes = ngenes))
  })
  genes <- bind_rows(genes)
  ngenes <- gs_in_c %>% filter(start >= max(cs_in_c$end) & end <= max(gs_in_c$end)) %>% nrow()
  genes <- genes %>% add_row(i = nrow(cs_in_c), ngenes = ngenes)
  genes$cum_sum <- cumsum(genes$ngenes)
  cs_in_c <- bind_cols(cs_in_c, genes)
  
  return(cs_in_c)
})

gb_ctfs <- bind_rows(gb_ctfs)

e_starts <- gb_ctfs %>% group_by(chr) %>% filter( i == 1)
e_ends <- gb_ctfs %>% group_by(chr) %>% filter( i == max(i))
g_betw <- gb_ctfs %>% group_by(chr) %>% filter( i > 1 & i < max(i))
g_betw$chr <- factor(g_betw$chr, levels = as.character(c(seq(1:22), "X")))

ceros <- g_betw %>% filter(ngenes == 0) %>% group_by(chr) %>% tally()
ceros$empty = "Empty"
noceros <- g_betw %>% filter(ngenes != 0) %>% group_by(chr) %>% tally()
noceros$empty = "Non Empty"

ctcfs_counts <- bind_rows(ceros, noceros)

png(filename = "../figures/ctcf_count_empty.png", width = 1200, height = 600)
ggplot(ctcfs_counts, aes(x = chr, y = n)) +
  geom_bar(
    aes(color = empty, fill = empty),
    stat = "identity", position = position_stack()
  ) +
  theme_few(base_size = 20) +
  theme(legend.title = element_blank()) +
  ylab("CTCF binding sites") +
  xlab("Chromosome") +
  scale_fill_tableau() +
  scale_color_tableau() 
dev.off()

noceros <- g_betw %>% filter(ngenes != 0) 
png(filename = "../figures/ctcf_genes_between.png", width = 1200, height = 600)
ggplot(noceros, aes(x = chr, y = ngenes, fill = chr)) +
  geom_boxplot() +
  theme_few(base_size = 20) +
  theme(legend.position = "none") +
  ylab("Genes between CTCF binding sites") +
  xlab("Chromosome") +
  scale_fill_manual(values = getPalette(23))
dev.off()

##### Cumulative sum to remove ctcfs with no genes
unique_ctcfs <- gb_ctfs %>% distinct(chr, cum_sum, .keep_all = TRUE)
unique_ctcfs$id <- paste0("CTCF", stringi::stri_pad_left(rownames(unique_ctcfs), 5, 0))
unique_ctcfs <- unique_ctcfs %>% select(id, chr, start, end, i, ngenes, cum_sum)
write_tsv(unique_ctcfs, path = "../data/significant_ctcfs.tsv")
