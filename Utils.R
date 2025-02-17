suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("GGally"))
suppressPackageStartupMessages(library("ggrepel"))
suppressPackageStartupMessages(library("tidyverse"))

get_sig_targets_rows <- function(results_df, targets, padj_cutoff, l2fc_cutoff) {
  # Given dataframe From Desegz results, return filtered dataframe
  # filter for targets in result_at and
  sig_targets <- filter(results_df,
                        padj <= padj_cutoff & (log2FoldChange <= -(l2fc_cutoff) | log2FoldChange >= l2fc_cutoff))
  if (!is.null(targets)){
    sig_targets <- filter(sig_targets, gene_id %in% targets$ensembl_id)
    return(nrow(sig_targets))
  }else{
    return(0)
  }
}


get_gene_ids <- function(target_genes) {
  # creates a list
  genes <- list()
  for (tg in target_genes){
    gene <- append(gene, tg$ensembl_id)
  }
  return(unlist(gene))
}


replicates<- function(condition,metadata){
  replicates <- metadata %>%
    filter(condition == condition) %>% 
    select(replicate) %>% 
    arrange(replicate)
  return(replicates)
}

replicates_tpm<-function(df,tpm_df){
  dataframe_tpm<-tpm_df%>%
    select(row.names(df))
  return(dataframe_tpm)
}

create_dds_object <- function(counts,
                              metadata,
                              conditions ,
                              reference_condition,
                              filter_level){
  column_data <- metadata %>%
    filter(condition %in% conditions) %>%
    select(condition, sample_name, replicate, sample_id)
  column_data$condition <- droplevels(column_data$condition)
  column_datascondition <- relevel(column_dataScondition, ref = reference_condition)
  counts_1 < select(counts, row.names(column_data))
  if (filter_level = "rowsum"){
    counts_filtered <- counts_1 %>% filter(rowSums(.) > 0)
    }
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_filtered),colData = column_data, design = ~ condition)
  
  return(dds)
}

get_results <- function(dds, independent_filtering, gene_id_map = ensembl_id_map){
  if (independent_filtering) {
    results <- results(dds)
  } else {
    results <- results(dds, independentFiltering=FALSE)
  }
  df <- merge(data.frame(results),gene_id_map,by.x = "row.names", by.y ="ensembl_gene_id")
  df <- df %>% rename(gene_id = Row.names, gene_name = external_gene_name)
  df <- df %>%
    mutate(gene_name = if_else(gene_name == '', paste(gene_id, gene_biotype, sep = '_'), gene_name))
  return(df)
}

get_all_significant_genes <- function(results_df, targets, padj_cutoff, log_cutoff) {
  significant_genes<-results_df%>%
    filter(padj <= padj_cutoff & (log2FoldChange <= -(log_cutoff) | log2FoldChange >= log_cutoff))
  if (!is.null(targets)){
    significant_genes <- filter(significant_genes, !gene_id %in% targets$ensembl_id)
  }
  significant_genes <- significant_genes %>% arrange(padj)
  return(significant_genes)
}

get_protein_coding_df <- function(df, ensembl_map, df_id_col = NULL, map_id_col = "ensembl_gene_id") {
  # Assumes biomart export with `gene_biotype` as available column
  # @type can be c("ensembl_gene_id") or transcript id
  # Returns vector of protein coding only ensembl ids
  protein_coding <- filter(ensembl_map, gene_biotype == "protein_coding")
  protein_coding_ids <- protein_coding[, map_id_col]
  if (is.null(df_id_col)) {
    df <- filter(df, rownames(df) %in% protein_coding_ids)
  } else {
    df <- filter(df, df[, df_id_col] %in% protein_coding_ids)
  }
  return(df)
}


get_condition_average <- function(df, metadata) {
  dataframe <- data.frame(matrix(nrow = nrow(df)), row.names = row.names(df))
  conditions<-levels(metadata$condition)
  for (cond in conditions ) {
    dataframe[[cond]] <- rowMeans(df[, rownames(metadata[metadata$condition == cond, ]), drop = FALSE], na.rm = TRUE)
  }
  return(dataframe)
}


get_condition_annotations <- function(condition,
                                      sample_metadata,
                                      target_metadata,
                                      effector_metadata) {
  # Assumes three metadata tables to be combined
  # Returns subset of metadata values for that condition
  metadata <- select(sample_metadata, -c(targeting_entity_name, dbd_effector_entity_name))
  metadata <- metadata[which(metadata$condition == condition), ]
  metadata$targeting_entity_id <- as.list(strsplit(metadata$targeting_entity_id, ", "))
  combined_meta <- unnest(metadata, cols = targeting_entity_id) %>%
    left_join(target_metadata, by = c("sstid", "targeting_entity_id"))
  
  combined_meta <- combined_meta %>%
    left_join(effector_metadata, by = c("sstid", "dbd_effector_entity_id"))
  # you can add in more columns you want according to the existing data.
  cond_annot <- combined_meta %>%
    select(c("experiment_condition",
             "condition",
             "experiment_name"
    ))
  cond_annot <- cond_annot %>% distinct()
  return(cond_annot)
}


get_condition_plot_title <- function(cond_1_annot, cond_2_annot=NULL) {
  
  # Extract common variables from cond_1_annot for simplicity
  exp_condition <- unique(cond_1_annot$experiment_condition)
  target_gene <- unique(cond_1_annot$target_gene)
  if (unique(cond_1_annot$condition) == unique(cond_2_annot$condition)) {
    # Compare same condition across different experiments
    title <- paste0(gsub("_", " ", cond_1_annot$condition), " - ", cond_1_annot$experiment_name)
  }else {
    title <- gsub("_", " ", cond_1_annot$condition)
  }
  return(title)
}

# Plotting functions
plot_dis <- function(dds, title, filename) {
  jpeg(filename)
  plotDispEsts(dds, main = title)
  dev.off()
}

plot_ma <- function(df, 
                    targets, 
                    title = NULL, 
                    padj, 
                    log, 
                    top_n,
                    protein_coding_only, 
                    x_limits, 
                    y_limits,
                    cond_name,
                    ref_name) {
  if (protein_coding_only) {
    df <- get_protein_coding_df(df, ensembl_map = ensembl_id_map)
  }
  
  tg <- filter(df, gene_id %in% targets$ensembl_id)
  sig_genes <- get_all_significant_genes(df, targets, padj_cutoff = padj, log_cutoff = log)
  if (top_n == 0) {
    top_sig_genes <- sig_genes
  } else { 
    top_sig_genes <- sig_genes[1:top_n, ]}
  down_genes <- filter(sig_genes, log2FoldChange <= -log)
  up_genes <- filter(sig_genes, log2FoldChange >= log)
  
  plot <- ggplot(df, aes(x = log(baseMean + 1, 10), y = log2FoldChange, label = gene_name)) +
    geom_point(alpha = 0.5, cex = 0.3, color = "blue")
  
  plot <- plot + 
    theme_bw() +
    theme(
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10)
    ) +
    labs(x = "log10 mean expression") +
    labs(y = "log2 fold change")
  
  return(plot)
}

plot_volcano <- function(df, 
                         targets, 
                         title = NULL, 
                         padj, 
                         log, 
                         top_n, 
                         protein_coding_only,
                         x_limits, 
                         y_limits,
                         cond_name,
                         ref_name) {
  plot_caption=NULL
  include_target_genes=TRUE
  # Assumes dataframe from DESeq2 results
  data_f<-subset(df,!is.na(padj))
  # getting the floating point just to handle any NAN's or minimum value
  minimum_value<-1e-323
  data_f$padj[data_f$padj == 0] <- minimum_value
  if (protein_coding_only) {
    data_f <- get_protein_coding_df(data_f, ensembl_map = ensembl_id_map, df_id_col = "gene_id")
  }
  #tg <- filter(df, gene_id == target_gene_id)
  t_df <- filter(data_f, gene_id %in% targets$ensembl_id)
  sig_genes <- get_all_significant_genes(data_f, targets, padj_cutoff = padj, log_cutoff = log)
  if (top_n == 0) {
    top_sig_genes <- sig_genes
  } else { 
    top_sig_genes <- sig_genes[1:top_n, ]}
  down_genes <- filter(top_sig_genes, log2FoldChange <= -log)
  up_genes <- filter(top_sig_genes, log2FoldChange >= log)
  genes_removed_null_padj <- data_f %>% filter(padj>0)
  xlim_max <- max(abs(data_f$log2FoldChange), na.rm = TRUE)
  
  num_de_genes <- nrow(sig_genes) + get_sig_targets_rows(data_f, targets, padj_cutoff = padj, l2fc_cutoff = log)
  for (i in 1:nrow(t_df)) {
    if (!is.na(t_df[i, ]$padj) && !is.na(t_df[i, ]$log2FoldChange)) {
      if (t_df[i, ]$padj <= padj && abs(t_df[i, ]$log2FoldChange) >= log) {
        num_de_genes <- num_de_genes + 1
      }
    }
  }
  
  if (all(x_limits == c(Inf, Inf))) { x_limits = NULL }
  if (all(y_limits == c(Inf, Inf))) { y_limits = NULL }
  
  # not considering the na values
  df_non_na<-data_f
  df_non_na$log<-(-log(data_f$padj,10))
  dataf_non_na<-subset(df_non_na,!is.na(log))
  x_limits<-c(min(data_f$log2FoldChange),max(data_f$log2FoldChange))
  y_limits<-c(min(dataf_non_na$log),max(dataf_non_na$log))
  maxi<-max(abs(x_limits))
  x_limits_1<-c(-maxi,maxi)
  ra=sum(abs(y_limits))*1.5/sum(abs(x_limits_1))
  plot <- ggplot(dataf_non_na, aes(x = log2FoldChange, y = -log(padj, 10), label = gene_name,na.rm=TRUE))+xlim(x_limits_1)+ylim(y_limits)+
    geom_point(alpha = 1, color = "#72B7BC")
  
  tg_up_bound = 10
  plot <- plot + 
    geom_text_repel(data = top_sig_genes[1:tg_up_bound, ], max.overlaps = 20, size = 3)
  
  if (include_target_genes) {
    gene<-list(get_list_supp_fun(targets,cond_name,ref_name,df))
    filtered_df <- t_df[t_df$gene_name %in% gene[[1]], ]
    # to not print caption when the gene ist is null for the conditions you don't want to print target genes
    if(is.null(gene[[1]])){
      plot_caption=NULL
    }else{
      plot_caption <- NULL
    }
    plot <- plot + 
      geom_point(data = filtered_df, shape=18,cex=5,color="pink")+
      geom_text_repel(data =  filtered_df) 
  }
  plot <- plot + 
    theme_bw() +
    labs(title = title, subtitle = paste0("DE genes: ", num_de_genes)) +
    labs(x = "log2 fold change") +
    labs(y = "-log10 padj")+ theme(panel.border = element_rect(color = "black", fill = NA, size = 2))+
    labs(caption = plot_caption)
  
  
  return(plot)
}

get_list_supp_fun<-function(targets,cond_name,ref_name,df){
  # this is a function which returns a list of target genes to be plotted
  df_1<-targets
  rownames(df_1)<-NULL
  gc<-df_1%>%
    select(gene_name,condition)%>%
    mutate(gene_name=unlist(gene_name))%>%
    unnest(condition)
  cd <- gc %>%
    group_by(gene_name) %>%
    summarize(condition = list(condition), .groups = 'drop')
  # Convert to named list (dictionary)
  cdl <- setNames(cd$condition, cd$gene_name)
  keys<-names(cdl)
  g<-c()
  cond_name<-as.character(cond_name)
  fl<-list()
  # iterating through the dictionary 
  for(i in 1:length(cdl)){
    # if the conditions have no values, the conditions will be printed 
    if(any(cdl[[i]]=="NULL")){
      fl<-c(fl,names(cdl[i]))
    }else{
      gene_name_1 <- keys[i]
      
      # Check if the condition or reference name is present in `g`
      if (any(c(cond_name, ref_name) %in% g)) {
        gene_cond_new <- gc %>%
          filter(gene_name == gene_name_1) %>%
          mutate(new = condition %in% c(cond_name, ref_name))
        
        # Collect the gene names where `new` is TRUE
        fl <- c(fl, unname(gene_cond_new$gene_name[gene_cond_new$new]))
      }
    }
  }
  # returning the list 
  return(unlist(fl))
}
