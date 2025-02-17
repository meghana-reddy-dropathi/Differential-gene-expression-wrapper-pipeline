#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("patchwork"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("yaml"))
suppressPackageStartupMessages(library("data.table"))
source("./Utils.R")
options <-list(make_option("--run", type = "character", default = NULL,
                           help = "Path to the config file."))
# Loading all the parameters of the config file
print("Loading the required data variables")
option_parser <- OptionParser(option_list = options)
options <- parse_args(option_parser)
config_options <- read_yaml(options$run_config)
date_of_run <- config_options$date_of_run
id <- config_options$id 
output_dir <- config_options$output_dir 
metadata <- config_options$path_to_metadata
metadata_path_for_tar <- config_options$metadata_path_for_tar 
metadata_path_for_dbd <- config_options$metadata_path_for_dbd
path_to_counts <- config_options$path_to_counts 
path_to_tpm <- config_options$path_to_tpm
ensembl_id_path <- config_options$ensembl_id_path 
target_genes_list <- config_options$target_genes
counts_filter_level <- config_options$deseq2$counts_filter_level 
independent_filter <- config_options$deseq2$independent_filter 
d_padj <- as.numeric(config_options$deseq2$padj_cutoff)
d_l2fc <- as.numeric(config_options$deseq2$log_fold_change_cutoff)
path_to_condition_pairs <- config_options$path_to_condition_pairs 
include_protein_coding <- config_options$protein_coding_only 
collapse_tech_reps <- config_options$collapse_tech_reps
# setting the initial x axis and y axis limits
x_limits <- c(Inf, Inf)
y_limits <- c(Inf, Inf)
#Creating all the output directories
print("Creating the required output directories")
if (!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}else{
  unlink(output_dir, recursive = TRUE)
  print("deleting the pre-exisiting directory")
}
de_unfiltered_output <- paste0(output_dir, "/diff_exp_unfiltered")
dir.create(de_unfiltered_output, recursive = TRUE)
ind_output <- paste0(output_dir, "/plots/individual")
dir.create(ind_output, recursive = TRUE)
dispests_output <- paste0(output_dir, "/plots/qc_dispests")
dir.create(dispests_output, recursive = TRUE)
diff_genes_output <- paste(output_dir, "/diff_exp-genes")
dir.create(diff_genes_output, recursive = TRUE)
# Creating the gene targets dataframe
gene_targets <- as.data.frame(do.call(rbind, target_genes_list))
# Loading data
metadata <- read.delim(metadata, header = TRUE, sep = ",", row.names = "ngs_id")
# Set condition column as factors, required for DESeg2
metadata$condition <- as.factor(metadata$condition)
targeting_metadata <- read.delim(metadata_path_for_tar, header = TRUE, sep = ",")
dbd_effector_metadata <- read.delim(metadata_path_for_dbd, header = TRUE, sep = ",")
print("reading the counts")
counts <- read.delim(path_to_counts, header = TRUE, sep = "\t", row.names = 1)
counts <- subset(counts, select = -transcript_id.s.)
counts <- counts %>%
  mutate(across(.col = where(is.double), .fns = as.integer))
print("reading tpm")
tpm <- read.delim(path_to_tpm, header = TRUE, sep = "\t", row.names = 1)
tpm <- subset(tpm, select = -transcript_id.s.)
print("reading biomart map")
ensembl_id_map <- read.delim(ensembl_id_path, sep = "\t")
print("reading condition pairs")
condition_pairs <- read.delim(path_to_condition_pairs, header = TRUE, sep = ",")
# Running DESeq2
average_tpm <- get_condition_average(df=tpm, metadata=metadata)
average_tpm<-average_tpm[,-1]
# Instantiate summary data frame
summary_data <- data.frame(condition = character(),
                           reference_condition= character(),
                           number_de_genes = integer(),
                           de_results_file = character())
for (row in 1:nrow(condition_pairs)){
  cond <- condition_pairs[row, "condition"]
  ref_cond <- condition_pairs[row, "reference_condition"]
  print(paste0("Running DESeq2 for: ", cond))
  dds <- create_dds_object(counts = counts,
                           metadata = metadata,
                           conditions = c(cond, ref_cond),
                           reference_condition = ref_cond,
                           filter_level = counts_filter_level
                           )
  if (collapse_tech_reps == TRUE){
    dds = collapseReplicates(dds, groupby = dds$sample_id, renameCols=TRUE)
  }else{
    print("No technical replicates were found")
  }
  dds <- DESeq(dds)
  results_df <- get_results(dds, independent_filtering = TRUE)
  results_file_path <- file.path(de_unfiltered_output, paste0("de_genes_", cond, "_unfil.txt"))
  write.table(results_df, results_file_path, sep = '\t', row.names = FALSE)
  de_genes <- get_all_significant_genes(results_df,
                                        targets = NULL,
                                        padj_cutoff = d_padj,
                                        log_cutoff = d_l2fc)
  de_genes_file_path <- file.path(diff_genes_output, paste0("de_genes_", cond,".txt"))
  write.table(de_genes, de_genes_file_path, sep = "\t", row.names = FALSE)
  summary_df <- summary_data%>%
    add_row(tibble_row(condition = cond,
                       reference_condition = ref_cond,
                       number_de_genes = nrow(de_genes),
                       de_results_file = paste0("de_genes_", ".txt")))
  x_limits <- sapply(x_limits, function(limit){
    if (limit == "Inf"){
      if (limit == x_limits[1]){
        return(-Inf)
      }else{
        return(Inf)
      }
    }
    as.numeric(limit)
  }, USE.NAMES = FALSE)
  y_limits <- sapply(y_limits, function(limit){
    if (limit == "Inf"){
      if (limit == y_limits[1]){
        return(-Inf)
      }else{
        return(Inf)
      }
    }
    as.numeric(limit)
  }, USE.NAMES = FALSE)
  print("Plotting various plots")
  condition_annot <- get_condition_annotations(cond,
                                               sample_metadata = metadata,
                                               target_metadata = targeting_metadata,
                                               effector_metadata = dbd_effector_metadata)
  ref_condition_annot <- get_condition_annotations(ref_cond,
                                                   sample_metadata = metadata,
                                                   target_metadata = targeting_metadata,
                                                   effector_metadata = dbd_effector_metadata)
  
  experiment_subtitle <- condition_annot$experiment_name
  condition_title <- get_condition_plot_title(condition_annot, ref_condition_annot)
  reference_condition_title <- get_condition_plot_title(ref_condition_annot, condition_annot)
  plot_title <- paste0(condition_title, "plot")
  plot_filename <- paste0(cond,"plot")
  unique_condition_zfp <- unique(condition_annot$zfp_id)
  cond_name <- unique(condition_annot$condition)
  ref_name <- unique(ref_condition_annot$condition)
  
  # plotting the dispersion data
  plot_dis(dds, title = plot_title, filename = file.path(dispects_output,paste0(plot_filename, "_dispect_plot.jpg")))
  save_plot <- function(plot_function,plot_args, file_suffix, protein_coding_only=FALSE, width=10, height=12){
    plot_args$protein_coding_only <- protein_coding_only
    p <- do.call(plot_function,plot_args)
    if(protein_coding_only){
      file_suffix <- paste0(file_suffix,"_pro_cod")
    }
    ggsave(filename = file.path(ind_output, paste0(plot_filename, file_suffix,".jpg")),
           plot = p,
           width = width,
           height = height,
           units = "cm")
  }
  plot_args <- list(df= results_df,
                    targets= gene_targets,
                    title = stringr::str_wrap(plot_title),
                    padj = d_padj,
                    log = d_l2fc,
                    protein_coding_only=FALSE,
                    top_n=0,
                    x_limits = x_limits,
                    y_limits = y_limits,
                    cond_name,
                    ref_name)
  save_plot(plot_volcano,plot_args,"_volcano_plot")
  save_plot(plot_ma, plot_args,"ma_plot")
  if(include_protein_coding){
    save_plot(plot_volcano,plot_args,"_volcano_plot",TRUE)
    save_plot(plot_ma,plot_args, "ma_plot",TRUE)
  }
}
summary_file_path = file.path(output_dir, paste0(tolower(id),"de_genes_summary.txt"))
write.table(summary_df, summary_file_path, sep="\t", row.names= FALSE)
print(paste0("DESeq2 output directory is in the following file path:",output_dir))
                           