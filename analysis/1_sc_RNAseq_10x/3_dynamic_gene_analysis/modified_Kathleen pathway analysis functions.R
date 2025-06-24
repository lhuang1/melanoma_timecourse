# These are the functions used in pathway analysis
library("GSA")
library("qvalue")
source(file = "/singerlab/kathleen/function scripts/gene name conversions.R")


#####################################################################################
# function to make .gmt file from .txt files
# input and output must have full path
# input is the folder we are drawing paths from
# output is the full path and name for the new gmt file
# .txt files must be of the format below, genes on their own line
#name (name of this pathway)
#source (where we got this pathway from, paper, online resource, etc)
#genes


make_gmt <- function(input_folder, output_file) {
  signature_files <- list.files(path = input_folder, full.names = TRUE)
  gmt_output <- c()
  
  for (file in signature_files) {
    signature <- scan(file, what = character(), quiet = TRUE, sep = "\n")
    # use this if you put spaces in your source names, no spaces anywhere
    #signature[2] <- gsub(" ", "_", signature[2])
    signature <- signature[signature != ""]
    
    gmt_output <- c(gmt_output, paste(signature, collapse = "\t"))
  }
  cat(paste(gmt_output, collapse = "\n"), "\n", file = output_file)
}



########################################################################################################

# this function processes de_genes against many annotated lists
pathway_annotated_loop <- function(de_genes, background, annotated_lists) {
  # work with the de and background genes because theyre the same for every list
  original_de_genes_size <- length(de_genes)
  
  #make sure we don't have duplicates
  de_genes <- de_genes[!duplicated(de_genes)]
  background <- background[!duplicated(background)]
  
  # filter the genes to only have ones in the background
  de_genes <- intersect(de_genes, background)
  
  # calculate the numbers we need for the hypergeometric test
  population_size <- length(background)
  sample_size <- length(de_genes)
  
  # make a data frame to collect our results
  data_results <- data.frame(term=character(), term_size=numeric(), terms_in_bkgd=numeric(), frac_term_in_bkgd=double(), odds_ratio=double(), sample_size=double(), num_in_term=numeric(), frac_in_term=double(), p_val=double(), website=character(), genes=character(), stringsAsFactors = FALSE)
  
  for (i in 1:length(annotated_lists[["geneset.names"]])) {
    annotated <- annotated_lists$genesets[[i]]
    ann_list_size <- length(annotated)
    
    #make sure annoatated list genes are unique
    annotated <- annotated[!duplicated(annotated)]
    
    #limit the annotated and de_genes list to genes that are in the background
    annotated <- intersect(annotated, background)
    genes <- intersect(de_genes, annotated)
    
    num_in_term <- length(genes)
    
    # calculate the terms for the hypergeometric test and results table
    frac_in_term <- num_in_term / sample_size
    frac_term_in_bkgd <- length(annotated) / ann_list_size
    ann_in_background <- length(annotated)
    
    # we want a filter because with less than 3 terms overlapping it doesnt mean much
    if (num_in_term < 3) {
      p_val <- 1
      odds_ratio <- 0
      
    } else {
      
      # calculate the p_val
      population_fail <- population_size - ann_in_background
      p_val <- 1 - phyper(num_in_term - 1, ann_in_background, population_fail, sample_size)
      
      #Calculate the odds ratio
      num_out_term <- sample_size - num_in_term
      double_fail <- population_fail - num_out_term
      # odds_ratio <- (num_in_term * double_fail) / (num_out_term * (ann_list_size - num_in_term))
      odds_ratio <- (num_in_term * double_fail) / (num_out_term * (ann_in_background - num_in_term))
      
    }
    
    website <- annotated_lists$geneset.descriptions[[i]]
    term <- annotated_lists$geneset.names[[i]]
    genes <- paste(genes, collapse = " ")
    
    data_results <- rbind(data_results, data.frame(term, ann_list_size, ann_in_background, frac_term_in_bkgd, odds_ratio, sample_size, num_in_term, frac_in_term, p_val, website, genes, stringsAsFactors = FALSE))
    
  }
  
  # sort the results table
  data_results <- data_results[order(data_results$p_val),]
  data_results$fdr <- tryCatch({qvalue(p = data_results$p_val, fdr.level = cheat_sheet_q_val)$qvalues}, error = function(e) {return(1)})
  
  # create a data structure to return
  
  all_info <- list(results_data = data_results, original_de_genes_size=original_de_genes_size, de_genes_in_background = sample_size, background = population_size)
  
  return(all_info)
}


##########################################################################################################################

default_gene_sets <- list( "KeggReactome" = "c2.kegg.reactome.v7.0.symbols.gmt",
                           "BiocartaPid" = "c2.pid.biocarta.v7.0.symbols.gmt",
                           "Perturbations" = "c2.cgp.v7.0.symbols.gmt",
                           "GoTerms" = "c5.all.v7.0.symbols.gmt",
                           "c7Immune" = "c7.all.v7.0.symbols.gmt",
                           "Literature" = "literature.symbols.gmt",
                           "Hallmark" = "h.all.v7.0.symbols.gmt"
)

#this function processes de_genes against all the annotated lists I downloaded from msigdb
standard_pathway_analysis <- function(de_genes, de_genes_name, background, background_name, dir_results, set_lists = default_gene_sets) {
  
  if (length(de_genes) > 2) {
    
    cheat_sheet <- data.frame(term=character(), term_size=numeric(), terms_in_bkgd=numeric(), frac_term_in_bkgd=double(), odds_ratio=double(), num_in_term=numeric(), frac_in_term=double(), p_val=double(), website=character(), genes=character(), fdr=double(), set_list=character(), stringsAsFactors = FALSE)
    
    for (GeneSet in names(set_lists)) {
      
      # this wrapper should suppress all the junk that prints to the console
      sink("/dev/null")
      ann_lists <- GSA.read.gmt(paste("/singerlab/linglin/b16/0_data/gene_lists/GMT/", set_lists[[GeneSet]], sep = ""))
      sink()
      
      # convert to upper case (added by Linglin)
      ann_lists$genesets <- lapply(ann_lists$genesets, toupper)
      de_genes <- toupper(de_genes)
      background <- toupper(background)
      
      set_results <- pathway_annotated_loop(de_genes = de_genes, background = background, annotated_lists = ann_lists)
      file_name <- paste(dir_results, de_genes_name, " ", GeneSet, " ", background_name, ".csv", sep = "")
      
      write(paste(GeneSet, ", ", set_results$background, " genes in background, ", set_results$original_de_genes_size, " genes entered, ", set_results$de_genes_in_background, " genes found in background\n"), file_name)
      data_out <- set_results$results_data
      write.table(data_out, file = file_name, append = TRUE, quote = TRUE, row.names = FALSE, sep = ",")
      
      # cheat sheet agreggation
      short_list <- data_out[data_out$fdr <= cheat_sheet_q_val,]
      if (nrow(short_list) > 0 ) {
        short_list$set_list <- GeneSet
        cheat_sheet <- rbind(cheat_sheet, short_list)
      }
      
    }
    
    # write cheat sheet
    write.table(cheat_sheet, file = paste(dir_results, de_genes_name, " ", background_name, " cheat sheet.csv", sep = ""), quote = TRUE, row.names = FALSE, sep = ",")
  }
  
}