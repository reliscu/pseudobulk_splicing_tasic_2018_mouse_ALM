library(dplyr)
library(tidyr)
library(tibble)
library(qvalue)
library(data.table)

source("/mnt/lareaulab/reliscu/code/fisher_test.R")

get_module_enrichments <- function(network_dir, ctype_gene_list, mod_def="PosFDR") {
    # Traverse networks to get cell type enrichments for each module
    networks <- list.dirs(file.path(getwd(), network_dir), full.names=TRUE, recursive=FALSE)
    networks <- networks[lengths(lapply(networks, list.files)) > 0]

    enrichments_list <- lapply(seq_along(networks), function(i) {
        kME_path <- list.files(networks[i])[grep("kME", list.files(networks[i]))]
        kME <- fread(file.path(networks[i], kME_path), data.table=FALSE)
        mod_assignment_col <- grep(mod_def, colnames(kME))
        mod_genes <- tapply(kME$Gene, kME[,mod_assignment_col], list)
        
        if (length(mod_genes) > 0) {
            all_genes <- kME$Gene

            # For each module: calculate enrichment for DE genes from each cell type
            mod_enrichments_list <- lapply(mod_genes, function(mod) {
                lapply(unlist(lapply(ctype_gene_list, function(set) {
                    fisher_test(set, mod, all=all_genes)
                })), c)
            })
            
            # Save the network the module came from
            network_id <- sapply(strsplit(networks[i], "/"), function(x) x[length(x)])
            mod_enrichments_df <- reshape2::melt(mod_enrichments_list)
            colnames(mod_enrichments_df) <- c("Pval", "Cell_type", "Module")
            
            # Save path to module eigengenes table for downstream analyses
            ME_path <- list.files(networks[i])[grep("eigengene", list.files(networks[i]))]
            
            data.frame(
                Network=network_id,
                kME_path=file.path(networks[i], kME_path),
                ME_path=file.path(networks[i], ME_path),
                mod_enrichments_df
            )
        }
    })
    enrichments_df <- do.call(rbind, enrichments_list)
    enrichments_df$Qval <- qvalue(enrichments_df$Pval)$qvalue

    enrichments_df
}

prep_DE_genes <- function(res_list, lfc_threshold=2, pairwise=FALSE, unique=TRUE) {

    pval_threshold <- .05

    if (pairwise) {
        ctypes <- unique(sapply(strsplit(names(res_list), "_vs_"), "[", 1))
        ctype_genes <- lapply(ctypes, function(target) {
            # Subset to pairwise tests with target cell type
            ctype_res_list <- res_list[grep(paste0("^", target), names(res_list))]
            # For each pairwise test, return genes that meet p-value threshold:
            ctype_genes_list <- lapply(ctype_res_list, function(df) {
                mask <- df$adj.P.Val < pval_threshold
                df[mask, 1]
            }) 
            # Restrict to genes that were identified in EVERY pairwise test
            Reduce(intersect, ctype_genes_list)
        })
        names(ctype_genes) <- ctypes

    } else {
        ctype_genes <- lapply(res_list, function(df) {
            mask <- (df['adj.P.Val'] < pval_threshold) & (abs(df['logFC']) > lfc_threshold)
            df[mask, 1]
        })
        names(ctype_genes) <- names(res_list) 
    }

    if (unique) {
        all_genes <- unlist(ctype_genes)
        duplicates <- unique(names(table(all_genes)[table(all_genes) > 1]))
        ctype_genes <- lapply(ctype_genes, function(x) x[!(x %in% duplicates)])
    }
    
    ctype_genes[lengths(ctype_genes) > 0]

}