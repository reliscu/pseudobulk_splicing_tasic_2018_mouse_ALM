source("/mnt/lareaulab/reliscu/code/fisher_test.R")
source("/mnt/lareaulab/reliscu/projects/NSF_GRFP/analyses/pseudobulk_test/tasic_2018/mouse_ALM/code/enrichment_fxns.R")

library(ggplot2)

source("/mnt/lareaulab/reliscu/code/ggplot_theme.R")

theme_set(default_theme())

get_top_corr_mods <- function(network_dir, pseudobulk_legend, top_qval_mods_df, ctype_genes_list, mod_def) {
    top_qval_mods_df <- top_qval_mods_df[top_qval_mods_df$Cell_type %in% names(ctype_genes_list),]
    ctypes <- unique(top_qval_mods_df$Cell_type) 
    top_corr_mods_list <- vector(mode="list", length=length(ctypes))
    total_cells_per_sample <- colSums(pseudobulk_legend[,-c(1, 2)])
    
    for (i in seq_along(ctypes)) {
        # Get # cells in each sample
        mask <- pseudobulk_legend$Cell.type == ctypes[i]
        n_cells_per_sample <- colSums(pseudobulk_legend[mask, -c(1, 2)])
        frac_per_sample <- n_cells_per_sample / total_cells_per_sample
        
        ctype_genes <- ctype_genes_list[[which(names(ctype_genes_list) %in% ctypes[i])]]

        if (var(frac_per_sample) > 0) {
            # Get most enriched cell type module
            old_mod <- top_qval_mods_df$Module[i]
            ME_df <- fread(top_qval_mods_df$ME_path[i], data.table=FALSE)
            ME_vec <- ME_df[,grep(paste0("^", old_mod, "$"), colnames(ME_df))]
            old_corr <- cor(frac_per_sample, ME_vec)

            # Traverse networks to find module most correlated to cell type abundance
            networks <- list.dirs(file.path(getwd(), network_dir), full.names=TRUE, recursive=FALSE)
            networks <- networks[lengths(lapply(networks, list.files)) > 0]

            ME_corrs_list <- lapply(seq_along(networks), function(j) {
                kME_path <- list.files(networks[j])[grep("kME", list.files(networks[j]))]
                ME_path <- list.files(networks[j])[grep("eigengene", list.files(networks[j]))]
                ME_df <- fread(file.path(networks[j], ME_path), data.table=FALSE)
                ME_corrs <- apply(ME_df[,-1, drop=FALSE], 2, function(ME) {
                    cor(ME, frac_per_sample)
                })
                new_mod <- names(which.max(ME_corrs))
                new_corr <- ME_corrs[which.max(ME_corrs)]
                
                # Save the network the new module came from
                network_id <- sapply(strsplit(networks[j], "/"), function(x) x[length(x)])
                
                data.frame(
                    Cell_type=ctypes[i],
                    Pseudobulk_SD=round(sd(n_cells_per_sample), 2),
                    Cor=new_corr,
                    Old_cor=old_corr,
                    Pval=NA,
                    Old_pval=top_qval_mods_df$Pval[i],
                    Module_genes=NA,
                    Old_module_genes=NA,
                    DE_genes=paste(ctype_genes[1:15], collapse=", "),
                    Module=new_mod,
                    Old_module=old_mod,
                    Network=network_id,
                    Old_network=top_qval_mods_df$Network[i],
                    ME_path=file.path(networks[j], ME_path),
                    Old_ME_path=top_qval_mods_df$ME_path[i],
                    kME_path=file.path(networks[j], kME_path),
                    Old_kME_path=top_qval_mods_df$kME_path[i]
                )
            })
            ME_corrs <- do.call(rbind, ME_corrs_list)
            
            top_corr_mods_list[[i]] <- ME_corrs %>%
                arrange(Network, Module) %>%
                slice_max(Cor, with_ties=F)
        }
    }

    top_corr_mods_df <- do.call(rbind, top_corr_mods_list)

    # Get top corr module enrichment for cell type genes
    enrich_pvals <- get_top_corr_mods_enrichments(top_corr_mods_df, ctype_genes_list, mod_def)
    top_corr_mods_df$Pval <- enrich_pvals

    # Get module genes 
    new_mod_genes <- lapply(1:nrow(top_corr_mods_df), function(i) {
        paste(
            get_mod_genes(top_corr_mods_df$kME_path[i], top_corr_mods_df$Module[i], mod_def, n_genes=15), 
            collapse=", "
        )
    }) 
    old_mod_genes <- lapply(1:nrow(top_corr_mods_df), function(i) {
        paste(
            get_mod_genes(top_corr_mods_df$Old_kME_path[i], top_corr_mods_df$Old_module[i], mod_def, n_genes=15), 
            collapse=", "
        )
    }) 
    top_corr_mods_df$Module_genes <- unlist(new_mod_genes) 
    top_corr_mods_df$Old_module_genes <- unlist(old_mod_genes)

    # Finally, add cell type DE genes


    top_corr_mods_df %>%
        arrange(Cor)
}

get_top_corr_mods_enrichments <- function(top_corr_mods_df, ctype_genes_list, mod_def="PosFDR") {
    enrich_pvals <- c()
    for (i in 1:nrow(top_corr_mods_df)) {
        # Get genes for most correlated module
        kME <- fread(top_corr_mods_df$kME_path[i], data.table=FALSE)
        mod_assignment_col <- grep(mod_def, colnames(kME))
        mod_genes <- kME$Gene[kME[,mod_assignment_col] %in% top_corr_mods_df$Module[i]]
        if (length(mod_genes) > 0) {
            all_genes <- kME$Gene
            ctype_genes <- ctype_genes_list[[which(names(ctype_genes_list) %in% top_corr_mods_df$Cell_type[i])]]     
            if (length(ctype_genes) > 0) {
                enrich_pvals[i] <- fisher_test(ctype_genes, mod_genes, all=all_genes)
            } else {
                enrich_pvals[i] <- "No DE genes"
            }
        } else {
            enrich_pvals[i] <- "No kME genes"
        }
    }
    enrich_pvals
}

get_mod_genes <- function(kME_path, module, mod_def="PosFDR", n_genes=10) {
    kME <- fread(kME_path, data.table=FALSE)
    kME <- kME[order(-kME[,paste0("kME", module)]),]
    mod_assignment_col <- grep(mod_def, colnames(kME))
    mod_genes <- kME$Gene[kME[,mod_assignment_col] %in% module]
    mod_genes[1:n_genes]
}

get_top_corr_mod_stats <- function(top_corr_mods_df) {
    ctypes <- top_corr_mods_df$Cell_type

    top_corr_mods_df$Mod_stats_path <- sapply(top_corr_mods_df$ME_path, function(x) gsub("eigengenes", "statistics", x))

    mod_stats <- fread(top_corr_mods_df$Mod_stats_path[1], data.table=FALSE)
    col_idx <- which(colnames(mod_stats) %in% c("PC1VE", "MeanExpr", "Specificity", "Homogeneity"))
    col_idx <- sort(c(col_idx, grep("^Unique", colnames(mod_stats))))

    stats <- c()
    for (i in seq_along(ctypes)) {
        mod <- top_corr_mods_df$Module[i]
        mod_stats <- fread(top_corr_mods_df$Mod_stats_path[i], data.table=FALSE)
        colnames(mod_stats)[grep("^Unique", colnames(mod_stats))] <- "Unique_members" 
        stats <- rbind(stats, mod_stats[mod_stats$Module == mod, col_idx])
    }

    top_corr_mods_df <- cbind(top_corr_mods_df, stats)
    top_corr_mods_df
}

plot_ctype_abundance_vs_top_qval_ME <- function(pseudobulk_legend, top_qval_mods_df) {
    options(repr.plot.width=6, repr.plot.height=6)

    total_cells_per_sample <- colSums(pseudobulk_legend[,-c(1, 2)])
    ctypes <- top_qval_mods_df$Cell_type

    for (i in seq_along(ctypes)) {
        # Get cell type proportion in each sample
        mask <- pseudobulk_legend$Cell.type == ctypes[i]
        n_cells_per_sample <- colSums(pseudobulk_legend[mask, -c(1, 2)])
        frac_per_sample <- n_cells_per_sample / total_cells_per_sample

        # Get eigengene from most enriched cell type module
        mod <- top_qval_mods_df$Module[i]
        ME_df <- fread(top_qval_mods_df$ME_path[i], data.table=FALSE)
        ME_vec <- ME_df[,grep(paste0("^", mod, "$"), colnames(ME_df))]

        df <- data.frame(Frac=frac_per_sample, ME=ME_vec)

        subtitle <- paste(
            top_qval_mods_df$Module[i], top_qval_mods_df$Network_short[i], "\n",
            "Cor:", round(cor(frac_per_sample, ME_vec), 2), "\n",
            "Qval =", formatC(top_qval_mods_df$Qval[i], format="e", digits=1) 
        )

        # Plot cell type proportion vs. module eigengene
        
        print(
            ggplot(df, aes(x=Frac*100, y=ME)) +
                geom_point() +
                theme(
                    plot.title=element_text(hjust=0.5),
                    plot.subtitle=element_text(hjust=0.5),
                    axis.title.x=element_text(size=14),
                    axis.title.y=element_text(size=14)
                ) +
                labs(
                    title=ctypes[i], 
                    subtitle=subtitle
                ) +
                xlab("% cells per sample") +
                ylab("Module eigengene")
        )
    }
    
}

plot_ctype_abundance_vs_top_corr_ME <- function(pseudobulk_legend, top_corr_mods_df) {
    options(repr.plot.width=6, repr.plot.height=6)

    top_corr_mods_df$Network_short <- gsub("Bicor-None_", "", gsub("_merge_ME_0.9_20151", "", top_corr_mods_df$Network))
 
    total_cells_per_sample <- colSums(pseudobulk_legend[,-c(1, 2)])
    ctypes <- top_corr_mods_df$Cell_type

    for (i in seq_along(ctypes)) {
        # Get cell type proportion in each sample
        mask <- pseudobulk_legend$Cell.type == ctypes[i]
        n_cells_per_sample <- colSums(pseudobulk_legend[mask, -c(1, 2)])
        frac_per_sample <- n_cells_per_sample / total_cells_per_sample

        # Get eigengene from most enriched cell type module
        mod <- top_corr_mods_df$Module[i]
        ME_df <- fread(top_corr_mods_df$ME_path[i], data.table=FALSE)
        ME_vec <- ME_df[,grep(paste0("^", mod, "$"), colnames(ME_df))]

        df <- data.frame(Frac=frac_per_sample, ME=ME_vec)

        subtitle <- paste(
            top_corr_mods_df$Module[i], top_corr_mods_df$Network_short[i], "\n",
            "Cor:", round(cor(frac_per_sample, ME_vec), 2), "\n",
            "Old cor:", round(top_corr_mods_df$Old_cor[i], 2) 
        )

        # Plot cell type proportion vs. module eigengene
        
        print(
            ggplot(df, aes(x=Frac*100, y=ME)) +
                geom_point() +
                theme(
                    plot.title=element_text(hjust=0.5),
                    plot.subtitle=element_text(hjust=0.5),
                    axis.title.x=element_text(size=14),
                    axis.title.y=element_text(size=14)
                ) +
                labs(
                    title=ctypes[i], 
                    subtitle=subtitle
                ) +
                xlab("% cells per sample") +
                ylab("Module eigengene")
        )
    }
    
}

plot_ctype_abundance <- function(pseudobulk_legend, top_qval_mods_df) {
    options(repr.plot.width=15, repr.plot.height=6)

    ctypes <- unique(pseudobulk_legend$Cell.type)

    for (i in seq_along(ctypes)) {
        # Get cell type proportion in each sample
        mask <- pseudobulk_legend$Cell.type == ctypes[i]
        n_cells_per_sample <- colSums(pseudobulk_legend[mask, -c(1, 2)])
        # frac_per_sample <- n_cells_per_sample/nrow(pseudobulk_legend) 

        n_cells_SD <- round(sd(n_cells_per_sample), 2)
        subtitle <- paste("SD =", n_cells_SD)

        # Plot cell type proportion vs. module eigengene
        df <- data.frame(Sample=1:length(n_cells_per_sample), No.cells=sort(n_cells_per_sample))

        options(repr.plot.width=15, repr.plot.height=6)

        print(
            ggplot(df, aes(x=Sample, y=No.cells)) +
                geom_point() +
                geom_line() +
                theme(
                    plot.title=element_text(hjust=0.5),
                    plot.subtitle=element_text(hjust=0.5),
                    axis.title.x=element_text(size=14),
                    axis.title.y=element_text(size=14)
                ) +
                labs(
                    title=ctypes[i], 
                    subtitle=subtitle
                ) +
                xlab("Sample") +
                ylab("# cells per sample")
        )
    }
    
}

plot_ME_corr_vs_stats <- function(cell_meta, top_corr_mods_df) {
    # Get original # cells per cell type
    top_corr_mods_df <- merge(
        top_corr_mods_df,
        data.frame(sort(table(cell_meta$cell_subclass))), 
        by.x="Cell_type", by.y="Var1"
    )
    corr_col <- grep("^Cor", colnames(top_corr_mods_df))
    top_corr_mods_df <- top_corr_mods_df[!is.na(top_corr_mods_df[,corr_col]),] 
    pal <- colorRampPalette(c("navy", "cyan", "yellow", "red"))(200)
    color_vec <- pal[as.numeric(cut(top_corr_mods_df[,corr_col], breaks=200))]
    plot(top_corr_mods_df$Freq, top_corr_mods_df$Pseudobulk_SD, col=color_vec, pch=16)
    plot(top_corr_mods_df$Freq, top_corr_mods_df[,corr_col])
}