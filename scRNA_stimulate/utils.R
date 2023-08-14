

# Colors for PBMCs
pbmc_colors = c("B" = "#66C2A5", 
              "DC" = "#FC8D62",
              "HSC" = "#8DA0CB",
              "MK" = "#E78AC3", 
              "Mono_CD14" = "#A6D854",
              "Mono_CD16" = "#f2ec72",
              "NK" = "#62AAEA", 
              "T_CD4" = "#D1C656",
              "T_CD8" = "#968763")

# Colors for pancreas
celltype.colors = c('alpha'="#ed2bb1",
                    'beta'="#239eb3",
                    'gamma'="#d1bfec",
                    'delta'= "#FF6347",
                    'stellate'="#11e38c",
                    'immune'="#812050",
                    'ductal'="#b2d27a",
                    'endothelial'="#4e2da6",
                    'acinar'="#f6bb86",
                    'schwann'="#115d52",
                    'epsilon'="#a1def0",
                    'mast'="#8fec2f")

# Colors for pancreas query donors (Baron et al., 2016)
querydonor.colors = c('human1' = '#b9dbf0',
                      'human2' = '#77a1ba',
                      'human3' = '#6c7ca8',
                      'human4' = '#364261',
                      'mouse1' = '#e68c8c',
                      'mouse2' = '#b35757')

simple.colors = c(   'CD141+DC'='#f2bd80',
                    'Central memory T cells '='#1d6d1f',
                    'Circulating NK '='#8c3ba0',
                    'Conventional dendritic cells(CD1C DC)'='#6533ed',
                    'Cytotoxicity CD8T'='#83e3f0',
                    'Effector memory T cells'='#fd5917',
                    'Exhausted CD8+ T (Tex) cells'='#4f8c9d',
                    'ILCs'='#eb1fcb',
                    'Liver-resident NK (lrNK) cell '='#f5cdaf',
                    'Lymphoid-B'='#9698dc',
                    'M1'='#20f53d',
                    'Mast'='#f3d426',
                    'Mono'='#f6932e',
                    'Myeloid-derived suppressor cells'='#caf243',
                    'NK'='#38b5fc',
                    'TAM-like'='#c82565',
                    'Th0'='#d6061a',
                    'Th1'='#e36f6f',
                    'Treg'='#1dfee1')

# Custom ordering to match original author publication ordering of states
group.ordering = c("HSC_MPP", "Pre pro B cell", 'pro-B cell', 'pre-B cell', 'B cell',
            'ILC precursor', 'Early lymphoid/T', 'NK', 'Neut-myeloid prog.',
            'pDC precursor','DC precursor', 'DC1', 'DC2', 'Monocyte precursor', 'Monocyte', 
            'Mono-Mac', 'Kupffer Cell', 'VCAM1+ EI macro.', 'MEMP', 'Mast cell',
            'Megakaryocyte', 'Early Erythroid', 'Mid Erythroid', 'Late Erythroid',
            'Endothelial cell', 'Fibroblast', 'Hepatocyte') 


#' Basic function to plot cells, colored and faceted by metadata variables
#' 
#' @param metadata metadata, with UMAP labels in UMAP1 and UMAP2 slots
#' @param title Plot title
#' @param color.by metadata column name for phenotype labels
#' @param facet.by metadata column name for faceting
#' @param color.mapping custom color mapping
#' @param show.legend Show cell type legend

plotBasic = function(umap_labels,                # metadata, with UMAP labels in UMAP1 and UMAP2 slots
                        title = 'Query',         # Plot title
                        color.by = 'cell_type',  # metadata column name for coloring
                        facet.by = NULL,         # (optional) metadata column name for faceting
                        color.mapping = NULL,    # custom color mapping
                        legend.position = 'right') {  # Show cell type legend
    
    p = umap_labels %>%
            dplyr::sample_frac(1L) %>% # permute rows randomly
            ggplot(aes(x = UMAP1, y = UMAP2)) + 
            geom_point_rast(aes(col = get(color.by)), size = 0.3, stroke = 0.2, shape = 16)
        if (!is.null(color.mapping)) { p = p + scale_color_manual(values = color.mapping) }
    
    # Default formatting
    p = p + theme_bw() +
            labs(title = title, color = color.by) + 
            theme(plot.title = element_text(hjust = 0.5)) +
            theme(legend.position=legend.position) +
            theme(legend.text = element_text(size=8), legend.title=element_text(size=12)) + 
            guides(colour = guide_legend(override.aes = list(size = 4))) + guides(alpha = 'none')

    if(!is.null(facet.by)) {
        p = p + facet_wrap(~get(facet.by)) +
                theme(strip.text.x = element_text(size = 12)) }    
    return(p)
}
