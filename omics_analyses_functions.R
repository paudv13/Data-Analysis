library(stringr)
library(phyloseq)
library(ggpubr)
library(ggplot2)
library(vegan)
library(knitr)

# add pseudocount to 0 values for log transform (pseudocount is equal to half the minimum non-zero value on a per-sample basis)
add_pseudocount <- function(count_table, margin = 2) {# add pseudocount to 0 values for log transform (pseudocount is equal to half the minimum non-zero value on a per-sample basis)
  samplemins <- apply(count_table, FUN = function(x){min(x[x !=0])/2}, MARGIN = margin)
  
  for (j in 1:ncol(count_table)) {
    currentrow <- count_table[,j]
    currentrow[currentrow == 0] <- samplemins[j]
    count_table[, j] <- currentrow
  }
  
  return(count_table)
}

metaphlan4_to_phyloseq <- function(metadata, bug_list, tree = T, abund = 0.001, prev = 0.1, duplicated_rownames_error = F){
  
  #Set the prevalence and the abundance to filter
  filter_prevalence_abundance_dataframe_no_loop <- function(unfilter_dataframe, abundance = abund, prevalence = prev){
    rows_to_keep=c()
    number_of_samples=dim(unfilter_dataframe)[2]
    colsums_vector = colSums(unfilter_dataframe)
    
    for (i in 1:dim(unfilter_dataframe)[1]) {
      row_vector = unfilter_dataframe[i,]
      relabun_row_vector = row_vector/colsums_vector
      # num_over_abundance = sum((relabun_row_vector > abundance) == TRUE) # old way that doesn't work anymore
      num_over_abundance = as.numeric(sum((relabun_row_vector > abundance), na.rm = TRUE))
      
      if (num_over_abundance/number_of_samples > prevalence) {
        rows_to_keep = c(rows_to_keep, i)
      }
    }
    filtered_dataframe <- unfilter_dataframe[rows_to_keep,]
    return(list(filtered_dataframe, rows_to_keep))
  }
  
  
  if (tree == T) {
    ######Isolate species from buglist
    just_species <- bug_list[grepl("t__",rownames(bug_list)),]
    
    just_species[is.na(just_species)] <- 0
    just_species <- just_species[!colSums(just_species) == 0] # remove samples with 0 total abundance
    
    
    
    ###### Filter by abundance and prevalence
    
    filtered <- filter_prevalence_abundance_dataframe_no_loop(just_species, abundance = abund, prevalence = prev)
    just_species_filtered <- filtered[[1]] #apply to species df
    just_species_filtered <- just_species_filtered[!colSums(just_species_filtered) == 0] # remove samples with 0 total abundance after filtering
    
    
    rm(filtered)
    
    #####Create taxa table
    
    # rm(taxa_table)
    for (i in row.names(just_species_filtered)) {
      tax_list <- str_split(i,"\\|", simplify = T)
      kingdom <- str_split(tax_list[1], pattern = "k__", simplify = T)[2]
      phylum<- str_split(tax_list[2], pattern = "p__", simplify = T)[2]
      class<- str_split(tax_list[3], pattern = "c__", simplify = T)[2]
      order<- str_split(tax_list[4], pattern = "o__", simplify = T)[2]
      family<- str_split(tax_list[5], pattern = "f__", simplify = T)[2]
      genus <- str_split(tax_list[6], pattern = "g__", simplify = T)[2]
      specie <- str_split(tax_list[7], pattern = "s__", simplify = T)[2]
      SGBID <- str_split(tax_list[8], pattern = "SGB", simplify = T)[2]
      
      if (is.na(SGBID)) {
        SGBID <- str_split(tax_list[8], pattern = "EUK", simplify = T)[2]
      }
      
      rm(tax_list)
      
      if (!exists("taxa_table")){
        taxa_table <- data.frame(kingdom, phylum, class, order, family, genus, specie, SGBID, stringsAsFactors = F)
      }
      else if (exists("taxa_table")){
        taxa_table <- rbind(taxa_table, c(kingdom, phylum, class, order, family, genus, specie, SGBID))
      }
    }
    
    row.names(taxa_table) <- row.names(just_species_filtered)
    
    taxa_table$SGBID <- str_replace_all(taxa_table$SGBID, pattern = "_group", replacement = "")
    
    
    ####Import metaphlan tree
      tree <- read_tree('/mnt/synology/GERARD/DATA/mpa_vJan21_CHOCOPhlAnSGB_202103.nwk')
      just_species_filtered <- just_species_filtered[rownames(taxa_table)[taxa_table$SGBID %in% tree$tip.label],]
      taxa_table <- taxa_table[rownames(taxa_table)[taxa_table$SGBID %in% tree$tip.label],]

      if (duplicated_rownames_error) {
        #Fix the only species with duplicated SGBID
        taxa_table["k__Bacteria|p__Firmicutes|c__Negativicutes|o__Acidaminococcales|f__Acidaminococcaceae|g__Acidaminococcus|s__Acidaminococcus_timonensis|t__SGB5741", 8] <- "SGB5741"
        taxa_table["k__Eukaryota|p__Eukaryota_unclassified|c__Eukaryota_unclassified|o__Eukaryota_unclassified|f__Hexamitidae|g__Giardia|s__Giardia_intestinalis|t__EUK5741", 8] <- "EUK5741"
        
      }
      #continue with the normal script
      rownames(taxa_table) <- taxa_table$SGBID
      
      rownames(just_species_filtered) <- taxa_table$SGBID
      
      tree <- phy_tree(tree)
      
      taxtab <- tax_table(taxa_table)
      taxa_names(taxtab) <- rownames(taxa_table)
      phy_obj <- phyloseq(otu_table(just_species_filtered, taxa_are_rows = T), sample_data(metadata), taxtab, tree)
      
      colnames(phy_obj@tax_table) <- c("kingdom", "phylum", "class", "order", "family", "genus", "specie", "SGBID") #set taxa_table colnames to be self-explanatory
      
    
  } else {
    ######Isolate species from buglist
    just_species <- bug_list[grepl("s__",rownames(bug_list)),]
    just_species <- just_species[!grepl("t__",rownames(just_species)),]
    
    just_species[is.na(just_species)] <- 0
    just_species <- just_species[!colSums(just_species) == 0] # remove samples with 0 total abundance
    
    
    
    ###### Filter by abundance and prevalence
    
    filtered <- filter_prevalence_abundance_dataframe_no_loop(just_species, abundance = abund, prevalence = prev)
    just_species_filtered <- filtered[[1]] #apply to species df
    just_species_filtered <- just_species_filtered[!colSums(just_species_filtered) == 0] # remove samples with 0 total abundance after filtering
    
    
    rm(filtered)
    
    #####Create taxa table
    
    # rm(taxa_table)
    for (i in row.names(just_species_filtered)) {
      tax_list <- str_split(i,"\\|", simplify = T)
      kingdom <- str_split(tax_list[1], pattern = "k__", simplify = T)[2]
      phylum<- str_split(tax_list[2], pattern = "p__", simplify = T)[2]
      class<- str_split(tax_list[3], pattern = "c__", simplify = T)[2]
      order<- str_split(tax_list[4], pattern = "o__", simplify = T)[2]
      family<- str_split(tax_list[5], pattern = "f__", simplify = T)[2]
      genus <- str_split(tax_list[6], pattern = "g__", simplify = T)[2]
      specie <- str_split(tax_list[7], pattern = "s__", simplify = T)[2]
      
      rm(tax_list)
      
      if (!exists("taxa_table")){
        taxa_table <- data.frame(kingdom, phylum, class, order, family, genus, specie, stringsAsFactors = F)
      }
      else if (exists("taxa_table")){
        taxa_table <- rbind(taxa_table, c(kingdom, phylum, class, order, family, genus, specie))
      }
    }
    
    row.names(taxa_table) <- row.names(just_species_filtered)
    
    #####Create phyloseq_object
    taxtab <- tax_table(taxa_table)
    taxa_names(taxtab) <- rownames(taxa_table)
    phy_obj <- phyloseq(otu_table(just_species_filtered, taxa_are_rows = T), sample_data(metadata), taxtab)
    
    colnames(phy_obj@tax_table) <- c("kingdom", "phylum", "class", "order", "family", "genus", "specie") #set taxa_table colnames to be self-explanatory
    
  }
  
  return(phy_obj)
  
}

humann3_to_phyloseq <- function(metadata, bug_list, tree = T, abund = 0.001, prev = 0.1){
  
  ######Isolate species from buglist
  to_look = "s__"
  
  for (i in row.names(bug_list)) {
    if (grepl(to_look,i)) {
      if (!exists("just_species")){
        just_species <- data.frame(bug_list[i,])
      } 
      else if (exists("just_species")){
        just_species <- rbind(just_species, bug_list[i,])
      }
    }
  }
  just_species[is.na(just_species)] <- 0
  just_species <- just_species[!colSums(just_species) == 0] # remove samples with 0 total abundance
  
  
  
  ###### Filter by abundance and prevalence
  
  #Set the prevalence and the abundance to filter
  filter_prevalence_abundance_dataframe <- function(unfilter_dataframe, abundance = abund, prevalence = prev){
    rows_to_keep=c()
    number_of_samples=dim(unfilter_dataframe)[2]
    colsums_vector = colSums(unfilter_dataframe)
    for (i in 1:dim(unfilter_dataframe)[1]) {
      row_vector = unfilter_dataframe[i,]
      relabun_row_vector = row_vector/colsums_vector
      # num_over_abundance = sum((relabun_row_vector > abundance) == TRUE) # old way that doesn't work anymore
      num_over_abundance = as.numeric(sum((relabun_row_vector > abundance), na.rm = TRUE))
      
      if (num_over_abundance/number_of_samples > prevalence) {
        rows_to_keep = c(rows_to_keep, i)
      }
    }
    filtered_dataframe <- unfilter_dataframe[rows_to_keep,]
    return(list(filtered_dataframe, rows_to_keep))
  }
  
  filtered <- filter_prevalence_abundance_dataframe(just_species)
  just_species_filtered <- filtered[[1]] #apply to species df
  just_species_filtered <- just_species_filtered[!colSums(just_species_filtered) == 0] # remove samples with 0 total abundance after filtering 
  
  
  rm(filtered)
  
  #####Create taxa table
  
  # rm(taxa_table)
  for (i in row.names(just_species_filtered)) {
    tax_list <- str_split(i,"\\|", simplify = T)
    kingdom <- str_split(tax_list[1], pattern = "k__", simplify = T)[2]
    phylum<- str_split(tax_list[2], pattern = "p__", simplify = T)[2]
    class<- str_split(tax_list[3], pattern = "c__", simplify = T)[2]
    order<- str_split(tax_list[4], pattern = "o__", simplify = T)[2]
    family<- str_split(tax_list[5], pattern = "f__", simplify = T)[2]
    genus <- str_split(tax_list[6], pattern = "g__", simplify = T)[2]
    specie <- str_split(tax_list[7], pattern = "s__", simplify = T)[2]
    
    rm(tax_list)
    
    if (!exists("taxa_table")){
      taxa_table <- data.frame(kingdom, phylum, class, order, family, genus, specie, stringsAsFactors = F)
    } 
    else if (exists("taxa_table")){
      taxa_table <- rbind(taxa_table, c(kingdom, phylum, class, order, family, genus, specie))
    }
  }
  
  row.names(taxa_table) <- row.names(just_species_filtered)
  
  ####Import metaphlan tree
  
  if (tree == T) {
  
    tree <- read_tree('/mnt/synology/GERARD/DATA/mpa_v30_CHOCOPhlAn_201901_species_tree.nwk')
    
    tree$tip.label <- sub("GCA_[0-9]+\\|", "", tree$tip.label)
    tree <- phy_tree(tree)
    
    taxtab <- tax_table(taxa_table)
    taxa_names(taxtab) <- rownames(taxa_table)
    phy_obj <- phyloseq(otu_table(just_species_filtered, taxa_are_rows = T), sample_data(metadata), taxtab, tree)
  
    
  } else {
  
  #####Create phyloseq_object
  taxtab <- tax_table(taxa_table)
  taxa_names(taxtab) <- rownames(taxa_table)
  phy_obj <- phyloseq(otu_table(just_species_filtered, taxa_are_rows = T), sample_data(metadata), taxtab)
  }
  
  colnames(phy_obj@tax_table) <- c("kingdom", "phylum", "class", "order", "family", "genus", "specie") #set taxa_table colnames to be self-explanatory
  
  return(phy_obj)
}

beta_diversity <- function(phy_obj, colorby, shapeby = NULL){ #colorby is the name of the metddata variable for which we will colour the plots
  phylo_ord <- ordinate(phy_obj, method = "PCoA", distance = 'bray') 
  p1 <- plot_ordination(phy_obj, phylo_ord, color = colorby, shape = shapeby )
  p1 <- p1 + theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    # geom_text(aes(label=colnames(otu_table(phy_obj)), hjust= 0, vjust = 2)) +
    ggtitle("Beta Diversity - Bray-Curtis Distance") +
    theme(plot.title = element_text(hjust = 0.5))
  
  d <- UniFrac(phy_obj, weighted = F)
  phylo_ord <- ordinate(phy_obj, method = "PCoA", distance = d) 
  p2 <- plot_ordination(phy_obj, phylo_ord, color = colorby, shape = shapeby )
  p2 <- p2 + theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    # geom_text(aes(label=colnames(otu_table(phy_obj)), hjust= 0, vjust = 2)) +
    ggtitle("Beta Diversity - UnWeighted UniFrac") +
    theme(plot.title = element_text(hjust = 0.5))
  
  d <- UniFrac(phy_obj, weighted = T)
  phylo_ord <- ordinate(phy_obj, method = "PCoA", distance = d) 
  p3 <- plot_ordination(phy_obj, phylo_ord, color = colorby, shape = shapeby )
  p3 <- p3 + theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    # geom_text(aes(label=colnames(otu_table(phy_obj)), hjust= 0, vjust = 2)) +
    ggtitle("Beta Diversity - Weighted UniFrac") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list( bray = p1, unweighted_UniFrac = p2, weighted_UniFrac = p3))
}

beta_bray <- function(phy_obj, colorby, shapeby = NULL){ #colorby is the name of the metddata variable for which we will colour the plots
  phylo_ord <- ordinate(phy_obj, method = "PCoA", distance = 'bray') 
  p1 <- plot_ordination(phy_obj, phylo_ord, color = colorby, shape = shapeby )
  p1 <- p1 + theme_classic() +
    scale_color_brewer(palette = "Dark2") +
    # geom_text(aes(label=colnames(otu_table(phy_obj)), hjust= 0, vjust = 2)) +
    ggtitle("Beta Diversity - Bray-Curtis Distance") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p1)
}

adonis_from_physeq <- function(phy_obj, variables, distance = "bray", parallel = "1", bywhat = "margin"){#variables <- string of variables names list sepparated by "," or by "+"
  #distance options = "bray" for Bray-Curtis, "wu" for weighted UniFrac, "uu" for unweighted UniFrac
  # bywhat = "margin" for marginal analysis, where order of factors does not matter, "terms" when order of factors matters
  trasposed_otutable <- as.data.frame(t(otu_table(phy_obj)))
  metadata <- data.frame(sample_data(phy_obj))
  
  variables <- str_replace_all(variables, ",", "+")
  
  if (distance == "bray") {
    formulaa <- as.formula(paste("trasposed_otutable", variables, sep = " ~ "))
  } else if (distance == "wu"){
    d <- UniFrac(phy_obj, weighted = T)
    formulaa <- as.formula(paste("d", variables, sep = " ~ "))
  } else if (distance == "uu"){
    d <- UniFrac(phy_obj, weighted = F)
    formulaa <- as.formula(paste("d", variables, sep = " ~ "))
  }
  
  set.seed(123)
  adonis_res <- adonis2(formulaa, data = metadata, parallel = parallel, by = bywhat, na.action = na.omit)
  
  adonis_res$variable = rownames(adonis_res)
  
  return(adonis_res)
}

simpson_from_physeq <- function(phy_obj, variable, palette = "jco"){ #Name of the variable to study
  simpson <- diversity(otu_table(phy_obj), index = "simpson", MARGIN = 2) #Margin = 2 for samples Margin = 1 for OTUs. diversity function from vegan
  diversity_data <- data.frame(simpson)
  metadata <- data.frame(sample_data(phy_obj))
  metadata$sam_name <- row.names(metadata)
  # metadata$Diet <- as.factor(metadata$Diet)
  diversity_data$sam_name <- row.names(diversity_data)
  diversity_df <- merge(diversity_data,metadata, by = "sam_name")
  colnames(diversity_df)[2] <- "simpson"
  diversity_df[,variable] <- as.factor(diversity_df[,variable])
  
  psimpsoncohort <- ggboxplot(diversity_df, 
                              x = variable, 
                              y = "simpson",
                              color = variable,
                              # color = "dodgerblue2",
                              # fill = "grey",
                              palette = palette,
                              # palette = c("tomato2", "dodgerblue2", "orange2"),
                              add = "jitter") +
    # ) +
    rremove("xlab") +
    theme(text = element_text(family = "sans"))
  
  diversity_df[,variable] <- as.factor(diversity_df[,variable])
  lev <- levels(diversity_df[,variable]) # get the variables
  
  
  # make a pairwise list that we want to compare.
  L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
  
  
  
  p2simpsoncohort <- psimpsoncohort + stat_compare_means(
    comparisons = L.pairs,
    method = "wilcox.test",
    # label = "p.signif",
    # symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), symbols = c("****", "***", "**", "*", "ns")),
  )
  return(list(index = diversity_df[,c(1,2)], plot = p2simpsoncohort, plot2 = psimpsoncohort))
}


shannon_from_physeq <- function(phy_obj, variable, palette = "jco", add = "jitter"){ #Name of the variable to study
  Shannon <- diversity(otu_table(phy_obj), index = "shannon", MARGIN = 2) #Margin = 2 for samples Margin = 1 for OTUs. diversity function from vegan
  diversity_data <- data.frame(Shannon)
  metadata <- data.frame(sample_data(phy_obj))
  metadata$sam_name <- row.names(metadata)
  # metadata$Diet <- as.factor(metadata$Diet)
  diversity_data$sam_name <- row.names(diversity_data)
  diversity_df <- merge(diversity_data,metadata, by = "sam_name")
  diversity_df[,variable] <- as.factor(diversity_df[,variable])
  
  pshannoncohort <- ggboxplot(diversity_df, 
                              x = variable, 
                              y = "Shannon",
                              color = variable,
                              # color = "dodgerblue2",
                              # fill = "grey",
                              palette = palette,
                              # palette = c("tomato2", "dodgerblue2", "orange2"),
                              add = add,
                              add.params = list(size = 0.5)) +
                              # ) +
    rremove("xlab") +
    theme(text = element_text(family = "sans"))
  
  diversity_df[,variable] <- as.factor(diversity_df[,variable])
  lev <- levels(diversity_df[,variable]) # get the variables
  
  
  # make a pairwise list that we want to compare.
  L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
 
  
  
  p2shannoncohort <- pshannoncohort + stat_compare_means(
    comparisons = L.pairs,
    method = "wilcox.test",
    # label = "p.signif",
    # symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1), symbols = c("****", "***", "**", "*", "ns")),
  )
  return(list(index = diversity_df[,c(1,2)], plot = p2shannoncohort, plot2 = pshannoncohort))
}


chao1_from_physeq_relab <- function(phy_obj, variable, palette = "jco", add = "jitter"){ #Name of the variable to study
  # chao1_df <- estimate_richness(otu_table(phy_obj)*100000, measures = c("Chao1")) #doesnt work in all the cases
  integer_otu_table <- otu_table(data.frame(lapply(as.data.frame(otu_table(phy_obj))*100000, as.integer), row.names = rownames(otu_table(phy_obj))), taxa_are_rows = T)
  chao1_df <- estimate_richness(integer_otu_table, measures = c("Chao1"))
  
  chao1_df$sam_name <- row.names(chao1_df)
  
  diversity_df <- cbind(chao1_df, data.frame(sample_data(phy_obj)))
  diversity_df[,variable] <- as.factor(diversity_df[,variable])
  
  pchao <- ggboxplot(diversity_df, 
                     x = variable, 
                     y = "Chao1",
                     color = variable,
                     # color = "dodgerblue2",
                     # fill = "grey",
                     palette = palette,
                     # palette = c("tomato2", "dodgerblue2", "orange2"),
                     add = add,
                     add.params = list(size = 0.5)) +
                    # ) +
    rremove("xlab") +
    theme(text = element_text(family = "sans"))
  
  # ggtitle("Chao1 index depending on disease")
  
  
  
  diversity_df[,variable] <- as.factor(diversity_df[,variable])
  lev <- levels(diversity_df[,variable]) # get the variables
  
  # make a pairwise list that we want to compare.
  L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
  
  # pval <- list(
  #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  #   symbols = c("****", "***", "**", "*", "n.s")
  # )
  
  p2chao <- pchao + stat_compare_means(
    comparisons = L.pairs,
    method = "wilcox.test",
    # label = "p.signif",
    # symnum.args = list(
    #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    #   symbols = c("", "", "", "", "")
    #   # symbols = c("****", "***", "**", "*", "n.s")
    # )
  )
  
  return(list(index = diversity_df[,c(1,2)], plot = p2chao, plot2 = pchao))
}

chao1_from_physeq_abs <- function(phy_obj, variable, palette = "jco", add = "jitter"){ #Name of the variable to study
  # chao1_df <- estimate_richness(otu_table(phy_obj)*100000, measures = c("Chao1")) #doesnt work in all the cases
  integer_otu_table <- otu_table(phy_obj)
  chao1_df <- estimate_richness(integer_otu_table, measures = c("Chao1"))
  
  chao1_df$sam_name <- row.names(chao1_df)
  
  diversity_df <- cbind(chao1_df, data.frame(sample_data(phy_obj)))
  diversity_df[,variable] <- as.factor(diversity_df[,variable])
  
  pchao <- ggboxplot(diversity_df, 
                     x = variable, 
                     y = "Chao1",
                     color = variable,
                     # color = "dodgerblue2",
                     # fill = "grey",
                     palette = palette,
                     # palette = c("tomato2", "dodgerblue2", "orange2"),
                     add = add,
                     add.params = list(size = 0.5)) +
    # ) +
    rremove("xlab") +
    theme(text = element_text(family = "sans"))
  
  # ggtitle("Chao1 index depending on disease")
  
  
  
  diversity_df[,variable] <- as.factor(diversity_df[,variable])
  lev <- levels(diversity_df[,variable]) # get the variables
  
  # make a pairwise list that we want to compare.
  L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
  
  # pval <- list(
  #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  #   symbols = c("****", "***", "**", "*", "n.s")
  # )
  
  p2chao <- pchao + stat_compare_means(
    comparisons = L.pairs,
    method = "wilcox.test",
    # label = "p.signif",
    # symnum.args = list(
    #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    #   symbols = c("", "", "", "", "")
    #   # symbols = c("****", "***", "**", "*", "n.s")
    # )
  )
  
  return(list(index = diversity_df[,c(1,2)], plot = p2chao, plot2 = pchao))
}

dysbiosis_score <- function(phy_obj, phy_controls, group_var, weighted_unifrac = T, addhealthies = F, boxplo = "*"){

  
  reference <- UniFrac(phy_controls, weighted = weighted_unifrac)
  
  # reference <- vegdist(t(healthies_abundance), method = "bray")
  reference <- as.data.frame(as.matrix(reference))
  reference_medians <- as.data.frame(apply(reference, 2, median))
  colnames(reference_medians) <- "Dysbiosis_score"
  reference_medians$cohort <- "Healthy"
  
  reflist <- colnames(reference)
  
  group_dists <- UniFrac(merge_phyloseq(phy_obj, phy_controls), weighted = weighted_unifrac)
  
  group_dists <- as.data.frame(as.matrix(group_dists))
  
  group_dists <- group_dists[,! colnames(group_dists) %in% reflist]
  group_dists <- group_dists[rownames(group_dists) %in% reflist,]
  
  group_medians <- as.data.frame(apply(group_dists, 2, median))
  colnames(group_medians) <- "Dysbiosis_score"
  group_medians[rownames(group_medians), "cohort"] <- sample_data(phy_obj)[rownames(group_medians), group_var] %>% c()
  group_medians <- group_medians[order(group_medians$cohort),]
  
  densplot <- ggplot(rbind(reference_medians, group_medians), aes(Dysbiosis_score, color = cohort, fill = cohort)) + 
    geom_density(alpha = 0.5) + 
    theme_classic()  +
    xlim(c(0.2, 0.8)) +
    ylab("Density") +
    xlab("Dysbiosis Score")
  
  
  if (addhealthies == T) {
    group_medians <- rbind(reference_medians, group_medians)
  }
  
  
  boxpl <- ggboxplot(group_medians, 
                     x = "cohort", 
                     y = "Dysbiosis_score",
                     color = "cohort",
                     # fill = "grey", 
                     # palette = "Dark2",
                     add = "jitter") +
    rremove("xlab") +
    theme(text = element_text(family = "sans"))
  
  
  group_medians[,"cohort"] <- as.factor(group_medians[,"cohort"])
  lev <- levels(group_medians[,"cohort"]) # get the variables
  
  # make a pairwise list that we want to compare.
  L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
  
  
  if (boxplo == "none"){
    return(return(list(values = group_medians, density_plot = densplot, boxplot = boxpl)))
  } else if (boxplo == "pval"){
    boxpl2 <- boxpl + stat_compare_means(
      comparisons = L.pairs,
      method = "wilcox.test",
      # label = "p.signif",
      # symnum.args = list(
      #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1),
      #   # symbols = c("", "", "", "", "")
      #   symbols = c("****", "***", "**", "*")
      # ), hide.ns = T
    )
  } else if (boxplo == "*") {
    boxpl2 <- boxpl + stat_compare_means(
      comparisons = L.pairs,
      method = "wilcox.test",
      label = "p.signif",
      symnum.args = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.1),
        # symbols = c("", "", "", "", "")
        symbols = c("****", "***", "**", "*")
      ), hide.ns = T
    )
  }
  
  return(return(list(values = group_medians, density_plot = densplot, boxplot = boxpl2)))
  
}

correlation_from_physeq <- function(phy_obj, var1, var2){ # var1 and var2 is a string containing the name of the metadata variables
  # to test correlation for
  
  #this function first performs a shapiro test to see if data distribution is normal, and then chooses spearman or pearson test
  #depending on the result
  
  library(ggpubr)
  
  phy_metadata <- data.frame(sample_data(phy_obj))
  phy_metadata <- phy_metadata[!is.na(phy_metadata[,var1]) & !is.na(phy_metadata[,var2]),]
  phy_metadata[,var1]<- as.numeric(phy_metadata[,var1])
  phy_metadata[,var2] <- as.numeric(phy_metadata[,var2])
  
  shapiro1 <- shapiro.test(phy_metadata[,var1])$p.value
  shapiro2 <- shapiro.test(phy_metadata[,var2])$p.value
  # phy_metadata[,c("hbi", "BMI", "disease", "patientID")] %>% View
  if (shapiro1 <= 0.05 | shapiro2 <= 0.05) {
    correlation <- cor.test(phy_metadata[,var1], phy_metadata[,var2], method = "spearman")
    
    correlation_plot <- ggscatter(phy_metadata, x = var1, y = var2,
                            add = "reg.line", conf.int = T,
                            cor.coef = T, cor.method = "spearman") + 
                  ggtitle(paste("Spearman's correlation plot between", var1, "and", var2, sep = " "))
  } else {
    correlation <- cor.test(phy_metadata[,var1], phy_metadata[,var2], method = "pearson")
    correlation_plot <- ggscatter(phy_metadata, x = var1, y = var2,
                            add = "reg.line", conf.int = T,
                            cor.coef = T, cor.method = "pearson") +
      ggtitle(paste("Pearson's correlation plot between", var1, "and", var2, sep = " "))
  }    
  
  return(list(correlation_test = correlation, correlation_plot = correlation_plot))
}

corrplot_from_physeq <- function(phy_obj, interest_metadata, is_taxa = T, corrmethod = "spearman", qval = 0.05, adjust_method = "BH") { #interest_metadata is a vector with the colnamaes of the metadata variables to correlate
  library(Hmisc)
  library(corrplot)
  current_metadata <- data.frame(sample_data(phy_obj))
  current_metadata <- current_metadata[,interest_metadata]
  current_metadata <- as.matrix(current_metadata)
  current_species <- as.matrix(data.frame(otu_table(phy_obj)))
  
  a <- rcorr(current_metadata, t(current_species), type = corrmethod)
  a$r[lower.tri(a$r)] <- NA
  a$P[lower.tri(a$P)] <- NA
  a$n[lower.tri(a$n)] <- NA
  
  melted_pvals <- reshape2::melt(a$P)
  melted_pvals <- melted_pvals[! is.na(melted_pvals$value),]
  melted_pvals <- melted_pvals[melted_pvals$Var1 %in% interest_metadata,]
  melted_pvals <- melted_pvals[!(melted_pvals$Var2 %in% interest_metadata),]
  melted_pvals$FDR <- p.adjust(melted_pvals$value, method = adjust_method)
  
  melted_Rsquared <- reshape2::melt(a$r)
  melted_Rsquared <- melted_Rsquared[! is.na(melted_Rsquared$value),]
  melted_Rsquared <- melted_Rsquared[melted_Rsquared$Var1 %in% interest_metadata,]
  melted_Rsquared <- melted_Rsquared[!(melted_Rsquared$Var2 %in% interest_metadata),]
  melted_pvals$R <- melted_Rsquared$value
  #create a column with R coefficient only if correlation is significant
  melted_pvals$r_if_sig <- melted_pvals$R
  melted_pvals$r_if_sig[melted_pvals$FDR > qval] <- NA
  
  
  rm(melted_Rsquared)
  
  #Subset only species that have at least one significant correlation in variable1
  significantvar1 <- unique(melted_pvals$Var1[melted_pvals$FDR <= qval]) %>% as.vector
  to_plot <- melted_pvals[melted_pvals$Var1 %in% significantvar1,]
  
  #Subset only species that have at least one significant correlation in variable2
  significantvar2 <- unique(melted_pvals$Var2[melted_pvals$FDR <= qval]) %>% as.vector
  to_plot <- to_plot[to_plot$Var2 %in% significantvar2,]
  
  if (is_taxa) {
    #Replace SGB ID for taxonomy name in Var2
    taxonomy <- tax_table(phy_obj) %>% data.frame #create auxiliar taxonomy df
    to_plot$Var2 <- taxonomy[to_plot$Var2 %>% as.vector() ,7]
  }
  
  #Plot the corr heatmap
  plot_result <- ggplot(data = to_plot, aes(x = Var1, y = Var2, fill = R, label=round(r_if_sig,2))) +
    geom_tile() +
    labs(x = NULL, y = NULL, fill = paste(corrmethod,"correlation", sep = "\n")) +
    # map a red, white and blue color scale to correspond to -1:1 sequential gradient 
    scale_fill_gradient2(mid="#FBFEF9",high="#0C6291",low="#A63446", limits=c(-1,1)) +
    geom_text() +
    theme_classic() +
    # remove excess space on x and y axes
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    # change global font to roboto
    theme(text=element_text(family="Arial"))
  
  
  return(plot_result)
}

boxplot_from_physeq <- function(phy_obj, varx, vary, palette = "jco"){
  phy_metadata <- data.frame(sample_data(phy_obj))
  phy_metadata <- phy_metadata[!is.na(phy_metadata[,vary]) & !is.na(phy_metadata[,varx]),]
  phy_metadata[,vary]<- as.numeric(phy_metadata[,vary])
  # phy_metadata[,var2] <- as.numeric(phy_metadata[,var2])
  
  shapiro1 <- shapiro.test(phy_metadata[,vary])$p.value
  # shapiro2 <- shapiro.test(phy_metadata[,var2])$p.value
  # phy_metadata[,c("hbi", "BMI", "disease", "patientID")] %>% View
  if (shapiro1 <= 0.05) {
    method = "wilcox.test"
  } else {
    method = "t.test"
  }    
  plot1 <- ggboxplot(phy_metadata, 
                     x = varx, 
                     y = vary,
                     color = varx,
                     # fill = "grey", 
                     palette = palette,
                     add = "jitter") +
    rremove("xlab") +
    theme(text = element_text(family = "sans"))+
    ggtitle(paste(method, "'s boxplot between", vary, " and ", varx, sep = ""))
  
  # ggtitle("Chao1 index depending on disease")
  
  
  
  phy_metadata[,varx] <- as.factor(phy_metadata[,varx])
  lev <- levels(phy_metadata[,varx]) # get the variables
  
  # make a pairwise list that we want to compare.
  L.pairs <- combn(seq_along(lev), 2, simplify = FALSE, FUN = function(i) lev[i])
  
  # pval <- list(
  #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
  #   symbols = c("****", "***", "**", "*", "n.s")
  # )
  
  plot2 <- plot1 + stat_compare_means(
    comparisons = L.pairs,
    method = method,
    # label = "p.signif",
    # symnum.args = list(
    #   cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    #   symbols = c("", "", "", "", "")
    #   # symbols = c("****", "***", "**", "*", "n.s")
    # )
  )
  return(plot2)
}

plot_taxa_abundances <- function(phy_obj, facet_condition, taxonomic_rank = "phylum", plotlegend = T){
  #extract data from physeq object
  phy_subset <- tax_glom(phy_obj, taxrank = taxonomic_rank)
  current_abundance <- data.frame(otu_table(phy_subset))
  current_metadata <- data.frame(sample_data(phy_subset))
  
  #rename taxa
  rownames(current_abundance) <- tax_table(phy_subset)[rownames(current_abundance), taxonomic_rank]
  
  #convert to long df format for plotting purposes
  library(reshape2)
  current_abundance <- data.frame(t(current_abundance))
  current_abundance$sample <- rownames(current_abundance)
  current_abundance$condition <- current_metadata[rownames(current_abundance), facet_condition]
  current_abundance_long <- melt(current_abundance, id.vars = c("sample", "condition"))
  colnames(current_abundance_long) <- c("sample", "condition", "taxa", "abundance")
  
  #set the palette to use depending on number of taxa
  library(RColorBrewer)
  getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
  final_palette <- getPalette(length(unique(current_abundance_long$taxa)))
  
  #plot
  p <- ggplot(current_abundance_long, aes(x  = sample, fill = taxa, y = abundance)) +
    geom_bar(stat = "identity", colour = "black") +
    theme_classic() +
    # theme(axis.text.x = element_text(angle = 90),
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          # axis.line.x = element_blank(),
          axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
          # legend.text = element_text(size = 12, face = "bold", colour = "black"), 
          axis.text.y = element_text(colour = "black", size = 8, face = "bold")) + 
    # scale_y_continuous(expand = c(0,0)) + 
    labs(x = "", y = "Relative Abundance (%)", fill = "Taxa") +
    scale_fill_manual(values = final_palette) +
    scale_y_continuous(expand = c(0,0)) +
    facet_grid(~condition, scales = "free", space = "free")
  
  if (plotlegend == F) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

filter_prevalence_abundance_dataframe <- function(unfilter_dataframe, abundance = 0.001, prevalence = 0.1){
  rows_to_keep=c()
  number_of_samples=dim(unfilter_dataframe)[2]
  colsums_vector = colSums(unfilter_dataframe)
  for (i in 1:dim(unfilter_dataframe)[1]) {
    row_vector = unfilter_dataframe[i,]
    relabun_row_vector = row_vector/colsums_vector
    num_over_abundance = sum((relabun_row_vector > abundance) == TRUE)
    if (num_over_abundance/number_of_samples > prevalence) {
      rows_to_keep = c(rows_to_keep, i)
    }
  }
  filtered_dataframe <- unfilter_dataframe[rows_to_keep,]
  return(list(filtered_dataframe, rows_to_keep))
}

filter_prevalence_abundance_dataframe_no_loop <- function(unfilter_dataframe, abundance = 0.001, prevalence = 0.1){ #samples in columns and species in rows
  unfiltered_data <- t(unfilter_dataframe)
  
  total_samples <- nrow(unfiltered_data)
  min_samples <- total_samples * prevalence
  
  
  # Filter by abundance using zero as value for NAs
  data_zeros <- unfiltered_data
  data_zeros[is.na(data_zeros)] <- 0
  
  colsums_vector = colSums(data_zeros)
  abundance <- colsums_vector * abundance
  abundance <- replicate(abundance, n = total_samples) %>% as.matrix() %>% t()
  
  filtered_data <-
    unfiltered_data[, 
                    colSums(data_zeros > abundance) > min_samples,
                    drop = FALSE]
  total_filtered_features <-
    ncol(unfiltered_data) - ncol(filtered_data)
  filtered_feature_names <-
    setdiff(colnames(unfiltered_data), colnames(filtered_data))
  
  return(list(t(filtered_data), filtered_feature_names))
}

rcorr_padjust <- function(x, method = "BH", ...) {
  stopifnot(class(x) == "rcorr")
  x$P[upper.tri(x$P)] <- p.adjust(x$P[upper.tri(x$P)], method = method)
  x$P[lower.tri(x$P)] <- p.adjust(x$P[lower.tri(x$P)], method = method)
  return(x)}

biplot <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                    obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                    ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                    alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                    varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, loadings.n = 10,
                    ...) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  library(ggrepel)
  
  pcloadings <- pcobj$rotation[,1:2]
  pcloadings <- rowSums(abs(pcloadings))
  pcloadings <- as.data.frame(pcloadings[order(pcloadings, decreasing = T)])
  
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), 
               FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v[rownames(pcloadings)[1:loadings.n],], aes(x = 0, y = 0, 
                                                                                xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                                                       "picas")), color = muted("red"))
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v[rownames(pcloadings)[1:loadings.n],], aes(label = varname, 
                                                                             x = xvar, y = yvar, angle = angle, hjust = hjust), 
                       color = "darkred", size = varname.size)
  }
  return(g)
}

biplot2 <- function (plot.data, loadings.data = NULL, colour = NULL, size = NULL, 
                     linetype = NULL, alpha = NULL, fill = NULL, shape = NULL, 
                     label = FALSE, label.label = "rownames", label.colour = colour, 
                     label.alpha = NULL, label.size = NULL, label.angle = NULL, 
                     label.family = NULL, label.fontface = NULL, label.lineheight = NULL, 
                     label.hjust = NULL, label.vjust = NULL, label.repel = FALSE, 
                     label.position = "identity", loadings = FALSE, loadings.arrow = grid::arrow(length = grid::unit(8, 
                                                                                                                     "points")), loadings.colour = "#FF0000", loadings.label = FALSE, 
                     loadings.label.label = "rownames", loadings.label.colour = "#FF0000", 
                     loadings.label.alpha = NULL, loadings.label.size = NULL, 
                     loadings.label.angle = NULL, loadings.label.family = NULL, 
                     loadings.label.fontface = NULL, loadings.label.lineheight = NULL, 
                     loadings.label.hjust = NULL, loadings.label.vjust = NULL, 
                     loadings.label.repel = FALSE, label.show.legend = NA, frame = FALSE, 
                     frame.type = NULL, frame.colour = colour, frame.level = 0.95, 
                     frame.alpha = 0.2, xlim = c(NA, NA), ylim = c(NA, NA), log = "", 
                     main = NULL, xlab = NULL, ylab = NULL, asp = NULL, loadings.n = length(plot.data), pcobj = NULL, ...) 
{
  pcloadings <- pcobj$rotation[,1:2]
  pcloadings <- rowSums(abs(pcloadings))
  pcloadings <- as.data.frame(pcloadings[order(pcloadings, decreasing = T)])
  
  
  plot.columns <- colnames(plot.data)
  mapping <- ggplot2::aes_string(x = plot.columns[1L], y = plot.columns[2L])
  if (is.logical(shape) && !shape && missing(label)) {
    label <- TRUE
  }
  p <- ggplot2::ggplot(data = plot.data, mapping = mapping)
  if (!is.logical(shape) || shape) {
    p <- p + ggfortify:::geom_factory(ggplot2::geom_point, plot.data, 
                                      colour = colour, size = size, linetype = linetype, 
                                      alpha = alpha, fill = fill, shape = shape)
  }
  p <- ggfortify:::plot_label(p = p, data = plot.data, label = label, 
                              label.label = label.label, label.colour = label.colour, 
                              label.alpha = label.alpha, label.size = label.size, 
                              label.angle = label.angle, label.family = label.family, 
                              label.fontface = label.fontface, label.lineheight = label.lineheight, 
                              label.hjust = label.hjust, label.vjust = label.vjust, 
                              label.repel = label.repel, label.show.legend = label.show.legend, 
                              label.position = label.position)
  if (loadings.label && !loadings) {
    loadings <- TRUE
  }
  if (loadings && !is.null(loadings.data)) {
    scaler <- min(max(abs(plot.data[, 1L]))/max(abs(loadings.data[, 
                                                                  1L])), max(abs(plot.data[, 2L]))/max(abs(loadings.data[, 
                                                                                                                         2L])))
    loadings.columns <- colnames(loadings.data)
    loadings.mapping <- ggplot2::aes_string(x = 0, y = 0, 
                                            xend = loadings.columns[1L], yend = loadings.columns[2L])
    loadings.data[, 1L:2L] <- loadings.data[, 1L:2L] * scaler * 
      0.8
    p <- p + geom_segment(data = loadings.data[rownames(pcloadings)[1:loadings.n],], mapping = loadings.mapping, 
                          arrow = loadings.arrow, colour = loadings.colour)
    p <- ggfortify:::plot_label(p = p, data = loadings.data[rownames(pcloadings)[1:loadings.n],], label = loadings.label, 
                                label.label = loadings.label.label, label.colour = loadings.label.colour, 
                                label.alpha = loadings.label.alpha, label.size = loadings.label.size, 
                                label.angle = loadings.label.angle, label.family = loadings.label.family, 
                                label.fontface = loadings.label.fontface, label.lineheight = loadings.label.lineheight, 
                                label.hjust = loadings.label.hjust, label.vjust = loadings.label.vjust, 
                                label.repel = loadings.label.repel, label.show.legend = label.show.legend, 
                                label.position = label.position)
  }
  if (missing(frame) && !is.null(frame.type)) {
    frame <- TRUE
  }
  . <- NULL
  if (frame) {
    if (is.null(frame.type) || frame.type == "convex") {
      if (is.null(frame.colour) || !(frame.colour %in% 
                                     colnames(plot.data))) {
        hulls <- plot.data[grDevices::chull(plot.data[, 
                                                      1L:2L]), ]
      }
      else {
        hulls <- plot.data %>% dplyr::group_by(.data[[frame.colour]]) %>% 
          dplyr::do(.[grDevices::chull(.[, 1L:2L]), 
          ])
      }
      mapping <- aes_string(colour = frame.colour, fill = frame.colour)
      p <- p + ggplot2::geom_polygon(data = hulls, mapping = mapping, 
                                     alpha = frame.alpha)
    }
    else if (frame.type %in% c("t", "norm", "euclid")) {
      mapping <- aes_string(colour = frame.colour, fill = frame.colour)
      p <- p + ggplot2::stat_ellipse(mapping = mapping, 
                                     level = frame.level, type = frame.type, geom = "polygon", 
                                     alpha = frame.alpha)
    }
    else {
      stop("frame.type must be convex, t, norm or euclid")
    }
  }
  p <- ggfortify:::post_autoplot(p = p, xlim = xlim, ylim = ylim, log = log, 
                                 main = main, xlab = xlab, ylab = ylab, asp = asp)
  return(p)
}



autoplot_pca <- function (object, data = NULL, scale = 1, x = 1, y = 2, variance_percentage = TRUE, 
                          ...) 
{
  plot.data <- ggplot2::fortify(object, data = data)
  plot.data$rownames <- rownames(plot.data)
  if (ggfortify:::is_derived_from(object, "prcomp")) {
    ve <- object$sdev^2/sum(object$sdev^2)
    PC <- paste0("PC", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    loadings.column <- "rotation"
    lam <- object$sdev[c(x, y)]
    lam <- lam * sqrt(nrow(plot.data))
  }
  else if (ggfortify:::is_derived_from(object, "princomp")) {
    ve <- object$sdev^2/sum(object$sdev^2)
    PC <- paste0("Comp.", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    loadings.column <- "loadings"
    lam <- object$sdev[c(x, y)]
    lam <- lam * sqrt(nrow(plot.data))
  }
  else if (ggfortify:::is_derived_from(object, "factanal")) {
    if (is.null(attr(object, "covariance"))) {
      p <- nrow(object$loading)
      ve <- colSums(object$loading^2)/p
    }
    else ve <- NULL
    PC <- paste0("Factor", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    scale <- 0
    loadings.column <- "loadings"
  }
  else if (ggfortify:::is_derived_from(object, "lfda")) {
    ve <- NULL
    PC <- paste0("PC", c(x, y))
    x.column <- PC[1]
    y.column <- PC[2]
    scale <- 0
    loadings.column <- NULL
  }
  else {
    stop(paste0("Unsupported class for autoplot.pca_common: ", 
                class(object)))
  }
  if (scale != 0) {
    lam <- lam^scale
    plot.data[, c(x.column, y.column)] <- t(t(plot.data[, 
                                                        c(x.column, y.column)])/lam)
  }
  plot.columns <- unique(c(x.column, y.column, colnames(plot.data)))
  plot.data <- plot.data[, plot.columns]
  if (!is.null(loadings.column)) {
    loadings.data <- as.data.frame(object[[loadings.column]][, 
    ])
    loadings.data$rownames <- rownames(loadings.data)
    loadings.columns <- unique(c(x.column, y.column, colnames(loadings.data)))
    loadings.data <- loadings.data[, loadings.columns]
  }
  else {
    loadings.data <- NULL
  }
  if (is.null(ve) | !variance_percentage) {
    labs <- PC
  }
  else {
    ve <- ve[c(x, y)]
    labs <- paste0(PC, " (", round(ve * 100, 2), "%)")
  }
  xlab <- labs[1]
  ylab <- labs[2]
  p <- biplot2(plot.data = plot.data, loadings.data = loadings.data, 
               xlab = xlab, ylab = ylab, pcobj = object, ...)
  return(p)
}