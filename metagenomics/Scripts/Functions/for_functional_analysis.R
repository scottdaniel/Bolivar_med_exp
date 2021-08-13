#for sunbeam sbx_gene_cluster extension results

read_gene_aln_results <- function(base_dir, s_seq, pattern = "_1.txt") {

  tibble(FileName = list.files(
    base_dir, pattern=paste0("*", pattern))) %>%
    group_by(FileName) %>%
    do(read_diamond_file(file.path(base_dir, .$FileName), stringsAsFactors = F)) %>%
    ungroup() %>%
    select(-FileName) %>%

    right_join(select(s_seq, SampleID), by="SampleID") %>%
    filter(!is.na(geneID)) %>%
    pivot_wider(id_cols = SampleID,
                names_from = c(geneID,taxon),
                names_sep = "__",
                values_from = count,
                values_fill = list(count = 0)) %>%
    pivot_longer(-SampleID, names_to = c("geneID","taxon"), names_sep = "__", values_to = "count") %>%


    mutate(database = basename(base_dir))

}

read_gene_aln_results_filter_taxon <- function(base_dir, s_seq, desired_taxon) {

  tibble(FileName = list.files(
    base_dir, pattern="*_1.txt")) %>%
    group_by(FileName) %>%
    do(read.delim(file.path(base_dir, .$FileName), stringsAsFactors = F)) %>%
    ungroup() %>%
    mutate(SampleID = sub("_1.txt", "", FileName, perl=T)) %>%
    select(-FileName) %>%

    filter(taxon %in% desired_taxon) %>%
    right_join(select(s_seq, SampleID), by="SampleID") %>%

    pivot_wider(id_cols = SampleID,
                names_from = c(geneID,taxon),
                names_sep = "__",
                values_from = count,
                values_fill = list(count = 0)) %>%
    pivot_longer(-SampleID, names_to = c("geneID","taxon"), names_sep = "__", values_to = "count") %>%


    mutate(database = basename(base_dir))

}

read_gene_results_no_taxon_v1 <- function(base_dir, s_seq) {

  tibble(FileName = list.files(
    base_dir, pattern="*_1.txt")) %>%
    group_by(FileName) %>%
    do(read.delim(file.path(base_dir, .$FileName), stringsAsFactors = F)) %>%
    ungroup() %>%
    mutate(SampleID = sub("_1.txt", "", FileName, perl=T)) %>%
    select(-FileName) %>%

    group_by(geneID, SampleID) %>%
    summarize(sum_by_kegg_term = sum(count)) %>%
    right_join(select(s_seq, SampleID), by="SampleID") %>%
    filter(!is.na(geneID)) %>%
    pivot_wider(id_cols = SampleID,
                names_from = geneID,
                values_from = sum_by_kegg_term,
                values_fill = list(sum_by_kegg_term = 0)) %>%
    pivot_longer(-SampleID, names_to = "geneID", values_to = "sum_by_kegg_term") %>%
    mutate(database = basename(base_dir))

}

read_gene_results_no_taxon_v2 <- function(base_dir, s_seq, pattern = "_1.txt") {

  tibble(FileName = list.files(
    base_dir, pattern=paste0("*", pattern))) %>%
    group_by(FileName) %>%
    do(read_diamond_file(file.path(base_dir, .$FileName), stringsAsFactors = F)) %>%
    ungroup() %>%
    select(-FileName) %>%

    group_by(geneID, SampleID) %>%
    summarize(sum_by_kegg_term = sum(count)) %>%
    right_join(select(s_seq, SampleID), by="SampleID") %>%
    filter(!is.na(geneID)) %>%
    pivot_wider(id_cols = SampleID,
                names_from = geneID,
                values_from = sum_by_kegg_term,
                values_fill = list(sum_by_kegg_term = 0)) %>%
    pivot_longer(-SampleID, names_to = "geneID", values_to = "sum_by_kegg_term") %>%
    mutate(database = basename(base_dir))

}

#read a single file of diamond results
read_diamond_file <- function(path_to_file, pattern = "_1.txt", ...) {

  if (!str_detect(basename(path_to_file), pattern)) {
    warning(paste("The pattern", pattern, "did not match the file name!",
                        "Your SampleID's will be the same as the filenames"))
  }

  my_df <- read.delim(path_to_file, ...) %>%
    mutate(SampleID = str_replace(basename(path_to_file), pattern = pattern, ""))

  my_df

}

#for the old pathways.py script
read_ko_table <- function (filepath, sample_suffix=".fastq.") {
  ko <- read_tsv(filepath)
  colnames(ko) <- sub(sample_suffix, "", colnames(ko), fixed = TRUE)
  ko
}

run_lm_kegg <- function(props_toTest, s_toTest, form1, p_cutoff) {
  props_toTest[,s_toTest$SampleID] %>%
    as.data.frame() %>%
    rownames_to_column(var = "Pathway") %>%
    gather(key = "SampleID", value = "props", -Pathway) %>%
    mutate(props = props+1E-6) %>%
    merge(s_toTest, by="SampleID") %>%
    mutate(props_logit = log(props/(1-props))) %>%
    group_by(Pathway) %>%
    do(tidy_lm(lm(as.formula(form1), data=., na.action=na.omit))) %>%
    setNames(c("Pathway","term","Estimate","Std.Error","t.value","p.value")) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value < p_cutoff)
}

#fixing zeros properly from JungJin on a long formate proportion table
fix_zeros <- function(props_toTest_long) {
  props_toTest_long %>%
    mutate(props = props + min(filter(props_toTest_long, props>0)$props) / 10) %>%
    mutate(props_logit = log(props/(1-props)))
}

#when you already have in long format
run_lm_on_props <- function(factor_to_test, props_toTest, s_toTest, form1, p_cutoff) {

  props_toTest %>%
    left_join(s_toTest, by = "SampleID") %>%
    group_by(!!sym(factor_to_test)) %>%
    do(tidy_lm(lm(as.formula(form1), data=., na.action=na.omit))) %>%
    setNames(c(factor_to_test,"term","Estimate","Std.Error","t.value","p.value")) %>%
    ungroup() %>%
    filter(term != '(Intercept)') %>%
    group_by(term) %>%
    mutate(fdr = p.adjust(p.value, method="BH")) %>%
    ungroup() %>%
    filter(p.value < p_cutoff)
}

read_table_if_no_rds <- function(name_of_table, rds_file_name, base_dir, s_seq, pattern = "_1.txt") {

  if (!file.exists(here("metagenomics", "Data",rds_file_name))) {
    name_of_table <- read_gene_results_no_taxon_v2(base_dir, s_seq, pattern = "_1.txt")
    write_rds(name_of_table, here("metagenomics", "Data",rds_file_name), "xz", compression = 9L)
    return(name_of_table)
  }else{
    read_rds(here("metagenomics", "Data",rds_file_name))
  }

}

read_table_if_no_rds_taxon_included <- function(name_of_table, rds_file_name, base_dir, s_seq) {

  if (!file.exists(here("metagenomics", "Data",rds_file_name))) {
    name_of_table <- read_gene_aln_results(base_dir, s_seq)
    write_rds(name_of_table, here("metagenomics", "Data",rds_file_name), "xz", compression = 9L)
    return(name_of_table)
  }else{
    read_rds(here("metagenomics", "Data",rds_file_name))
  }

}

tidy_permanova <- function(anov){
  data.frame(Term = rownames(anov$aov.tab), anov$aov.tab, row.names = NULL) %>%
    rename(p.value = Pr..F.)
}
permanova_test <- function(dist_matrix, s_toTest, form1, perm, strata=NULL){
  set.seed(42)
  if (!grepl("~", form1)) {
    form1 <- paste0("dist_matrix ~ ", form1)
  }
  dist_matrix <- dist_subset(dist_matrix, s_toTest$SampleID)
  form1 <- as.formula(form1)
  if(is.null(strata)) {
    tidy_permanova(adonis(form1, data=s_toTest, permutations=perm))
  } else {
    tidy_permanova(adonis(form1, data=s_toTest, permutations=perm, strata=s_toTest[,strata]))
  }
}

permanova_posthoc <- function(dist_matrix, s_toTest, form1, perm, strata=NULL, group_label, p_cutoff=0.05){
  if (!grepl("~", form1)) {
    form1 <- paste0("dist_matrix ~ ", form1)
  }
  a_ixn <- permanova_test(dist_matrix, s_toTest, form1, perm, strata) %>%
    mutate(comparison = "all") %>%
    mutate(group1 = NA) %>%
    mutate(group2 = NA)

  combs <- combn(as.character(unique(s_toTest[[group_label]])), 2)
  num_tests <- dim(combs)[2]

  # do post hoc tests
  if (filter(a_ixn, Term == group_label)$p.value < p_cutoff) {
    for (i in 1:num_tests){
      s_temp <- filter(s_toTest, .data[[group_label]] %in% combs[,i])
      #dist_toTest = dist_subset(dist_matrix, s_temp$SampleID)
      a_ixn <- rbind(a_ixn,
                     permanova_test(dist_matrix, s_temp, form1, perm, strata) %>%
                       mutate(comparison = paste(combs[,i], collapse=' - ')) %>%
                       mutate(group1 = combs[,i][1]) %>%
                       mutate(group2 = combs[,i][2])
      )
    }
  }
  a_ixn
}

