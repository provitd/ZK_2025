# Aim: Define function relative to pathway analysis

library(msigdbr)
library(clusterProfiler)
library(dplyr)

get_terms_msigdbr <- function(x, species = "Mus musculus", gene = "ensembl_gene"){
  db <- msigdbr(species = species,
                category = ifelse(is.na(x$cat), NA, x$cat),
                subcategory = ifelse(is.na(x$subcat), NA, x$subcat )) %>%
    dplyr::select(gene, x$id, x$description)
  dict <- list()
  dict[['TERM2GENE']] <- db %>% dplyr::select(x$id, gene) %>%
    unique() %>% rename("TERM" = x$id, "GENE" = gene)
  if (x$id != x$description){
    dict[['TERM2NAME']] <- db %>% dplyr::select(x$id, x$description) %>%
      unique() %>% rename("TERM" = x$id, "DESCRIPTION" = x$description)
  } else {
    dict[['TERM2NAME']] <- NA
  }
  return(dict)
}


terms <- list()
terms[["Hallmark"]] = list(
  name = "Hallmark",
  cat = "H",
  id = "gs_name",
  description = "gs_name"
)
terms[["REACTOME"]] <- list(
  name = "REACTOME",
  cat = "C2",
  subcat = "CP:REACTOME",
  id = "gs_exact_source",
  description = "gs_name"
)
terms[["KEGG"]] <- list(
  name = "KEGG",
  cat = "C2",
  subcat = "CP:KEGG",
  id = "gs_exact_source",
  description = "gs_name"
)
terms[["GO:BP"]] <- list(
  name = "GO:BP",
  cat = "C5",
  subcat = "GO:BP",
  id = "gs_exact_source",
  description = "gs_name"
)
terms[["GO:MF"]] <- list(
  name = "GO:MF",
  cat = "C5",
  subcat = "GO:MF",
  id = "gs_exact_source",
  description = "gs_name"
)
terms[["GO:CC"]] <- list(
  name = "GO:CC",
  cat = "C5",
  subcat = "GO:CC",
  id = "gs_exact_source",
  description = "gs_name"
)
terms[["All"]] <- list(
  name = "All",
  cat = NA,
  subcat = NA,
  id = "gs_name",
  description = "gs_exact_source"
)


# Compute terms
terms <- lapply(terms, function(x) {
  terms[[x$name]]$dict <- get_terms_msigdbr(x, gene = "gene_symbol")
  return(terms[[x$name]])
})




pipeline_enrichment <- function(opt, species = "Mus musculus") {
  write(paste0('Running enrichment analysis ', opt$name, '.'), stderr())
  
  TERM2GENE = opt$term$dict[["TERM2GENE"]]
  TERM2NAME = opt$term$dict[["TERM2NAME"]]
  
  if (opt$method == 'ORA') {
    x <- do.call(clusterProfiler::enricher,
                 modifyList(opt$enricher, list(TERM2GENE = TERM2GENE,
                                               TERM2NAME = TERM2NAME)
                 ))
  } else { # GSEA
    x <- do.call(clusterProfiler::GSEA,
                 modifyList(opt$GSEA, list(TERM2GENE = TERM2GENE,
                                           TERM2NAME = TERM2NAME)
                 ))
  }
  
  is.noenrich <- function(res) {
    if (!is.null(res)){
      if (dim(res@result)[1] != 0){return(FALSE)}
      else {return(TRUE)}
    } else {return(TRUE)}
  }
  
  if (is.noenrich(x)){
    warning(paste('No enrichment result found for', opt$name, 'comparison.'))
    return(x)
  }
  
  x@organism = species
  if (opt$method == 'ORA') {
    x@ontology = opt$term$name
  } else {
    x@setType = opt$term$name
  }
  
  # dim(x) outputs the dim of x@result
  # filtered by the p.adjust according to pvalueCutoff
  if (dim(x)[1] > 0){
    # Calculate pairwise termsim
    optPT <- opt$pairwise_termsim
    optPT$x <- x
    x <- do.call(enrichplot::pairwise_termsim, optPT)
  }
  
  return(x)
}

