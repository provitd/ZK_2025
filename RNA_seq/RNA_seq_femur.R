# Aim: DESeq analysis of ZK vs OIL in Femur
#      followed by GSEA analysis

library(glue)
library(here)
source(here("utils/fct_pathway.R"))

sample_type <- "fem"
dir <- "RNA_seq/"

data <- read.csv(here(glue(dir, "data/GSE302870_S20193_alldata.tsv")),
                 row.names = 1, sep='\t')
metadata <- structure(list(Code = c("MNGC394", "MNGC396", "MNGC431", "MNGC432",
                                    "MNGC433", "MNGC434", "MNGC435"),
                           Name = c("s9_fem_109", "s11_fem_122", "s10_fem_116_new",
                                   "s12_fem_135_new", "s13_fem_82_new",
                                   "s14_fem_83_new", "s15_fem_90_new"),
                           condition = c("oil", "oil", "oil", "oil",
                                    "zk", "zk", "zk")),
                      row.names = c(NA, 7L), class = "data.frame")
counts <- data[metadata$Name]
gene_mapping <- data["Gene.name"]

# DESeq analysis

dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = metadata,
                                      design = ~ condition)
dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)
vsd <- DESeq2::vst(dds)
DESeq2::plotPCA(vsd)
DESeq2::plotDispEsts(dds)

deseq_obj <- DESeq2::results(dds, alpha = 0.05,
                     cooksCutoff = T, independentFiltering = T,
                     contrast = c("condition","zk","oil"))

DESeq2::plotMA(deseq_obj)

deseq_res <- deseq_obj %>% as.data.frame() %>%
  base::merge(., gene_mapping, by.x = "row.names", by.y = "row.names",
              all.x = TRUE, all.y = FALSE) %>%
  tibble::column_to_rownames("Row.names") %>%
  base::merge(., counts, by = "row.names", all.x = TRUE) %>%
  dplyr::rename("Ensembl Gene ID" = "Row.names")

dir.create(here(dir, "table"))
write.table(deseq_res, here(dir, "table", glue("{sample_type}_deseq.tsv")),
            sep='\t', quote = FALSE, na = '', dec = ',', row.names = FALSE)

# Remove duplicated genes so that we can perform a pathway analysis using Symbol
gene_to_remove <- deseq_res$Symbol[duplicated(deseq_res$Symbol)]
dim(deseq_res)
deseq_res <- deseq_res %>% dplyr::filter(!Symbol %in% gene_to_remove)
dim(deseq_res)

# Pathway analysis

# Compute ordered geneList
geneList <- deseq_res %>%
  tibble::column_to_rownames("Symbol") %>%
  dplyr::mutate(metric = log2FoldChange) %>%
  dplyr::arrange(desc(metric)) %>%
  dplyr::select(metric) %>%
  tibble::rownames_to_column() %>%
  tibble::deframe()

# Compute enrichment
OPT <- list()
OPT$comparison <- list()
for (term in terms) { # terms is defined in script/fct_pathway.R file
  OPT$comparison[[glue("GSEA {term$name}")]] <- list(
    name = glue("GSEA {term$name}"),
    term = terms[[term$name]],
    method = "GSEA",
    GSEA = list(geneList = geneList),
    pairwise_termsim = list(method = "JC")
  )
}

OPT$comparison <- lapply(OPT$comparison, function(x){
  x$gsea_res <- pipeline_enrichment(x) # pipeline_enrichment is defined in script/fct_pathway.R file
  return(x)
})

# Save results as tsv
dir.create(here(dir, "table/pathway"))
lapply(OPT$comparison, function(x){
  file <- gsub('\\W', '_', x$name, perl=TRUE)
  directory <- here(dir, "table/pathway")
  file <- here(directory, glue("{sample_type}_{file}.tsv"))
  write.table(x$gsea_res, file,
              sep='\t', quote = FALSE, na = '', dec = ',', row.names = FALSE)
})
