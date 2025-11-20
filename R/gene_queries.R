# ============================================================================
# Gene Query Functions - New Naming Convention (snake_case)
# ============================================================================

#' Query Comparative Genomics Data for a Gene
#'
#' Retrieves comparative genomics data showing evolutionary conservation and
#' ortholog information for a specified gene across different species.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol (e.g., "THBS2").
#'
#' @return A data frame containing comparative genomics data including ortholog
#'   relationships, conservation scores, and cross-species functional annotations.
#'
#' @examples
#' \dontrun{
#' comp_data <- get_gene_comparative_genomics("THBS2")
#' comp_data <- get_gene_comparative_genomics("TP53")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen compGenomicsQuery
#' @export
get_gene_comparative_genomics <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::compGenomicsQuery(ensembl_id)
  return(res)
}

#' Query DepMap Essentiality Data for a Gene
#'
#' Retrieves DepMap (Dependency Map) essentiality data showing which cell lines
#' are dependent on the specified gene for survival.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#'
#' @return A data frame containing DepMap essentiality scores across different
#'   cell lines.
#'
#' @examples
#' \dontrun{
#' depmap_data <- get_gene_depmap("THBS2")
#' depmap_data <- get_gene_depmap("MYC")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen depMapQuery
#' @export
get_gene_depmap <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::depMapQuery(ensembl_id)
  return(res)
}

#' Query Gene Ontology Data for a Gene
#'
#' Retrieves Gene Ontology (GO) annotations for a specified gene, including
#' biological processes, molecular functions, and cellular components.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#'
#' @return A data frame containing GO annotations including GO terms, evidence codes,
#'   and functional classifications.
#'
#' @examples
#' \dontrun{
#' go_data <- get_gene_ontology("THBS2")
#' go_data <- get_gene_ontology("TP53")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen geneOntologyQuery
#' @export
get_gene_ontology <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::geneOntologyQuery(ensembl_id)
  return(as.data.frame(res))
}

#' Query Genetic Constraint Data for a Gene
#'
#' Retrieves genetic constraint data showing how tolerant a gene is to
#' loss-of-function variants.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#'
#' @return A data frame containing constraint metrics including pLI and LOEUF scores.
#'
#' @examples
#' \dontrun{
#' constraint_data <- get_gene_constraint("THBS2")
#' constraint_data <- get_gene_constraint("MYC")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen geneticConstraintQuery
#' @export
get_gene_constraint <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::geneticConstraintQuery(ensembl_id)
  return(as.data.frame(res))
}

#' Query Cancer Hallmarks Data for a Gene
#'
#' Retrieves cancer hallmarks data showing which hallmarks of cancer are
#' associated with the specified gene.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#'
#' @return A data frame containing cancer hallmarks annotations.
#'
#' @examples
#' \dontrun{
#' hallmarks_data <- get_gene_hallmarks("THBS2")
#' hallmarks_data <- get_gene_hallmarks("MYC")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen hallmarksQuery
#' @export
get_gene_hallmarks <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::hallmarksQuery(ensembl_id)
  return(as.data.frame(res))
}

#' Query Protein Interactions Data for a Gene
#'
#' Retrieves protein-protein interaction data for the specified gene from
#' various interaction databases.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param source_database Character string. Source database for interactions
#'   (e.g., "intact", "string"). Default: "intact".
#' @param index Integer. Page index for pagination (default: 0).
#' @param size Integer. Number of records to retrieve (default: 10).
#'
#' @return A data frame containing protein interaction data.
#'
#' @examples
#' \dontrun{
#' interactions <- get_gene_interactions("THBS2")
#' interactions <- get_gene_interactions("THBS2", source_database = "string", size = 50)
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen interactionsQuery
#' @export
get_gene_interactions <- function(gene_id, source_database = "intact", index = 0, size = 10) {
  validate_gene_id(gene_id)
  if (!is.character(source_database) || length(source_database) != 1) {
    stop("'source_database' must be a single character string", call. = FALSE)
  }
  validate_nonnegative_integer(index, "index")
  validate_positive_integer(size, "size")

  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::interactionsQuery(
    ensembl_id,
    sourceDatabase = source_database,
    index = index,
    size = size
  )
  return(as.data.frame(res))
}

#' Query Known Drugs Data for a Gene
#'
#' Retrieves information about drugs that target the specified gene.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param size Integer. Number of records to retrieve (default: 10).
#' @param cursor Character or NULL. Cursor for pagination (default: NULL).
#'
#' @return A data frame containing drug information.
#'
#' @examples
#' \dontrun{
#' drugs_data <- get_gene_drugs("THBS2")
#' drugs_data <- get_gene_drugs("EGFR", size = 50)
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen knownDrugsGeneQuery
#' @export
get_gene_drugs <- function(gene_id, size = 10, cursor = NULL) {
  validate_gene_id(gene_id)
  validate_positive_integer(size, "size")

  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::knownDrugsGeneQuery(ensembl_id, cursor = cursor, size = size)
  return(as.data.frame(res))
}

#' Query Mouse Phenotypes Data for a Gene
#'
#' Retrieves mouse phenotype data for orthologs of the specified human gene.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#'
#' @return A data frame containing mouse phenotype data.
#'
#' @examples
#' \dontrun{
#' mouse_data <- get_gene_mouse_phenotypes("THBS2")
#' mouse_data <- get_gene_mouse_phenotypes("HOX11")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen mousePhenotypesQuery
#' @export
get_gene_mouse_phenotypes <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::mousePhenotypesQuery(ensembl_id)
  return(as.data.frame(res))
}

#' Query Pathways Data for a Gene
#'
#' Retrieves biological pathway information for the specified gene.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#'
#' @return A data frame containing pathway annotations.
#'
#' @examples
#' \dontrun{
#' pathways_data <- get_gene_pathways("THBS2")
#' pathways_data <- get_gene_pathways("AKT1")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen pathwaysQuery
#' @export
get_gene_pathways <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::pathwaysQuery(ensembl_id)
  return(as.data.frame(res))
}

#' Query Pharmacogenomics Data for a Gene
#'
#' Retrieves pharmacogenomics data showing how genetic variants in the specified
#' gene affect drug response.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#'
#' @return A data frame containing pharmacogenomics associations.
#'
#' @examples
#' \dontrun{
#' pgx_data <- get_gene_pharmacogenomics("THBS2")
#' pgx_data <- get_gene_pharmacogenomics("CYP2D6")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen pharmacogenomicsGeneQuery
#' @export
get_gene_pharmacogenomics <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::pharmacogenomicsGeneQuery(ensembl_id)
  return(as.data.frame(res))
}

#' Query Safety Liabilities Data for a Gene
#'
#' Retrieves safety liabilities data indicating potential adverse effects
#' associated with targeting the specified gene.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#'
#' @return A data frame containing safety information.
#'
#' @examples
#' \dontrun{
#' safety_data <- get_gene_safety("THBS2")
#' safety_data <- get_gene_safety("EGFR")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen safetyQuery
#' @export
get_gene_safety <- function(gene_id) {
  validate_gene_id(gene_id)
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::safetyQuery(ensembl_id)
  return(as.data.frame(res))
}
