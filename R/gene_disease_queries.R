# ============================================================================
# Gene-Disease Association Query Functions - New Naming Convention (snake_case)
# ============================================================================

#' Query ChEMBL Evidence for Gene-Disease Associations
#'
#' Retrieves ChEMBL evidence data for associations between a specified gene and disease.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param disease Character string. A disease name that will be resolved to EFO/MONDO ID.
#' @param size Integer. Number of records to retrieve (default: 10).
#' @param cursor Character or NULL. Cursor for pagination (default: NULL).
#'
#' @return A data frame containing ChEMBL evidence data.
#'
#' @examples
#' \dontrun{
#' chembl_data <- get_assoc_chembl("PSEN1", "Alzheimer disease")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen chemblQuery
#' @export
get_assoc_chembl <- function(gene_id, disease, size = 10, cursor = NULL) {
  validate_gene_id(gene_id)
  validate_disease_name(disease)
  validate_positive_integer(size, "size")

  disease_id <- get_trait_id(disease, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::chemblQuery(ensembl_id, disease_id, cursor = cursor, size = size)
  return(as.data.frame(res))
}

#' Query Europe PMC Evidence for Gene-Disease Associations
#'
#' Retrieves Europe PMC literature evidence for associations between gene and disease.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param disease Character string. A disease name.
#' @param size Integer. Number of records to retrieve (default: 50).
#' @param cursor Character or NULL. Cursor for pagination (default: NULL).
#'
#' @return A data frame containing literature evidence.
#'
#' @examples
#' \dontrun{
#' lit_data <- get_assoc_literature("PSEN1", "Alzheimer disease")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen europePMCQuery
#' @export
get_assoc_literature <- function(gene_id, disease, size = 50, cursor = NULL) {
  validate_gene_id(gene_id)
  validate_disease_name(disease)
  validate_positive_integer(size, "size")

  disease_id <- get_trait_id(disease, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::europePMCQuery(ensembl_id, disease_id, cursor = cursor, size = size)
  return(as.data.frame(res))
}

#' Query Gene Burden Evidence for Gene-Disease Associations
#'
#' Retrieves gene burden evidence showing how genetic variants contribute to disease risk.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param disease Character string. A disease name.
#' @param size Integer. Number of records to retrieve (default: 3500).
#'
#' @return A data frame containing gene burden evidence.
#'
#' @examples
#' \dontrun{
#' burden_data <- get_assoc_gene_burden("SORL1", "Alzheimer disease")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen geneBurdenQuery
#' @export
get_assoc_gene_burden <- function(gene_id, disease, size = 3500) {
  validate_gene_id(gene_id)
  validate_disease_name(disease)
  validate_positive_integer(size, "size")

  disease_id <- get_trait_id(disease, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::geneBurdenQuery(ensembl_id, disease_id, size = size)
  return(as.data.frame(res))
}

#' Query Genomics England Evidence for Gene-Disease Associations
#'
#' Retrieves evidence from Genomics England PanelApp.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param disease Character string. A disease name.
#' @param size Integer. Number of records to retrieve (default: 3500).
#'
#' @return A data frame containing Genomics England evidence.
#'
#' @examples
#' \dontrun{
#' ge_data <- get_assoc_genomics_england("SORL1", "Alzheimer disease")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen genomicsEnglandQuery
#' @export
get_assoc_genomics_england <- function(gene_id, disease, size = 3500) {
  validate_gene_id(gene_id)
  validate_disease_name(disease)
  validate_positive_integer(size, "size")

  disease_id <- get_trait_id(disease, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::genomicsEnglandQuery(ensembl_id, disease_id, size = size)
  return(as.data.frame(res))
}

#' Query GWAS Credible Sets for Target-Disease Associations
#'
#' Retrieves GWAS credible sets data showing genetic variants likely causal for associations.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param disease Character string. A disease name.
#' @param size Integer. Number of records to retrieve (default: 500).
#'
#' @return A data frame containing GWAS credible sets data.
#'
#' @examples
#' \dontrun{
#' gwas_data <- get_assoc_gwas_credible_sets("SORL1", "Alzheimer disease")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen gwasCredibleSetsQuery
#' @export
get_assoc_gwas_credible_sets <- function(gene_id, disease, size = 500) {
  validate_gene_id(gene_id)
  validate_disease_name(disease)
  validate_positive_integer(size, "size")

  disease_id <- get_trait_id(disease, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::gwasCredibleSetsQuery(ensembl_id, disease_id, size = size)
  return(as.data.frame(res))
}

#' Query Orphanet Evidence for Gene-Disease Associations
#'
#' Retrieves evidence from Orphanet rare disease database.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param disease Character string. A disease name.
#' @param size Integer. Number of records to retrieve (default: 3500).
#'
#' @return A data frame containing Orphanet evidence.
#'
#' @examples
#' \dontrun{
#' orphan_data <- get_assoc_orphanet("THBS2", "Alzheimer disease")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen orphanetQuery
#' @export
get_assoc_orphanet <- function(gene_id, disease, size = 3500) {
  validate_gene_id(gene_id)
  validate_disease_name(disease)
  validate_positive_integer(size, "size")

  disease_id <- get_trait_id(disease, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::orphanetQuery(ensembl_id, disease_id, size = size)
  return(as.data.frame(res))
}

#' Query UniProt Literature Evidence for Gene-Disease Associations
#'
#' Retrieves UniProt literature evidence linking genes to diseases.
#'
#' @param gene_id Character string. A valid human HGNC gene symbol.
#' @param disease Character string. A disease name.
#' @param size Integer. Number of records to retrieve (default: 3500).
#'
#' @return A data frame containing UniProt literature evidence.
#'
#' @examples
#' \dontrun{
#' uniprot_data <- get_assoc_uniprot("APOE", "Alzheimer disease")
#' }
#'
#' @importFrom BioEnricher GeneMapper
#' @importFrom otargen uniprotLiteratureQuery
#' @export
get_assoc_uniprot <- function(gene_id, disease, size = 3500) {
  validate_gene_id(gene_id)
  validate_disease_name(disease)
  validate_positive_integer(size, "size")

  disease_id <- get_trait_id(disease, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  ensembl_id <- gene_to_ensembl(gene_id)
  res <- otargen::uniprotLiteratureQuery(ensembl_id, efoId = disease_id, size = size)
  return(as.data.frame(res))
}
