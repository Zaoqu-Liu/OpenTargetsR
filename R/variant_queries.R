# ============================================================================
# Variant Query Functions - New Naming Convention (snake_case)
# ============================================================================

#' Query Pharmacogenomics Data for a Genetic Variant
#' @param variant_id Character string. Variant ID (format: chr_pos_ref_alt).
#' @return A data frame containing pharmacogenomics data.
#' @examples
#' \dontrun{
#' pgx_data <- get_variant_pharmacogenomics("12_111446804_T_C")
#' }
#' @importFrom otargen pharmacogenomicsVariantQuery
#' @export
get_variant_pharmacogenomics <- function(variant_id) {
  if (missing(variant_id) || is.null(variant_id) || !is.character(variant_id) ||
    length(variant_id) != 1 || nchar(variant_id) == 0) {
    stop("'variant_id' must be a single non-empty character string", call. = FALSE)
  }

  if (!grepl("^\\d+_\\d+_[ATCG]+_[ATCG]+$", variant_id)) {
    warning("variant_id should be in format 'chr_pos_ref_alt' (e.g., '12_111446804_T_C')")
  }

  res <- otargen::pharmacogenomicsVariantQuery(variant_id)
  return(as.data.frame(res))
}

#' Query QTL Credible Sets Data for a Genetic Variant
#' @param variant_id Character string. Variant ID (format: chr_pos_ref_alt).
#' @param size Integer. Number of records to retrieve (default: 500).
#' @param index Integer. Page index for pagination (default: 0).
#' @return A data frame containing QTL credible sets data.
#' @examples
#' \dontrun{
#' qtl_data <- get_variant_qtl_credible_sets("19_10352442_G_C")
#' }
#' @importFrom otargen qtlCredibleSetsQuery
#' @export
get_variant_qtl_credible_sets <- function(variant_id, size = 500, index = 0) {
  if (missing(variant_id) || is.null(variant_id) || !is.character(variant_id) ||
    length(variant_id) != 1 || nchar(variant_id) == 0) {
    stop("'variant_id' must be a single non-empty character string", call. = FALSE)
  }
  validate_positive_integer(size, "size")
  validate_nonnegative_integer(index, "index")

  if (!grepl("^\\d+_\\d+_[ATCG]+_[ATCG]+$", variant_id)) {
    warning("variant_id should be in format 'chr_pos_ref_alt' (e.g., '19_10352442_G_C')")
  }

  res <- otargen::qtlCredibleSetsQuery(variantId = variant_id, size = size, index = index)
  return(as.data.frame(res))
}

#' Query UniProt Variants Data for a Genetic Variant
#' @param variant_id Character string. Variant ID (format: chr_pos_ref_alt).
#' @return A data frame containing UniProt variants data.
#' @examples
#' \dontrun{
#' uniprot_data <- get_variant_uniprot("4_1804392_G_A")
#' }
#' @importFrom otargen uniProtVariantsQuery
#' @export
get_variant_uniprot <- function(variant_id) {
  if (missing(variant_id) || is.null(variant_id) || !is.character(variant_id) ||
    length(variant_id) != 1 || nchar(variant_id) == 0) {
    stop("'variant_id' must be a single non-empty character string", call. = FALSE)
  }

  if (!grepl("^\\d+_\\d+_[ATCG]+_[ATCG]+$", variant_id)) {
    warning("variant_id should be in format 'chr_pos_ref_alt' (e.g., '4_1804392_G_A')")
  }

  res <- otargen::uniProtVariantsQuery(variantId = variant_id)
  return(as.data.frame(res))
}

#' Query Variant Effect Predictor Data for a Genetic Variant
#' @param variant_id Character string. Variant ID (format: chr_pos_ref_alt).
#' @return A data frame containing VEP annotation data.
#' @examples
#' \dontrun{
#' vep_data <- get_variant_vep("4_1804392_G_A")
#' }
#' @importFrom otargen variantEffectPredictorQuery
#' @export
get_variant_vep <- function(variant_id) {
  if (missing(variant_id) || is.null(variant_id) || !is.character(variant_id) ||
    length(variant_id) != 1 || nchar(variant_id) == 0) {
    stop("'variant_id' must be a single non-empty character string", call. = FALSE)
  }

  if (!grepl("^\\d+_\\d+_[ATCG]+_[ATCG]+$", variant_id)) {
    warning("variant_id should be in format 'chr_pos_ref_alt' (e.g., '4_1804392_G_A')")
  }

  res <- otargen::variantEffectPredictorQuery(variantId = variant_id)
  return(as.data.frame(res))
}

#' Query Variant Effect Data for a Genetic Variant
#' @param variant_id Character string. Variant ID (format: chr_pos_ref_alt).
#' @return A data frame containing variant effect data.
#' @examples
#' \dontrun{
#' effect_data <- get_variant_effect("4_1804392_G_A")
#' }
#' @importFrom otargen variantEffectQuery
#' @export
get_variant_effect <- function(variant_id) {
  if (missing(variant_id) || is.null(variant_id) || !is.character(variant_id) ||
    length(variant_id) != 1 || nchar(variant_id) == 0) {
    stop("'variant_id' must be a single non-empty character string", call. = FALSE)
  }

  if (!grepl("^\\d+_\\d+_[ATCG]+_[ATCG]+$", variant_id)) {
    warning("variant_id should be in format 'chr_pos_ref_alt' (e.g., '4_1804392_G_A')")
  }

  res <- otargen::variantEffectQuery(variantId = variant_id)
  return(as.data.frame(res))
}
