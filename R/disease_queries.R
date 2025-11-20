# ============================================================================
# Disease Query Functions - New Naming Convention (snake_case)
# ============================================================================

#' Query Shared Trait Studies Data for Diseases
#'
#' Retrieves shared trait studies data for a specified disease, showing studies
#' and genetic associations that are shared between different diseases or traits.
#'
#' @param disease Character string. A disease name that will be resolved to EFO/MONDO ID.
#' @param size Integer. Number of records to retrieve (default: 500).
#' @param index Integer. Page index for pagination (default: 0).
#'
#' @return A data frame containing shared trait studies data.
#'
#' @examples
#' \dontrun{
#' shared_studies <- get_disease_shared_studies("Alzheimer disease")
#' shared_studies <- get_disease_shared_studies("Alzheimer disease", size = 1000)
#' }
#'
#' @importFrom otargen sharedTraitStudiesQuery
#' @export
get_disease_shared_studies <- function(disease, size = 500, index = 0) {
  validate_disease_name(disease)
  validate_positive_integer(size, "size")
  validate_nonnegative_integer(index, "index")

  disease_id <- get_trait_id(disease, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  res <- otargen::sharedTraitStudiesQuery(disease_id)
  return(as.data.frame(res))
}
