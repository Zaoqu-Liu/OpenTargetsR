# ============================================================================
# Drug Query Functions - New Naming Convention (snake_case)
# ============================================================================

#' Query Adverse Events Data for a Drug
#' @param drug_name Character string. A drug name.
#' @param index Integer. Page index (default: 0).
#' @param size Integer. Number of records (default: 10).
#' @return A data frame containing adverse events data.
#' @examples
#' \dontrun{
#' adverse_events <- get_drug_adverse_events("aspirin")
#' }
#' @importFrom otargen adverseEventsQuery
#' @export
get_drug_adverse_events <- function(drug_name, index = 0, size = 10) {
  if (missing(drug_name) || is.null(drug_name) || !is.character(drug_name) ||
    length(drug_name) != 1 || nchar(drug_name) == 0) {
    stop("'drug_name' must be a single non-empty character string", call. = FALSE)
  }
  validate_nonnegative_integer(index, "index")
  validate_positive_integer(size, "size")

  chemblid <- get_chembl_id(drug_name, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  res <- otargen::adverseEventsQuery(chemblid, index = index, size = size)
  return(as.data.frame(res))
}

#' Query Drug Indications Data
#' @param drug_name Character string. A drug name.
#' @return A data frame containing drug indications data.
#' @examples
#' \dontrun{
#' indications <- get_drug_indications("aspirin")
#' }
#' @importFrom otargen indicationsQuery
#' @export
get_drug_indications <- function(drug_name) {
  if (missing(drug_name) || is.null(drug_name) || !is.character(drug_name) ||
    length(drug_name) != 1 || nchar(drug_name) == 0) {
    stop("'drug_name' must be a single non-empty character string", call. = FALSE)
  }

  chemblid <- get_chembl_id(drug_name, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  res <- otargen::indicationsQuery(chemblid)
  return(as.data.frame(res))
}

#' Query Known Drugs Data by Drug Name
#' @param drug_name Character string. A drug name.
#' @param size Integer. Number of records (default: 10).
#' @param cursor Character or NULL. Cursor for pagination (default: NULL).
#' @return A data frame containing related drugs data.
#' @examples
#' \dontrun{
#' related_drugs <- get_drug_related_drugs("aspirin")
#' }
#' @importFrom otargen knownDrugsChemblQuery
#' @export
get_drug_related_drugs <- function(drug_name, size = 10, cursor = NULL) {
  if (missing(drug_name) || is.null(drug_name) || !is.character(drug_name) ||
    length(drug_name) != 1 || nchar(drug_name) == 0) {
    stop("'drug_name' must be a single non-empty character string", call. = FALSE)
  }
  validate_positive_integer(size, "size")

  chemblid <- get_chembl_id(drug_name, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  res <- otargen::knownDrugsChemblQuery(chemblid, cursor = cursor, size = size)
  return(as.data.frame(res))
}

#' Query Mechanisms of Action Data for a Drug
#' @param drug_name Character string. A drug name.
#' @return A data frame containing mechanisms of action data.
#' @examples
#' \dontrun{
#' moa_data <- get_drug_mechanisms("aspirin")
#' }
#' @importFrom otargen mechanismsOfActionQuery
#' @export
get_drug_mechanisms <- function(drug_name) {
  if (missing(drug_name) || is.null(drug_name) || !is.character(drug_name) ||
    length(drug_name) != 1 || nchar(drug_name) == 0) {
    stop("'drug_name' must be a single non-empty character string", call. = FALSE)
  }

  chemblid <- get_chembl_id(drug_name, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  res <- otargen::mechanismsOfActionQuery(chemblid)
  return(as.data.frame(res))
}

#' Query Pharmacogenomics Data for a Drug
#' @param drug_name Character string. A drug name.
#' @return A data frame containing pharmacogenomics associations.
#' @examples
#' \dontrun{
#' pgx_data <- get_drug_pharmacogenomics("aspirin")
#' }
#' @importFrom otargen pharmacogenomicsChemblQuery
#' @export
get_drug_pharmacogenomics <- function(drug_name) {
  if (missing(drug_name) || is.null(drug_name) || !is.character(drug_name) ||
    length(drug_name) != 1 || nchar(drug_name) == 0) {
    stop("'drug_name' must be a single non-empty character string", call. = FALSE)
  }

  chemblid <- get_chembl_id(drug_name, exact = FALSE, keep_all = FALSE)[[1]][[2]]
  res <- otargen::pharmacogenomicsChemblQuery(chemblid)
  return(as.data.frame(res))
}
