#' Resolve a drug name to ChEMBL molecule ID(s)
#'
#' This function queries the ChEMBL web service to resolve a drug or molecule
#' common name to its corresponding ChEMBL molecule identifier(s). It supports
#' both exact and fuzzy matching, and can return either the best match or all
#' matching results.
#'
#' @param term Character string. The drug or molecule common name to search for.
#' @param exact Logical. If TRUE, performs exact matching on preferred name;
#'   if FALSE (default), performs fuzzy search.
#' @param keep_all Logical. If FALSE (default), returns only the single best
#'   hit; if TRUE, returns all matching hits.
#'
#' @return A named list where each element corresponds to a ChEMBL molecule.
#'   Each element contains:
#'   \itemize{
#'     \item \code{pref_name}: The preferred name of the molecule in ChEMBL
#'     \item \code{chembl_id}: The ChEMBL molecule identifier
#'   }
#'   The list is named with the ChEMBL IDs for easy access.
#'
#' @details The function connects to the ChEMBL REST API to search for molecules.
#'   For exact searches, it queries the \code{pref_name__iexact} endpoint.
#'   For fuzzy searches, it uses the general search endpoint with a limit of
#'   500 results. Results are ranked by score when available, with higher
#'   scores indicating better matches.
#'
#' @note This function requires an internet connection to access the ChEMBL API.
#'   The search is case-insensitive for both exact and fuzzy matching.
#'
#' @examples
#' \dontrun{
#' # Get the ChEMBL ID for aspirin
#' get_chembl_id("aspirin")
#'
#' # Get all ChEMBL IDs for paracetamol
#' results <- get_chembl_id("paracetamol", keep_all = TRUE)
#'
#' # Exact search for a specific drug name
#' get_chembl_id("ASPIRIN", exact = TRUE)
#' }
#'
#' @importFrom httr GET timeout stop_for_status content
#' @importFrom jsonlite fromJSON
#' @importFrom tibble tibble
#' @importFrom dplyr distinct arrange desc slice
#' @importFrom purrr pmap
#' @importFrom rlang set_names
#' @importFrom utils URLencode
#' @export
get_chembl_id <- function(term, exact = FALSE, keep_all = FALSE) {
  # Silence R CMD check NOTEs about dplyr variables
  chembl_id <- score <- NULL

  # Input validation
  if (missing(term) || is.null(term) || !is.character(term) || length(term) != 1 || nchar(term) == 0) {
    stop("'term' must be a single non-empty character string", call. = FALSE)
  }
  if (!is.logical(exact) || length(exact) != 1) {
    stop("'exact' must be a single logical value (TRUE or FALSE)", call. = FALSE)
  }
  if (!is.logical(keep_all) || length(keep_all) != 1) {
    stop("'keep_all' must be a single logical value (TRUE or FALSE)", call. = FALSE)
  }

  base <- "https://www.ebi.ac.uk/chembl/api/data"

  # Build query URL
  url <- if (exact) {
    sprintf(
      "%s/molecule.json?pref_name__iexact=%s",
      base, utils::URLencode(term, reserved = TRUE)
    )
  } else {
    sprintf(
      "%s/molecule/search.json?q=%s&limit=500",
      base, utils::URLencode(term, reserved = TRUE)
    )
  }

  # Query ChEMBL API
  res <- tryCatch(
    {
      httr::GET(url, httr::timeout(15))
    },
    error = function(e) {
      stop(
        "Failed to connect to ChEMBL API. Please check your internet connection.\n",
        "Error: ", e$message,
        call. = FALSE
      )
    }
  )

  httr::stop_for_status(res)
  dat <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))

  # Extract hits
  hits <- if (!is.null(dat$molecules)) dat$molecules else dat
  if (length(hits) == 0 || nrow(hits) == 0) {
    stop(
      "No ChEMBL match found for '", term, "'.\n",
      "Suggestions:\n",
      "  - Check spelling\n",
      "  - Try a different name or synonym\n",
      "  - Use exact = FALSE for broader search",
      call. = FALSE
    )
  }

  # Create result table
  tbl <- tibble::tibble(
    pref_name = hits$pref_name,
    chembl_id = hits$molecule_chembl_id,
    score = if ("score" %in% names(hits)) hits$score else NA_real_
  ) |>
    dplyr::distinct(chembl_id, .keep_all = TRUE)

  # Sort by score if available
  if (!all(is.na(tbl$score))) {
    tbl <- dplyr::arrange(tbl, dplyr::desc(score))
  }

  # Keep only best hit if requested
  if (!keep_all) {
    tbl <- dplyr::slice(tbl, 1)
  }

  # Return as named list
  purrr::pmap(tbl, ~ list(pref_name = ..1, chembl_id = ..2)) |>
    rlang::set_names(tbl$chembl_id)
}

#' Resolve a disease/trait name to EFO/MONDO identifier(s)
#'
#' This function queries the EBI Ontology Lookup Service (OLS4) to resolve
#' a disease or trait name to its corresponding EFO (Experimental Factor Ontology)
#' or MONDO (Monarch Disease Ontology) identifier(s). These identifiers are
#' required for many Open Targets API queries.
#'
#' @param term Character string. The disease or trait name to search for.
#' @param exact Logical. If TRUE, uses exact label matching; if FALSE (default),
#'   performs fuzzy search.
#' @param keep_all Logical. If FALSE (default), returns only the single best
#'   hit; if TRUE, returns all cleaned matching hits.
#'
#' @return A named list where each element corresponds to an ontology term.
#'   Each element contains:
#'   \itemize{
#'     \item \code{label}: The preferred label of the disease/trait
#'     \item \code{obo_id}: The ontology identifier (EFO_xxxxxxx or MONDO_xxxxxxx)
#'     \item \code{iri}: The full IRI (Internationalized Resource Identifier)
#'   }
#'   The list is named with the ontology IDs for easy access.
#'
#' @details The function queries both EFO and MONDO ontologies simultaneously
#'   through the OLS4 API. Results are filtered to include only EFO_ and MONDO_
#'   prefixed terms, deduplicated, and ranked by relevance score when available.
#'   Higher scores indicate better matches to the search term.
#'
#' @note This function requires an internet connection to access the OLS4 API.
#'   The search covers both ontologies comprehensively and handles term
#'   normalization automatically.
#'
#' @examples
#' \dontrun{
#' # Get the best matching identifier for familial hypercholesterolemia
#' get_trait_id("familial hypercholesterolemia")
#'
#' # Get identifier for Alzheimer disease
#' get_trait_id("Alzheimer disease")
#'
#' # Get all matching identifiers
#' results <- get_trait_id("Parkinson disease", exact = TRUE, keep_all = TRUE)
#' }
#'
#' @importFrom httr GET timeout stop_for_status content
#' @importFrom jsonlite fromJSON
#' @importFrom tibble tibble
#' @importFrom dplyr filter distinct arrange desc slice
#' @importFrom stringr str_replace str_detect
#' @importFrom purrr pmap
#' @importFrom rlang set_names
#' @export
get_trait_id <- function(term, exact = FALSE, keep_all = FALSE) {
  # Silence R CMD check NOTEs about dplyr variables
  obo_id <- score <- NULL

  # Input validation
  if (missing(term) || is.null(term) || !is.character(term) || length(term) != 1 || nchar(term) == 0) {
    stop("'term' must be a single non-empty character string", call. = FALSE)
  }
  if (!is.logical(exact) || length(exact) != 1) {
    stop("'exact' must be a single logical value (TRUE or FALSE)", call. = FALSE)
  }
  if (!is.logical(keep_all) || length(keep_all) != 1) {
    stop("'keep_all' must be a single logical value (TRUE or FALSE)", call. = FALSE)
  }

  # Query OLS4 API
  res <- tryCatch(
    {
      httr::GET(
        "https://www.ebi.ac.uk/ols4/api/search",
        query = list(
          q = term,
          ontology = "efo,mondo",
          type = "class",
          rows = 500,
          exact = tolower(exact)
        ),
        httr::timeout(15)
      )
    },
    error = function(e) {
      stop(
        "Failed to connect to OLS4 API. Please check your internet connection.\n",
        "Error: ", e$message,
        call. = FALSE
      )
    }
  )

  httr::stop_for_status(res)
  docs <- jsonlite::fromJSON(httr::content(res, "text", encoding = "UTF-8"))$response$docs

  if (length(docs) == 0 || nrow(docs) == 0) {
    stop(
      "No EFO/MONDO match found for '", term, "'.\n",
      "Suggestions:\n",
      "  - Check spelling\n",
      "  - Try a different disease name or synonym\n",
      "  - Use exact = FALSE for broader search",
      call. = FALSE
    )
  }

  # Clean and filter results
  tbl <- tibble::tibble(
    label = docs$label,
    obo_id = stringr::str_replace(docs$obo_id, ":", "_"),
    score = if ("score" %in% names(docs)) docs$score else NA_real_,
    iri = docs$iri
  ) |>
    dplyr::filter(stringr::str_detect(obo_id, "^(EFO|MONDO)_")) |>
    dplyr::distinct(obo_id, .keep_all = TRUE)

  # Rank by score if available
  if (!all(is.na(tbl$score))) {
    tbl <- dplyr::arrange(tbl, dplyr::desc(score))
  }

  # Keep only best hit if requested
  if (!keep_all) {
    tbl <- dplyr::slice(tbl, 1)
  }

  if (nrow(tbl) == 0) {
    stop(
      "No EFO/MONDO match found for '", term, "' after filtering.\n",
      "The term may not be in the EFO or MONDO ontologies.",
      call. = FALSE
    )
  }

  # Return as named list
  purrr::pmap(tbl, ~ list(label = ..1, obo_id = ..2, iri = ..4)) |>
    rlang::set_names(tbl$obo_id)
}

# Internal helper function for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

# Internal helper function for input validation
validate_gene_id <- function(gene_id, arg_name = "gene_id") {
  if (missing(gene_id)) {
    stop(
      "Argument '", arg_name, "' is missing with no default.\n",
      "Please provide a gene symbol (e.g., 'TP53', 'APOE')",
      call. = FALSE
    )
  }
  if (is.null(gene_id) || !is.character(gene_id) ||
    length(gene_id) != 1 || nchar(gene_id) == 0) {
    stop(
      "'", arg_name, "' must be a single non-empty character string (HGNC gene symbol)",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_disease_name <- function(disease, arg_name = "disease") {
  if (missing(disease) || is.null(disease) || !is.character(disease) ||
    length(disease) != 1 || nchar(disease) == 0) {
    stop(
      "'", arg_name, "' must be a single non-empty character string (disease name)",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_positive_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || x <= 0 || x != as.integer(x)) {
    stop("'", arg_name, "' must be a single positive integer", call. = FALSE)
  }
  invisible(TRUE)
}

validate_nonnegative_integer <- function(x, arg_name) {
  if (!is.numeric(x) || length(x) != 1 || x < 0 || x != as.integer(x)) {
    stop("'", arg_name, "' must be a single non-negative integer", call. = FALSE)
  }
  invisible(TRUE)
}

# Convert gene symbol to Ensembl ID with error handling
gene_to_ensembl <- function(gene_id) {
  gene_res <- tryCatch(
    {
      BioEnricher::GeneMapper(gene_id, id.type.to = "ensembl_id")
    },
    error = function(e) {
      stop(
        "Failed to map gene symbol '", gene_id, "' to Ensembl ID.\n",
        "Please check:\n",
        "  - Gene symbol spelling (case-sensitive)\n",
        "  - Use official HGNC symbols (e.g., 'TP53', not 'p53')\n",
        "  - Your internet connection\n",
        "  - Try using a gene alias or previous symbol\n",
        "\nOriginal error: ", e$message,
        call. = FALSE
      )
    }
  )

  ensembl_id <- gene_res$geneID_res$to_geneID[1]

  if (is.na(ensembl_id)) {
    stop(
      "Could not map gene symbol '", gene_id, "' to Ensembl ID.\n",
      "Please check:\n",
      "  - Gene symbol spelling (case-sensitive)\n",
      "  - Use official HGNC symbols (e.g., 'TP53', not 'p53')\n",
      "  - Try using a gene alias or previous symbol",
      call. = FALSE
    )
  }

  ensembl_id
}
