#' OpenTargetsR: Enhanced R Interface for Open Targets Platform
#'
#' An enhanced R interface for querying the Open Targets Platform, providing
#' user-friendly wrapper functions around the otargen package. This package
#' simplifies data retrieval by accepting human-readable gene symbols and
#' disease names, automatically converting them to the required identifiers
#' for API queries.
#'
#' @details
#' The OpenTargetsR package provides comprehensive access to Open Targets data
#' through intuitive function interfaces. Key features include:
#'
#' \itemize{
#'   \item Automatic ID resolution for genes and diseases
#'   \item Consistent function naming and parameter conventions
#'   \item Comprehensive documentation with examples
#'   \item Support for all major Open Targets data types
#' }
#'
#' @section Gene-Disease Queries:
#' Functions for querying associations between genes and diseases:
#' \itemize{
#'   \item \code{\link{ChemblGeneDiseaseQuery}}
#'   \item \code{\link{EuropePMCGeneDiseaseQuery}}
#'   \item \code{\link{GeneBurdenGeneDiseaseQuery}}
#'   \item \code{\link{GenomicsEnglandGeneDiseaseQuery}}
#'   \item \code{\link{GWASCredibleSetsTargetDiseaseQuery}}
#'   \item \code{\link{OrphanetGeneDiseaseQuery}}
#'   \item \code{\link{UniProtLiteratureGeneDiseaseQuery}}
#' }
#'
#' @section Gene Queries:
#' Functions for querying gene-specific information:
#' \itemize{
#'   \item \code{\link{ComparativeGenomicsGeneQuery}}
#'   \item \code{\link{DepMapEssentialityGeneQuery}}
#'   \item \code{\link{GeneOntologyGeneQuery}}
#'   \item \code{\link{GeneticConstraintGeneQuery}}
#'   \item \code{\link{HallmarkGeneQuery}}
#'   \item \code{\link{InteractionsGeneQuery}}
#'   \item \code{\link{KnownDrugsGeneQuery}}
#'   \item \code{\link{MousePhenotypesGeneQuery}}
#'   \item \code{\link{PathwaysGeneQuery}}
#'   \item \code{\link{GenePharmacogenomicsQuery}}
#'   \item \code{\link{SafetyLiabilitiesGeneQuery}}
#' }
#'
#' @section Drug Queries:
#' Functions for querying drug-related information:
#' \itemize{
#'   \item \code{\link{DrugAdverseEventsQuery}}
#'   \item \code{\link{DrugIndicationsQuery}}
#'   \item \code{\link{DrugKnownDrugsQuery}}
#'   \item \code{\link{DrugMechanismsOfActionQuery}}
#'   \item \code{\link{DrugPharmacogenomicsQuery}}
#' }
#'
#' @section Variant Queries:
#' Functions for querying genetic variant information:
#' \itemize{
#'   \item \code{\link{VariantPharmacogenomicsQuery}}
#'   \item \code{\link{QTLCredibleSetsVariantQuery}}
#'   \item \code{\link{UniProtVariantsVariantQuery}}
#'   \item \code{\link{VariantEffectPredictorVariantQuery}}
#'   \item \code{\link{VariantEffectVariantQuery}}
#' }
#'
#' @section Utility Functions:
#' Helper functions for ID resolution:
#' \itemize{
#'   \item \code{\link{get_chembl_id}}
#'   \item \code{\link{get_trait_id}}
#' }
#'
#' @author Zaoqu Liu \email{liuzaoqu@163.com}
#' @keywords package
"_PACKAGE"
