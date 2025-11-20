# OpenTargetsR

An enhanced R interface for querying the Open Targets Platform, providing user-friendly wrapper functions around the otargen package.

## Overview

OpenTargetsR simplifies data retrieval from the Open Targets Platform by accepting human-readable gene symbols and disease names, automatically converting them to the required identifiers for API queries. This package provides a more intuitive interface compared to the original otargen package while maintaining full functionality.

## Features

- **User-friendly input**: Use HGNC gene symbols and disease names instead of complex identifiers
- **Automatic ID resolution**: Seamlessly converts gene symbols to Ensembl IDs and disease names to EFO/MONDO IDs
- **Comprehensive coverage**: Access all major Open Targets data types including:
  - Gene-disease associations
  - Drug information and mechanisms
  - Genetic variant effects
  - Pathway and ontology data
  - Safety and pharmacogenomics information
- **Consistent interface**: All functions follow consistent naming and parameter conventions
- **Rich documentation**: Extensive examples and parameter descriptions for all functions

## Installation

```r
# Install from GitHub (once published)
# devtools::install_github("liuzaoqu/OpenTargetsR")

# For now, install locally after downloading
# devtools::install_local("/path/to/OpenTargetsR")
```

## Quick Start

```r
library(OpenTargetsR)

# Get ChEMBL evidence for PSEN1 and Alzheimer disease
chembl_data <- ChemblGeneDiseaseQuery("PSEN1", "Alzheimer disease")

# Get adverse events for aspirin
adverse_events <- DrugAdverseEventsQuery("aspirin")

# Get gene ontology annotations for THBS2
go_data <- GeneOntologyGeneQuery("THBS2")

# Get variant effect predictions
vep_data <- VariantEffectPredictorVariantQuery("4_1804392_G_A")
```

## Function Categories

### Gene-Disease Queries
- `ChemblGeneDiseaseQuery()` - ChEMBL evidence data
- `EuropePMCGeneDiseaseQuery()` - Literature evidence
- `GeneBurdenGeneDiseaseQuery()` - Gene burden analysis
- `GenomicsEnglandGeneDiseaseQuery()` - Clinical genomics evidence
- `GWASCredibleSetsTargetDiseaseQuery()` - GWAS credible sets
- `OrphanetGeneDiseaseQuery()` - Rare disease evidence
- `UniProtLiteratureGeneDiseaseQuery()` - Protein literature evidence

### Gene-Only Queries
- `ComparativeGenomicsGeneQuery()` - Evolutionary conservation
- `DepMapEssentialityGeneQuery()` - Gene essentiality scores
- `GeneOntologyGeneQuery()` - GO annotations
- `GeneticConstraintGeneQuery()` - Constraint metrics
- `HallmarkGeneQuery()` - Cancer hallmarks
- `InteractionsGeneQuery()` - Protein interactions
- `KnownDrugsGeneQuery()` - Targeting drugs
- `MousePhenotypesGeneQuery()` - Mouse model data
- `PathwaysGeneQuery()` - Pathway annotations
- `GenePharmacogenomicsQuery()` - Pharmacogenomics
- `SafetyLiabilitiesGeneQuery()` - Safety profiles

### Drug Queries
- `DrugAdverseEventsQuery()` - Adverse events
- `DrugIndicationsQuery()` - Approved indications
- `DrugKnownDrugsQuery()` - Related compounds
- `DrugMechanismsOfActionQuery()` - Mechanisms of action
- `DrugPharmacogenomicsQuery()` - Genetic effects

### Variant Queries
- `VariantPharmacogenomicsQuery()` - PGx associations
- `QTLCredibleSetsVariantQuery()` - QTL data
- `UniProtVariantsVariantQuery()` - Protein consequences
- `VariantEffectPredictorVariantQuery()` - VEP annotations
- `VariantEffectVariantQuery()` - General effects

### Disease Queries
- `SharedTraitStudiesDiseasesQuery()` - Cross-trait genetics

### Utility Functions
- `get_chembl_id()` - Resolve drug names to ChEMBL IDs
- `get_trait_id()` - Resolve disease names to EFO/MONDO IDs

## Examples

### Gene-Disease Analysis
```r
# Comprehensive analysis of APOE and Alzheimer disease
gene <- "APOE"
disease <- "Alzheimer disease"

# Get multiple evidence types
chembl_evidence <- ChemblGeneDiseaseQuery(gene, disease)
literature_evidence <- EuropePMCGeneDiseaseQuery(gene, disease)
gwas_evidence <- GWASCredibleSetsTargetDiseaseQuery(gene, disease)
uniprot_evidence <- UniProtLiteratureGeneDiseaseQuery(gene, disease)
```

### Drug Safety Analysis
```r
# Analyze aspirin safety profile
drug <- "aspirin"

adverse_events <- DrugAdverseEventsQuery(drug)
mechanisms <- DrugMechanismsOfActionQuery(drug)
indications <- DrugIndicationsQuery(drug)
pharmacogenomics <- DrugPharmacogenomicsQuery(drug)
```

### Gene Function Analysis
```r
# Comprehensive gene function analysis
gene <- "TP53"

go_annotations <- GeneOntologyGeneQuery(gene)
pathways <- PathwaysGeneQuery(gene)
interactions <- InteractionsGeneQuery(gene, size = 100)
hallmarks <- HallmarkGeneQuery(gene)
constraint <- GeneticConstraintGeneQuery(gene)
safety <- SafetyLiabilitiesGeneQuery(gene)
```

## Dependencies

- otargen: Core Open Targets API interface
- BioEnricher: Gene ID mapping
- dplyr, tibble, tidyr: Data manipulation
- httr, jsonlite: Web API access
- stringr, purrr, rlang: Utility functions

## Author

Zaoqu Liu (liuzaoqu@163.com)

## License

MIT License

## Contributing

Please report issues and suggestions on the GitHub repository.