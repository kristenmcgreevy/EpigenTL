# EpigenTL
EpigenTL (Epigenetic Transfer Learning) provides tools to integrate diverse biological datasets with transfer learning (TL), with a focus on but not limited to, epigenetic DNA methylation (DNAm) data. Designed to facilitate cross-tissue and cross-dataset predictions, this package is a powerful tool for researchers looking to leverage external data sources to enhance their analyses, biomarker development, and cross tissue predictions in a penalized regression framework. While initially tailored for methylation data, EpigenTL's algorithms are broadly applicable, offering a flexible framework for any data type where leveraging auxiliary datasets can provide added insight.

This package also provides algorithms to estimate blood DNAm biomarkers from saliva DNAm profiles in humans, and skin from mammals, foregoing the need for invasive blood draws. 

### Key Features:
**Saliva to Blood DNAm Biomarker Prediction:** Utilizes saliva DNAm profiles to predict blood DNAm biomarkers, supporting both "C" (CpG-only) and "C+S" (CpG + Saliva DNAm Biomarkers) algorithms.

**Skin to Blood DNAm Biomarker Prediction in Mammals:** Utilizes skin DNAm profiles to predict blood DNAm biomarkers, supporting "C" (CpG-only) algorithms for 17 of the DNAm Biomarkers.

**Flexible Transfer Learning with Lasso Regression:** Implements transfer learning with Lasso regression to enhance algorithmic prediction accuracy. 

**Flexible Levels of Auxiliary Information:** Our functions accommodate various levels of auxiliary dataset information, including whether the auxiliary datasets are known to be informative for the target (Oracle, Oracle 1df) and whether their relevance needs estimated (Estimate A0), ensuring flexibility and methods to best capture information in your data.

**Dynamic Information Estimation:** Uniquely capable of estimating the informativeness of auxiliary datasets, with user controls on how much data is considered, optimizing the integration process without prior knowledge of dataset relevance.

**Optimized Performance:** Offers options for lambda calculations used in Lasso ("Constant" or "CV" for cross-validation) and coefficient thresholding strategies for fine-tuning model performance.

**Versatile Coefficient Aggregations:** Offers multiple coefficient aggregation methods, allowing users to choose the best fit for their analysis needs.


### Getting Started:
Ensure the `glmnet` and `dplyr` packages are installed in R. Download EpigenTL's files and source the functions to begin utilizing its capabilities.

Detailed documentation and examples are provided to help users quickly familiarize themselves with the package's capabilities and apply them to their research projects. 

EpigenTL is not just a tool for epigenetic research; it's a comprehensive solution for any scientific investigation looking to enrich analyses with auxiliary data insights. Whether you're exploring methylation patterns across tissues or integrating disparate biological datasets, EpigenTL provides the algorithms and tools to push your research forward.
