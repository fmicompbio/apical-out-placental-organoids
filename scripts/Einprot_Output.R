
# script to recreate the einprot output from Calvi et al., 2026
#
# Charlotte Soneson and Hans-Rudolf Hotz
# 04-Aug-2026


# download the data stored at https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b6e6d9caa31f41d594634ed7b5fe6502
# in particular: 
#   A5_260722_JS_4658_ExperimentSetupOverview_PAF_einprot (default BGS).txt
#   A5_260722_JS_4658_Report_einprot_Pivot_PG.all (Pivot).tsv
#   Spectronaut_4758_PG.Quantity_MinProbGlobal_center.median.shared_limma_einprot0.10.2_sampleAnnot.txt

## Load packages
suppressPackageStartupMessages({
    library(einprot) # einprot version: 0.10.2
})

## ------------------------------------------------------------------------- ##
## Input/output file paths
## ------------------------------------------------------------------------- ##

## Output directory where all output will be stored - will be created if needed
outputDir = "wo_S07_S08_S11_S12_CScore7"

if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
}

## Path to Spectronaut setup.txt file 
spectronautLogFile <- "Calvi_et_al_2026/A5_260722_JS_4658_ExperimentSetupOverview_PAF_einprot (default BGS).txt"

## Path to Spectronaut output file ('Report_Pivot_PG.tsv' or 'Report.tsv')
spectronautFile <- "Calvi_et_al_2026/A5_260722_JS_4658_Report_einprot_Pivot_PG.all (Pivot).tsv"

## File type ("pg_pivot", "long_format")
spectronautFileType <- "pg_pivot"

## Desired output level ("pg")
outLevel <- "pg"

## ------------------------------------------------------------------------- ##
## Information about the experiment. The experiment ID will be used to retrieve
## additional information.
## ------------------------------------------------------------------------- ##

experimentId <- 4758
analysisDetails <- ""

cysAlkylation <- "CAM-Cys"
sampleIs <- "digested"
enzymes <- "Lys-C + Trypsin"

species = "Human"

## ------------------------------------------------------------------------- ##
## Spectronaut intensity parameters
## ------------------------------------------------------------------------- ##

## The column to use for the main assay
aName <- "PG.Quantity"

## Regular expression defining columns in PG pivot output that represent
## intensities
## Typically ".PG.Quantity$" or ".PEP.Quantity$"
## Note that when analyzing peptide-level data, you may want to set
## addInteractiveVolcanos and/or addHeatmaps to FALSE to speed up the analysis
## and avoid excessively large output files
iColPattern <- ".PG.Quantity$"

## Assay(s) to use for exported values and bar plots. Should be the name of
## one or more assays with absolute abundance values, if available. If not,
## can be set to NULL, in which case the main assay will be used.
assaysForExport <- NULL

## Whether to include abundance values (using the assaysForExport defined
## above) in the result tables.
addAbundanceValues <- TRUE

## ------------------------------------------------------------------------- ##
## Sample annotations.
## ------------------------------------------------------------------------- ##

## read from tab-delimited text file 
sampleAnnotFile <- "Calvi_et_al_2026/Spectronaut_4758_PG.Quantity_MinProbGlobal_center.median.shared_limma_einprot0.10.2_sampleAnnot.txt"
sampleAnnot <- read.delim(sampleAnnotFile, header = TRUE, sep = "\t")

## 'includeOnlySamples' specifies a subset of samples to include in the
## analysis. The provided value will be matched with the sample names.
## Similarly, 'excludeSamples' specifies a subset of samples to exclude.
## Max one of these can be specified.
includeOnlySamples <- ""
# excludeSamples <- ""
excludeSamples <- c("A5_260722_JS_4758_ApicalOUT_Z149_S07", "A5_260722_JS_4758_ApicalIN_Z149_S08", "A5_260722_JS_4758_ApicalOUT_Z143_S11", "A5_260722_JS_4758_ApicalIN_Z143_S12")

## ------------------------------------------------------------------------- ##
## Comparisons
## ------------------------------------------------------------------------- ##
## Two or more groups can be merged to create a new group that will be used
## for comparisons. Each merged group should be an entry in the named list
## below. E.g.
## mergeGroups <- list(new_group_name = c("old_group_1", "old_group_2"))
## If not merging should be done, set mergeGroups <- list()
mergeGroups <- list()

## Specification of control group and comparisons to make in the statistical
## analysis.
## If allPairwiseComparisons is TRUE, each group will be compared to
## each other group, including the control group.
## If allPairwiseComparisons is set to FALSE, each non-control group will
## be compared to the control group.
## If comparisons is a non-empty list, only comparisons specified here will
## be performed. E.g.:
## comparisons <- list(c("group1", "group2"), c("group2", "group3"))
## To create a list of comparisons of each group vs its complement:
## comparisons <- lapply(unique(sampleAnnot$group), function(g) {
##   c("complement", g)
## })
comparisons <- list(c("ApicalOUT", "ApicalIN"),
                    c("Ctrl", "ApicalIN"),
                    c("Ctrl", "ApicalOUT")
                    )
ctrlGroup <- ""
allPairwiseComparisons <- FALSE

## ------------------------------------------------------------------------- ##
## Analysis parameters - preprocessing
## ------------------------------------------------------------------------- ##

## Filter functions
## filtersDF defines filtering functions that are applied after reading the
## long format text file (before reshaping into an intensity matrix). Each
## function should take a data.frame as input and return a filtered data.frame.
## filtersSE defines filtering functions that are applied to the
## SummarizedExperiment object. Each function should take a SummarizedExperiment
## object as input and return a logical vector of the same length as the number
## of proteins in the object, where TRUE means that the protein should be
## retained.
## To turn off, set filtersSE <- list() and filtersDF <- list()
filtersDF <- list()
filtersSE <- list(
    min7_Score = function(sce) {
        if ("PG.Cscore" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$PG.Cscore >= 7 & !is.na(rowData(sce)$PG.Cscore)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    min2_Peptides = function(sce) {
        if ("PG.NrOfStrippedSequencesIdentified.Experiment.wide" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$PG.NrOfStrippedSequencesIdentified.Experiment.wide >= 2 &
                !is.na(rowData(sce)$PG.NrOfStrippedSequencesIdentified.Experiment.wide)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Contaminants = function(sce) {
        if ("PG.ProteinGroups" %in% colnames(rowData(sce))) {
            kp <- !grepl("^contam_", rowData(sce)$PG.ProteinGroups) &
                !is.na(rowData(sce)$PG.ProteinGroups)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Decoys = function(sce) {
        if ("PG.ProteinGroups" %in% colnames(rowData(sce))) {
            kp <- !grepl("_Decoy$", rowData(sce)$PG.ProteinGroups) &
                !is.na(rowData(sce)$PG.ProteinGroups)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    }
)

## Imputation method (either 'impSeqRob', 'MinProb', 'MinProbGlobal' or 'custom')
imputeMethod <- "MinProbGlobal"

## Additional imputation arguments
## For example, to run MinProb imputation row-wise with imputation around
## the 0.4 quantile, set imputeMethod = "MinProb" and
## imputeArgs = list(q = 0.4, MARGIN = 1L)
imputeArgs <- list()

## Normalization method (set to 'none' to skip normalization)
## Options: 'center.mean', 'center.median', 'diff.median', 'vsn',
## 'center.mean.shared', 'center.median.shared'
## See ?MsCoreUtils::normalizeMethods() [https://rdrr.io/bioc/MsCoreUtils/man/normalize.html]
## for more details on the different methods.
normMethod <- "center.median.shared"

## Vector of spike-in features used to estimate normalization factors
spikeFeatures = NULL

## ------------------------------------------------------------------------- ##
## Analysis parameters - testing
## ------------------------------------------------------------------------- ##

## Statistical test (one of "limma", "ttest", "proDA" or "none")
stattest <- "limma"

## Whether to do a single model fit and extract results via contrasts
## (if FALSE, will fit a separate model for each comparison)
singleFit <- FALSE

## Whether to explicitly subtract the "baseline"/reference value in each
## batch. Should only be used when the same technical sample is repeated
## across batches. baselineGroup corresponds to the group ID representing the
## technical sample.
subtractBaseline <- FALSE
baselineGroup <- ""

## Minimum number of valid (non-imputed) values required to compare two groups
minNbrValidValues <- 4

## Minimal fold change to consider in the statistical tests (only used if
## stattest = "limma"). The null hypothesis of the tests will be that the
## absolute log2FC is smaller than or equal to this value.
## Set this to a value below which you are definitely not interested in log2FC.
minlFC <- 0

## Whether to use the SAM statistic+permutations (similar to Perseus) for
## determining significance for the t-test. Setting this to FALSE will speed
## up the analysis for large data sets.
samSignificance <- TRUE

## ------------------------------------------------------------------------- ##
## Analysis parameters - plotting
## ------------------------------------------------------------------------- ##

## Parameters for labeling the features in the volcano plots
## For limma, features passing `volcanoAdjPvalThr` and `volcanoLog2FCThr`
## will be colored in the plots, and the top `volcanoMaxFeatures` of
## these (furthest away from the origin) will be labeled. Features in
## `volcanoFeaturesToLabel` will always be labeled in the plots.
## For the t-test, volcanoAdjPvalThr and volcanoS0 will be used (via the
## SAM test) to determine which features are significant, mimicking the
## Perseus approach.
## volcanoMaxFeatures will also be used to determine the maximal number of top
## features to show in the waterfall plot (selected as the ones with the
## largest log-fold changes).
## Set volcanoLabelSign to 'pos' or 'neg' (rather than 'both') to only label
## features with positive/negative logFC in the volcano plots.
## If `addInteractiveVolcanos` is TRUE, additional interactive volcano plots
## will be added to the report.
## `interactiveDisplayColumns` should be a named vector with the columns
## that should be shown when hovering over a point in the interactive
## volcano plots. The `interactiveGroupColumn` groups features such that when
## you hover over a point in the interactive heatmaps, you see all the
## other points with the same value of this column highlighted.
## If `addHeatmaps` is TRUE, heatmaps will be generated. For large data sets,
## one may want to set this to FALSE to speed up the analysis.
volcanoAdjPvalThr <- 0.05
volcanoLog2FCThr <- 1
volcanoMaxFeatures <- 40
volcanoLabelSign <- "both"
volcanoS0 <- 0.1
volcanoFeaturesToLabel <- c("")
addInteractiveVolcanos <- TRUE
interactiveDisplayColumns <- c(Label = "einprotLabel",
                               pid = "pid",
                               Gene = "einprotGene",
                               Accession = "einprotProtein",
                               Description = "PG.Protein.Descriptions",
                               Organism = "PG.Organisms",
                               adjP = "adj.P.Val",
                               logFC = "logFC",
                               Peptides = "PG.NrOfStrippedSequencesIdentified.Experiment.wide",
                               Score = "PG.Cscore")
interactiveGroupColumn <- NULL
addHeatmaps <- TRUE

## ------------------------------------------------------------------------- ##
## Analysis parameters - gene sets/complexes
## ------------------------------------------------------------------------- ##

## Feature collections to test (a subset of c('complexes', 'GO', 'pathways'))
includeFeatureCollections <- c('complexes', 'GO', 'pathways')

## Smallest size (=number of features overlapping with data set) required
## for a feature set to keep it for testing
minSizeToKeepSet <- 2

## Custom complexes. Should be a named (or empty) list
## E.g., list(my_complex = c("Gene1", "Gene2"))
customComplexes <- list()

## Significance threshold for complexes (to generate individual plots)
complexFDRThr <- 0.1

## Which complexes to test. Either "current" or "all".
## If "current", only complexes defined for the current species are used.
## If "all", all complexes (mapped orthologs) are used.
complexSpecies <- "current"

## Path to the complex database. If set to NULL, the latest complex database
## contained in einprot will be used (this only contains the complexes from
## the public databases).
complexDbPath <- NULL

## The maximum number of complexes to generate plots for
maxNbrComplexesToPlot <- 10

## The maximum similarity allowed between a complex C and complexes for which
## volcano plots have already been generated in order for a volcano plot for
## C to also be generated. A value of 1 means that all complexes are plotted,
## smaller values exclude some complexes.
maxComplexSimilarity <- 1

## STRING version and path to STRING files
## If stringDir is set to "", einprot will download the required files to a
## temporary directory
stringVersion <- "12.0"
stringDir <- ""

## ------------------------------------------------------------------------- ##
## Other parameters
## ------------------------------------------------------------------------- ##

## Feature identifiers
## The ID column will be used as the primary identifier for the features, and
## can not contain duplicated entries.
## The label column will be used to label points in plots.
## The gene ID column will be used to match against GO terms and annotated
## complexes.
## The protein ID column will be interpreted as (possibly semi-colon separated)
## UniProt IDs and will be used to auto-generate URLs for lookup.
## The STRING ID column will be used to match against STRING identifiers.
## Each of these variables can be either the name of a column from the input
## data file (after applying make.names() to the column names), or a function
## that takes a data frame as input (corresponding to the input file) and
## returns a character vector.
## For pg_pivot input:
idCol <- function(df) combineIds(df, combineCols = c("PG.Genes",
                                                     "PG.ProteinGroups"))
labelCol <- function(df) combineIds(df, combineCols = c("PG.Genes",
                                                        "PG.ProteinGroups"))
geneIdCol <- function(df) getFirstId(df, colName = "PG.Genes")
proteinIdCol <- "PG.ProteinGroups"
stringIdCol <- function(df) getFirstId(df, colName = "PG.ProteinGroups")

## Extra feature annotation columns to extract. Can be a named list, an
## empty list or NULL. Each list entry should be either the name of a
## column or a function that takes a data frame and returns a character
## vector, like the annotation functions above.
extraFeatureCols <- NULL

## Extra columns to include in the link table (can be a regular expression)
linkTableColumns <- c("^PG.ProteinDescriptions$", "^PG.Organisms$",
                      "^PG.FastaFiles$", "^PG.Protein.Existence$",
                      "^PG.NrOfStrippedSequencesIdentified.Experiment.wide$")

## Random seed for any non-deterministic calculations
seed <- 123

## Report title
reportTitle <- "Spectronaut data processing"

## Report author
reportAuthor <- ""

## ------------------------------------------------------------------------- ##
## Generate the base name of the output based on the parameters above
## ------------------------------------------------------------------------- ##

outputBaseName <- paste0("Spectronaut_", experimentId, "_", aName, "_",
                         imputeMethod, "_", normMethod, "_", stattest,
                         "_einprot", packageVersion("einprot"))

## ------------------------------------------------------------------------- ##
## Run analysis, using parameters as specified above.
## ------------------------------------------------------------------------- ##

## Generate experimentInfo and species information
experimentInfo <- c(
    list("Experiment ID" = experimentId
         #"Species" = species
    ),
    list("Analysis details" = analysisDetails,
         "Sample annotation file" = sampleAnnotFile)
)
if (!is.null(cysAlkylation)) {
    experimentInfo$cysAlkylation <- cysAlkylation
}
if (!is.null(sampleIs)) {
    experimentInfo$`Sample is` <- sampleIs
}
if (!is.null(enzymes)) {
    experimentInfo$Enzymes <- enzymes
}
#species <- experimentInfo$Species
#experimentInfo$Species <- NULL


## Run analysis
html <- runSpectronautAnalysis(
        templateRmd = system.file("extdata/process_basic_template.Rmd",
                                  package = "einprot"),
        outputDir = outputDir,
        outputBaseName = outputBaseName,
        reportTitle = reportTitle,
        reportAuthor = reportAuthor,
        forceOverwrite = TRUE,
        experimentInfo = experimentInfo,
        species = species,
        spectronautFile = spectronautFile,
        spectronautFileType = spectronautFileType,
        outLevel = outLevel,
        spectronautLogFile = spectronautLogFile,
        aName = aName,
        idCol = idCol,
        labelCol = labelCol,
        geneIdCol = geneIdCol,
        proteinIdCol = proteinIdCol,
        stringIdCol = stringIdCol,
        extraFeatureCols = extraFeatureCols,
        iColPattern = iColPattern,
        sampleAnnot = sampleAnnot,
        includeOnlySamples = includeOnlySamples,
        excludeSamples = excludeSamples,
        filtersDF = filtersDF,
        filtersSE = filtersSE,
        imputeMethod = imputeMethod,
        imputeArgs = imputeArgs,
        assaysForExport = assaysForExport,
        addAbundanceValues = addAbundanceValues,
        addHeatmaps = addHeatmaps,
        comparisons = comparisons,
        ctrlGroup = ctrlGroup,
        allPairwiseComparisons = allPairwiseComparisons,
        singleFit = singleFit,
        subtractBaseline = subtractBaseline,
        baselineGroup = baselineGroup,
        normMethod = normMethod,
        spikeFeatures = spikeFeatures,
        stattest = stattest,
        minNbrValidValues = minNbrValidValues,
        minlFC = minlFC,
        samSignificance = samSignificance,
        #nperm = nperm,
        volcanoAdjPvalThr = volcanoAdjPvalThr,
        volcanoLog2FCThr = volcanoLog2FCThr,
        volcanoMaxFeatures = volcanoMaxFeatures,
        volcanoLabelSign = volcanoLabelSign,
        volcanoS0 = volcanoS0,
        volcanoFeaturesToLabel = volcanoFeaturesToLabel,
        addInteractiveVolcanos = addInteractiveVolcanos,
        interactiveDisplayColumns = interactiveDisplayColumns,
        interactiveGroupColumn = interactiveGroupColumn,
        complexFDRThr = complexFDRThr,
        maxNbrComplexesToPlot = maxNbrComplexesToPlot,
        maxComplexSimilarity = maxComplexSimilarity,
        seed = seed,
        includeFeatureCollections = includeFeatureCollections,
        minSizeToKeepSet = minSizeToKeepSet,
        customComplexes = customComplexes,
        complexSpecies = complexSpecies,
        complexDbPath = complexDbPath,
        stringVersion = stringVersion,
        stringDir = stringDir,
        mergeGroups = mergeGroups,
        linkTableColumns = linkTableColumns,
        #customYml = customYml,
        doRender = TRUE
)


