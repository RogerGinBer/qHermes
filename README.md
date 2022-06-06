# qHermes

qHermes is the quantitative adaptation of the HERMES strategy for LC-MS signal characterization presented [here](https://www.nature.com/articles/s41592-021-01307-z) and implemented in the `RHermes` package at [RogerGinBer/RHermes](https://github.com/RogerGinBer/RHermes). qHermes uses `xcms` as a peak detection/alignment/grouping framework and annotates the signals either using a peak-free strategy (EIC-based) or a feature-based one. 


Four quantification functions are currently implemented: 
- `fastSOI` (EIC quantification, formula-adduct DB annotation), 
- `fastSOIFromList` (EIC quantification, HERMES-derived SOI list), 
- `annotateFeatures` (XCMS quantification, formula-adduct DB annotation) and 
- `annotateFeaturesFromList` (XCMS quantification, HERMES-derived SOI list).

Two plotting functions are implemented to check each feature quantification, whether be from an `XCMSnExp` object or a `SummarizedExperiment` object:

- `plotFeatures`
- `plotSummarizedExp`

## Installation

```
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("RogerGinBer/qHermes")
```

## Usage

## Issues

Bugs, suggestions and feature requests are welcome at https://github.com/RogerGinBer/qHermes/issues

The package is still in an **unstable/development stage**, so function names and parameters can (and will) change.
