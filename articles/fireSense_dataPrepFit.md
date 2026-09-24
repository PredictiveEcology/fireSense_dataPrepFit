---
title: "fireSense_dataPrepFit Manual"
subtitle: "v.1.2.0.9010"
date: "Last updated: 2026-09-24"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
    theme: sandstone
    number_sections: false
    df_print: paged
    keep_md: yes
editor_options:
  chunk_output_type: console
bibliography: citations/references_fireSense_dataPrepFit.bib
link-citations: true
always_allow_html: true
pkgdown:
  as_is: true
---

# fireSense_dataPrepFit Module

<!-- the following are text references used in captions for LaTeX compatibility -->

(ref:fireSense-dataPrepFit) *fireSense_dataPrepFit*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut, cre], Eliot McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut], Alex M Chubaty <achubaty@for-cast.ca> [ctb] <!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

Prepares the data needed to fit `fireSense_IgnitionFit`, `fireSense_EscapeFit` and `fireSense_SpreadFit`.

### Module summary

The module combines historical fire records, vegetation, land cover and climate into the covariate tables used to fit the three `fireSense` processes [@Marchal:2017a; @Marchal:2017b; @Marchal:2019]: ignition, escape (an ignition that grows beyond one pixel) and spread.

Vegetation, land cover and stand age are snapshots, one per year in `P(sim)$dataYears` (default 2000, 2010, 2020).
Each year in `P(sim)$fireYears` uses the snapshot at or before it, so no fire year may precede the first data year, and every data year needs at least one fire year.
Climate is annual and is matched to the year of each fire.
Ignition and escape share climate variables (`sim$climateVariablesForFire$ignition`); spread can use others (`$spread`).

By default the module builds its own inputs:

- land cover (`rstLCCs`) with `fireSenseUtils::makeFireSenseLCC()` and stand age with `LandR::prepInputsStandAgeMap()`, per data year;
- `cohortDatas` and `pixelGroupMaps`, by running `Biomass_borealDataPrep` (and `Biomass_speciesData`, if it is in the project) in a nested `simInitAndSpades()` per data year. This is cached. Missing modules are downloaded to this module's `submodules` folder;
- fire polygons from the latest NBAC release, and ignition points (lightning- and natural-caused only) from the current NFDB release.

`historicalClimateRasters` has no default and must be supplied, e.g., by `canClimateData`.
It must cover every fire year; the module stops otherwise.

##### Previous SpreadFit results

At `init`, the module looks in a Google Drive ledger (`P(sim)$spreadFitGoogleDriveFolder`, `P(sim)$spreadFitFilename`) for SpreadFit results that overlap `sim$studyArea`.
If there are any, `sppEquiv`, the fuel classes, the non-forest groups and the spread climate variables are taken from that fit, and fuel classes are not estimated.
`sim$climateVariables` must then be supplied.

##### Fuel classes

If `P(sim)$estimateFuelClasses` is `TRUE`, and the user supplied none of `nonForestedLCCGroups`, `fuelClassTable` or a modified `FuelClass` column in `sppEquiv`, fuel classes are estimated with `fireSenseUtils::assessFuelClasses()`.

Non-forest land cover classes are grouped by their historical propensity to burn, estimated with a GLM, into two groups.

For forest, each species starts as its own fuel class with an *a priori* class, by default an interpretation of the [FBP fuel types](https://cwfis.cfs.nrcan.gc.ca/background/fueltypes/c1) (`LandR::sppEquivalencies_CA$FuelClass`).
A GLM relates the biomass of each species, in pixels where it is present, to whether the pixel burned.
Species are then merged until `P(sim)$targetFuelClasses` (default 5) remain:

1. Abundance is the proportion of forested pixels where the species is more than 10% of the biomass.
2. Rare species (abundance < 5%) are merged with another species of the same *a priori* class, least abundant first. This can leave fewer classes than the target.
3. If there are still too many classes, the next least abundant species are merged, but only with species of the same *a priori* class and the same coefficient sign (or a non-significant one).

::: {.example #of semi-automated fuel classes}
7 species were initially present on the landscape; 2 must be combined to achieve the target 5
in this case, Betu_pap and Popu_tre were combined due to their low abundance (above10PctRelB)
Pice_gla and Pice_eng were combined due as Abie_las was the most abundant. Pice_mar and Pinu_con
were not combined as there were no other species with similar FuelClass values


| species  | coef     | sign     | FuelClass  | Abundance | newFuelClass |
|----------|----------|----------|------------|-----------|--------------|
| Pice_mar | -0.01320 | negative | BlkSprc    | 0.066     | Pice_mar     |
| Pinu_con | 0.00154  | positive | LdJkPine   | 0.625     | Pinu_con     |
| Betu_pap | -0.00170 | negative | PopBrch    | 0.051     | Bt_pa.Pp_tr  |
| Pice_eng | -0.01910 | negative | SprcFrLrch | 0.095     | Pc_en.Pc_gl  |
| Popu_tre | -0.02100 | negative | PopBrch    | 0.223     | Bt_pa.Pp_tr  |
| Pice_gla | -0.00920 | negative | SprcFrLrch | 0.352     | Pc_en.Pc_gl  |
| Abie_las | -0.00630 | negative | SprcFrLrch | 0.560     | Abie_las     |

:::

##### youngAge, missingLCCgroup and treed wetland

Pixels disturbed within `P(sim)$cutoffForYoungAge` years (default 15) are `youngAge`, a fuel class of its own.
Time since disturbance comes from the stand age maps and the fire polygons (`firePolysForAge`), per data year.
Non-forest pixels can also be `youngAge`, unless `P(sim)$nonForestCanBeYoungAge` is `FALSE`.

`sim$nonForestedLCCGroups` takes precedence over `P(sim)$forestedLCC`: a land cover class listed in a non-forest group is a categorical fuel even if it is also in `forestedLCC`.
This lets treed wetland, for example, be simulated as forest by LandR but be a land cover fuel in `fireSense`.
A flammable class in neither gets no fuel covariate.

Pixels with forested land cover but no cohorts (usually because there was no species data) have no biomass, so they are assigned to one of the non-forest groups, `sim$missingLCCgroup`.
It is estimated along with the fuel classes; otherwise it defaults to the first group.

### Usage

The module is normally run with `canClimateData` (for `historicalClimateRasters` and `climateVariables`) ahead of the three fit modules.
`fireSense_SpreadFit` and `fireSense_EscapeFit` preparation both need the ignition preparation, so keep `fireSense_IgnitionFit` in `P(sim)$whichModulesToPrepare`.


``` r
out <- SpaDES.project::setupProject(
  paths = list(projectPath = "~/fireSenseFit"),
  modules = c("PredictiveEcology/canClimateData@development", # needs its own parameters
              "PredictiveEcology/fireSense_dataPrepFit@development"),
  times = list(start = 2020, end = 2020),
  studyArea = mySA, # buffered to limit edge effects
  params = list(fireSense_dataPrepFit = list(
    .useCache = c(".inputObjects", "dataPrepBuild", "prepIgnitionFitData",
                  "prepEscapeFitData", "prepSpreadFitData")))
)
sim <- do.call(SpaDES.core::simInitAndSpades, out)
```

### Module inputs and parameters

Table \@ref(tab:moduleInputs-fireSense-dataPrepFit) shows the full list of module inputs.

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleInputs-fireSense-dataPrepFit)(\#tab:moduleInputs-fireSense-dataPrepFit)List of (ref:fireSense_dataPrepFit) input objects and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
   <th style="text-align:left;"> sourceURL </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> climateVariables </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> Climate variable definitions for `climateData::prepClimateLayers` (canClimateData). Unless supplied, `.inputObjects` builds them from `climateVariablesForFire`. Declared as an input because `.inputObjects` sets it: a cached `.inputObjects` restores only declared inputs, so without this a cache hit returned no `climateVariables` and canClimateData failed. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> climateVariablesForFire </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List with elements `ignition` and `spread`, each a character vector of climate variable names, with or without underscores (e.g., `CMD_sm` or `CMDsm`). IgnitionFit uses all of `ignition`; SpreadFit uses `spread`. Default: `ignition = c('CMD', 'cumMDC', 'CMD_sm', 'CMD_sp')`, `spread = 'auto'`: the `ignition` variable that best separates the study area's worst fire years (see `spreadClimateSelection`). Unless supplied, `climateVariables` is built from these. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cohortDatas </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of `cohortData` data.tables, one per `dataYears`, named `year&lt;year&gt;`. If not supplied, built by running Biomass_borealDataPrep for each data year. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFirePoints </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of spatial points, one per fire year, named `year&lt;year&gt;`; each point is the ignition location of one fire in `firePolys`. The default is the polygon centroids. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitAdditionalColNames </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> Names of the ledger (`spreadFitPreRun`) columns to read. The default is `fireSenseUtils::spreadFitAdditionalColNamesTxt`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> firePolys </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of fire polygons, one per year in `fireYears`, named `year&lt;year&gt;`. The default is the latest NBAC release. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> firePolysForAge </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of annual fire polygons used for time since disturbance; as `firePolys`, but starting `cutoffForYoungAge` years before the first of `fireYears`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> historicalFireRaster </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Optional raster of fire year, 1985-2020. If supplied it replaces `firePolysForAge` for time since disturbance. Only downloaded when `useRasterizedFireForSpread = TRUE`. </td>
   <td style="text-align:left;"> https://opendata.nfis.org/downloads/forest_change/CA_Forest_Fire_1985-2020.zip </td>
  </tr>
  <tr>
   <td style="text-align:left;"> historicalClimateRasters </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of SpatRasters of historical climate, named by climate variable, with layers named `year&lt;year&gt;`. Must be supplied, and must cover `fireYears`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionFirePoints </td>
   <td style="text-align:left;"> SpatVector </td>
   <td style="text-align:left;"> Points of annual ignitions, of every fire size, with columns `YEAR` and `SIZE_HA`. The default is the lightning- and natural-caused fires of the current NFDB release. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> missingLCCgroup </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> The `nonForestedLCCGroups` name given to forested pixels that are absent from `cohortData`. The default is the first name; it is replaced if fuel classes are estimated. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForestedLCCGroups </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> Named list of non-forest land cover classes, e.g. `list(wetland = c(19, 23, 32))`. Each group becomes a fuel covariate. The default is one group, `nf`, of every class that is neither forested nor non-flammable; it is replaced if fuel classes are estimated. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pixelGroupMaps </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of `pixelGroupMap` SpatRasters matching `cohortDatas`, named `year&lt;year&gt;`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> propFlammables </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of SpatRasters of the proportion of flammable land cover in a pixel, one per `dataYears`. Built with `rstLCCs` when that is not supplied. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Template raster for `studyArea`. The default is 240 m, from SCANFI land cover. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatchLarge </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Optional larger template. If supplied it defines `studyArea_biomassParam`; if not, `.inputObjects` sets it to `rasterToMatch` for Biomass_speciesData. Declared as an input because `.inputObjects` sets it (a cached `.inputObjects` restores only declared inputs). </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch_biomassParam </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Template raster for `studyArea_biomassParam`, passed to Biomass_borealDataPrep. Expected to cover at least `rasterToMatch` (formerly `rasterToMatchLarge`). </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstLCCs </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of land cover SpatRasters, one per `dataYears`, named `year&lt;year&gt;`, on `rasterToMatch_biomassParam`. The default is from `fireSenseUtils::makeFireSenseLCC`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table of LandR species equivalencies. The default is from `LandR::speciesInStudyArea`. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> standAgeMaps </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of stand age SpatRasters, one per `dataYears`, named `year&lt;year&gt;`; used to create `cohortDatas` and time since disturbance. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFirePolys </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> Not needed from the user: `firePolys` in the CRS of `rasterToMatch`, declared as an input because a later event modifies it. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea </td>
   <td style="text-align:left;"> SpatVector </td>
   <td style="text-align:left;"> Study area for all data. Should be buffered to limit edge effects on fire spread. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea_biomassParam </td>
   <td style="text-align:left;"> SpatVector </td>
   <td style="text-align:left;"> study area passed to Biomass_borealDataPrep for vegetation calibration </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Parameters are in Table \@ref(tab:moduleParams-fireSense-dataPrepFit).

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleParams-fireSense-dataPrepFit)(\#tab:moduleParams-fireSense-dataPrepFit)List of (ref:fireSense-dataPrepFit) parameters and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> paramName </th>
   <th style="text-align:left;"> paramClass </th>
   <th style="text-align:left;"> default </th>
   <th style="text-align:left;"> min </th>
   <th style="text-align:left;"> max </th>
   <th style="text-align:left;"> paramDesc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> areaMultiplier </td>
   <td style="text-align:left;"> numeric,.... </td>
   <td style="text-align:left;"> ::, fire.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Size of the unburned buffer sampled around each fire: a scalar (buffer area is `areaMultiplier fireSize`) or a quoted function of `fireSize`. See `?fireSenseUtils::bufferToArea`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bufferForFireRaster </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1000 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Buffer distance within which separate patches of burned pixels count as one fire. Only used when `useRasterizedFireForSpread = TRUE`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cutoffForYoungAge </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 15 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Age at and below which pixels are considered 'young' (`young &lt;- age &lt;= cutoffForYoungAge`) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> dataYears </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 1985, 19.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Two or more increasing years for which vegetation, land cover and stand age are built (`cohortDatas`, `rstLCCs`, `standAgeMaps`, ...). Each fire year uses the data year at or before it, so no `fireYears` may precede the first, and every data year needs at least one fire year before the next data year. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> estimateFuelClasses </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Estimate fuel classes with `fireSenseUtils::assessFuelClasses`? Skipped if the user supplies `nonForestedLCCGroups`, `fuelClassTable` or a `FuelClass` column that differs from LandR's, or if a previous SpreadFit exists for the study area. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireYears </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 1985, 19.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Years of fire records to use for fitting. None may precede the first of `dataYears`, and `historicalClimateRasters` must cover all of them. The default runs from 1985, the first SCANFI V2 year, to the latest year with historical climate for every tile (`climateData::latestHistoricalYear()`); climate is the last of the inputs to reach a year. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammabilityThreshold </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0.1 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> Minimum proportion of flammable fine-resolution land cover for a `rasterToMatch` pixel to be flammable. Only used when `rstLCCs` is built here. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> forestedLCC </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 81, 210,.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Forested land cover classes - these differ from non-forest because the biomass and composition of fuels are taken into account by fireSense, while non-forest classes are treated categorically </td>
  </tr>
  <tr>
   <td style="text-align:left;"> igAggFactor </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Aggregation factor (number of `rasterToMatch` cells per side) for the ignition and escape covariates. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fuelClassCol </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> FuelClass </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> the column in `sppEquiv` that defines unique fuel classes. A column named `FuelClass` exists in the `LandR::sppEquivalencies_CA` and will be used by default. To change the `FuelClass` classifications, add a column to that table, or to `sim$sppEquiv` and then modify this `fuelClassCol` parameter </td>
  </tr>
  <tr>
   <td style="text-align:left;"> minBufferSize </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 5000 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Minimum number of cells in each fire's burned-plus-buffer sample, applied after `areaMultiplier`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonflammableLCC </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 0, 20, 3.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Non-flammable classes in `rstLCCs`; the default is water, snow/ice, rock and barren land. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForestCanBeYoungAge </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> if TRUE, burned non-forest will be treated as `youngAge`. Recommended to be TRUE as burned forest is often classified as non-forest </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquivCol </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> LandR </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> column name in `sppEquiv` object that defines unique species in `cohortData` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitGoogleDriveFolder </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> https://.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> URL of the Google Drive folder holding the ledger of previous SpreadFit results (`spreadFitFilename`), read with `reproducible::CacheGeo`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitFilename </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> fireSens.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Name of the ledger file in `spreadFitGoogleDriveFolder`: study area polygons with their fitted SpreadFit parameters. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> targetFuelClasses </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> the target number of unique fuel classes when using semi-automated approach </td>
  </tr>
  <tr>
   <td style="text-align:left;"> useRasterizedFireForSpread </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Use `historicalFireRaster` in place of fire polygons for spread? Not currently supported: `TRUE` stops with an error when preparing SpreadFit. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> whichModulesToPrepare </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> fireSens.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Which fireSense fit modules to prep? defaults to all 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> `studyArea` name used in file names and cache tags; `NULL` derives it from `sim$studyArea`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is intended for data-type modules, where stochasticity and time are not relevant </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCacheArgs </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list(.ca.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Extra `reproducible::Cache()` arguments, by event. A cached event's digest covers this module's code but not the package functions it calls, so `dataPrepBuild` passes those in `.cacheExtra`: a changed function then re-runs the event. </td>
  </tr>
</tbody>
</table>

### Events

All events run once, at `start(sim)`.

##### init

Checks Google Drive for a previous SpreadFit (see above). It cannot be cached, because that would freeze the answer.

##### dataPrepBuild

Runs straight after `init`, before other modules' `init`.
Puts `rstLCCs` on `rasterToMatch`, and builds, per data year, `flammableRTMs`, `landcoverDTs` and `nonForest_timeSinceDisturbances`.
Fuel classes are estimated here.
Can be cached: add `"dataPrepBuild"` to `.useCache`.

##### prepIgnitionFitData

Fuel, climate and lightning rasters are aggregated by `P(sim)$igAggFactor`, and ignitions are counted per coarse pixel and year.
`sim$fireSense_ignitionCovariates` has one row per coarse pixel and year.
Fitted probabilities only apply at that resolution, so the template `sim$ignitionFitRTM` is also output, with attributes `nonNAs` (the number of rows) and `meanForestB`.
No ignition formula is built: fireSense_IgnitionFit fits with xgboost, which does not use one.

##### prepEscapeFitData

Adds an `escapes` column to the ignition covariates: the ignitions that grew beyond one `rasterToMatch` pixel (not one coarse pixel).
Builds `sim$fireSense_escapeFormula` if it was not supplied. Needs `prepIgnitionFitData`.

##### prepSpreadFitData

Spread uses fires of every cause, where ignition uses lightning- and natural-caused ones only, and only fires larger than one pixel.
Each fire needs an ignition point in a flammable pixel inside its polygon (`sim$spreadFirePoints`; default the polygon centroid); `fireSenseUtils::harmonizeFireData()` enforces this.
Each fire polygon is buffered, with the buffer size set by `P(sim)$areaMultiplier` and `P(sim)$minBufferSize`, and only burned and buffer pixels are kept (`sim$fireBufferedListDT`).
The module stops if buffer pixels are fewer than five times the burned pixels.

To keep the objects small for the optimizer, the covariates are split in two:
`sim$fireSense_annualSpreadFitCovariates`, one table of climate per fire year, and `sim$fireSense_nonAnnualSpreadFitCovariates`, one table of fuels per data year.
Fuel columns that are all zero are dropped.
Builds `sim$fireSense_spreadFormula` if it was not supplied.
`P(sim)$useRasterizedFireForSpread = TRUE` is not currently supported.

##### plotAndMessage, cleanUp

`plotAndMessage` is a placeholder; `cleanUp` frees memory.

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-fireSense-dataPrepFit)).

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:moduleOutputs-fireSense-dataPrepFit)(\#tab:moduleOutputs-fireSense-dataPrepFit)List of (ref:fireSense-dataPrepFit) outputs and their description.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> objectName </th>
   <th style="text-align:left;"> objectClass </th>
   <th style="text-align:left;"> desc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> climateVariablesForFire </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> As the input. If a previous SpreadFit exists, `spread` becomes the climate variables of that fit, and they are added to `ignition`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFitPreRun </td>
   <td style="text-align:left;"> data.frame </td>
   <td style="text-align:left;"> Ledger rows of previous SpreadFit results that overlap `studyArea`, from `CacheGeo`: a geometry column (convert with `sf::st_as_sf`), plus `polygonID` and the columns in `spreadFitAdditionalColNames`. `NULL` if there is no previous fit. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaWithSpreadParams </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> Same as `spreadFitPreRun`; not created if there is no previous fit. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppColorVect </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> Named vector of hex colours, one per species. Only created if a previous SpreadFit exists, from the `sppEquiv` stored with it. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppNameVector </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> Sorted species names (`sppEquivCol`) from the `sppEquiv` stored with a previous SpreadFit; only created if one exists. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadClimateSelection </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> With `climateVariablesForFire$spread = 'auto'`: each candidate's AUC for separating the worst quarter of fire years (by area burned), its Spearman correlation, and which was chosen. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> climateVariables </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> Climate variable definitions, as used by `climateData::prepClimateLayers` (canClimateData). Unless supplied, built from `climateVariablesForFire` for `fireYears` (and projected years unless canClimateData's `climateGCM` is 'NRV'). If a previous SpreadFit exists, the variables of that fit are added. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireBufferedListDT </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list of data.tables with fire id, `pixelID`, and buffer status </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fuelClassTable </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Fuel class assigned to each tree species; only created if fuel classes are estimated. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFirePolys </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of annual fire polygons used for SpreadFit: larger than one pixel and matched to `spreadFirePoints`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_annualSpreadFitCovariates </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of data.tables, one per fire year, of `pixelID` and the spread climate covariates in the fire buffers. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_escapeCovariates </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> ignition covariates with added column of escapes </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_escapeFormula </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> formula for escape, using fuel classes and landcover, as character </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_ignitionCovariates </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> table of aggregated ignition covariates with annual ignitions </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_nonAnnualSpreadFitCovariates </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of data.tables, one per `dataYears`, of `pixelID` and the fuel covariates in the fire buffers. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_spreadFormula </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> formula for spread, using climate and vegetation covariates, as character </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionFirePoints </td>
   <td style="text-align:left;"> SpatVector </td>
   <td style="text-align:left;"> The input, in the CRS of `rasterToMatch` and clipped to `studyArea`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionFitRTM </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Template raster with the resolution and extent of `fireSense_ignitionCovariates`. Attributes: `nonNAs`, the number of rows in that table, and `meanForestB`, the mean forest biomass per pixel. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> landcoverDTs </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of data.tables, one per `dataYears`, of `pixelID` and a 0/1 column per non-forest group, for flammable pixels. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> lightningMaps </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A 4-layer SpatRaster of lightning: lightningDays, lightningDensity, positiveCG, positiveCGdensity </td>
  </tr>
  <tr>
   <td style="text-align:left;"> missingLCCgroup </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> As the input, or the estimated group if fuel classes are estimated. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForestedLCCGroups </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> As the input, or the estimated groups if fuel classes are estimated. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForest_timeSinceDisturbances </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of SpatRasters, one per `dataYears`, of years since disturbance in flammable pixels, forested or not. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstLCCs </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> The input, on `rasterToMatch`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableRTMs </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of binary SpatRasters of flammable land cover on `rasterToMatch`, one per `dataYears`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> sppEquiv table potentially modified with new or overwritten fuel class </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFirePoints </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of ignition points, one per fire year, for fires larger than one pixel, harmonized with `spreadFirePolys`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> propFlammable </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Last element of `propFlammables`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> standAgeMap </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Last element of `standAgeMaps`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstLCC_RTM </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Last element of the output `rstLCCs`, i.e., on `rasterToMatch`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstLCC </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Last element of the input `rstLCCs`, i.e., on `rasterToMatch_biomassParam`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableRTM </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Last element of `flammableRTMs`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> landcoverDT </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Last element of `landcoverDTs`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForest_timeSinceDisturbance </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Last element of `nonForest_timeSinceDisturbances`. </td>
  </tr>
</tbody>
</table>

### Links to other modules

This module links to [fireSense_IgnitionFit](https://github.com/PredictiveEcology/fireSense_IgnitionFit), [fireSense_EscapeFit](https://github.com/PredictiveEcology/fireSense_EscapeFit), and [fireSense_SpreadFit](https://github.com/PredictiveEcology/fireSense_SpreadFit)

### Getting help

Contact the authors for help. 

<https://github.com/PredictiveEcology/fireSense_dataPrepFit/issues>

## References

<!-- autogenerated from bibligraphy -->
