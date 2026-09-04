---
title: "fireSense_dataPrepFit Manual"
subtitle: "v.1.1.1"
date: "Last updated: 2025-04-08"
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
---

# fireSense_dataPrepFit Module

<!-- the following are text references used in captions for LaTeX compatibility -->

(ref:fireSense-dataPrepFit) *fireSense_dataPrepFit*



[![made-with-Markdown](figures/markdownBadge.png)](https://commonmark.org)

<!-- if knitting to pdf remember to add the pandoc_args: ["--extract-media", "."] option to yml in order to get the badge images -->

#### Authors:

Ian Eddy <ian.eddy@nrcan-rncan.gc.ca> [aut, cre], Eliot McIntire <eliot.mcintire@nrcan-rncan.gc.ca> [aut], Alex M Chubaty <achubaty@for-cast.ca> [ctb] <!-- ideally separate authors with new lines, '\n' not working -->

## Module Overview

Prepare data required by `fireSense_IginitionFit`, `fireSense_EscapeFit`, and `fireSense_SpreadFit` modules.

### Module summary

This module uses historical fire, vegetation, and climate raster data to generate data used to fit the three processes captured by `fireSense` [@Marchal:2017a; @Marchal:2017b; @Marchal:2019]: fire ignition, fire escape (the probability that an ignition results in a fire greater than the cell resolution), and fire spread.
The landscape data is collected in two snapshots: 2001 and 2011, and the vegetation composition and biomass in these years is assumed to have influenced fires during two periods, from 2002-2011 and 2012-2021, respectively. 
Unlike the vegetation data, the climate data is assumed to be seasonal or annual, and is matched to the year each fire occurred. 
The same climate variable(s) must be used to model fire ignition and escape, but can vary for spread.

##### Fuel classes:

This module can optionally estimate fuel classes for both forested and non-forested areas.
The non-forested mechanism is simpler: it combines landcover classes based on having a similar historical propensity to burn as estimated using a GLM, reducing an arbitrary amount of classes to a target amount (by default, two).
For example, bryoids and non-treed wetland might be combined in one class, with shrubland and grassland in another.

The estimation of forested fuel classes is more complex. 
The module will initially treat each species as a separate fuel class and then attempt to reduce the number of fuel classes to a target amount.
It uses a GLM to estimate the relationship between the biomass of each species, where present in pixel, to the likelihood of the likelihood of the pixel burning.
For each species, the model fitting data is subset by removing pixels where the species is not present. 
Additionally, each species will have an *a priori* fuel class.
The default fuel classes are an interpretation of the [FBP classes](https://cwfis.cfs.nrcan.gc.ca/background/fueltypes/c1) and can be viewed via `LandR::sppEquivalencies_CA$FuelClass`.

Individual tree species are grouped into fuel classes according to the following logic: first, the abundance of each species is assessed, calculated as the proportion of the landscape where a given species represents more than 10% of the biomass in a pixel.
"Rare" species (those with abundance \< 5%) are combined if another species is present with the same fuel class.
If there are two or more options, they are grouped in order of abundance, from least to most.
It is possible the resulting number of fuel classes will be fewer than the target amount if multiple species were "rare".
If the resulting number of unique fuel classes is greater than the target amount, the next least abundant species will be combined.
However, for non-rare species to be combined, they must also have the same coefficient sign, or a non-significant sign.
If there are more than one species with a similar fuel class and coefficient sign, they are grouped based on abundance, from least to most.
By default, the module will aim for five forested fuel classes, two non-forested fuel classes, and `youngAge` (described below), for a total of eight.

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

##### Special cases: youngAge, missingLCC, and treed wetland

`fireSense` uses a `youngAge` class to distinguish areas that have recently burned (by default within the past 15 years), and updates the landscape snapshots accordingly.
For example, if a conifer pixel burned in 2009, subsequent observations of the same pixel are converted to `youngAge`, with prior observations unchanged.
Additionally, the 2011-2011 landscape snapshot is adjusted to ensure areas are correctly classified if they burned prior to 2001 within the period relevant to the `youngAge` class.
By default, both forested and non-forested classes will convert to `youngAge` following fire, however a user can disable this behaviour in the latter class via the parameter `nonForestCanBeYoungAge`.

A related consideration is whether treed wetland (a land cover class in the NTEMS land cover product used by default) should be treated as a variant of forested landcover or non-forest cover.
This decision determines whether the biomass of the tree species in treed wetland is included as a covariate in the fireSense models, the same as non-wetland forest species, or whether treed wetland is instead considered as a land cover class, similar to shrubland and grassland, for example.
If wetland tree species are distinct from non-wetland species, then the decision may have less consequence.
The object `nonForestedLCCGroups` takes precedence when determining fuel classes; if the land cover value for non-treed wetland is included here, then it will be counted as non-forest regardless of whether it appears the parameter `forestedLCCclasses`.
This approach allows treed wetland to be simulated in forest simulation models (i.e., with LandR) without it being explicitly treated as forest in `fireSense`.
If treed wetland is absent from both `nonForestedLCCGroups` and `forestedLCCclasses`, then it is considered non-flammable land cover, as with any other land cover class.

If `fireSense_dataPrepFit` is to estimate fuel classes, it will also do so for a special class termed `missingLCC`.
This class encompasses pixels that have a forested land cover class but are absent in the objects `cohortData2011` and `cohortData2001`.
The most frequent reason for the exclusion of forest pixels in LandR is that there was no tree species data for the pixel.
Therefore, with no available species biomass, these pixels are treated as non-forest land cover, and must be assigned to one of the groups in `nonForestedLCCGroups`.

Finally, for fuel classes that are modelled with biomass instead of presence/absence, the log of biomass is used, after converting all instances of zero biomass for a fuelclass to one log below the minimum observed level (by default, 100 $g/m2).

### Module inputs and parameters

Table \@ref(tab:moduleInputs-fireSense-dataPrepFit) shows the full list of module inputs.

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
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
   <td style="text-align:left;"> climateVariablesForFire </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> A list detailing which climate variables in `sim$historicalClimateRasters` to use for which fire processes (ignition and spread). If the list is length one, both processes will use the same variables. The default is to use 'MDC'. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cohortData2001 </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table that defines the cohorts by pixelGroup in 2001 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> cohortData2011 </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> Table that defines the cohorts by pixelGroup in 2011 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFirePoints </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> named list of spatial points for each fire year with each point denoting an ignition location. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> firePolys </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> List of sf polygon objects representing annual fire polygons.List must be named with followign convention: `year&lt;numeric year&gt;` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> firePolysForAge </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list of fire polygons used to classify `timeSinceDisturbance` in nonforest LCC </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> historicalFireRaster </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> a raster with values representing fire year 1985-2020 </td>
   <td style="text-align:left;"> https://opendata.nfis.org/downloads/forest_change/CA_Forest_Fire_1985-2020.zip </td>
  </tr>
  <tr>
   <td style="text-align:left;"> historicalClimateRasters </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> length-one list of containing a raster stack of historical climate list named after the variable and raster layers named as `year&lt;numeric year&gt;` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionFirePoints </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> list of sf polygon objects representing annual ignition locations. This includes all fires regardless of size. It should have the same CRS as `sim$rasterToMatch` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> missingLCCgroup </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> if a pixel is forested but is absent from `cohortData`, it will be grouped in this class. It can be estimated if `P(sim)$estimateFuelClasses` is TRUE. If supplied, it must be one of the names in `sim$nonForestedLCCGroups` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForestedLCCGroups </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> a named list of non-forested landcover groups, e.g. list('wetland' = c(19, 23, 32)) These will become fuel covariates, and the groups will be estimated if `P(sim)$estimateFuelClasses` is TRUE </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pixelGroupMap2001 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> defines the `pixelGroups` for cohortData table in 2001 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> pixelGroupMap2011 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> defines the `pixelGroups` for cohortData table in 2011 </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rasterToMatch </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> template raster for study area. Assumes some buffering of core area to limit edge effect of fire. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstLCC2001 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Raster of land cover - will use Biomass_borealDataPrep to generate if missing. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> rstLCC2011 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> Raster of land cover - will use Biomass_borealDataPrep to generate if missing. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> sppEquiv </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> table of LandR species equivalencies </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> standAgeMap2001 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> map of stand age in 2001 used to create `cohortData2001` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> standAgeMap2011 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> map of stand age in 2011 used to create `cohortData2011` </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyArea </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> study area that determines spatial boundaries of all data. Should be buffered to accomodate edge effects </td>
   <td style="text-align:left;"> NA </td>
  </tr>
  <tr>
   <td style="text-align:left;"> studyAreaReporting </td>
   <td style="text-align:left;"> sf </td>
   <td style="text-align:left;"> (optional) study area used for reporting purposes, specifically whether fires inside the studyAreaReporting polygon are being removed for also falling partially outside the studyArea polygon, indicating the buffered studyArea shoudl be expanded. </td>
   <td style="text-align:left;"> NA </td>
  </tr>
</tbody>
</table>

Provide a summary of user-visible parameters (Table \@ref(tab:moduleParams-fireSense-dataPrepFit))

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
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
   <td style="text-align:left;"> Either a scalar that will buffer `areaMultiplier fireSize` or a quoted function of `fireSize`. See `?fireSenseUtils::bufferToArea`. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> bufferForFireRaster </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 1000 </td>
   <td style="text-align:left;"> 0 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> The distance that determine whether separate patches of burned pixels originated from the same fire. Only relevant when `useRasterizedFireForSpread = TRUE`. This param is separate from `minBufferSize`, which is used to determine the minimum sample of burned and unburned pixels to include in each fire. </td>
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
   <td style="text-align:left;"> estimateFuelClasses </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> estimate fuel classes from combination of data and P(sim)$fuelClassCol? </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireYears </td>
   <td style="text-align:left;"> integer </td>
   <td style="text-align:left;"> 2002, 20.... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> A numeric vector indicating which years should be extracted from the fire databases to use for fitting. Should not include years prior to 2002, to ensure correct intialization from data. </td>
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
   <td style="text-align:left;"> 40 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> aggregation factor for rasters during ignition prep. </td>
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
   <td style="text-align:left;"> Minimum number of cells in buffer and nonbuffer. This is imposed after the multiplier on the `bufferToArea` fn </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonflammableLCC </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 20, 31, .... </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> non-flammable LCC in rstLCC layers - defaulting to water, snow/ice, rock, barren land. </td>
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
   <td style="text-align:left;"> targetFuelClasses </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:left;"> the target number of unique fuel classes when using semi-automated approach </td>
  </tr>
  <tr>
   <td style="text-align:left;"> useCentroids </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> TRUE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should fire ignitions start at the `sim$firePolygons` centroids or at the ignition points in `sim$firePoints`? </td>
  </tr>
  <tr>
   <td style="text-align:left;"> useRasterizedFireForSpread </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should rasterized fire be used in place of a vectorized fire dataset? This method attributes burned pixels to specific fires, only examines the latest fire in a pixel, and may be subject to temporal error. is therefore more appropriate in areas with low rates of fire, or where the NFDB dataset may be incomplete (e.g., northern Ontario). </td>
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
   <td style="text-align:left;"> .plotInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Describes the simulation time at which the first plot event should occur. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .plotInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Describes the simulation time interval between plot events. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInitialTime </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Describes the simulation time at which the first save event should occur. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .saveInterval </td>
   <td style="text-align:left;"> numeric </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> This describes the simulation time interval between save events. </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .studyAreaName </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> `studyArea` name that will be appended to file-backed rasters </td>
  </tr>
  <tr>
   <td style="text-align:left;"> .useCache </td>
   <td style="text-align:left;"> logical </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Should this entire module be run with caching activated? This is intended for data-type modules, where stochasticity and time are not relevant </td>
  </tr>
</tbody>
</table>

### Events

##### init

During this event, objects are created and prepared if they are used by all three fire processes, such as the `flammableMap`, `landcoverDT`, and `nonForest_timeSinceDisturbance` for each fire period.
This is also when fuel class estimation occurs.

##### prepIgnitionFitData

This event prepares ignition data.
The pixels of the individual landscape snapshots of 2002-2011 and 2012-2020 are converted to voxels by aggregating the fuel class rasters according to the parameter `P(sim)$igAggFactor` and calculating the mean cover or biomass within.
Ignitions are summed within voxels for every year, such that the data will ultimately contain one row for each voxel and year, with the number of ignitions, the mean annual climate within that voxel, and the mean amount of each fuel class (whether biomass or cover) as additional columns.
Because the probability estimates of the fitted model can only be interpreted at the resolution of the fitted data, a template form of the aggregated raster is output as `sim$ignitionFitRTM`, and contains the attribute `nonNAs`, the number of pixels in the aggregated raster with non-NA data.  

##### prepEscapeFitData

This rather simple event adds an "escapes" column to the data table prepared by `prepIgnitionFitData`, where an escape is an ignition that resulted in a fire greater than the resolution of the `flammableMap` raster (which is *not* the area of the voxel).
It also builds the formula object if it hasn't been supplied, and performs some minor data quality control (e.g., ensuring that fire escapes do not outnumber ignitions in any voxel).

##### prepSpreadFitData

The data used to fit the fire spread model differs from that of ignition and escape in several ways.
First, the Fire Spread data includes fires with anthropogenic origins, whereas by default the Ignition dataset only includes those of natural origin.
Second, it utilizes a spatial subset of the landscape by buffering the fire polygons and extracting the climate and landscape data within.
This massively reduces the spread dataset by removing pixels that are of no interest to the model calibration, i.e., due to lack of observed fires.
The size of the buffer is dependent upon the parameter `P(sim)$areaMultiplier`. 

Second, each fire must contain a corresponding point of ignition that is located inside the fire polygon in a flammable pixel.
This object is called `spreadFirePoints` to differentiate it from the similar object used by the Ignition process.
The latter is different in that it contains fires of all sizes (i.e., smaller than the resolution of `sim$rasterToMatch`), does not have a corresponding polygon object, and is limited to fires of natural origin.
If `spreadFirePoints` is not supplied, the default behaviour is to calculate it from the centroids of the Fire Spread data and select the nearest flammable pixel inside the fire polygon in the event the centroid does not satisfy both of these criteria.

Lastly, the `youngAge` category is calculated for each year individually.
For example, by default if in 2001 a pixel was estimated as 12 years removed from the last disturbance in 2001, if the same pixel fell inside a fire buffer from 2006, it would no longer be youngAge, as 12 + (2006-2001) would exceed the parameter P(sim)$youngeAgeThreshold.
Similar adjustments are made for pixels that were initially non-youngAge, burned in an earlier year, and then re-appear in the buffer of later fires.
By contrast, the Ignition Fire data does not annually resolve the `youngAge class. 

As in `prepareIgnitionFit` and `prepareEscapeFit`, the model formula object will be created if it hasn't been supplied by a user.
Lastly, to reduce object size to a bare minimum for the sake of the lengthy optimization procedure employed in spreadFit, the data is organized into the annually-varying and non-annual components.
The former include the climate variables and `youngAge` class, and are organized into a list of `data.tables`, named by year, containing the aforementioned covariates, whether the pixel was burned or in the buffer, the fire ID to which it belongs, and the year of the fire.
The non-annual components are sorted into a list of two data.tables, one for each fire period, the names of which include every year (e.g `sim$nonAnnual_SpreadFitCovariates$Year2001_Year2002_Year2003...`).  

### Module outputs

Description of the module outputs (Table \@ref(tab:moduleOutputs-fireSense-dataPrepFit)).

<table class="table" style="color: black; margin-left: auto; margin-right: auto;">
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
   <td style="text-align:left;"> fireBufferedListDT </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list of data.tables with fire id, `pixelID`, and buffer status </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fuelClassTable </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> table with assigned fuel class of each tree species, after running assessFuelClasses </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFirePolys </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list of sf polygon objects representing annual fires </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_annualSpreadFitCovariates </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list of tables with climate covariates, `youngAge`, burn status, `polyID`, and `pixelID` </td>
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
   <td style="text-align:left;"> fireSense_ignitionFormula </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> formula for ignition, using climate and vegetation covariates, as character </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_nonAnnualSpreadFitCovariates </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> list of two tables with vegetation covariates, burn status, polyID, and `pixelID` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> fireSense_spreadFormula </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> formula for spread, using climate and vegetation covariates, as character </td>
  </tr>
  <tr>
   <td style="text-align:left;"> ignitionFitRTM </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> A (template) raster with information with regards to the spatial resolution and geographical extent of `fireSense_ignitionCovariates`. Used to pass this information onto `fireSense_ignitionFitted` Needs to have number of non-NA cells as attribute (`attributes(ignitionFitRTM)$nonNAs`). </td>
  </tr>
  <tr>
   <td style="text-align:left;"> landcoverDT2001 </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> data.table with `pixelID` and relevant landcover classes for flammable pixels in 2001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> landcoverDT2011 </td>
   <td style="text-align:left;"> data.table </td>
   <td style="text-align:left;"> data.table with `pixelID` and relevant landcover classes for flammable pixels in 2011 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> missingLCCgroup </td>
   <td style="text-align:left;"> character </td>
   <td style="text-align:left;"> if estimating fuel classes, the nonforest class to assign forested pixels absent from `sim$cohortData` </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForestedLCCgroups </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> a named list of non-forested landcover groups forming distinct fuel classes e.g. list('wetland' = c(19, 23, 32)) </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForest_timeSinceDisturbance2001 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> time since burn for non-forested pixels in 2001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> nonForest_timeSinceDisturbance2011 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> time since burn for non-forested pixels in 2011 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableRTM2001 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> binary raster of flammable landcover for 2001 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> flammableRTM2011 </td>
   <td style="text-align:left;"> SpatRaster </td>
   <td style="text-align:left;"> binary raster of flammable landcover for 2011 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> spreadFirePoints </td>
   <td style="text-align:left;"> list </td>
   <td style="text-align:left;"> Named list of `sf` polygon objects representing annual fire centroids. This only includes fires that escaped (e.g. `size &gt; res(flammableRTM)`. </td>
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
