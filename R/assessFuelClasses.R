# for a future example when moved to package
# tempDF <- data.table(species =c ("Pice_mar", "Pinu_con", "Popu_tre", "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"),
#                      coef = c(-0.0132, 0.00154, -0.021, -0.0017, -0.0191, -0.0092, -0.0063),
#                      sign = c("negative", "positive", "negative", "negative", "negative", "negative", "negative"),
#                      FuelClass = c("BlkSprc", "LdJkPine", "PopBrch", "PopBrch", "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"),
#                      above10PctRelB = c(0.066, 0.625, 0.223, 0.051, 0.095, 0.352, 0.560))
# tempDF2 <- data.table(species =c ("Pice_mar", "Pinu_con", "Popu_tre", "Betu_pap", "Pice_eng", "Pice_gla", "Abie_las"),
#                       coef = c(-0.0132, 0.00154, 0.021, -0.0017, -0.0191, -0.0092, -0.0063),
#                       sign = c("negative", "positive", "negative", "positive", "negative", "negative", "negative"),
#                       FuelClass = c("BlkSprc", "LdJkPine", "PopBrch", "PopBrch", "SprcFrLrch", "SprcFrLrch", "SprcFrLrch"),
#                       above10PctRelB = c(0.066, 0.625, 0.223, 0.051, 0.095, 0.352, 0.560))


assessFuelClasses <- function(landscape, fuelCol, thresholds = c(0.05, 0.1, 0.15),
                              sppEquiv = sim$sppEquiv, sppEquivCol) {

  landscape$B_MgHa <- landscape$B/100
  makeGLM <- function(species, landscape) {
    landscape <- landscape[speciesCode == species,]
    out <- glm(data = landscape, formula = burned ~ B_MgHa, family = binomial(link = "logit"))
  }

  actualSpecies <- unique(landscape[!is.na(B),]$speciesCode)
  fuelGLMs <- lapply(actualSpecies, makeGLM, landscape = landscape)
  names(fuelGLMs) <- actualSpecies
  coeffs <- sapply(fuelGLMs, coefficients)
  pvalues <- sapply(fuelGLMs, FUN = function(x){summary(x)$coefficients[2,4]})
  speciesStats <- data.table(species = names(fuelGLMs), coef = coeffs[2,], pvalue = pvalues)


  #diagnostic plots
  Bmeans <- landscape[speciesCode %in% actualSpecies, .(meanB = as.integer(mean(B_MgHa)),
                                                        maxB = as.integer(max(B_MgHa))),
                      .(speciesCode)]
  Bmeans <- melt(Bmeans, id.vars = "speciesCode", value.name = "B_MgHa")
  preds <- lapply(actualSpecies, function(spec, GLM = fuelGLMs,
                                          toJoin = Bmeans) {

    GLM <- GLM[[spec]]
    df <- data.table(B_MgHa = seq(0, max(toJoin$B_MgHa), 2),
                     speciesCode = spec)
    #make sure mean and max are in pred frame, for plotting
    toJoin <- toJoin[speciesCode %in% spec, .(B_MgHa, speciesCode)]
    df <- rbind(toJoin, df, fill = TRUE)
    df[, burnprob := predict(GLM, newdata = df)]
    return(df)
  })
  preds <- rbindlist(preds)
  predsOfNote <- preds[Bmeans, on = c("speciesCode", "B_MgHa")]
  #TODO: decide whether this is helpful
  # ggplot(preds, aes(x = B_MgHa, y = burnprob, col = speciesCode)) +
  #          geom_line() +
  #          labs(x = "B (Mg/ha)", y = "burn pred") +
  #   geom_point(data = predsOfNote, aes(x = B_MgHa, y = burnprob, shape = variable))

  ####TODO: figure out what to do next - probably plot these predictions by biomass.
  forest <- landscape[!is.na(B),]
  fuels <-  sppEquiv[, .SD, .SDcols = c(sppEquivCol, fuelCol)]

  setnames(fuels, old = sppEquivCol, new = "speciesCode")
  forest <- fuels[forest, on = c("speciesCode")]
  forest[, propB := B/totalBiomass]
  forestPix <- nrow(forest[, .N, .(cell, year)]) #this is the only way because cells change between years
  selectCols <- c("speciesCode", fuelCol)
  above10PctRelB <- forest[propB > 0.1, .(above10PctRelB = .N/forestPix), by = c("speciesCode", fuelCol)]
  speciesStats[, sig := pvalue < 0.001]
  speciesStats[, sign := ifelse(coef > 0, "positive", "negative")]
  speciesStats[sig == FALSE, sign := "neutral"]
  speciesStats <- speciesStats[, .(species, coef, pvalue, sign)][above10PctRelB, on = c("species" = "speciesCode")]
  #   1. start by identifying species that have <5% of pixels with >10% biomass for merging
  #   2. Only species that are named in same FBPS class can be merged (we need to look at that now)... e.g., Pinu_con and Pinu_ban
  #   3. Ensure these are fitting within "negatives", "positives", and "neutrals"
  #   So, totally unrelated ones can't be merged, but otherwise allows merging of e.g. pines
  #   These would only be merged if they fulfill Criterion 1 above. Otherwise merged.
  #   if there is a small amount of some species left over after all merging, then put it in one of:
  #   "other positive", "other negative",
  #   If the species is allowed to be merged based on the Potential Fuel Class column,
  #   Revisit this classification if there are >5 classes, starting with higher thresholds for "Criterion 1"

    # Function to assign NewFuelClass

  out <- combine_fuel_classes(df = speciesStats)
  return(out)
}

combine_fuel_classes <- function(df, thresholds = c(0.05, 0.1, 0.15)) {


  # Initialize the new assignedFuelClass column
  df$assignedFuelClass <- df$species
  if (length(unique(df$species) > 5)) {
    df[, alreadyMerged := FALSE]

    for (multiplier in c(1, 2, 4)) { #this will increase thresholds
      toMerge <- 0
      count <- 1
      neededMerge <- length(unique(df$species)) - 5

      #for simplicity, remove those that won't merge
      guaranteedUnique <- df[, .N, .(FuelClass)]
      guaranteedUnique <- guaranteedUnique[N == 1]
      guaranteedUnique <- df[FuelClass %in% guaranteedUnique$FuelClass,]
      #TODO: Ignore if it is unique but under some threshold of above10PctRelB?

      possMerge <- df[!species %in% guaranteedUnique$species]
      while (count < length(thresholds) | toMerge < neededMerge) {
        #start by identifying fuelClasses below the biomass threshold
        BelowThresh <- possMerge[above10PctRelB <= thresholds[count] * multiplier,]
        toMerge <- nrow(BelowThresh)
        count <- count + 1
      }

      #now, iteratively join BelowThresh speciies to ones in possMerge
      #track which ones are merged and edit assigned fuel class

      #this resets the matching at each threshold - not sure if necessary
      matched <- possMerge[0]

      for (i in 1:length(BelowThresh$species)) {
        #this will iterate over one that is now joined
        toMerge <- BelowThresh[i,]
        if (toMerge$alreadyMerged == FALSE) {
          possMatches <- possMerge[FuelClass == toMerge$FuelClass &
                                     c(sign == toMerge$sign | sign == "neutral")
                                   & species != toMerge$species]
          #cannot merge with itself, otherwise fuelclass must be identical, and sign same or neutral (non-sig)
          if (nrow(possMatches) > 0) {
            if (nrow(possMatches) > 1) {
              #take nearest coefficient
              possMatches[, coef2 := abs(coef - toMerge$coef)]
              setkey(possMatches, coef2)
              possMatches <- possMatches[1,]
              possMatches[, coef2 := NULL]
            }
            matchedTo <- possMatches #guaranteed length 1
            #adjust for future species matching
            possMerge[species == matchedTo$species, assignedFuelClass := FuelClass]
            #adjust the matched species fuel class
            toMerge[, assignedFuelClass := matchedTo$FuelClass]
            #track them inside the for loop
            matched <- rbind(matched, toMerge)
            #avoid
            matched[, alreadyMerged := TRUE]
          }
        }
      }

      #get the ones that weren't matched, combine with matched and uniques
      unmerged <- possMerge[!species %in% matched$species]
      allSpecies <- rbind(guaranteedUnique, unmerged, matched)

      # Iterate through each unique FuelClass
      #correct for different signs being assigned to same FuelClass
      safetyCatch <- allSpecies[, .N, .(assignedFuelClass, sign)]
      duplicated <- safetyCatch[duplicated(safetyCatch$assignedFuelClass)]$assignedFuelClass
      if (length(duplicated) > 0) {
        allSpecies[assignedFuelClass %in% duplicated,
                   assignedFuelClass := paste0(assignedFuelClass, "_", sign)]
      }
      #the above new class name is only relevant if god forbid 4 species are in one group with 2 each opposing signs

      #keep name the same if there is only one
      allSpecies[, N := .N, .(assignedFuelClass)]
      allSpecies[N == 1, assignedFuelClass := species] # catches groups that started with multiple

      setkey(allSpecies, species)
      setkey(df, species)
      if (!identical(allSpecies$species, df$species)) {
        stop("species do not match - possible error in assessFuelClasses")
      }
      df <- allSpecies
      if (length(unique(df$assignedFuelClass)) <= 5) {
        break
      }
    }
    if (length(unique(df$assignedFuelClass)) > 5){
      #some combinations may never work depending on e.g GLM sign,
      warning("could not adequately resolve fuel classes - consider editing FuelClass")
    }
    df[, alreadyMerged := NULL]
  } else {
    df[, assignedFuelClass := species]
  }

  return(df)
}


