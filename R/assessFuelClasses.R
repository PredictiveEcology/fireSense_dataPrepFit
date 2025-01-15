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

combine_fuel_classes <- function(df) {
  # Sort by FuelClass, sign, and above10 in ascending order for consistent processing
  df <- df[order(df$FuelClass, df$sign, df$above10), ]

  # Initialize the new assignedFuelClass column
  df$assignedFuelClass <- df$FuelClass

  # Group counter to number new FuelClasses only when needed
  group_counters <- list()

  # Iterate through each unique FuelClass
  for (fuel_class in unique(df$FuelClass)) {
    # Filter species within the same FuelClass
    subset_df <- df[df$FuelClass == fuel_class, ]

    # If there's only one species in the FuelClass, use the species name
    if (nrow(subset_df) == 1) {
      df$assignedFuelClass[df$species == subset_df$species] <- subset_df$species
      next
    }

    # Process positive and negative signs separately
    for (sign_group in c("positive", "negative")) {
      species_subset <- subset_df[subset_df$sign %in% c(sign_group, "neutral"), ]

      # While there are ungrouped species in the subset
      while (nrow(species_subset) > 0) {
        # Pick the species with the smallest above10 value
        seed_species <- species_subset[which.min(species_subset$above10), ]

        # Calculate the absolute difference in coefficients
        species_subset$coef_diff <- abs(species_subset$coef - seed_species$coef)

        # Exclude the seed_species itself to find the closest match
        species_subset <- species_subset[order(species_subset$coef_diff), ]

        if (nrow(species_subset) > 1) {
          # Merge with the closest match
          closest_species <- species_subset[2, ]  # Second row is the closest match

          # Determine if numbering is necessary
          if (!fuel_class %in% names(group_counters)) {
            group_counters[[fuel_class]] <- 1
          } else {
            group_counters[[fuel_class]] <- group_counters[[fuel_class]] + 1
          }
          new_class <- paste0(fuel_class, "_", group_counters[[fuel_class]])

          # Assign the new FuelClass to both species
          df$assignedFuelClass[df$species %in% c(seed_species$species, closest_species$species)] <- new_class

          # Remove these species from the subset
          species_subset <- species_subset[!species_subset$species %in% c(seed_species$species, closest_species$species), ]
        } else {
          # Assign the original FuelClass for unmerged species
          df$assignedFuelClass[df$species == seed_species$species] <- fuel_class

          # Remove the seed_species from the subset
          species_subset <- species_subset[species_subset$species != seed_species$species, ]
        }
      }
    }
  }

 df[, N := .N, .(assignedFuelClass)]
 df[N == 1, assignedFuelClass := species] # catches groups that started with multiple

  return(df)
}


