Sys.setenv(LANG = "en")
setwd("D:/.IMPERIAL/.70062_Research_Project/Datasets")

library(dplyr)
library(stringr)
library(sf)
library(ggplot2)
library(glmmTMB)
library(ggeffects)
library(sjPlot)

# ==============================================================================
# PART 1: GEOGRAPHIC FRAMEWORK & SPATIAL CHECKPOINT (CACHED)
# ==============================================================================
# This section stitches together England, Wales, Scotland, and Northern Ireland 
# into one unified UK map. Because shapefiles are massive, we use an if/else 
# statement to only build this once. If the .rds file exists, it skips the math.
# ==============================================================================
cat("\n--- CHECKING FOR SPATIAL MASTER ---\n")
spatial_rds <- "Boundaries/UK_Spatial_Master.rds"

if (file.exists(spatial_rds)) {
  uk_spatial_master <- readRDS(spatial_rds)
  cat("Loaded existing UK_Spatial_Master.rds (Skipping geometry build)\n")
} else {
  cat("Building Spatial Master from scratch... this may take a minute.\n")
  
  # A. Load and format population data for England & Wales
  ew_pop <- read.csv("Boundaries/Wards/Ward_pop_estimate.csv", sep=";", check.names = FALSE) %>%
    dplyr::mutate(Total_Pop = as.numeric(stringr::str_replace_all(Total, " ", ""))) %>%
    dplyr::select(WD23CD = `Ward 2023 Code`, Total_Pop)
  
  ew_boundaries <- read.csv("Boundaries/Wards/WD_MAY_2023_UK_BGC_OnlyEngland&Wales.csv", sep=";", check.names = FALSE)
  
  # Calculate baseline metrics: Area and Population Density per Ward
  ew_master <- ew_boundaries %>%
    dplyr::left_join(ew_pop, by = "WD23CD") %>%
    dplyr::filter(!is.na(Total_Pop)) %>%
    dplyr::mutate(
      ID = dplyr::row_number(),
      Area_SqKm = Shape__Area / 1000000,
      Pop_Density_SqKm = Total_Pop / Area_SqKm
    ) %>%
    dplyr::select(ID, Zone_Code = WD23CD, Zone_Name = WD23NM, Total_Pop, Area_SqKm, Pop_Density_SqKm)
  
  # B. Load Scotland and Northern Ireland equivalent census areas
  scot_data <- read.csv("Boundaries/SG_IntermediateZoneBdry_2022/SG_IntermediateZoneBdry_2022_MHW_ExportTable.csv", sep = ";", check.names = FALSE) %>%
    dplyr::select(Zone_Code = IZCode, Zone_Name = IZName, Total_Pop = TotPop2022, Area_SqKm = StdAreaKM2, Pop_Density_SqKm = `Population_density_km²`)
  
  ni_data <- read.csv("Boundaries/geography-sdz2021-esri-shapefile/SDZ2021_ExportTable.csv", sep = ";", check.names = FALSE) %>%
    dplyr::select(Zone_Code = SDZ2021_cd, Zone_Name = SDZ2021_nm, Total_Pop = All_usual_residents, Area_SqKm = Area_SqKm, Pop_Density_SqKm = `Population_density_km²`)
  
  uk_master <- dplyr::bind_rows(ew_master, scot_data, ni_data) %>%
    dplyr::mutate(ID = dplyr::row_number()) %>%
    dplyr::select(ID, Zone_Code, Zone_Name, Total_Pop, Area_SqKm, Pop_Density_SqKm)
  
  # C. Load the physical polygon geometries (Shapefiles)
  ew_sf <- sf::st_read("Boundaries/Wards/WD_MAY_2023_UK_BGC.shp", quiet = TRUE) %>%
    sf::st_transform(27700) %>% # 27700 is the British National Grid (meters)
    dplyr::select(Zone_Code = WD23CD) %>%
    dplyr::filter(grepl("^[EW]", Zone_Code))
  
  scot_sf <- sf::st_read("Boundaries/SG_IntermediateZoneBdry_2022/SG_IntermediateZoneBdry_2022_MHW.shp", quiet = TRUE) %>%
    sf::st_transform(27700) %>%
    dplyr::select(Zone_Code = IZCode)
  
  ni_sf <- sf::st_read("Boundaries/geography-sdz2021-esri-shapefile/SDZ2021.shp", quiet = TRUE) %>%
    sf::st_transform(27700) %>%
    dplyr::select(Zone_Code = SDZ2021_cd)
  
  # D. Bind all geometries and scale the population density for statistical modeling
  uk_spatial_master <- rbind(ew_sf, scot_sf, ni_sf) %>%
    dplyr::left_join(uk_master, by = "Zone_Code") %>%
    dplyr::mutate(Pop_Density_Z = scale(Pop_Density_SqKm)[,1]) %>%
    dplyr::filter(!is.na(Total_Pop))
  
  # Save the checkpoint so we never have to run this block again
  saveRDS(uk_spatial_master, spatial_rds)
  cat("Saved new UK_Spatial_Master.rds\n")
}

# ==============================================================================
# PART 2: BIOLOGICAL & CLIMATE DATA EXTRACTION
# ==============================================================================
# Here we map the literal X/Y coordinates of the beetle sightings and the Met 
# Office climate grids into our census polygons using spatial intersections.
# ==============================================================================
cat("\n--- PROCESSING BIOLOGICAL & CLIMATE DATA ---\n")

# Load points
temp_grid <- read.csv("ClimateData_MetOffice/Annual_Average_Temperature_Change.csv", sep=";")
grid_sf <- sf::st_as_sf(temp_grid, coords = c("square_lon", "square_lat"), crs = 4326) %>% 
  sf::st_transform(27700)

ca_pts <- read.csv("Chrysolina/Chrysolina.csv", sep=";") %>%
  sf::st_as_sf(coords = c("Chryso_lon", "Chryso_lat"), crs = 4326) %>% 
  sf::st_transform(27700)

pc_pts <- read.csv("Pyrochroa/Pyrochroa.csv", sep=";") %>%
  sf::st_as_sf(coords = c("Pyro_lon", "Pyro_lat"), crs = 4326) %>% 
  sf::st_transform(27700)

# Extract local temperatures: Assigns the nearest Met Office grid to each Ward
zone_centroids <- sf::st_centroid(uk_spatial_master)
zone_climate <- sf::st_join(zone_centroids, grid_sf, join = sf::st_nearest_feature) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(Zone_Code, Temp_Hist = baseline_1981_2000) %>% 
  dplyr::filter(!is.na(Temp_Hist)) %>%
  dplyr::mutate(Temp_Hist_Z = scale(Temp_Hist)[,1])

# Aggregate target species (C. americana) counts per Ward, per Year (1994-2021)
ca_zone <- sf::st_join(ca_pts, uk_spatial_master, join = sf::st_intersects) %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(Year = Chryso_year) %>%
  dplyr::filter(Year >= 1994 & Year <= 2021 & !is.na(Zone_Code)) %>%
  dplyr::group_by(Zone_Code, Year) %>%
  dplyr::summarise(Count_Ca = dplyr::n(), .groups = "drop")

# Aggregate calibration species (P. coccinea) to measure observer bias
pc_zone <- sf::st_join(pc_pts, uk_spatial_master, join = sf::st_intersects) %>%
  sf::st_drop_geometry() %>%
  dplyr::mutate(Year = Pyro_year) %>%
  dplyr::filter(Year >= 1994 & Year <= 2021 & !is.na(Zone_Code)) %>%
  dplyr::group_by(Zone_Code, Year) %>%
  dplyr::summarise(Count_Pc = dplyr::n(), .groups = "drop")


# ==============================================================================
# PART 3: MASTER PANEL & GIS SUMMARIES (AND MEMORY CLEANUP)
# ==============================================================================
# We must create a "zero-inflated" panel. If a ward wasn't visited in 1998, R 
# doesn't know it exists. Expand_grid forces R to acknowledge the zeroes (absences).
# ==============================================================================
cat("\n--- BUILDING MASTER PANELS & SUMMARIES ---\n")

active_zones <- unique(c(ca_zone$Zone_Code, pc_zone$Zone_Code))

# A. Longitudinal Panel for Statistical Models (Ward + Year)
master_panel <- tidyr::expand_grid(Zone_Code = active_zones, Year = 1994:2021) %>%
  dplyr::left_join(ca_zone, by = c("Zone_Code", "Year")) %>%
  dplyr::left_join(pc_zone, by = c("Zone_Code", "Year")) %>%
  dplyr::left_join(zone_climate, by = "Zone_Code") %>%
  dplyr::left_join(sf::st_drop_geometry(uk_spatial_master)[, c("Zone_Code", "Pop_Density_Z")], by = "Zone_Code") %>%
  dplyr::mutate(
    # Clean NAs into true ecological zeroes
    Count_Ca = ifelse(is.na(Count_Ca), 0, Count_Ca),
    Count_Pc = ifelse(is.na(Count_Pc), 0, Count_Pc),
    # Create Binary Presence/Absence for the Hurdle Model
    Presence_Ca = ifelse(Count_Ca > 0, 1, 0),
    Presence_Pc = ifelse(Count_Pc > 0, 1, 0),
    # Center the year so the model doesn't mathematically calibrate to Year 0 AD
    Year_Centered = Year - 1994
  ) %>%
  dplyr::filter(!is.na(Temp_Hist_Z) & !is.na(Pop_Density_Z))

# B. Zone Summaries for Final Map Export (Collapsed across all years)
ca_totals <- ca_zone %>%
  dplyr::group_by(Zone_Code) %>%
  dplyr::summarise(Total_Ca = sum(Count_Ca), .groups = "drop")

zone_gis_summary <- sf::st_drop_geometry(uk_spatial_master) %>%
  dplyr::left_join(ca_totals, by = "Zone_Code") %>%
  dplyr::mutate(
    Total_Ca = ifelse(is.na(Total_Ca), 0, Total_Ca),
    Density_Coefficient = Total_Ca / Pop_Density_SqKm # Observer bias metric
  ) %>%
  dplyr::select(Zone_Code, Total_Ca, Density_Coefficient)

# STREAMLINING: Free up RAM before running heavy statistics
# Removed zone_climate, ca_zone, and pc_zone from the kill list so Part 6 doesn't crash!
rm(ew_boundaries, ew_master, ew_pop, ew_sf, ni_data, ni_sf, scot_data, scot_sf, 
   ca_pts, pc_pts, temp_grid, grid_sf)
gc() # Garbage collection


# ==============================================================================
# PART 4: STATISTICAL ANALYSIS (THERMAL / PRESENCE MODEL)
# ==============================================================================
cat("\n--- CHECKING THERMAL MODEL ---\n")
thermal_rds <- "thermal_model.rds"

if (file.exists(thermal_rds)) {
  thermal_model <- readRDS(thermal_rds)
  cat("Loaded existing thermal_model.rds\n")
} else {
  cat("Calculating thermal model...\n")
  # A binomial GLM to determine the strict climate gateway for initial establishment
  thermal_model <- glm(
    Presence_Ca ~ Temp_Hist_Z + Presence_Pc, 
    family = binomial(link = "logit"), 
    data = master_panel
  )
  saveRDS(thermal_model, thermal_rds)
}
print(summary(thermal_model))


# ==============================================================================
# PART 5: ABUNDANCE MODEL & INTERACTION PLOT
# ==============================================================================
cat("\n--- CHECKING ABUNDANCE MODEL & PREDICTIONS ---\n")
abundance_rds <- "abundance_model_tmb.rds"
predictions_rds <- "abundance_predictions_tmb.rds"

if (file.exists(abundance_rds)) {
  abundance_model_tmb <- readRDS(abundance_rds)
  cat("Loaded existing abundance_model_tmb.rds\n")
} else {
  cat("Calculating abundance model (this is heavy, please wait)...\n")
  # Negative Binomial (nbinom2) handles overdispersed ecological counts.
  # (1 | Zone_Code) is the random intercept controlling for spatial clustering.
  abundance_model_tmb <- glmmTMB::glmmTMB(
    Count_Ca ~ Year_Centered * Pop_Density_Z + Temp_Hist_Z + Count_Pc + (1 | Zone_Code),
    data = master_panel, family = nbinom2
  )
  saveRDS(abundance_model_tmb, abundance_rds)
}
print(summary(abundance_model_tmb))

# ==============================================================================
# H3: DENSITY-DEPENDENT ABUNDANCE PLOT (UPDATED)
# ==============================================================================

# Cache the plotting math with an explicit condition for zero observer bias
if (file.exists(predictions_rds)) {
  predictions_tmb <- readRDS(predictions_rds)
} else {
  predictions_tmb <- ggeffects::ggpredict(
    abundance_model_tmb, 
    terms = c("Year_Centered", "Pop_Density_Z [-1, 0, 2]"),
    condition = c(Count_Pc = 0) # Forces the plot to assume zero observer bias
  )
  saveRDS(predictions_tmb, predictions_rds)
}

# 1. Create a label dataframe with fixed coordinates for the bottom right
label_data <- as.data.frame(predictions_tmb) %>%
  dplyr::group_by(group) %>%
  dplyr::filter(x == min(x) | x == max(x)) %>%
  dplyr::arrange(x) %>%
  dplyr::summarize(
    slope = (log10(dplyr::last(predicted)) - log10(dplyr::first(predicted))) / (dplyr::last(x) - dplyr::first(x)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    group_letter = dplyr::case_when(
      group == "-1" ~ "(a)",
      group == "0"  ~ "(b)",
      group == "2"  ~ "(c)"
    ),
    label_text = sprintf("%s Slope: %.3f", group_letter, slope),
    # Set fixed X and Y coordinates to place the text in the bottom right corner
    # Y coordinates are spaced logarithmically to appear evenly distributed
    text_x = 2020, 
    text_y = dplyr::case_when(
      group == "-1" ~ 0.018,   # Top label
      group == "0"  ~ 0.006,   # Middle label
      group == "2"  ~ 0.002    # Bottom label
    )
  )

# 2. Build the updated plot
abundance_plot <- ggplot2::ggplot(predictions_tmb, ggplot2::aes(x = x + 1994, y = predicted, color = group)) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.1, color = NA) +
  
  # Add the generated text labels inside the plot space
  ggplot2::geom_text(
    data = label_data, 
    ggplot2::aes(x = text_x, y = text_y, label = label_text),
    hjust = 1,          # Right-aligns the text against the text_x coordinate
    show.legend = FALSE, 
    fontface = "bold",
    size = 4
  ) +
  
  # Updated legend labels with the specific SD and National Average notation
  ggplot2::scale_color_manual(
    values = c("#3182bd", "#756bb1", "#de2d26"), 
    labels = c("(a) Rural (-1 SD)", "(b) Suburban (0 SD, National Average)", "(c) Urban (+2 SD)")
  ) +
  ggplot2::scale_fill_manual(
    values = c("#3182bd", "#756bb1", "#de2d26"), 
    guide = "none"
  ) +
  
  # Restored standard X-axis limits since the text is now safely inside the plot
  ggplot2::scale_x_continuous(
    breaks = seq(1995, 2020, by = 5),
    limits = c(1994, 2023) 
  ) +
  
  ggplot2::scale_y_log10() + 
  
  # Updated Y-axis title
  ggplot2::labs(
    title = "Density-Dependent Expansion of C. americana (1994–2021)", 
    subtitle = "Predicted abundance index highlighting rates of population increase across urbanisation levels",
    x = "Year", 
    y = "Predicted Annual Abundance Index (Log10)", 
    color = "Census Area Classification:"
  ) +
  
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    legend.position = "bottom"
  )

print(abundance_plot)


# ==============================================================================
# PART 6: FUTURE CLIMATE PROJECTIONS
# ==============================================================================
cat("\n--- GENERATING FUTURE CLIMATE PREDICTIONS ---\n")

# Re-extract the original climate means to ensure standardized future math
temp_hist_vector <- sf::st_drop_geometry(uk_spatial_master) %>%
  dplyr::left_join(zone_climate, by = "Zone_Code") %>% dplyr::pull(Temp_Hist)

temp_mean <- mean(temp_hist_vector, na.rm = TRUE)
temp_sd   <- sd(temp_hist_vector, na.rm = TRUE)

# Base template for forecasting (Holds Year at 27 (2021) and Observer Bias at 0)
base_pred_data <- uk_spatial_master %>%
  sf::st_drop_geometry() %>%
  dplyr::left_join(zone_climate, by = "Zone_Code") %>%
  dplyr::mutate(Year_Centered = 27, Count_Pc = 0)

future_predictions_table <- base_pred_data %>% dplyr::select(Zone_Code)

# Predict Baseline
future_predictions_table$Pred_Abund_Baseline <- predict(abundance_model_tmb, newdata = base_pred_data, type = "response", re.form = NA)

# STREAMLINING: Loop through scenarios instead of repeating code blocks
scenarios <- c(1.5, 2.5, 4.0)
for (temp_inc in scenarios) {
  # sprintf("%.1f") forces R to keep the trailing .0 for the 4.0 scenario
  formatted_temp <- sprintf("%.1f", temp_inc) 
  col_name <- paste0("Pred_Abund_", gsub("\\.", "", formatted_temp), "C")
  
  # Standardize the artificially warmed temperature
  temp_scenario_data <- base_pred_data %>%
    dplyr::mutate(Temp_Hist_Z = ((Temp_Hist + temp_inc) - temp_mean) / temp_sd)
  
  # Predict and attach to table
  future_predictions_table[[col_name]] <- predict(abundance_model_tmb, newdata = temp_scenario_data, type = "response", re.form = NA)
}


# ==============================================================================
# PART 7: MASTER UNIFIED SPATIAL EXPORT
# ==============================================================================
cat("\n--- CHECKING GIS EXPORT ---\n")
map_rds <- "Boundaries/Future_Projections_Data.rds"

if (file.exists(map_rds)) {
  final_master_map <- readRDS(map_rds)
  cat("Loaded existing Future_Projections_Data.rds. (Delete file to regenerate GeoPackage)\n")
} else {
  cat("Building and exporting final GeoPackage to ArcGIS...\n")
  
  # Combine geometries, historical summaries, and future predictions
  final_master_map <- uk_spatial_master %>%
    dplyr::left_join(zone_gis_summary, by = "Zone_Code") %>%
    dplyr::left_join(future_predictions_table, by = "Zone_Code") %>%
    dplyr::select(Zone_Code, Zone_Name, Total_Pop, Pop_Density_SqKm, Pop_Density_Z, 
                  Total_Ca, Density_Coefficient, Pred_Abund_Baseline, 
                  Pred_Abund_15C, Pred_Abund_25C, Pred_Abund_40C)
  
  # Dummy scale locks ArcGIS symbology color ramps across multiple layers
  final_master_map$Dummy_Scale <- seq(0, 20, length.out = nrow(final_master_map))
  
  # Export map to the Outputs directory
  sf::st_write(final_master_map, "../Outputs/Chrysolina_Invasion_Map.gpkg", layer = "main", append = FALSE)
  saveRDS(final_master_map, map_rds)
  cat("Success! Combined spatial layer written to GeoPackage layer 'main' in the Outputs folder.\n")
}

# ==============================================================================
# PART 8: EXPORTING PUBLICATION-READY TABLES
# ==============================================================================
cat("\n--- EXPORTING TABLES ---\n")

# Export to the Outputs directory
sjPlot::tab_model(thermal_model, 
                  digits = 3, 
                  file = "../Outputs/Table_1_Thermal_Model.html",
                  title = "Table 1: Logistic Hurdle Model (Presence Probability)")

sjPlot::tab_model(abundance_model_tmb, 
                  digits = 3,                 
                  file = "../Outputs/Table_2_Abundance_Model.html",
                  title = "Table 2: Negative Binomial Count Model (Abundance Index)")

cat("\n--- PIPELINE COMPLETE! ---\n")
