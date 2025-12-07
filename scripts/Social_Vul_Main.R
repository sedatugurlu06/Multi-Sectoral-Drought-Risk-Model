# ==============================================================================
# PROJECT: DROUGHT RISK MODEL - SOCIAL VULNERABILITY COMPONENT
# AUTHOR: [SEDAT UGURLU]
# DATE: 2025-12-06
# ==============================================================================

# 1. AUTO-INSTALL & LOAD LIBRARIES
# ------------------------------------------------------------------------------
# Define the list of required packages
required_packages <- c(
  "sf",           # Spatial vector data
  "tidyverse",    # Data manipulation and visualization
  "leaflet",      # Interactive maps
  "giscoR",       # Official Eurostat boundary data
  "countrycode",  # Convert country names to 4codes
  "htmlwidgets"   # Save interactive maps as HTML
)

# Loop through the list: Install if missing, then load
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("Installing missing package:", pkg))
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  } else {
    message(paste("Loaded package:", pkg))
  }
}

# Set Working Directory (Change this to your actual path)
# setwd("D:/Thesis/Code/R_Language/Gender_Equality")

# ==============================================================================
# 2. CONFIGURATION BLOCK
# ==============================================================================
# This block controls the entire model. Change weights or parameters here.

CONFIG <- list(
  # 1. Gender Equality (Country Level)
  gender = list(
    name = "Gender Equality",
    file = "../data/cleaned_gender_final.csv",
    join_col = "ISO_CODE",
    value_col = "Gender_Score",
    level = "Country",
    weight = 0.15,
    s_curve = "decreasing",         # High Gender Eq = Low Vulnerability
    steepness = 0.2
  ),
  
  # 2. Rural Population (NUTS 2)
  rural = list(
    name = "Rural Population",
    file = "../data/clean_rural_pop_percentage.csv",
    join_col = "geo",
    value_col = "Avg_Rural_Pop_Pct",
    level = "NUTS2",
    weight = 0.30,
    s_curve = "increasing",         # High Rural % = High Vulnerability
    steepness = 0.2
  ),
  
  # 3. AROPE - Poverty (NUTS 2)
  arope = list(
    name = "Poverty Risk (AROPE)",
    file = "../data/clean_poverty_risk_rate.csv",
    join_col = "geo",
    value_col = "Avg_Poverty_Rate",
    level = "NUTS2",
    weight = 0.15,
    s_curve = "increasing",         # High Poverty = High Vulnerability
    steepness = 0.2
  ),
  
  # 4. Social Dependency (NUTS 3 or 2)
  dependency = list(
    name = "Social Dependency",
    file = "../data/clean_social_dependency.csv",
    join_col = "geo",
    value_col = "Avg_Dependency_Ratio",
    level = "NUTS3",                # Using NUTS 3 for higher precision
    weight = 0.15,
    s_curve = "increasing",         # High Dependency = High Vulnerability
    steepness = 0.2
  ),
  
  # 5. HDI (Country)
  hdi = list(
    name = "Human Development Index",
    file = "../data/clean_hdi_aquastat.csv",
    join_col = "ISO_CODE",   
    value_col = "Avg_HDI_National",
    level = "Country",                
    weight = 0.25,
    s_curve = "decreasing",         # High HDI = Low Vulnerability
    steepness = 0.2
  )
)

# ==============================================================================
# 3. DEFINE S-CURVE FUNCTIONS (VECTORIZED FIX)
# ==============================================================================

# Function 1: Increasing S-Curve (For "Bad" indicators like Poverty)
s_curve_increasing <- function(x, b, k) {
  # If midpoint is 0, avoid division by zero error; return input
  if(is.na(b) || b == 0) return(x)
  
  # The formula works for vectors automatically
  # If x is 0, b/x becomes Inf, resulting in 0 vulnerability (Correct)
  return(1 / (1 + (b / x)^k)) 
}

# Function 2: Decreasing S-Curve (For "Good" indicators like HDI)
s_curve_decreasing <- function(x, b, k) {
  # Handle empty midpoint
  if(is.na(b)) return(x)
  
  # The explicit "if (x == 0)" check caused the error.
  # We remove it because the formula naturally handles 0:
  # (0 / b)^k = 0  -->  1 / (1 + 0) = 1 (Max Vulnerability)
  return(1 / (1 + (x / b)^k))
}

# Wrapper function
normalize_data <- function(values, type, steepness) {
  # Calculate average (midpoint) ignoring NAs
  midpoint <- mean(values, na.rm = TRUE) 
  
  # Safety check: If no data exists, return NAs
  if(is.nan(midpoint)) return(rep(NA, length(values)))
  
  if (type == "increasing") {
    return(s_curve_increasing(values, midpoint, steepness))
  } else {
    return(s_curve_decreasing(values, midpoint, steepness))
  }
}

# ==============================================================================
# 4. LOAD GEOMETRY & FILTER STUDY AREA (FIXED)
# ==============================================================================

print("Loading Basin Shapefile...")

# 1. DEFINE STUDY AREA (Europe minus RU/TR)
# We first get the list of country codes (CNTR_ID) for Europe
# gisco_get_countries SUPPORTS 'region', so we use it here.
europe_geo <- gisco_get_countries(year = "2020", region = "Europe")

# Filter out Russia (RU) and Turkey (TR) from this list
allowed_countries_geo <- europe_geo %>%
  filter(!CNTR_ID %in% c("RU", "TR")) 

# Extract the list of codes (e.g., "AT", "DE", "FR"...)
allowed_codes <- allowed_countries_geo$CNTR_ID

# 2. LOAD YOUR BASINS
# Make sure this points to your actual file
 basin_shp <- st_read("../shp/hybas_eu_lev06_v1c.shp")

# --- MOCK DATA (Corrected) ---
if (!exists("basin_shp")) {
  message("Using Mock Basin Data for Europe")
  # FIX: Use 'country = allowed_codes' instead of 'region = "Europe"'
  basin_shp <- gisco_get_nuts(year = "2021", nuts_level = 2, country = allowed_codes) %>% 
    st_centroid() %>% st_buffer(60000) %>% 
    select(NUTS_ID) %>% rename(HYBAS_ID = NUTS_ID)
}
# ---------------------------------------------------------------

# 3. STANDARDIZE CRS
basin_shp <- st_transform(basin_shp, 3035)

# 4. FILTER BASINS
print("Filtering basins: Excluding Russia and Turkey...")

# Create the spatial mask from our allowed countries
europe_mask <- allowed_countries_geo %>% st_transform(3035)

# Apply the filter (Keep only basins inside the mask)
basin_shp <- st_filter(basin_shp, europe_mask)

print(paste("Basins remaining after filter:", nrow(basin_shp)))


# 5. LOAD ADMIN BOUNDARIES (For Data Joining)
print("Downloading Admin Boundaries for Data Join...")

# Country Level (Use the codes we defined)
countries_shp <- gisco_get_nuts(year = "2021", nuts_level = 0, country = allowed_codes) %>%
  st_transform(3035)

# NUTS 2 (Use the codes we defined)
nuts2_shp <- gisco_get_nuts(year = "2021", nuts_level = 2, country = allowed_codes) %>%
  st_transform(3035)

# NUTS 3 (Use the codes we defined)
nuts3_shp <- gisco_get_nuts(year = "2021", nuts_level = 3, country = allowed_codes) %>%
  st_transform(3035)

print("Geometry Loaded and Filtered.")
# ==============================================================================
# 5. PROCESSING FUNCTION (FINAL & ROBUST)
# ==============================================================================

process_indicator <- function(basin_sf, config_item, admin_sf) {
  
  message(paste(">>> Processing:", config_item$name, "..."))
  
  # A. LOAD DATA
  if (!file.exists(config_item$file)) {
    stop(paste("ERROR: File not found:", config_item$file))
  }
  df <- read_csv(config_item$file, show_col_types = FALSE)
  
  # B. CHECK COLUMNS
  # 1. Check Join Column (ID)
  if (!config_item$join_col %in% names(df)) {
    stop(paste0("ERROR in ", config_item$name, ": Join column '", config_item$join_col, 
                "' not found. Found: ", paste(names(df), collapse=", ")))
  }
  
  # 2. Check Value Column (Data)
  if (!config_item$value_col %in% names(df)) {
    stop(paste0("ERROR in ", config_item$name, ": Value column '", config_item$value_col, 
                "' not found. Found: ", paste(names(df), collapse=", ")))
  }
  
  # C. STANDARDIZE COLUMNS
  # We rename the specific columns we identified to standard names
  # This prevents deleting columns if naming overlaps occur
  df$GEO_ID <- df[[config_item$join_col]]
  df$RAW_VALUE <- df[[config_item$value_col]]
  
  # Select ONLY what we need (ID and Value) to avoid conflicts
  df <- df %>% select(GEO_ID, RAW_VALUE)
  
  # FORCE NUMERIC
  if(!is.numeric(df$RAW_VALUE)) {
    message("   Converting values to numeric...")
    df$RAW_VALUE <- as.numeric(df$RAW_VALUE)
  }
  
  # D. PREPARE ADMIN GEOMETRY & CODES
  if(config_item$level == "Country") {
    # If using cleaned_gender_final.csv, IDs are already Codes (e.g. "AT", "EL")
    # We just need to ensure Greece is consistent if not already fixed
    df$GEO_ID[df$GEO_ID == "GR"] <- "EL"
    df$GEO_ID[df$GEO_ID == "GB"] <- "UK"
    
    admin_sf <- admin_sf %>% rename(GEO_ID = NUTS_ID)
  } else {
    admin_sf <- admin_sf %>% rename(GEO_ID = NUTS_ID)
  }
  
  # Remove NAs
  df <- df %>% filter(!is.na(GEO_ID) & !is.na(RAW_VALUE))
  
  # E. JOIN DATA TO MAP
  map_layer <- admin_sf %>% inner_join(df, by = "GEO_ID")
  
  # Check Match
  if(nrow(map_layer) == 0) {
    stop(paste0("ERROR: No matching regions found for ", config_item$name, 
                ". Check if CSV IDs match NUTS/Country codes."))
  } else {
    message(paste("   Matched", nrow(map_layer), "regions."))
  }
  
  # F. SPATIAL INTERSECTION
  intersection <- st_intersection(basin_sf, map_layer)
  intersection$area_w <- as.numeric(st_area(intersection))
  
  # G. AGGREGATE
  basin_stats <- intersection %>%
    st_drop_geometry() %>%
    group_by(HYBAS_ID) %>%
    summarise(
      weighted_raw_val = sum(RAW_VALUE * area_w, na.rm=T) / sum(area_w, na.rm=T)
    )
  
  # H. NORMALIZE
  basin_stats$norm_score <- normalize_data(
    basin_stats$weighted_raw_val, 
    config_item$s_curve, 
    config_item$steepness
  )
  
  result <- basin_stats %>% select(HYBAS_ID, norm_score, weighted_raw_val)
  names(result) <- c("HYBAS_ID", paste0("VULN_", config_item$name), paste0("RAW_", config_item$name))
  
  return(result)
}

# ==============================================================================
# 6. EXECUTE MODEL (ROBUST VERSION)
# ==============================================================================

# Initialize master data with Basin IDs
master_data <- basin_shp %>% st_drop_geometry() %>% select(HYBAS_ID)

# Function to safely run the process and join
safe_process_and_join <- function(current_data, config_item, admin_shape) {
  message(paste(">>> Processing:", config_item$name, "..."))
  
  # 1. Check if file exists
  if (!file.exists(config_item$file)) {
    stop(paste("CRITICAL ERROR: File not found:", config_item$file, 
               "\nPlease check your working directory or file name."))
  }
  
  # 2. Process
  result <- process_indicator(basin_shp, config_item, admin_shape)
  
  # 3. Check if result is valid
  if (nrow(result) == 0) {
    warning(paste("WARNING: Result for", config_item$name, "is empty! Check your joins."))
  }
  
  # 4. Join to Master Data
  new_data <- left_join(current_data, result, by = "HYBAS_ID")
  
  message(paste(">>> Success:", config_item$name, "added."))
  return(new_data)
}

# --- RUN INDICATORS SEQUENTIALLY ---
# If any step fails, the script will stop HERE so you can fix it.

# 1. Gender Equality
master_data <- safe_process_and_join(master_data, CONFIG$gender, countries_shp)

# 2. Rural Population
master_data <- safe_process_and_join(master_data, CONFIG$rural, nuts2_shp)

# 3. Poverty Risk (AROPE)
master_data <- safe_process_and_join(master_data, CONFIG$arope, nuts2_shp)

# 4. Social Dependency
master_data <- safe_process_and_join(master_data, CONFIG$dependency, nuts3_shp)

# 5. HDI
master_data <- safe_process_and_join(master_data, CONFIG$hdi, countries_shp)

# ==============================================================================
# 7. CALCULATE FINAL SOCIAL VULNERABILITY & FILTER
# ==============================================================================

# Verify all columns exist before calculation
required_cols <- c("VULN_Gender Equality", "VULN_Rural Population", 
                   "VULN_Poverty Risk (AROPE)", "VULN_Social Dependency", 
                   "VULN_Human Development Index")

missing_cols <- required_cols[!required_cols %in% names(master_data)]

if (length(missing_cols) > 0) {
  stop(paste("CANNOT CALCULATE INDEX. Missing columns:", paste(missing_cols, collapse=", ")))
}

# 1. CALCULATE INDEX
master_data <- master_data %>%
  mutate(
    Social_Vulnerability_Index = 
      (`VULN_Gender Equality` * CONFIG$gender$weight) +
      (`VULN_Rural Population` * CONFIG$rural$weight) +
      (`VULN_Poverty Risk (AROPE)` * CONFIG$arope$weight) +
      (`VULN_Social Dependency` * CONFIG$dependency$weight) +
      (`VULN_Human Development Index` * CONFIG$hdi$weight)
  )

# 2. FILTER: REMOVE BASINS WITH 0 VULNERABILITY
# This removes basins that had missing data or resulted in 0
print(paste("Basins before filtering zeros:", nrow(master_data)))

master_data <- master_data %>%
  filter(Social_Vulnerability_Index > 0)

print(paste("Basins after filtering zeros:", nrow(master_data)))

# 3. ATTACH TO SHAPEFILE
# We use inner_join to ensure we only keep the shapes that remain in our filtered data
final_map_sf <- basin_shp %>%
  inner_join(master_data, by = "HYBAS_ID") %>%
  st_transform(4326) # Ready for Leaflet

print("Saving data for the App...")

# We save it into the 'data' folder (one level up)
saveRDS(final_map_sf, "../data/ready_to_map_data.rds")

print("DONE! You can now run app.R")
