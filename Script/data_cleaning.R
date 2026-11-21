suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(haven)
  library(janitor)
  library(data.table)
  library(tidycensus)
})

##### Distance Matrix ####
inflow_name <- "/Users/leon/Documents/SavingSIAM/Data/stateinflow1516.csv"
outflow_name <- "/Users/leon/Documents/SavingSIAM/Data/stateoutflow1516.csv"
inflow_raw <- read_csv(inflow_name, col_types = cols(.default = "c"))

clean_migration_flows <- inflow_raw %>%
  mutate(
    origin_fips = as.numeric(y1_statefips),
    dest_fips   = as.numeric(y2_statefips),
    returns_n1  = as.numeric(n1)
  ) %>%
  filter(!dest_fips %in% c(96, 97, 98) & dest_fips <= 56) %>%
  filter(!origin_fips %in% c(96, 97, 98) & origin_fips <= 56) %>%
  filter(origin_fips != dest_fips) %>%
  select(origin_fips, dest_fips, returns_n1) %>%
  arrange(origin_fips, dest_fips)

flow_matrix_wide <- clean_migration_flows %>%  
  pivot_wider(
    names_from = dest_fips,
    values_from = returns_n1,
    values_fill = list(returns_n1 = 0)
  ) %>%
  arrange(origin_fips) %>%
  select(origin_fips, num_range("", sort(unique(clean_migration_flows$dest_fips))))

write_matrix_ready <- flow_matrix_wide %>% select(-origin_fips)
write_csv(write_matrix_ready, "/Users/leon/Documents/SavingSIAM/Output/empirical_F_matrix.csv", col_names = FALSE)


#### USDA Natural Amenities scale ####

amenities_name <- "/Users/leon/Documents/SavingSIAM/Data/amenities/NATAMENF.csv"
census90_name <- "/Users/leon/Documents/SavingSIAM/Data/amenities/cencounts.csv"

census90_data <- read.csv(census90_name)
amenities_data <- fread(amenities_name,skip = 3, col.names = c(
  "fips_code", "fips_used", "state", "county_name", "division",
  "rural_urban_code", "influence_code", "jan_temp_mean", "jan_sun_mean",
  "july_temp_mean", "july_hum_mean", "topography_code", "water_area_pct",
  "water_area_log", "jan_temp_z", "jan_sun_z", "july_temp_z", "july_hum_z",
  "topog_z", "water_area_z", "nat_amenity_scale", "nat_amenity_rank"
))

census90_data$pop1990 <- as.integer(census90_data$pop1990)

census90 <- census90_data %>% 
  select(pop1990, name) %>%
  mutate(
    state = substr(name, 1, 2),
    county_name = str_sub(name, 4, -1) %>%
      str_remove_all("(?i) (COUNTY|PARISH|OF|THE)") %>%
      str_trim() %>%
      str_to_upper()
  ) %>%  
  group_by(state) %>% 
  arrange(desc(pop1990)) %>%
  slice(-1) %>%
  ungroup()

amenities <- right_join(amenities_data, census90 , by = c("state", "county_name"))

df_clean <- amenities %>%
  clean_names() %>%
  filter(!is.na(state)) %>%
  mutate(
    pop_1990 = as.numeric(pop1990),
    amenity_score = as.numeric(nat_amenity_scale)
  )  %>% 
  group_by(state) %>%
  summarize(
    simple_mean_amenity = mean(amenity_score, na.rm = TRUE),
    pop_weighted_amenity = sum(amenity_score * pop_1990, na.rm = TRUE) / sum(pop_1990, na.rm = TRUE)
  ) %>%
  ungroup()

state_crosswalk <- tibble(
  state = state.abb,
  state_name = state.name
)

final_export <- df_clean %>%
  left_join(state_crosswalk, by = "state") %>%
  mutate(state_name = case_when(
    state == "DC" ~ "District of Columbia",
    TRUE ~ state_name
  )) %>%
  filter(!is.na(state_name)) %>%
  arrange(state_name)

write_csv(final_export, "/Users/leon/Documents/SavingSIAM/Output/state_natural_amenities.csv")


#### NAEP (educ) ####

clean_naep <- function(file_path, subject_name) {
  read_csv(file_path, show_col_types = FALSE) %>%
    clean_names() %>%
    filter(national_school_lunch_program_eligibility_3_categories == "Not eligible") %>%
    mutate(average_scale_score = suppressWarnings(as.numeric(average_scale_score))) %>%
    filter(!is.na(average_scale_score)) %>%
    group_by(jurisdiction) %>%
    summarize(!!sym(paste0(subject_name, "_score")) := mean(average_scale_score, na.rm = TRUE)) %>%
    ungroup()
}

math_data <- clean_naep("/Users/leon/Documents/SavingSIAM/Data/amenities/NDECMath.csv", "math")
reading_data <- clean_naep("/Users/leon/Documents/SavingSIAM/Data/amenities/NDECReading.csv", "reading")

education_amenities <- math_data %>%
  inner_join(reading_data, by = "jurisdiction") %>%
  filter(jurisdiction %in% c(state.name, "District of Columbia")) %>%
  rename(state_name = jurisdiction) %>%
  mutate(
    math_z = (math_score - mean(math_score)) / sd(math_score),
    reading_z = (reading_score - mean(reading_score)) / sd(reading_score),
    edu_amenity_index = (math_z + reading_z) / 2
  ) %>%
  arrange(state_name)

write_csv(education_amenities, "/Users/leon/Documents/SavingSIAM/Output/education_amenities.csv")

#### Economic Freedom ####
efna_raw <- read_csv("/Users/leon/Documents/SavingSIAM/Data/amenities/EconFreedom.csv", name_repair = "unique_quiet", show_col_types = FALSE)

efna_clean <- efna_raw %>%
  select(
    state_name = `Overall subn Scores`,
    matches("^(19|20)\\d{2}(\\.\\.\\.([2-9]|[1-3]\\d|4[0-4]))?$")
  ) %>%
  rename_with(~ str_remove(., "\\.\\.\\..*"), -state_name) %>%
  filter(state_name %in% state.name) %>%
  mutate(across(-state_name, ~ suppressWarnings(as.numeric(.)))) %>%
  pivot_longer(
    cols = -state_name, 
    names_to = "year", 
    values_to = "score"
  ) %>%
  mutate(year = as.numeric(year)) %>%
  filter(year >= 2005 & year <= 2010) %>%
  group_by(state_name) %>%
  summarize(efna_score_raw = mean(score, na.rm = TRUE)) %>%
  ungroup()

dc_row <- tibble(state_name = "District of Columbia", efna_score_raw = NA_real_)

efna_final <- bind_rows(efna_clean, dc_row) %>%
  arrange(state_name) %>%
  mutate(efna_z = (efna_score_raw - mean(efna_score_raw, na.rm = TRUE)) / sd(efna_score_raw, na.rm = TRUE))

write_csv(efna_final, "/Users/leon/Documents/SavingSIAM/Output/freedom_amenities.csv")

#### FBI Crime Report ####
crime_raw <- read_csv("/Users/leon/Documents/SavingSIAM/Data/amenities/estimated_crimes_1979_2024.csv", show_col_types = FALSE)

crime_clean <- crime_raw %>%
  filter(year >= 2010 & year <= 2015) %>%
  filter(!is.na(state_abbr) & state_abbr != "") %>%
  group_by(state_abbr) %>%
  summarize(
    avg_violent_rate = mean((violent_crime / population) * 100000, na.rm = TRUE),
    avg_property_rate = mean((property_crime / population) * 100000, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  left_join(tibble(state_abbr = state.abb, state_name = state.name), by = "state_abbr") %>%
  mutate(state_name = case_when(
    state_abbr == "DC" ~ "District of Columbia",
    TRUE ~ state_name
  )) %>%
  filter(!is.na(state_name)) %>%
  select(state_name, avg_violent_rate, avg_property_rate) %>%
  arrange(state_name) %>%
  mutate(crime_z = (avg_violent_rate - mean(avg_violent_rate)) / sd(avg_violent_rate))

write_csv(crime_clean, "/Users/leon/Documents/SavingSIAM/Output/amenities_crime.csv")

#### CBP ####

cbp_raw <- read_csv("/Users/leon/Documents/SavingSIAM/Data/amenities/CBP2015.csv", show_col_types = FALSE)

cbp_alaska <- cbp_raw %>% 
  filter(NAME == "Alaska") %>%
  filter(LFO_LABEL == "All establishments" & EMPSZES_LABEL == "All establishments")
  

  
cbp_clean <- cbp_raw %>%
  filter(GEO_ID != "Geographic identifier code") %>%
  filter(LFO_LABEL == "All establishments" & EMPSZES_LABEL == "All establishments") %>%
  filter(NAME %in% c(state.name, "District of Columbia")) %>%
  select(state_name = NAME, naics = NAICS2012, estab = ESTAB) %>%
  mutate(estab = as.numeric(estab)) %>%
  filter(naics %in% c("00", "71", "72")) %>%
  pivot_wider(names_from = naics, values_from = estab, names_prefix = "naics_") %>%
  mutate(
    amenity_estab = naics_71 + naics_72,
    amenity_share = amenity_estab / naics_00,
    urban_z = (amenity_share - mean(amenity_share, na.rm = TRUE)) / sd(amenity_share, na.rm = TRUE)
  ) %>%
  arrange(state_name) %>%
  select(state_name, amenity_share, urban_z)

write_csv(cbp_clean, "/Users/leon/Documents/SavingSIAM/Output/amenities_variety.csv")

#### Density (Census) ####

options(tigris_use_cache = TRUE)

# census_api_key("4c0e5e7c8b558d086d1f735e135f70a7ea2d0420", install = TRUE, overwrite=TRUE)

get_state_density <- function(st_name) {
  tryCatch({
    tract_data <- get_decennial(
      geography = "tract",
      variables = "P001001",
      year = 2010,
      geometry = TRUE,
      state = st_name
    )
    Sys.sleep(2)
    tract_clean <- tract_data %>%
      mutate(
        land_area_m2 = as.numeric(sf::st_area(geometry)),
        land_area_km2 = land_area_m2 / 1e6
      ) %>%
      sf::st_drop_geometry() %>%
      rename(pop_tract = value) %>%
      filter(land_area_km2 > 0 & pop_tract > 0)
    return(tract_clean)
  }, error = function(e) {
    return(NULL)
  })
}

spatial_universe <- c(state.name, "District of Columbia")
all_tracts_clean <- map_dfr(spatial_universe, get_state_density)

final_density <- all_tracts_clean %>%
  separate(NAME, into = c("tract_name", "county_name", "state_name"), sep = ", ") %>%
  mutate(state_name = str_trim(state_name)) %>%
  group_by(state_name) %>%
  mutate(
    pop_state = sum(pop_tract, na.rm = TRUE),
    weight_tract = pop_tract / pop_state,
    density_tract = pop_tract / land_area_km2
  ) %>%
  summarize(
    experienced_density = sum(weight_tract * density_tract, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(urban_density_z = (experienced_density - mean(experienced_density)) / sd(experienced_density)) %>%
  arrange(state_name)

write_csv(final_density, "/Users/leon/Documents/SavingSIAM/Output/amenities_density.csv")

#### Damage (Hsiang, 2017)####

state_map <- tibble(GeoName = state.name, state_code = state.abb) %>%
  add_row(GeoName = "District of Columbia", state_code = "DC")

#agriculture share
agri_share_raw <- read_csv(
  "/Users/leon/Documents/SavingSIAM/Data/SAGDP2.csv",
  na           = c("(NA)", "(D)", "(L)", ""),  
  show_col_types = FALSE
)

agri_share <- agri_share_raw  %>%
  mutate(GeoFIPS = str_remove_all(GeoFIPS, '"')) %>%
  filter(
    !is.na(GeoName),
    GeoFIPS >= "01000", 
    GeoFIPS <= "56000"  
  )

year_cols <- as.character(1997:2025)

agri_long <- agri_share %>% filter(LineCode == 3) %>% select(GeoFIPS, GeoName, all_of(year_cols)) %>%
  pivot_longer(all_of(year_cols), names_to = "Year", values_to = "Agri_GDP_M")
total_long <- agri_share %>% filter(LineCode == 1) %>% select(GeoFIPS, GeoName, all_of(year_cols))%>%
  pivot_longer(all_of(year_cols), names_to = "Year", values_to = "Total_GDP_M")


share_long_2010s <- agri_long %>%
  left_join(total_long, by = c("GeoFIPS", "GeoName", "Year")) %>%
  left_join(state_map, by = c("GeoName")) %>%
  mutate(
    Year       = as.integer(Year),
    Agri_Share = Agri_GDP_M / Total_GDP_M   
  ) %>% 
  filter(Year %in% c("2010", "2011", "2012", "2013", "2014", "2015")) %>%
  group_by(GeoName) %>%
  reframe(agri_share = mean(Agri_Share), state = unique(state_code)) 


path_data <- "/Users/leon/Documents/SavingSIAM/Data/"
path_out  <- "/Users/leon/Documents/SavingSIAM/Output/data/"

VSL <- 9000000         
alpha_coastal <- 0.0
agri_share_def <- 0.015 

state_map <- tibble(GeoName = state.name, state_code = state.abb) %>%
  add_row(GeoName = "District of Columbia", state_code = "DC")

elec_raw <- read_csv(paste0(path_data, "elect_exp_states.csv"), skip = 3, show_col_types = FALSE) %>%
  select(GeoName, elec_2010 = `2010`) %>%
  mutate(elec_2010 = as.numeric(elec_2010))

inc_raw <- read_csv(paste0(path_data, "person_incom_states.csv"), skip = 3, show_col_types = FALSE) %>%
  select(GeoName, inc_2010 = `2010`) %>%
  mutate(inc_2010 = as.numeric(inc_2010))

energy_shares <- elec_raw %>%
  inner_join(inc_raw, by = "GeoName", relationship = "many-to-many") %>%
  mutate(energy_share = elec_2010 / inc_2010) %>%
  inner_join(state_map, by = "GeoName") %>%
  select(state_code, energy_share)

county_raw <- read_csv(paste0(path_data, "county_damages_by_sector.csv"), show_col_types = FALSE)

county_clean <- county_raw %>%
  rename(
    state          = `State Code`,
    pop            = `County Population (in 2012)`,
    income_pc      = `County Income (in 2012)`, 
    agri_raw       = `Agricultural Damage (%) (4 major crops)`,
    mortality_raw  = `Mortality (deaths per 100k)`,
    energy_raw     = `Energy Expenditures (%)`,
    labor_low      = `Labor Low-risk (%)`,
    labor_high     = `Labor High-risk (%)`,
    coastal_log    = `Coastal damage (log10(% county income))`,
    crime_prop     = `Property Crime (%)`,
    crime_viol     = `Violent Crime (%)`
  ) %>%
  mutate(
    total_income = pop * income_pc,
    coastal_pct = ifelse(is.na(coastal_log), 0, 10^coastal_log),
    mortality_usd_pc = (mortality_raw / 100000) * VSL,
    agri_raw = ifelse(is.na(agri_raw), 0, agri_raw)
  )

state_damages <- county_clean %>%
  left_join(state_map, by = c("state" = "state_code")) %>%
  group_by(GeoName) %>%
  summarise(
    mortality_usd_pc = weighted.mean(mortality_usd_pc, w = pop, na.rm = TRUE),
    crime_prop = weighted.mean(crime_prop, w = pop, na.rm = TRUE),
    crime_viol = weighted.mean(crime_viol, w = pop, na.rm = TRUE),
    income_pc  = weighted.mean(income_pc, w = pop, na.rm = TRUE),
    labor_tot  = weighted.mean((labor_low + labor_high)/2, w = total_income, na.rm = TRUE),
    agri_yield = weighted.mean(agri_raw, w = total_income, na.rm = TRUE),
    energy_exp = weighted.mean(energy_raw, w = total_income, na.rm = TRUE),
    coastal    = weighted.mean(coastal_pct, w = total_income, na.rm = TRUE),
    state = unique(state),
    .groups = "drop"
  )
ucr_data <- read_csv(paste0(path_data, "amenities/estimated_crimes_1979_2024.csv"), show_col_types = FALSE)

ucr_baseline_us <- ucr_data %>%
  filter(year >= 2010, year <= 2015, state_name == "United States Total") %>%
  summarise(
    violent_crime = sum(violent_crime, na.rm = TRUE),
    homicide = sum(homicide, na.rm = TRUE),
    rape = sum(coalesce(rape_revised, rape_legacy), na.rm = TRUE),
    robbery = sum(robbery, na.rm = TRUE),
    aggravated_assault = sum(aggravated_assault, na.rm = TRUE),
    property_crime = sum(property_crime, na.rm = TRUE),
    burglary = sum(burglary, na.rm = TRUE),
    larceny = sum(larceny, na.rm = TRUE),
    motor_vehicle_theft = sum(motor_vehicle_theft, na.rm = TRUE)
  )

# McCollister et al. (2010) Table 5 Unit Costs
cost_murder   <- 8982907
cost_rape     <- 240776
cost_robbery  <- 42310
cost_assault  <- 107020

cost_burglary <- 6462
cost_larceny  <- 3532
cost_motor    <- 10772

pct_murder  <- ucr_baseline_us$homicide / ucr_baseline_us$violent_crime
pct_rape    <- ucr_baseline_us$rape / ucr_baseline_us$violent_crime
pct_robbery <- ucr_baseline_us$robbery / ucr_baseline_us$violent_crime
pct_assault <- ucr_baseline_us$aggravated_assault / ucr_baseline_us$violent_crime

vsc_violent <- (pct_murder * cost_murder) + (pct_rape * cost_rape) + 
               (pct_robbery * cost_robbery) + (pct_assault * cost_assault)

pct_burglary <- ucr_baseline_us$burglary / ucr_baseline_us$property_crime
pct_larceny  <- ucr_baseline_us$larceny / ucr_baseline_us$property_crime
pct_motor    <- ucr_baseline_us$motor_vehicle_theft / ucr_baseline_us$property_crime

vsc_property <- (pct_burglary * cost_burglary) + (pct_larceny * cost_larceny) + 
                (pct_motor * cost_motor)
# -------------------------------------------------------------------------
###### Computation of the average earning #####

cps_raw <- read_csv(paste0(path_data, "cps.csv.gz"), show_col_types = FALSE)

cps_filtered <- cps_raw %>%
  filter(
    YEAR >= 2010,
    YEAR <= 2015,
    AGE >= 25,
    AGE <= 60,
    EMPSTAT %in% c(10, 12),
    QEARNWEE == 0 | is.na(QEARNWEE),
    !is.na(EARNWEEK2),
    EARNWEEK2 < 9999,        
    EARNWEEK2 < 999999       
  ) %>% 
  mutate(
    EARNWEEK_ADJ = ifelse(max(EARNWEEK2, na.rm = TRUE) > 20000, EARNWEEK2 / 100, EARNWEEK2),
    implicit_hourly_wage = EARNWEEK_ADJ / UHRSWORKORG
  ) %>%
  mutate(
    EARNWEEK_ADJ = ifelse(max(EARNWEEK2, na.rm = TRUE) > 20000, EARNWEEK2 / 100, EARNWEEK2),
    implicit_hourly_wage = EARNWEEK_ADJ / UHRSWORKORG
  ) %>%
  filter(
    EARNWEEK2 >= (5 * 3.625)
  )

bounds <- quantile(cps_filtered$EARNWEEK2, probs = c(0.01, 0.99), na.rm = TRUE)

cps_trimmed <- cps_filtered %>%
  filter(EARNWEEK2 >= bounds[1] & EARNWEEK2 <= bounds[2])

mean_weekly_earnings <- weighted.mean(cps_trimmed$EARNWEEK2, w = cps_trimmed$EARNWT, na.rm = TRUE)

Y_bar_US <- mean_weekly_earnings * 52

print(paste("Votre scalaire national Y_bar_US est :", round(Y_bar_US, 2), "USD"))

gamma <- 3.0 

final_shocks <- state_damages %>%
  left_join(energy_shares, by = c("state" = "state_code")) %>%
  left_join(share_long_2010s, by = "state") %>%
  left_join(crime_clean, by = c("GeoName.x" = "state_name")) %>%
  mutate(
    shock_wage = (labor_tot / 100) + ((agri_yield / 100) * agri_share) - ((coastal / 100) * (1 - alpha_coastal)),
    
    extra_violent_pc = (crime_viol / 100) * (avg_violent_rate / 100000),
    extra_property_pc = (crime_prop / 100) * (avg_property_rate / 100000),
    cost_crime_pc = (extra_violent_pc * 205782) + (extra_property_pc * 4799),
    
    total_amenity_damage_usd = cost_crime_pc + mortality_usd_pc,
    
    damage_model_units = total_amenity_damage_usd / Y_bar_US,
    
    shock_amenity = gamma * damage_model_units,
    
    shock_housing = ((energy_exp / 100) * energy_share) - ((coastal / 100) * alpha_coastal)
  )

df_wage <- final_shocks %>% select(GeoName.x, shock_wage)
write_csv(df_wage, paste0(path_out, "shock_wage_2080.csv"))

df_amenity <- final_shocks %>% select(GeoName.x, shock_amenity)
write_csv(df_amenity, paste0(path_out, "shock_amenity_2080.csv"))

df_housing <- final_shocks %>% select(GeoName.x, shock_housing)
write_csv(df_housing, paste0(path_out, "shock_housing_2080.csv"))

print("Data cleaning and processing completed successfully.") 

