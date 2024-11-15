

library(dplyr); library(ggplot2); library(readr); library(tidyr); library(lubridate);
library(strucchange); library(bfast); library(segmented)
library(sp); library(sf); library(countrycode); library(r2country);
library(MASS); library(ISOweek); library(stlplus); library(sfsmisc)
library(spdep)
library(metaforest)
library(caret)
library(parallel)
library(doParallel)
# Define function to complete weekly data for each year, which will be used later
# Define the function to calculate tile area, used for computing density
# The functions 'DTECTslope_ro' and 'DTECTslope_0_ro' for detecting robust linear regressions should be run first.

complete_adj_week <- function(data) {
  years <- unique(data$Year)
  list_data <- lapply(years, function(y) {
    year_data <- data[data$Year == y, ]
    start_week <- min(year_data$adj_Week)
    if (y == 2020 || y == 2021) {
      end_week <- 52  # Fix the end week as 52 for the years 2020 and 2021
    } else {
      end_week <- max(year_data$adj_Week) 
    } 
    full_weeks <- data.frame(adj_Week = start_week:end_week)
    full_data <- merge(full_weeks, year_data, by = "adj_Week", all.x = TRUE)
    full_data$Year <- y 
    return(full_data)
  })
  
  full_data <- do.call(rbind, list_data)
  return(full_data)
}

##################################################

calculate_tile_area <- function(lat, level) {
  earth_circumference_km <- 40075  # Earth's circumference at the equator in kilometers
  degrees_per_tile <- 360 / (2^level)  # Degrees of longitude per tile
  rad_per_deg <- pi / 180  # Radians per degree
  
  # Calculate the width of a tile in kilometers at the given latitude
  tile_width_km <- cos(lat * rad_per_deg) * earth_circumference_km / (2^level)
  tile_height_km <- degrees_per_tile / 360 * earth_circumference_km
  
  # Area of the tile in square kilometers
  tile_area_km <- tile_width_km * tile_height_km
  return(tile_area_km)
}

######################################################
####Detect trend for each tile
####The breakpoint for each tile was derived from the breakpoint detected for this country.

directory_path <- "*"
country_dirs <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

bp_data <- read_csv('*/country_bp_urban_rlm.csv')
record <- read_csv('*/country_name.csv')

# Define function to process data for each country directory
country_data <- function(country_dir) {
  country_name <- basename(country_dir)
  level <- record %>%
    filter(country == country_name) %>%
    pull(level)
  night_data <- file.path(country_dir, "lulc_night_update.csv")
  
  if(!file.exists(night_data )) {
    warning(paste("night.csv not found in", country_dir))
    return(NULL)
    } else {
      pop <- read_csv(night_data)
      pop$est[pop$est < 0] <- 0
    
    # Define workdays based on country
      if (country_name %in% c('Algeria', 'Egypt', 'Jordan', 'Qatar', 'Saudi Arabia')) {
      workday <- c('Sun', 'Mon', 'Tue', 'Wed', 'Thu')
      } else if (country_name == 'Nepal') {
      workday <- c('Sun', 'Mon', 'Tue', 'Wed', 'Thu', 'Fri')
      } else if (country_name == 'United Arab Emirates') {
      workday <- c('Mon', 'Tue', 'Wed', 'Thu')
      } else {
      workday <- c('Mon', 'Tue', 'Wed', 'Thu', 'Fri')  # Default case for other countries
      }
      
      pop_week <- pop %>% filter(week %in% workday)
      pop_week <-  pop_week %>%
        mutate(
          Year = isoyear(Date),
          Week = isoweek(Date),
          adj_Week = if_else(Year == 2020, Week - 1, Week),
          weekdate = Date - lubridate::wday(Date, week_start = 1) + 1   # Assuming the week starts on Monday
      )
      pop_week$weekdate <- as.Date(pop_week$weekdate)
    
      pop_week_grouped <- pop_week %>%
      dplyr::group_by(lat, lon, urban_category, weekdate, urban, tree, crop, water, adj_Week, Week, Year) %>%
      dplyr::summarise(ave = mean(est, na.rm = TRUE)) %>%
      ungroup()
    
      tile_area <- calculate_tile_area(pop_week_grouped$lat, level)
      pop_week_grouped$area <- tile_area
      pop_week_grouped$den <-  pop_week_grouped$ave/pop_week_grouped$area
    
      pop_week_list <- pop_week_grouped %>%
        group_by(lat, lon, area, urban_category, urban, tree, crop, water) %>%
        group_split()
    
      pop_week_ts_list <- lapply(pop_week_list, function(df) {
        if (any(df$ave == 0)| length(df$den) <= 104) {
          return(NULL)
          }
        df$location_id = paste(df$lat[1], df$lon[1], sep = "_")
        df$lat = df$lat[1]
        df$lon = df$lon[1]
        df$area = df$area[1]
        df$urban_category = df$urban_category[1]
        df$urban = df$urban[1]
        df$crop = df$crop[1]
        df$tree = df$tree[1]
        df$water = df$water[1]
        df$weekdate = df$weekdate
        df$Week = df$Week
        df$adj_Week = df$adj_Week
        return(df)
        })
      
      pop_week_ts_list <- Filter(Negate(is.null), pop_week_ts_list)
    
      pop_week_ts_list <- lapply(pop_week_ts_list, function(list_item) {
        if (is.null(list_item)) return(NULL) 
      
        completed_df = complete_adj_week(list_item) 
        completed_df$Week = if_else(completed_df$Year == 2020, completed_df$adj_Week + 1, completed_df$adj_Week)
        final_df <- completed_df  %>%
          arrange(Year, adj_Week, Week) %>%
          group_by(Year) %>%
          mutate(
            weekdate = ifelse(is.na(weekdate), ISOweek2date(sprintf("%d-W%02d-%d", Year, Week, 1)), weekdate),
            ave = ifelse(is.na(ave), NA, ave),
            den = ifelse(is.na(den), NA, den)
            ) %>%
          ungroup()
        final_df$weekdate = as.Date(final_df$weekdate)
        final_df$number <- as.numeric(as.factor((final_df$weekdate)))

        final_df$ts_data <- ts(final_df$den, frequency = 52, start = c(2020, final_df$adj_Week[1]))
        return(final_df)
        })
    
      pop_week_ts_list <- Filter(Negate(is.null), pop_week_ts_list)

      stl_results <- list() 
    
    #Use breakpoint detected for this country before
      bps <- bp_data %>%
        filter(country == country_name) %>%
        pull(bp_num)
      check_missing_subseries <- function(ts_data) {
        subseries <- split(ts_data, cycle(ts_data))
        all_missing <- sapply(subseries, function(sub) all(is.na(sub)))
        return(any(all_missing))
        }
      for (i in seq_along(pop_week_ts_list)) {
        print(paste("Processing index:", i))  # Print the current index
        df <- pop_week_ts_list[[i]]
        if (!is.null(df) && length(df$ts_data) > 2 * frequency(df$ts_data)) {
          if (check_missing_subseries(df$ts_data)) {
            print(paste("Skipping index:", i, "due to all values missing in at least one subseries"))
            stl_results[[i]] <- NULL
            next
            }
          df <- as.list(df)
          stlv <- stlplus(df$ts_data, s.window = "periodic")
          trend_vector <- as.vector(stlv[["data"]][["trend"]])
          season_vector <- as.vector(stlv[["data"]][["seasonal"]])
          time_vector <- as.vector(time(df$ts_data))
          df$trend <- trend_vector
          df$season <- season_vector
          df$deseason <- df$ts_data - df$season
          df$ts_time <- time_vector
        
          df_filtered <-  as.data.frame(df) %>% filter(!is.na(ts_data))
          if (!is.na(bps)) {
            bptest <- DTECTslope_ro(df_filtered, bps)
            } else {
              bptest <- DTECTslope_0_ro(df_filtered)
              }
        
          df$bp <- bptest
          stl_results[[i]] <- df
          
          } else {
            stl_results[[i]] <- NULL  # Or handle short series differently
          }
        }
    
      nonull_count <- sum(!sapply(stl_results, is.null))
      print(nonull_count)

      s1 <- numeric()
      s2 <- numeric()
      sall <- numeric()
      p1 <- numeric()
      p2 <- numeric()
      pall <- numeric()
      lat <- numeric()
      lon <- numeric()
      urban <- numeric()
      tree <- numeric()
      crop <- numeric()
      water <- numeric()
      
      for (result in stl_results) {
        if (!is.null(result)) {
          lat <- c(lat, result$lat[1])
          lon <- c(lon, result$lon[1])
          urban <- c(urban, result$urban[1])
          water <- c(water, result$water[1])
          tree <- c(tree, result$tree[1])
          crop <- c(crop, result$crop[1])
          bp <- result$bp
        
          s1 <- c(s1, ifelse(!is.null(bp$s1), bp$s1, NA))
          s2 <- c(s2, ifelse(!is.null(bp$s2), bp$s2, NA))
          sall <-  c(sall, ifelse(!is.null(bp$sall), bp$sall, NA))
          p1 <- c(p1, ifelse(!is.null(bp$p1), bp$p1, NA))
          p2 <- c(p2, ifelse(!is.null(bp$p2), bp$p2, NA))
          pall <- c(pall, ifelse(!is.null(bp$pall), bp$pall, NA))
        }
        }
      result_df <- data.frame(
      lat = lat,
      lon = lon,
      urban= urban,
      tree =tree,
      crop = crop,
      water = water,
      s1 = s1,
      s2 = s2,
      sall = sall,
      p1 = p1,
      p2 = p2,
      pall = pall
      )
    
      # Calculate Local Moran's I 
    result_sf <- st_as_sf(result_df, coords = c("lon", "lat"), crs = 4326)
    neighbors <- knn2nb(knearneigh(st_coordinates(result_sf), k = 8))
    weights <- nb2listw(neighbors, style = "W", zero.policy = TRUE)
    
    if (!is.na(result_sf$s1)) {
      local_moran <- localmoran(result_sf$s1, weights)
      result_df$local_moran_i <- local_moran[,1]  # Local Moran's I values
      result_df$local_moran_p_value <- local_moran[,5]
      classification_df <- attr(local_moran, "quadr")
      median_classification <- classification_df$median
      result_df$classification <- median_classification
      result_df <- result_df %>%
      mutate(classification_s1 = case_when(
        local_moran_p_value >= 0.05 ~ "No Significant",
        local_moran_p_value < 0.05 & (classification %in% c("High-Low", "Low-High")) ~ "Outliers",
        local_moran_p_value < 0.05 ~ classification
      ))
    
    local_moran2 <- localmoran(result_sf$s2, weights)
    result_df$local_moran_2 <- local_moran2[,1]  # Local Moran's I values
    result_df$local_moran_p_value2 <- local_moran2[,5]
    classification_df2 <- attr(local_moran2, "quadr")
    median_classification2 <- classification_df2$median
    result_df$classification2 <- median_classification2
    result_df <- result_df %>%
      mutate(classification_s2 = case_when(
        local_moran_p_value2 >= 0.05 ~ "No Significant",
        local_moran_p_value2 < 0.05 & (classification2 %in% c("High-Low", "Low-High")) ~ "Outliers",
        local_moran_p_value2 < 0.05 ~ classification2
      ))
    
    } else {
    local_moran <- localmoran(result_sf$sall, weights)
    result_df$local_moran_i <- local_moran[,1]  # Local Moran's I values
    result_df$local_moran_p_value <- local_moran[,5]
    classification_df <- attr(local_moran, "quadr")
    median_classification <- classification_df$median
    result_df$classification <- median_classification
    result_df <- result_df %>%
      mutate(classification_s1 = case_when(
        local_moran_p_value >= 0.05 ~ "No Significant",
        local_moran_p_value < 0.05 & (classification %in% c("High-Low", "Low-High")) ~ "Outliers",
        local_moran_p_value < 0.05 ~ classification
        ),
      classification_s2 = classification_s1
      )
    }
    write.csv(result_df, paste0(country_dir, "/", country_name, "_slopes_den.csv"))
    }
}


###################################################
####This function provides correlations between trend and built-up density across urban areas for each country in the list.

directory_path <- "*"

country_dirs <- list.dirs(directory_path, full.names = TRUE, recursive = FALSE)

country_data <- function(country_dir_list) {
  result_list <- list()
  for (country_dir in country_dir_list) {
    country_name <- basename(country_dir)
    grid_file <- file.path(country_dir, paste0(country_name, "_slopes_den.csv"))
    
    if (!file.exists(grid_file)) {
      warning(paste(country_name, "_slopes.csv not found in", country_dir))
      return(NULL)

      } else {
        s <- read_csv(grid_file)
        s <- s %>% filter(urban > 0.25)
      
        if (country_name %in% c('Poland', 'Turkey','Egypt', 'Italy','Nepal')) {
          
        s1 <- s %>% filter(pall < 0.05)
        cor_result1 <- cor.test(s1$sall, s1$urban, conf.level = 0.95, method = "spearman")
        corr1 <- cor_result1$estimate
        p1 <- cor_result1$p.value
        
        corr2 <-corr1
        p2 <- cor_result1$p.value
        n <- nrow(s)
        } else {

          s1 <- s %>% filter(p1 < 0.05)
          s2 <- s %>% filter(p2 < 0.05)
          cor_result1 <- cor.test(s1$s1, s1$urban, conf.level = 0.95, method = "spearman")
          cor_result2 <- cor.test(s2$s2, s2$urban, conf.level = 0.95, method = 'spearman')
          corr1 <- cor_result1$estimate
          p1 <- cor_result1$p.value
        
          corr2 <- cor_result2$estimate
          p2 <- cor_result2$p.value
          n <- nrow(s)
          }
        df <- data.frame(country = country_name, n=n, corr1 = corr1, corr2 = corr2, 
                       p1 = p1, p2=p2)
      result_list[[country_name]] <- df
      return(df)
      }
    }
  result_df <- do.call(rbind, result_list)
  return(result_df)
}

all_cor <- lapply(country_dirs, country_data)
all_cor2 <- do.call(rbind, all_cor)
write.csv(all_cor2, '*/urban_grid_corre_spearman.csv')


##################################################
#### Using meta-forest and meta-regression models to explore the factors impacting heterogeneity

name <- read_csv('*/country_name.csv')
all <- read_csv('*/urban_grid_corre_spearman.csv')
#Fisher-Z transformation of correlations to obtain effect sizes
trans_1 <- escalc(measure = "ZCOR",  
                  ri = corr1,        
                  ni = n,           
                  data = all
)
trans_2 <- escalc(measure = "ZCOR", 
                  ri = corr2,         
                  ni = n,           
                  data = all
)
link1 <- trans_1   %>%
  left_join(name , by ='country')
link2 <- trans_2   %>%
  left_join(name , by ='country')

# Select relevant columns for meta-forest analysis
data1 <- link1[, c("yi", "vi","hdi_fina", "hdi_change_fina", "gdp_grow", "gdp_indu_fina", "tourism", "protect_fina", "pmdif2", "greendifcore")]
data2 <- link2[, c("yi", "vi", "hdi_fina", "hdi_change_fina", "gdp_grow", "gdp_indu_fina", "tourism", "protect_fina", "pmdif2", "greendifcore")]

moderators <- c( "hdi_fina", "hdi_change_fina", "gdp_grow", "gdp_indu_fina", "tourism", "protect_fina", "pmdif2", "greendifcore")

# MetaForest modeling to check convergence
set.seed(8)
check_conv <- MetaForest(z1~.,
                         data = data1,
                         whichweights = "random",
                         num.trees = 20000)
plot(check_conv)

# MetaForest model with specified number of trees
set.seed(2603)
mf_rep  <- MetaForest(yi~.,
                      data = data1,
                      whichweights = "random",
                      num.trees = 10000)

# Perform recursive preselection
pre <- preselect(mf_rep, replications = 200, algorithm = "recursive")
plot(pre)
# Cross-validation setup
cv_folds <- trainControl(method = "cv", 10, allowParallel = TRUE)
# Set up a tuning grid for the three tuning parameters of MetaForest
tuning_grid <- expand.grid(whichweights = c("random", "fixed", "unif"),
                           mtry = 1:6,
                           min.node.size = 2:6)

# Select variables based on preselection
X <- dplyr::select(data1, "vi", preselect_vars(pre, cutoff = .95))

# Register parallel processing (using 31 cores)
registerDoParallel(31) 
# Train MetaForest model with cross-validation
set.seed(2582)
mf_cv <- train(y = data1$yi,
               x = X,
               method = ModelInfo_mf(),
               trControl = cv_folds,
               tuneGrid = tuning_grid,
               keep.inbag = TRUE,
               verbose=TRUE,
               num.trees = 10000)

# Extract best result based on RMSE and corresponding R-squared
mf_cv$results[which.min(mf_cv$results$RMSE), ]
r2_cv <- mf_cv$results$Rsquared[which.min(mf_cv$results$RMSE)]

# Extract final model
forest <- mf_cv$finalModel
r2_oob <- forest$forest$r.squared

# Plot convergence and variable importance
plot(forest)
importance.export <- varImp(mf_cv)$importance
imp <- VarImpPlot(forest)

pi <- PartialDependence(forest, vars = names(forest$forest$variable.importance)[order(forest$forest$variable.importance, decreasing = TRUE)][1:6], 
                        plot_int = TRUE, bw=F)

set.seed(2269)
mf_rep2  <- MetaForest(yi~.,
                      data = data2,
                      whichweights = "random",
                      num.trees = 10000)
pre2 <- preselect(mf_rep2, replications = 200, algorithm = "recursive")
plot(pre2)

cv_folds2 <- trainControl(method = "cv", 10, allowParallel = TRUE)

tuning_grid2 <- expand.grid(whichweights = c("random", "fixed", "unif"),
                           mtry = 1:6,
                           min.node.size = 2:6)

X2 <- dplyr::select(data2, "vi", preselect_vars(pre2, cutoff = .95))

registerDoParallel(31)
set.seed(2269)
mf_cv2 <- train(y = data2$yi,
               x = X2,
               method = ModelInfo_mf(),
               trControl = cv_folds2,
               tuneGrid = tuning_grid2,
               keep.inbag = TRUE,
               verbose=TRUE,
               num.trees = 10000)
mf_cv2$results[which.min(mf_cv$results$RMSE), ]
r2_cv2 <- mf_cv$results$Rsquared[which.min(mf_cv$results$RMSE)]
forest2 <- mf_cv$finalModel
r2_oob2 <- forest$forest$r.squared
plot(forest2)
imp2 <- VarImpPlot(forest2)
pi2 <- PartialDependence(forest2, vars = names(forest2$forest$variable.importance)[order(forest2$forest$variable.importance, decreasing = TRUE)][1:6], 
                        plot_int = TRUE, bw=F)


# Perform meta-analysis models with moderators selected from meta-forest
# Apply transformations to obvious nonlinear variables
model1 <- rma(yi , vi, mods = ~ hdi + hdi_grow + gdp_grow + I(gdp_grow^2)+ gdp_indu +
                protect + I(protect^2)+ greendif, data = link1 , method = 'ML')
summary(model1)
# Retain only significant variables for analysis
model11 <- rma(yi, vi, mods = ~ hdi + hdi_change + gdp_indu, data = link1, method = 'ML')
summary(model11)
# Compare full and pruned models
anova(model1, model11)

model2 <- rma(yi, vi, mods = ~ hdi + gdp_indu + pmdif2 + I(pmdif2^2), 
              data = link2, method = 'ML')
summary(model2)
model22 <- rma(yi, vi, mods = ~ hdi + pmdif2 + I(pmdif2^2), 
               data = link2, method = 'ML')
summary(model22)
anova(model2, model22)

