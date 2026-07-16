# Functions

# Subset to desired fixrate ----

# TESTING
# data = dat %>% 
#   filter(ID == 1806) %>% 
#   make_track(x_, y_, t_, crs = 32612, all_cols = TRUE)

# Arguments
#   - data: one individual's data as a track_xyt
#   - rate: desired fixrate, see ?seq.POSIXt
subset_fixrate <- function(data,
                           rate = "5 hours",
                           starts = c(1, 2, 3)) {
  # IMPORTANT: assume we can round times by simply setting minutes
  # and seconds to 00:00
  
  data$tr <- data$t_
  lubridate::minute(data$tr) <- 0
  lubridate::second(data$tr) <- 0
  
  # Sequences of desired fixes
  frs <- lapply(starts, function(s) {
    seq(from = data$tr[s], to = data$tr[nrow(data)], by = rate)
  })
  
  # See how many actually appear in the data
  fr_count <- sapply(starts, function(s) {
    x <- frs[[which(s == starts)]]
    return(sum(x %in% data$tr))
  })
  
  # Choose the one with the max
  # If there are ties, choose the first
  fr_which <- which(fr_count == max(fr_count))[1]
  fix_times <- frs[[fr_which]]
  
  # Keep just the data for these fixes
  res <- data[which(data$tr %in% fix_times), ]
  
  # Return
  return(res)
}

# Example
# dat %>%
#   filter(ID == 1806) %>%
#   make_track(x_, y_, t_, crs = 32612, all_cols = TRUE) %>%
#   subset_fixrate(rate = "1 hour") %>%
#   summarize_sampling_rate()

# Attach covariates ----
# For all these functions:
#   x -- a data.frame of an individual elk's steps (incl random steps)

# ... attach ID ----
# Just adds an elk's ID to the data
add_id <- function(x, id){
  x$ID <- id
  return(x)
}

# ... attach winter ----
# Just adds the winter to the data
add_winter <- function(x, winter){
  x$winter <- winter
  return(x)
}

# ... attach snow data ----
# Specific to a day
attach_snow <- function(x, verbose = FALSE) {
  if(verbose) {
    cat("\nElk ID:", x$ID[1], " Attaching snow...\n")
  }
  # Load terra package
  library(terra)
  # Give all rows a unique ID
  x$row_id_ <- 1:nrow(x)
  # Start time
  x$yr1 <- lubridate::year(x$t1_)
  x$yd1 <- stringr::str_pad(lubridate::yday(x$t1_), 
                            width = 3, side = "left", pad = "0")
  # End time
  x$yr2 <- lubridate::year(x$t2_)
  x$yd2 <- stringr::str_pad(lubridate::yday(x$t2_), 
                            width = 3, side = "left", pad = "0")
  
  # Combine into year and yday
  x$yyd1 <- paste(x$yr1, x$yd1, sep = "_")
  x$yyd2 <- paste(x$yr2, x$yd2, sep = "_")
  
  # Split into list for each day
  sp1 <- split(x, x$yyd1)
  sp2 <- split(x, x$yyd2)
  
  # Attach all start snow
  snow1 <- do.call(rbind,
                   lapply(sp1, function(s) {
                     
                     # Start raster file path
                     fn1 <- paste0("../elk_RT_data/geo/",
                                   "daymet_swe/", s$yyd1[1], "_swe.tif")
                     
                     # Day 366 doesn't exist in Daymet data
                     if (s$yd1[1] == 366) {
                       s$swe_start <- NA
                     } else {
                       # As long as day is not 366
                       try({
                         # Load raster
                         r <- terra::rast(fn1)
                         names(r) <- "swe_start"
                         # Attach
                         s <- amt::extract_covariates(s, r, where = "start")
                       })
                       if(!exists("swe_start", where = s)){
                         s$swe_start <- NA
                       }
                     }
                     
                     # Keep just row_id_ and swe_start
                     s <- s[, c("row_id_", "swe_start")]
                     # Return
                     return(s)
                   })
  )
  
  # Attach all end snow
  snow2 <- do.call(rbind, 
                   lapply(sp2, function(s) {
                     # End raster file path
                     fn2 <- paste0("../elk_RT_data/geo/", 
                                   "daymet_swe/", s$yyd2[2], "_swe.tif")
                     
                     # Day 366 doesn't exist in Daymet data
                     if (s$yd2[1] == 366) {
                       s$swe_end <- NA
                     } else {
                       try({
                         # Load raster
                         r <- terra::rast(fn2)
                         names(r) <- "swe_end"
                         # Attach
                         s <- amt::extract_covariates(s, r, where = "end")
                       })
                       if(!exists("swe_end", where = s)){
                         s$swe_end <- NA
                       }
                     }
                     # Keep just row_id_ and swe_end
                     s <- s[, c("row_id_", "swe_end")]
                     # Return
                     return(s)
                   })
  )
  
  # Join start and end to data
  x <- dplyr::left_join(x, snow1, by = "row_id_")
  x <- dplyr::left_join(x, snow2, by = "row_id_")
  
  # Remove row_id_
  x$row_id_ <- NULL
  
  # Remove date columns
  x$yr1 <- x$yr2 <- x$yd1 <- x$yd2 <- x$yyd1 <- x$yyd2 <- NULL
  
  # Return
  return(x)
}

# ... attach RAP data ----
# Specific to a year
attach_RAP <- function(x, verbose = FALSE) {
  if(verbose) {
    cat("\nElk ID:", x$ID[1], " Attaching RAP...\n")
  }
  # Load terra package
  library(terra)
  # Winter
  wint <- x$winter[1]
  
  # Load biomass rasters
  ann <- terra::rast(paste0("../elk_RT_data/geo/RAP/", 
                            "biomass_annuals_raster_stack.tif"))
  per <- terra::rast(paste0("../elk_RT_data/geo/RAP/", 
                            "biomass_perennials_raster_stack.tif"))
  
  # Load openness raster
  open <- terra::rast(paste0("../elk_RT_data/geo/RAP/", 
                             "openness_raster_stack.tif"))
  
  # Attach RAP
  # Raster names
  nm <- paste0("RAP", wint)
  # Sum biomass rasters
  b <- ann[[nm]] + per[[nm]]
  names(b) <- "bio"
  # Get openness raster
  o <- open[[nm]]
  names(o) <- "open"
  
  # Attach
  x <- amt::extract_covariates(x, b, where = "both")
  x <- amt::extract_covariates(x, o, where = "both")
  
  # Return
  return(x)
}

# ... attach DEM data ----
# Temporally static
attach_DEM <- function(x, verbose = FALSE) {
  if(verbose) {
    cat("\nElk ID:", x$ID[1], " Attaching DEM...\n")
  }
  # Load terra package
  library(terra)
  
  # Load elevation
  e <- terra::rast("../elk_RT_data/geo/DEM/dem_raster.tif")
  names(e) <- "elev"
  # Load roughness
  r <- terra::rast("../elk_RT_data/geo/DEM/roughness_raster.tif")
  names(r) <- "rough"
  # Load cosine aspect
  c <- terra::rast("../elk_RT_data/geo/DEM/cos_aspect_raster.tif")
  names(c) <- "cos_asp"
  # Load sine aspect
  s <- terra::rast("../elk_RT_data/geo/DEM/sin_aspect_raster.tif")
  names(s) <- "sin_asp"
  
  # Attach
  x <- amt::extract_covariates(x, e, where = "both")
  x <- amt::extract_covariates(x, r, where = "both")
  x <- amt::extract_covariates(x, c, where = "both")
  x <- amt::extract_covariates(x, s, where = "both")
  
  # Return
  return(x)
}

# ... attach hour transformations ----

# h can be a vector
# season should be length 1
cougar_act <- function(h, season) {
  if (!(season %in% c("early", "late"))) {
    stop("Argument 'season' must be 'early' or 'late'.")
  }
  if (length(season) != 1) {
    stop("Argument 'season' must be length 1.")
  }
  act <- read.csv("data/cougar_activity_table.csv")
  return(act[match(h, act$h), paste0("index_", season)])
}

wolf_act <- function(h, season) {
  if (!(season %in% c("early", "late"))) {
    stop("Argument 'season' must be 'early' or 'late'.")
  }
  if (length(season) != 1) {
    stop("Argument 'season' must be length 1.")
  }
  act <- read.csv("data/wolf_activity_table.csv")
  return(act[match(h, act$h), paste0("index_", season)])
}

# This is NOT vectorized
# This function assigns values in radians between -pi and pi for all
# hours of any given day. The difference between this and something
# like suncalc::getSunlightTimes()$altitude is that it ignores the
# fact that the sun gets higher on some days than others; i.e., this
# function always returns pi at solar noon, not something less than that
# on non-summer-solstice days.
sun_position <- function(time, x, y, tz = "US/Mountain") {
  date <- as.Date(time)
  
  # Suntimes for previous day, focal day, and next day
  st <- suncalc::getSunlightTimes(date = date + (-1:1),
                                  lat = y,
                                  lon = x) %>% 
    lubridate::with_tz(tz)
  
  # Probably a more efficient way to do this, but this is what came to me
  # first. Seems fine unless time cost is prohibitive
  # Previous nadir (smallest negative number)
  nadir_diff <- (st$nadir - time)
  nadir_neg <- which(nadir_diff < 0)
  #[1] takes care of exact ties
  nadir <- st$nadir[nadir_neg[which(nadir_neg == max(nadir_neg))[1]]]
  # Round to nearest minute
  nadir <- round(nadir, units = "mins")
  # Sunrise (closest after nadir = smallest positive number)
  sunrise_diff <- st$sunrise - nadir
  sunrise_pos <- which(sunrise_diff > 0)
  sunrise <- st$sunrise[sunrise_pos[which(sunrise_pos == min(sunrise_pos))[1]]]
  # Round to nearest minute
  sunrise <- round(sunrise, units = "mins")
  # Solar noon (closest after sunrise = smallest positive number)
  noon_diff <- st$solarNoon - sunrise
  noon_pos <- which(noon_diff > 0)
  noon <- st$solarNoon[noon_pos[which(noon_pos == min(noon_pos))[1]]]
  # Round to nearest minute
  noon <- round(noon, units = "mins")
  # Sunset (next after solar noon = smallest positive number)
  sunset_diff <- st$sunset - noon
  sunset_pos <- which(sunset_diff > 0)
  sunset <- st$sunset[sunset_pos[which(sunset_pos == min(sunset_pos))[1]]]
  # Round to nearest minute
  sunset <- round(sunset, units = "mins")
  # Next nadir(next after sunset = smallest positive number)
  nadir2_diff <- st$nadir - sunset
  nadir2_pos <- which(nadir_diff > 0)
  nadir2 <- st$nadir[nadir2_pos[which(nadir2_pos == min(nadir2_pos))[1]]]
  # Round to nearest minute
  nadir2 <- round(nadir2, units = "mins")
  
  # Evenly space intervals between solar events
  pi0 <- seq(pi, 0, length.out = 300)
  mpi0 <- seq(-pi, 0, length.out = 300)
  
  # Nadir to sunrise
  nd_sr <- data.frame(period = "nadir_to_sunrise",
                      time = seq(from = nadir, to = sunrise, by = "1 min"))
  nd_sr$radians <- seq(from = -pi,
                       to = 0,
                       length.out = nrow(nd_sr))
  # Sunrise to solar noon
  sr_sn <- data.frame(period = "sunrise_to_solarnoon",
                      time = seq(from = sunrise, to = noon, by = "1 min"))
  sr_sn$radians <- seq(from = 0,
                       to = pi,
                       length.out = nrow(sr_sn))
  # Solar noon to sunset
  sn_ss <- data.frame(period = "solarnoon_to_sunset",
                      time = seq(from = noon, to = sunset, by = "1 min"))
  sn_ss$radians <- seq(from = pi,
                       to = 0,
                       length.out = nrow(sn_ss))
  # Sunset to nadir
  ss_nd2 <- data.frame(period = "sunset_to_nadir",
                       time = seq(from = sunset, to = nadir2, by = "1 min"))
  ss_nd2$radians <- seq(from = 0, 
                        to = -pi,
                        length.out = nrow(ss_nd2))
  # Combine
  rads <- rbind(nd_sr, sr_sn, sn_ss, ss_nd2)
  
  # On some dates, there are two matches for a rounded time
  # Keep the first
  if (max(table(rads$time)) > 1) {
    res <- rads[which(!duplicated(rads$time)), ]
    return(res)
  } else {
    # Return
    return(rads)
  }
  
}

# # Examples
# library(ggplot2)
# library(patchwork)
# t <- lubridate::ymd_hms("2020-12-21 12:00:00", tz = "US/Mountain")
# (sp <- sun_position(time = t, x = -110.4089, y = 44.98101))
# sp$Period <- factor(sp$period,
#                    levels = c("nadir_to_sunrise",
#                               "solarnoon_to_sunset",
#                               "sunrise_to_solarnoon",
#                               "sunset_to_nadir"),
#                    labels = c("Nadir to Sunrise",
#                               "Solar Noon to Sunset",
#                               "Sunrise to Solar Noon",
#                               "Sunset to Nadir"))
# # Plot radians
# (a <- ggplot(sp, aes(x = time, y = radians, color = Period)) +
#     geom_line(linewidth = 1) +
#     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
#     # annotate(geom = "label", label = "Focal Time", x = t, y = 0) +
#     scale_x_datetime(name = "Time", date_labels = "%H:%M") +
#     scale_y_continuous(breaks = c(-pi, -pi/2, 0, pi/2, pi),
#                        labels = expression(-pi, -pi/2, 0, pi/2, pi)) +
#     theme_bw())
# # Plot cos(radians)
# (b <- ggplot(sp, aes(x = time, y = cos(radians), color = Period)) +
#     geom_line(linewidth = 1) +
#     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
#     # annotate(geom = "label", label = "Focal Time", x = t, y = 0) +
#     scale_x_datetime(name = "Time", date_labels = "%H:%M") +
#     theme_bw())
# # Plot sin(radians)
# (c <- ggplot(sp, aes(x = time, y = sin(radians), color = Period)) +
#     geom_line(linewidth = 1) +
#     geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
#     # annotate(geom = "label", label = "Focal Time", x = t, y = 0) +
#     scale_x_datetime(name = "Time", date_labels = "%H:%M") +
#     theme_bw())
# # Combine
# sol_time <- a + b + c +
#   plot_layout(nrow = 1, guides = "collect") +
#   plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave("fig/solar_time_explainer.tif",
#        plot = sol_time, device = ragg::agg_tiff, width = 7, height = 3.5,
#        units = "in", dpi = 300, compression = "lzw")

attach_hour <- function(x, verbose = FALSE) {
  if(verbose) {
    cat("\nElk ID:", x$ID[1], " Attaching hour...\n")
  }
  
  # Round start time to the nearest minute
  x$t1_round <- round.POSIXt(x$t1_, units = "mins")
  
  # Hours
  x$hour_start <- lubridate::hour(x$t1_)
  x$hour_end <- lubridate::hour(x$t2_)
  
  # Cougar activity
  x$cougar_act_start <- cougar_act(x$hour_start, x$season[1])
  x$cougar_act_end <- cougar_act(x$hour_end, x$season[1])
  
  # Wolf activity
  x$wolf_act_start <- wolf_act(x$hour_start, x$season[1])
  x$wolf_act_end <- wolf_act(x$hour_end, x$season[1])
  
  # Solar times
  # Only for start time, which is the same for all strata
  s <- unique(x[, c("step_id_", "x1_", "y1_", "t1_round")])
  # Geographic coordinates
  utm_start <- sf::st_sfc(sf::st_multipoint(cbind(s$x1_, s$y1_)), crs = 32612)
  geo_start <- sf::st_coordinates(
    sf::st_transform(utm_start, 4326))[, 1:2, 
                                       drop = FALSE]
  s <- cbind(s, geo_start)
  
  # This is not vectorized, so split first
  sl <- split(s, 1:nrow(s))
  
  sts <- lapply(sl, function(xx) {
    sp <- sun_position(time = xx$t1_round,
                       x = xx$X,
                       y = xx$Y)
    xxx <- xx %>% 
      select(-X, -Y) %>% 
      left_join(sp, by = c("t1_round" = "time"))
    return(xxx)
  }) %>% 
    bind_rows()
  
  # Join
  x <- left_join(x, sts, by = c("step_id_", "x1_", "y1_", "t1_round")) %>% 
    # rename radians to solar time
    rename(solar_time = radians)
  
  # Return
  return(x)
}

# ... attach predation risk ----
attach_risk <- function(x, verbose = TRUE) {
  if(verbose) {
    cat("\nElk ID:", x$ID[1], " Attaching risk...\n")
  }
  # Get winter and season
  w <- x$winter[1]
  s <- x$season[1]
  
  # Load cougar risk
  cr <- terra::rast(paste0("../elk_RT_data/geo/cougar_risk/cougar_risk_",
                           w, "_", s, ".tif"))
  names(cr) <- "cougar_risk"
  # Load wolf risk
  wr <- terra::rast(paste0("../elk_RT_data/geo/wolf_risk/wolf_risk_",
                           w, "_", s, ".tif"))
  names(wr) <- "wolf_risk"
  
  # Combine
  pred <- c(cr, wr)
  
  # Attach
  x <- amt::extract_covariates(x, pred, where = "both")
  
  # Return
  return(x)
}

# Scale and center ----
scale_dat <- function(dat, scale_df) {
  # If terms with "_orig" already exist, the data have already been scaled
  if (any(grepl("_orig", names(dat)))) {
    warning(paste("Columns named '*_orig' already exist.",
                  "Assuming data have already been scaled. Returning input."))
    return(dat)
  }
  
  for (i in 1:nrow(scale_df)) {
    term <- scale_df$term[i]
    mu <- scale_df$mean[i]
    sig <- scale_df$sd[i]
    # Copy unscaled data as "<term>_orig"
    orig_name <- paste0(term, "_orig")
    dat[[orig_name]] <- dat[[term]]
    # Scale and center
    dat[[term]] <- (dat[[term]] - mu)/sig
  }
  return(dat)
}

# Get solar_time ----
# x and y default to a location in Mammoth
# Vectorized over times
get_solar_time <- function(time, 
                           x = -110.697706, y = 44.975903, 
                           tz = "US/Mountain") {
  # Vectorized across times
  l <- lapply(time, function(t) {
    tt <- round.POSIXt(t, units = "mins")
    xx <- sun_position(time = tt, x = x, y = y, tz = tz)
    yy <- xx[which(xx$time == tt)[1], ]
    return(yy)
  })
  res <- do.call(rbind, l)
  names(res)[3] <- "solar_time"
  
  # Double-check to avoid hidden problems later
  if (nrow(res) != length(time)) {
    stop("Problem in get_solar_time(). ", 
         "Rows of result do not match length of input.")
  }
  return(res)
}

# Create prediction data ----
create_pred_dat <- function(scale_df,
                            mod = NULL,
                            bio_start         = scale_df["bio_start", "mean"],
                            bio_end           = scale_df["bio_end", "mean"],
                            open_start        = scale_df["open_start", "mean"],
                            open_end          = scale_df["open_end", "mean"],
                            rough_start       = scale_df["rough_start", "mean"],
                            rough_end         = scale_df["rough_end", "mean"],
                            wolf_risk_start   = scale_df["wolf_risk_start", "mean"],
                            wolf_risk_end     = scale_df["wolf_risk_end", "mean"],
                            cougar_risk_start = scale_df["cougar_risk_start", "mean"],
                            cougar_risk_end   = scale_df["cougar_risk_end", "mean"],
                            elev_start        = scale_df["elev_start", "mean"],
                            elev_end          = scale_df["elev_end", "mean"],
                            swe_start         = scale_df["swe_start", "mean"],
                            swe_end           = scale_df["swe_end", "mean"],
                            sl_               = 100,
                            ta_               = 0,
                            start_time        = lubridate::ymd_hms(
                              "2020-03-01 00:00:00", tz = "US/Mountain"),
                            duration          = lubridate::hours(1)
                            
) {
  # Base combination
  res <- expand.grid(bio_start         = bio_start,
                     bio_end           = bio_end, 
                     open_start        = open_start,
                     open_end          = open_end, 
                     rough_start       = rough_start,
                     rough_end         = rough_end,  
                     wolf_risk_start   = wolf_risk_start,   
                     wolf_risk_end     = wolf_risk_end,     
                     cougar_risk_start = cougar_risk_start,  
                     cougar_risk_end   = cougar_risk_end,    
                     elev_start        = elev_start,   
                     elev_end          = elev_end,     
                     swe_start         = swe_start,    
                     swe_end           = swe_end,      
                     sl_               = sl_,          
                     ta_               = ta_,          
                     start_time        = start_time)
  
  # If we passed the model, use one of the step_id_ from the model
  # If not, just fill in "2"
  if (is.null(mod)) {
    sid <- "2"
  } else {
    sid <- mod$model$model[["strata(step_id_)"]][1]
  }
  res$step_id_ <- sid
  
  # Functions of time
  res$end_time <- res$start_time + duration
  
  # Solar time
  # Using x, y of Mammoth (defaults for 'get_solar()')
  res$start_time_round <- round.POSIXt(res$start_time, units = "mins")
  # Only get the solar time once for each unique start time
  # (using a join next that could duplicate)
  ss <- get_solar_time(unique(res$start_time_round))
  
  res <- dplyr::left_join(res, ss, by = c("start_time_round" = "time"))
  
  # Season
  res$season <- ifelse(lubridate::month(start_time) %in% 11:12,
                       "early", "late")
  
  # Predators
  res$wolf_act_start <- wolf_act(lubridate::hour(res$start_time), 
                                 season = res$season[1])
  res$wolf_act_end <- wolf_act(lubridate::hour(res$end_time), 
                               season = res$season[1])
  res$cougar_act_start <- cougar_act(lubridate::hour(res$start_time), 
                                     season = res$season[1])
  res$cougar_act_end <- cougar_act(lubridate::hour(res$end_time), 
                                   season = res$season[1])
  
  # log(SWE)
  res$log_swe_start <- log(swe_start + 1)
  res$log_swe_end <- log(swe_end + 1)
  
  return(res)
}

# Predict bio coefficient ----
# ... this is the old version without boostrapping ----
predict_bio_coef <- function(newdata, m) {
  
  ## Predict beta_bio_end
  # Relevant coefficients
  b <- coef(m$model)[grep("bio", names(coef(m$model)), fixed = TRUE)]
  # Create formula
  f <- as.formula(paste0("~ 0 + ", paste(names(b), collapse = " + ")))
  # Create model matrix
  mm <- model.matrix(f, newdata)
  # Calculate prediction
  p <- mm %*% b
  # Divide out bio_end
  beta_bio_end <- (p/newdata$bio_end)[,1]
  
  # Combine in data.frame
  res <- newdata
  res$beta_bio_end <- beta_bio_end
  
  # Return
  return(res)
}


# Update step-length distribution ----
# ... this is the old version without bootstrapping ----
predict_sl <- function(newdata, m) {
  # Tentative distribution
  tent_sl <- m$sl_
  
  ## Predict beta_log_sl
  # Relevant coefficients
  b <- coef(m$model)[grep("log(sl_)", names(coef(m$model)), fixed = TRUE)]
  # Create formula
  f <- as.formula(paste0("~ 0 + ", paste(names(b), collapse = " + ")))
  # Create model matrix
  mm <- model.matrix(f, newdata)
  # Calculate prediction
  p <- mm %*% b
  # Divide out log(sl_)
  beta_log_sl <- (p/log(newdata$sl_))[,1]
  
  ## Update tentative distribution
  # Update shape
  shp <- tent_sl$params$shape + beta_log_sl
  # Not updating scale
  scl <- tent_sl$params$scale
  
  # Combine in data.frame
  res <- newdata
  res$shp <- shp
  res$scl <- scl
  res$mean_sl <- shp * scl
  
  # Return
  return(res)
}


# ... this is the new version with bootstrapping ----
# Note: the start date-time should probably not vary in 'newdata'.
# The function will subset 'model_df' to match the season of the date-time
# in 'newdata' and is expecting a single value.

# Note: this function is specific to the model formula used; it would
# not work correctly if the model formula were changed.

predict_sl <- function(model_df, niter, newdata, verbose = TRUE) {
  
  # Bootstrap iterations
  bb <- lapply(1:niter, function(i) {
    if (verbose) {
      cat("Iteration", i, "of", niter, "              \r")
    }
    # Subset model_df to just those seasons relevant to the time in newdata
    uni_seas <- unique(newdata$season)
    model_df <- model_df %>% 
      filter(season %in% uni_seas)
    
    # In each bootstrap iteration, return shape, scale, and mean sl for all 
    # individual-winter-seasons
    inds <- lapply(1:nrow(model_df), function(j) {
      # Get model
      mod <- model_df$issf[[j]]$model
      
      # Need dummy value of step_id_ to agree with levels in the data
      # This will not affect calculation
      newdata$step_id_ <- mod$model$`strata(step_id_)`[1]
      
      # Create model matrix for the beta-star for log(sl_)
      X <- as.matrix(data.frame(int = 1,
                 `cos(ta_)` = cos(newdata$ta_),
                 log_swe_start = newdata$log_swe_start,
                 `sin(solar_time)` = sin(newdata$solar_time),
                 `cos(solar_time)` = cos(newdata$solar_time)))
      
      # Order of coefficients for beta-star
      bs_coef <- c("log(sl_)", 
                   "log(sl_):cos(ta_)",
                   "log(sl_):log_swe_start",
                   "log(sl_):sin(solar_time)",
                   "log(sl_):cos(solar_time)")
      
      # Resample coefficients
      B <- mvtnorm::rmvnorm(n = 1, 
                            mean = coef(mod)[bs_coef], 
                            sigma = vcov(mod)[bs_coef, bs_coef])[1, ]
      
      # Calculate beta-star for log(sl)
      bs_lsl <- unname((X %*% B)[, 1])
      
      
      # sl_ is not included in the model formula, interactions, so 
      # its "beta-star" is assumed to be 0, but replicate so it matches the
      # length of the beta-star for log(sl)
      bs_sl <- rep(0, length(bs_lsl))
      
      # Update tentative distribution
      tent_sl <- sl_distr(model_df$issf[[j]])
      upd_sl <- update_gamma(tent_sl,
                             beta_sl = bs_sl,
                             beta_log_sl = bs_lsl)
      
      # Combine with model information and data
      res <- cbind(
        data.frame(ID = rep(model_df$ID[[j]], length(bs_lsl)),
                   winter = model_df$winter[[j]],
                   iter = i),
        newdata,
        shape = upd_sl$params$shape,
        scale = upd_sl$params$scale
      )
      res$mean_sl <- res$shape * res$scale
      
      return(res)
    }) %>% 
      # Combine all individuals
      bind_rows()
    
    # Return from bootstrap iteration
    return(inds)
    
  }) %>% 
    # Combine all bootstrap iterations
    bind_rows()
  
  # Return
  return(bb)
  
}
  

# Split ID_winter_season column ----
split_ID_winter_season <- function(data) {
  # Split
  s <- strsplit(data$ID_winter_season, split = "_")
  # Get IDs
  ID <- as.numeric(sapply(s, getElement, 1))
  # Get winters
  winter <- as.numeric(sapply(s, getElement, 2))
  # Get seasons
  season <- sapply(s, getElement, 3)
  # Manipulate data
  data <- data %>% 
    dplyr::select(-ID_winter_season) %>% 
    dplyr::mutate(ID = ID, 
                  winter = as.numeric(winter),
                  season = season) %>% 
    dplyr::select(ID, winter, season, dplyr::everything())
  # Return
  return(data)
}

# Bootstrap resampling of iSSA betas ----
bootstrap_betas <- function(data, col, iter = 2000, verbose = TRUE) {
  LL <- lapply(1:iter, function(i) {
    # Possibly report status
    if (verbose) {
      cat("\r", i, "of", iter, "     ")
    }
    # Resample betas for all individuals
    ll <- lapply(data[[col]], resamp_beta)
    # Add combo of ID, winter, and season as names
    names(ll) <- paste(data$ID, data$winter, data$season, sep = "_")
    # Combine list elements in data.frame
    ii <- dplyr::bind_rows(ll, .id = "ID_winter_season")
    # Label iteration
    ii$iter <- i
    # Return from single iteration
    return(ii)
  })
  
  res <- dplyr::bind_rows(LL)
  return(res)
}

resamp_beta <- function(model) {
  # Get beta estimates
  b <- coef(model$model)
  # Get variance-covariance matrix
  S <- vcov(model$model)
  # Sample new betas
  new <- mvtnorm::rmvnorm(n = 1, mean = b, sigma = S)[1,]
  # Build data.frame
  df <- data.frame(term = names(new), estimate = unname(new))
  # Return
  return(df)
}

boot_rss <- function(model_df, niter, x1, x2, verbose = TRUE) {
  # Check that x1 and x2 have exactly the same number of rows.
  # Note this is less strict than amt::log_rss() which always requires
  # x2 to have exactly 1 row.
  if (!identical(nrow(x1), nrow(x2))) {
    stop("'x1' and 'x2' must have exactly the same number of rows.")
  }
  
  # Bootstrap iterations
  bb <- lapply(1:niter, function(i) {
    if (verbose) {
      cat("Iteration", i, "of", niter, "              \r")
    }
    # In each bootstrap iteration, return log_rss(x1, x2) for all 
    # individual-winter-seasons
    inds <- lapply(1:nrow(model_df), function(j) {
      # Get model
      mod <- model_df$issf[[j]]$model
      
      # Subset x1 and x2 to correct season
      X1 <- x1[which(x1$season == model_df$season[[j]]), ]
      X2 <- x2[which(x2$season == model_df$season[[j]]), ]
      
      # Need dummy value of step_id_ to agree with levels in the data
      # This will not affect calculation
      X1$step_id_ <- mod$model$`strata(step_id_)`[1]
      X2$step_id_ <- mod$model$`strata(step_id_)`[1]
      
      # Create model matrices
      mm1 <- model.matrix(mod, X1)
      mm2 <- model.matrix(mod, X2)
      
      # Get matrix of differences (x1 - x2)
      delta_mm <- mm1 - mm2
      
      # Resample coefficients
      B <- mvtnorm::rmvnorm(n = 1, mean = coef(mod), sigma = vcov(mod))[1, ]
      
      # Calculate log-RSS
      log_rss <- unname((delta_mm %*% B)[, 1])
      
      # Combine identifiers with data and log-RSS
      res <- cbind(
        data.frame(ID = rep(model_df$ID[[j]], length(log_rss)),
                   winter = model_df$winter[[j]],
                   iter = i),
        X1,
        log_rss = log_rss
      )
      
      return(res)
    }) %>% 
      # Combine all individuals
      bind_rows()
    
    # Return from bootstrap iteration
    return(inds)
    
  }) %>% 
    # Combine all bootstrap iterations
    bind_rows()
  
  # Return
  return(bb)
  
}

# Plot risk taking ----

risk_plots <- function(x_dat, var, xlab, 
                       wolf_ylim = c(NA, NA),
                       cougar_ylim = c(NA, NA)){
  # Summarize bootstrap iterations
  x_summ <- x_dat %>% 
    group_by(model, {{ var }}) %>% 
    summarize(mean = mean(log_rss),
              lwr = quantile(log_rss, 0.025),
              upr = quantile(log_rss, 0.975))
  
  # ... ... wolf_mean
  
  x_wolf_mean_plot <- x_dat %>% 
    filter(model == "wolf_mean") %>% 
    ggplot() +
    geom_line(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_ribbon(data = filter(x_summ, model == "wolf_mean"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                linetype = "dashed",
                fill = NA,
                color = "#fcbe11", linewidth = 0.7) +
    geom_line(data = filter(x_summ, model == "wolf_mean"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", linewidth = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(xlab) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = wolf_ylim) +
    # theme(axis.text.x = element_blank()) +
    # ggtitle("Wolf")
    NULL
  
  # ... ... wolf_max
  
  x_wolf_max_plot <- x_dat %>% 
    filter(model == "wolf_max") %>% 
    ggplot() +
    geom_line(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_ribbon(data = filter(x_summ, model == "wolf_max"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                linetype = "dashed",
                fill = NA,
                color = "#fcbe11", linewidth = 0.7) +
    geom_line(data = filter(x_summ, model == "wolf_max"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", linewidth = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = wolf_ylim) +
    theme(axis.text.x = element_blank()) +
    NULL
  
  # ... ... wolf_min
  
  x_wolf_min_plot <- x_dat %>% 
    filter(model == "wolf_min") %>% 
    ggplot() +
    geom_line(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_ribbon(data = filter(x_summ, model == "wolf_min"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                linetype = "dashed",
                fill = NA,
                color = "#fcbe11", linewidth = 0.7) +
    geom_line(data = filter(x_summ, model == "wolf_min"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", linewidth = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = wolf_ylim) +
    theme(axis.text.x = element_blank()) +
    NULL
  
  # ... ... cougar_mean
  
  x_cougar_mean_plot <- x_dat %>% 
    filter(model == "cougar_mean") %>% 
    ggplot() +
    geom_line(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_ribbon(data = filter(x_summ, model == "cougar_mean"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                linetype = "dashed",
                fill = NA,
                color = "#fcbe11", linewidth = 0.7) +
    geom_line(data = filter(x_summ, model == "cougar_mean"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", linewidth = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(xlab) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = cougar_ylim) +
    # theme(axis.text.x = element_blank()) +
    # ggtitle("Cougar") +
    NULL
  
  # ... ... cougar_max
  
  x_cougar_max_plot <- x_dat %>% 
    filter(model == "cougar_max") %>% 
    ggplot() +
    geom_line(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_ribbon(data = filter(x_summ, model == "cougar_max"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                linetype = "dashed",
                fill = NA,
                color = "#fcbe11", linewidth = 0.7) +
    geom_line(data = filter(x_summ, model == "cougar_max"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", linewidth = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = cougar_ylim) +
    theme(axis.text.x = element_blank()) +
    NULL
  
  # ... ... cougar_min
  
  x_cougar_min_plot <- x_dat %>% 
    filter(model == "cougar_min") %>% 
    ggplot() +
    geom_line(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_ribbon(data = filter(x_summ, model == "cougar_min"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                linetype = "dashed",
                fill = NA,
                color = "#fcbe11", linewidth = 0.7) +
    geom_line(data = filter(x_summ, model == "cougar_min"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", linewidth = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = cougar_ylim) +
    theme(axis.text.x = element_blank()) +
    NULL
  
  # ... ... combine
  ## Row labels
  # Row 1
  row1_title <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "Riskiest Time"),
              size = 6, angle = 90) +
    theme_void()
  
  # Row 2
  row2_title <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "Safest Time"),
              size = 6, angle = 90) +
    theme_void()
  
  # Row 3
  row3_title <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "Daily Mean"),
              size = 6, angle = 90) +
    theme_void()
  
  row_labels <- wrap_elements(full = row1_title + 
                                row2_title +
                                row3_title + 
                                plot_layout(ncol = 1, nrow = 3,
                                            heights = c(0.32, 0.32, 0.36)))
  
  col1 <- wrap_elements(
    full = 
      x_wolf_max_plot + x_wolf_min_plot + x_wolf_mean_plot +
      plot_layout(ncol = 1, nrow = 3,
                  heights = c(0.32, 0.32, 0.36), tag_level = "keep") +
      plot_annotation(title = "Wolf",
                      tag_levels = list(c("A", "C", "E")),
                      tag_prefix = "(",
                      tag_suffix = ")",
                      theme = theme(
                        plot.title = element_text(size = 18, 
                                                  hjust = 0.7))) & 
      theme(plot.tag = element_text(face = "bold", size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12))
  )
  
  col2 <- wrap_elements(
    full = 
      x_cougar_max_plot + x_cougar_min_plot + x_cougar_mean_plot +
      plot_layout(ncol = 1, nrow = 3,
                  heights = c(0.32, 0.32, 0.36), tag_level = "keep") +
      plot_annotation(title = "Cougar",
                      tag_levels = list(c("B", "D", "F")),
                      tag_prefix = "(",
                      tag_suffix = ")",
                      theme = theme(
                        plot.title = element_text(size = 18, 
                                                  hjust = 0.7))) & 
      theme(plot.tag = element_text(face = "bold", size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12))
  )
  
  x_plot <- row_labels + col1 + col2 +
    plot_layout(nrow = 1, ncol = 3,
                widths = c(0.04, 0.48, 0.48))
  
  # Return
  return(x_plot)
}

risk_plots_discrete <- function(x_dat, var, xlab, 
                                wolf_ylim = c(NA, NA),
                                cougar_ylim = c(NA, NA)){
  # Summarize bootstrap iterations
  x_summ <- x_dat %>% 
    group_by(model, {{ var }}) %>% 
    summarize(mean = mean(log_rss),
              lwr = quantile(log_rss, 0.025),
              upr = quantile(log_rss, 0.975))
  
  # ... ... wolf_mean
  
  x_wolf_mean_plot <- x_dat %>% 
    filter(model == "wolf_mean") %>% 
    ggplot() +
    geom_point(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_errorbar(data = filter(x_summ, model == "wolf_mean"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                color = "#fcbe11", width = 0.5) +
    geom_point(data = filter(x_summ, model == "wolf_mean"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(xlab) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = wolf_ylim) +
    # theme(axis.text.x = element_blank()) +
    # ggtitle("Wolf")
    NULL
  
  # ... ... wolf_max
  
  x_wolf_max_plot <- x_dat %>% 
    filter(model == "wolf_max") %>% 
    ggplot() +
    geom_point(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_errorbar(data = filter(x_summ, model == "wolf_max"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                color = "#fcbe11", width = 0.5) +
    geom_point(data = filter(x_summ, model == "wolf_max"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = wolf_ylim) +
    theme(axis.text.x = element_blank()) +
    NULL
  
  # ... ... wolf_min
  
  x_wolf_min_plot <- x_dat %>% 
    filter(model == "wolf_min") %>% 
    ggplot() +
    geom_point(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_errorbar(data = filter(x_summ, model == "wolf_min"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                color = "#fcbe11", width = 0.5) +
    geom_point(data = filter(x_summ, model == "wolf_min"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = wolf_ylim) +
    theme(axis.text.x = element_blank()) +
    NULL
  
  # ... ... cougar_mean
  
  x_cougar_mean_plot <- x_dat %>% 
    filter(model == "cougar_mean") %>% 
    ggplot() +
    geom_point(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_errorbar(data = filter(x_summ, model == "cougar_mean"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                color = "#fcbe11", width = 0.5) +
    geom_point(data = filter(x_summ, model == "cougar_mean"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(xlab) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = cougar_ylim) +
    # theme(axis.text.x = element_blank()) +
    # ggtitle("Cougar") +
    NULL
  
  # ... ... cougar_max
  
  x_cougar_max_plot <- x_dat %>% 
    filter(model == "cougar_max") %>% 
    ggplot() +
    geom_point(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_errorbar(data = filter(x_summ, model == "cougar_max"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                color = "#fcbe11", width = 0.5) +
    geom_point(data = filter(x_summ, model == "cougar_max"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = cougar_ylim) +
    theme(axis.text.x = element_blank()) +
    NULL
  
  # ... ... cougar_min
  
  x_cougar_min_plot <- x_dat %>% 
    filter(model == "cougar_min") %>% 
    ggplot() +
    geom_point(aes(x = {{ var }}, y = log_rss, group = iter),
              alpha = 0.05) +
    geom_errorbar(data = filter(x_summ, model == "cougar_min"),
                aes(x = {{ var }}, ymin = lwr, ymax = upr),
                color = "#fcbe11", width = 0.5) +
    geom_point(data = filter(x_summ, model == "cougar_min"),
              aes(x = {{ var }}, y = mean),
              color = "#fcbe11", size = 1) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    xlab(NULL) +
    ylab("Risk Taking") +
    coord_cartesian(ylim = cougar_ylim) +
    theme(axis.text.x = element_blank()) +
    NULL
  
  # ... ... combine
  ## Row labels
  # Row 1
  row1_title <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "Riskiest Time"),
              size = 6, angle = 90) +
    theme_void()
  
  # Row 2
  row2_title <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "Safest Time"),
              size = 6, angle = 90) +
    theme_void()
  
  # Row 3
  row3_title <- ggplot() +
    geom_text(aes(x = 0, y = 0, label = "Daily Mean"),
              size = 6, angle = 90) +
    theme_void()
  
  row_labels <- wrap_elements(full = row1_title + 
                                row2_title +
                                row3_title + 
                                plot_layout(ncol = 1, nrow = 3,
                                            heights = c(0.32, 0.32, 0.36)))
  
  col1 <- wrap_elements(
    full = 
      x_wolf_max_plot + x_wolf_min_plot + x_wolf_mean_plot +
      plot_layout(ncol = 1, nrow = 3,
                  heights = c(0.32, 0.32, 0.36), tag_level = "keep") +
      plot_annotation(title = "Wolf",
                      tag_levels = list(c("A", "C", "E")),
                      tag_prefix = "(",
                      tag_suffix = ")",
                      theme = theme(
                        plot.title = element_text(size = 18, 
                                                  hjust = 0.7))) & 
      theme(plot.tag = element_text(face = "bold", size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12))
  )
  
  col2 <- wrap_elements(
    full = 
      x_cougar_max_plot + x_cougar_min_plot + x_cougar_mean_plot +
      plot_layout(ncol = 1, nrow = 3,
                  heights = c(0.32, 0.32, 0.36), tag_level = "keep") +
      plot_annotation(title = "Cougar",
                      tag_levels = list(c("B", "D", "F")),
                      tag_prefix = "(",
                      tag_suffix = ")",
                      theme = theme(
                        plot.title = element_text(size = 18, 
                                                  hjust = 0.7))) & 
      theme(plot.tag = element_text(face = "bold", size = 14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12))
  )
  
  x_plot <- row_labels + col1 + col2 +
    plot_layout(nrow = 1, ncol = 3,
                widths = c(0.04, 0.48, 0.48))
  
  # Return
  return(x_plot)
}

# Summarize bootstrap "significance" ----

# Get difference in log-RSS for min and max for each iteration
#   -- Tricky for the ones that are non-monotonic
#   -- Rather than using min/max of 'var', find the value of var with
#       min/max of log-RSS (but sorted by 'var').
#   -- Argument 'mono' decides this; mono = FALSE uses this option, mono = TRUE
#       uses min/max of 'var'. (short for monotonic)
iter_diff <- function(x, var, mono = FALSE) {
  # Get data.frame
  df <- x[[var]]
  # Split by model
  l <- split(df, df$model)
  # Get effect by model
  e <- lapply(l, function(ll) {
    if (mono) {
      # Get min of var
      mn <- ll[which(ll[[var]] == min(ll[[var]], na.rm = TRUE)), ]
      # Get max of var
      mx <- ll[which(ll[[var]] == max(ll[[var]], na.rm = TRUE)), ]
      # Difference
      diff <- mx$log_rss - mn$log_rss
    } else {
      # Get min of log-RSS
      mn <- ll[which(ll[["log_rss"]] == min(ll[["log_rss"]], na.rm = TRUE)), ]
      # Get max of log-RSS
      mx <- ll[which(ll[["log_rss"]] == max(ll[["log_rss"]], na.rm = TRUE)), ]
      # Combine and sort by 'var'
      cmb <- rbind(mn, mx)
      cmb <- cmb[order(cmb[[var]]), ]
      # Difference
      diff <- cmb$log_rss[2] - cmb$log_rss[1]
    }
    return(diff)
  })
  # Coerce to data.frame
  edf <- as.data.frame(e)
  return(edf)
}

# Calculate a CI for the difference from a data.frame of model effects
diff_ci <- function(x, CI = 0.95) {
  # Get alpha from CI
  alpha <- 1 - CI
  # Get bounds of CI from alpha
  qlwr <- alpha/2
  qupr <- 1 - qlwr
  # Get empirical CI for the difference
  cols <- x
  cols$iter <- NULL
  lwr <- apply(cols, 2, quantile, qlwr)
  mean <- apply(cols, 2, mean)
  upr <- apply(cols, 2, quantile, qupr)
  # Combine
  xx <- rbind(lwr, mean, upr)
  return(xx)
}

# Calculate proportion of iterations with an effect in the direction of the mean
# (Empirical estimate of Pr(effect))
diff_prop <- function(x) {
  # Get rid of iteration column
  cols <- x
  cols$iter <- NULL
  # Calculate proportion
  prop <- apply(cols, 2, function(xx) {
    if (mean(xx) > 0) {
      return(sum(xx > 0)/length(xx))
    }
    
    if (mean(xx) < 0) {
      return(sum(xx < 0)/length(xx))
    }
  })
  # Return
  return(prop)
}