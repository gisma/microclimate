## -----------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(units)

#' Calculate sensible and latent heat fluxes from micrometeorological data
#'
#' This function computes the sensible heat flux (H) using the bulk aerodynamic method,
#' and estimates latent heat flux (LE) as the residual of the surface energy balance.
#'
#' Sensible heat flux is calculated as:
#'     H = ρ * cp * (T_surface - T_air) / ra
#' where:
#'   - ρ is air density [kg/m³]
#'   - cp is specific heat capacity of air [J/kg/K]
#'   - ra is aerodynamic resistance [s/m]
#'   - T_surface is approximated by air temperature at 2 m (Ta_2m)
#'   - T_air is approximated by air temperature at 10 m (Ta_10m)
#'
#' This definition follows standard micrometeorological literature:
#'  - Monteith & Unsworth (2013). *Principles of Environmental Physics*, 4th ed.
#'  - Foken (2008). *Micrometeorology*, Springer.
#'  - Campbell & Norman (1998). *An Introduction to Environmental Biophysics*.
#'
#' Note: T_surface is often approximated by near-surface air temperature (2 m)
#' when radiometric surface temperature or soil skin temperature is not available.
#' This ensures that positive H means upward flux from surface to atmosphere (daytime).
#'
#' LE is derived as the residual from:
#'     LE = Rn - G - H
#'
calc_fluxes <- function(data,
                        z1 = 2, z2 = 10,
                        k = 0.41,
                        rho_air = set_units(1.225, "kg/m^3"),
                        cp_air  = set_units(1005, "J/kg/K"),
                        lambda  = 2.45e6,             # latent heat of vaporization [J/kg]
                        ra_max  = 500,               # max aerodynamic resistance [s/m]
                        WS_min  = 0.2,               # min wind speed [m/s]
                        H_max = 300,                 # max absolute sensible heat flux [W/m²]
                        dT_thresh = 0.1,             # minimum relevant temperature gradient [K]
                        filter_unrealistic_H = TRUE  # remove H if unrealistically large
) {
  
  # Temperature gradient (approx. ΔT = T_surface - T_air)
  data$delta_T <- data$Ta_2m - data$Ta_10m  # positive → surface warmer than air → upward flux
  
  # Mean wind speed at reference height
  data$WS_mean <- rowMeans(data[c("WS_2m", "WS_10m")], na.rm = TRUE)
  
  # Aerodynamic resistance ra
  data$ra <- with(data, ifelse(
    WS_mean < WS_min | is.na(WS_mean), NA,
    log(z2 / z1) / (k * WS_mean)
  ))
  data$ra[data$ra > ra_max] <- NA  # filter unphysical values
  
  # Sensible heat flux H [W/m²]
  data$H <- with(data, ifelse(
    is.na(ra), NA,
    drop_units(rho_air * cp_air * delta_T / ra)
  ))
  
  # Optionally filter out unrealistically high H under weak gradients
  if (filter_unrealistic_H) {
    filter_idx <- which(abs(data$H) > H_max & abs(data$delta_T) < dT_thresh)
    if (length(filter_idx) > 0) {
      data$H[filter_idx] <- NA
    }
  }
  
  # Latent heat flux LE (residual closure)
  data$LE <- with(data, ifelse(is.na(H), NA, Rn - G - H))
  
  # Daily evapotranspiration [mm/day]
  data$ET_mm_day <- with(data, ifelse(is.na(LE), NA, (LE / lambda) * 86400))
  
  return(data)
}

## -----------------------------------------------------------------------------
# Load energy balance data
data_file <- "data/energie_bil_wiese.csv"
energy <- read_csv(data_file) %>%
  mutate(datetime = dmy_hm(datetime))

## -----------------------------------------------------------------------------
# Rename columns and derive additional variables
## -------------------------------------------------------------------------------------------------------------------------------------------------
# Rename columns, correct G sign (soil heat flux), and add derived variables
energy <- read_csv("data/energie_bil_wiese.csv") %>%
  mutate(datetime = dmy_hm(datetime)) %>%
  rename(
    Rn     = rad_bil,
    G_raw  = heatflux_soil,  # raw original, to negate
    Ta_2m  = Ta_2m,
    Ta_10m = Ta_10m,
    WS_2m  = Windspeed_2m,
    WS_10m = Windspeed_10m
  ) %>%
  mutate(
    G           = -G_raw,  # ⚠️ invert sign: positive = flux into soil (day), negative = out of soil (night)
    month_num   = month(datetime),
    month_label = case_when(
      month_num == 6  ~ "June",
      month_num == 11 ~ "November",
      TRUE            ~ as.character(month(datetime, label = TRUE))
    )
  )

## -----------------------------------------------------------------------------
# Compute energy fluxes using custom function
energy <- calc_fluxes(energy)

## -----------------------------------------------------------------------------
# Visualization: diagnostic plots for key variables
plot_diagnostics_by_month <- function(df, month_name, output_pdf = NULL) {
  df_month <- df[df$month_label == month_name, ]
  
  if (!is.null(output_pdf)) {
    pdf(output_pdf, width = 10, height = 6)
  }
  
  plot(df_month$datetime, df_month$delta_T, type = "l", col = "darkblue", lwd = 1.5,
       xlab = "Time", ylab = expression(Delta * "T (K)"),
       main = paste("Temperature Gradient ΔT –", month_name))
  abline(h = 0, col = "grey80", lty = 2)
  
  plot(df_month$datetime, df_month$WS_mean, type = "l", col = "darkgreen", lwd = 1.5,
       xlab = "Time", ylab = "Wind Speed (m/s)",
       main = paste("Mean Wind Speed –", month_name))
  abline(h = 0, col = "grey80", lty = 2)
  
  plot(df_month$datetime, df_month$ra, type = "l", col = "orange", lwd = 1.5,
       xlab = "Time", ylab = "ra (s/m)",
       main = paste("Aerodynamic Resistance –", month_name))
  abline(h = 0, col = "grey80", lty = 2)
  
  plot(df_month$datetime, df_month$H, type = "l", col = "red", lwd = 1.5,
       xlab = "Time", ylab = "H (W/m²)",
       main = paste("Sensible Heat Flux –", month_name))
  abline(h = 0, col = "grey80", lty = 2)
  
  plot(df_month$datetime, df_month$LE, type = "l", col = "blue", lwd = 1.5,
       xlab = "Time", ylab = "LE (W/m²)",
       main = paste("Latent Heat Flux –", month_name))
  abline(h = 0, col = "grey80", lty = 2)
  
  if (!is.null(output_pdf)) {
    dev.off()
  }
}

pdf("data/diagnostic_plots.pdf", width = 10, height = 6)
plot_diagnostics_by_month(energy, "June")
plot_diagnostics_by_month(energy, "November")
dev.off()

## -----------------------------------------------------------------------------
# Convert to long format for ggplot
energy_long <- energy %>%
  select(datetime, month_label, Rn, G, H, LE) %>%
  pivot_longer(cols = c(Rn, G, H, LE), names_to = "flux", values_to = "value")

plot_fluxes <- function(df, title) {
  ggplot(df, aes(x = datetime, y = value, color = flux)) +
    geom_line() +
    labs(title = title, x = "Time", y = "Flux (W/m²)") +
    theme_minimal()
}

print(plot_fluxes(filter(energy_long, month_label == "June"), "Energy Fluxes in June"))
print(plot_fluxes(filter(energy_long, month_label == "November"), "Energy Fluxes in November"))

## -----------------------------------------------------------------------------
# Monthly average ET (mm/day)
monthly_et <- energy %>%
  filter(!is.na(ET_mm_day)) %>%
  group_by(month_label) %>%
  summarise(mean_ET = mean(ET_mm_day, na.rm = TRUE), .groups = "drop")

ggplot(monthly_et, aes(x = month_label, y = mean_ET, fill = month_label)) +
  geom_col(show.legend = FALSE) +
  labs(
    title = "Mean Daily Evapotranspiration (ET) from Residual Method",
    x = "Month",
    y = "Evapotranspiration (mm/day)"
  ) +
  scale_fill_manual(values = c("June" = "tomato", "November" = "turquoise3")) +
  theme_minimal()

## -----------------------------------------------------------------------------
# Export processed data
write_csv(energy, "data/processed_energy_fluxes.csv")

# Ergebnis speichern
write_csv(energy, "data/processed_energy_fluxes.csv")
print(plot_fluxes(filter(energy_long, month_label == "June"), "Energy Fluxes in June"))
print(plot_fluxes(filter(energy_long, month_label == "November"), "Energy Fluxes in November"))
plot_diagnostics_by_month(energy, "June")
plot_diagnostics_by_month(energy, "November")

