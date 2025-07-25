# micropoint_prepare.R

library(tidyverse)
library(lubridate)
library(zoo)
library(microclimf)

# === 1. Daten einlesen ===
filepath <- "data/302_caldernklimawiese_complete/302_KlimaturmCaldernWiese_complete.csv"
rawdata <- read_csv(filepath, na = c("NULL", "", "NA", "NaN", "-9999"), show_col_types = FALSE)

# Stelle sicher, dass datetime korrekt geparst ist
rawdata <- rawdata %>%
  filter(!is.na(datetime)) %>%
  mutate(datetime = parse_date_time(datetime, orders = c("ymd HMS", "ymd HM"), tz = "UTC")) %>%
  filter(!is.na(datetime)) %>%
  arrange(datetime)

# === 2. Lücken durch lineare Interpolation füllen ===

# Erzeuge komplette Zeitreihe im 5-Minuten-Abstand
time_seq <- tibble(datetime = seq(from = min(rawdata$datetime),
                                  to = max(rawdata$datetime),
                                  by = "5 min"))

# Join mit Originaldaten
data_full <- time_seq %>%
  left_join(rawdata, by = "datetime") %>%
  arrange(datetime)

# Zähle wie viele NAs in numerischen Spalten interpoliert werden müssen
numeric_cols <- names(data_full)[sapply(data_full, is.numeric)]
na_before <- sum(is.na(data_full[numeric_cols]))

# Interpolation durchführen
data_interp <- data_full %>%
  mutate(across(all_of(numeric_cols), ~ na.approx(., x = datetime, na.rm = FALSE)))

# Nach Interpolation: Zähle neue NAs (sollten deutlich weniger oder null sein)
na_after <- sum(is.na(data_interp[numeric_cols]))
filled_values <- na_before - na_after

cat(sprintf("Interpolierte Werte: %d\n", filled_values))

# === 3. Tagesblöcke extrahieren von 00:00 bis 23:55 ===
data_interp <- data_interp %>%
  mutate(date = as.Date(datetime),
         hour = hour(datetime),
         minute = minute(datetime)) %>%
  filter(!(is.na(date) | is.na(hour) | is.na(minute)))

# Gruppiere pro Tag und speichere nur Tage mit exakt 288 Zeitpunkten
daily_blocks <- data_interp %>%
  filter(hour * 60 + minute %% 1440 <= 1435) %>%
  group_by(date) %>%
  filter(n() == 288) %>%
  group_split()

cat(sprintf("Gefundene vollständige Tagesblöcke: %d\n", length(daily_blocks)))

# === 4. Ersten Block ausgeben ===
cdata<- daily_blocks[[1]]
cdata <- as.data.frame(cdata)

# Neue Spalten ergänzen

rpdata <- data.frame(
  obs_time  =  as.POSIXct(cdata$datetime),
  temp      = as.double(cdata$Ta_2m),
  relhum    = as.double(cdata$Huma_2m),
  pres      = as.double(101300),
  swdown    = as.double(cdata$rad_sw_in),
  difrad    = as.double(cdata$rad_sw_in * 0.3),
  lwdown    = as.double(cdata$rad_lw_in),
  windspeed = as.double(cdata$Windspeed_2m),
  winddir   = as.double(cdata$Wind_direction_3m),
  precip    = as.double(0)
)
# ===== Vergleichs-Modellierung starten =====



# ===== Parameter definieren =====
vegparams_user <- forestparams
vegparams_user$canht <- as.double(0.6)
vegparams_user$pai   <- as.double(1.2)

soilparams_user <- groundparams
soilparams_user$soildepth <- as.double(0.3)

paii_user <- PAIgeometry(PAI = as.double(1.2), skew = as.double(0.3), spread = as.double(0.7), n = as.double(20))

# ===== Modelllauf: Benutzerdaten =====
cat(sprintf("[%s] Starte Modelllauf mit BENUTZERDATEN\n", Sys.time()), file = logfile, append = TRUE)
mout_user <- microclimf::runpointmodel(
  weather = rpdata ,
  reqhgt = 1,
  vegp = vegparams_user,
  groundp = soilparams_user,
  lat = as.double(50.8),
  long = as.double(8.8)
)
# Lade Standard-Testdaten aus micropoint
data("climdata", package = "micropoint")

# Optional: umbenennen zur Klarheit
climdata_default <- climdata
# ===== Modelllauf: Standarddaten aus Package =====
cat(sprintf("[%s] Starte Modelllauf mit STANDARDDATEN\n", Sys.time()), file = logfile, append = TRUE)
mout_default <- runpointmodel(
  climdata = climdata_default,
  reqhgt = 1,
  vegp = forestparams,
  paii = NA,
  groundp = groundparams,
  lat = 50.8,
  long = 8.8
)

# ===== Vergleichs-Plot =====
plot(climdata_user$obs_time, mout_user$tair, type = "l", col = "darkgreen",
     ylim = range(c(mout_user$tair, mout_default$tair), na.rm = TRUE),
     ylab = "Air Temperature (°C)", xlab = "Time", main = "Micropoint Modellvergleich")
lines(climdata_default$obs_time, mout_default$tair, col = "red", lty = 2)
legend("topright", legend = c("Benutzerdaten", "Standarddaten"),
       col = c("darkgreen", "red"), lty = c(1, 2), bty = "n")

# ===== Abschluss =====
cat(sprintf("[%s] Vergleichs-Modellierung abgeschlossen\n", Sys.time()), file = logfile, append = TRUE)
