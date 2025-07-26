# Load necessary libraries
library(ncdf4)
library(ggplot2)
library(terra)
close(nc)
#system("cdo selname,T,WindDir,WindSpd /media/creu/54D84A9D5B30CD28/proj-data/envimet-hasenkopf-simulation/core_plan_v2/core_plan_v2_output/NetCDF/core_plan_v2_2025-07-21_00.00.00.nc ~/out_plan.nc")

#system("cdo selname,T,WindDir,WindSpd /media/creu/54D84A9D5B30CD28/proj-data/envimet-hasenkopf-simulation/core_orig_v2/NetCDF/core_orig_v2_2025-07-21_00.00.00.nc ~/out_orig.nc")
printout= FALSE

# Path to the NetCDF file
file_path <- "/home/creu/out_orig.nc"

# Open the NetCDF file
nc <- nc_open(file_path)
# Extract temperature, altitude, and time variables
temperature <- ncvar_get(nc, "T")   # Assuming "temperature" is the variable name
#wd <- ncvar_get(nc, "WindDir")   # Assuming "temperature" is the variable name
wind_speed <- ncvar_get(nc, "WindSpd")   # Assuming "temperature" is the variable name
altitude <- ncvar_get(nc, "GridsK")         # Assuming "altitude" is the variable name
time <- ncvar_get(nc, "Time")         # Assuming "altitude" is the variable name
# Check the dimensions of temperature (should return: [lat, lon, altitude, time])
dim(wind_speed)

# Initialize a matrix to store the inversion altitudes and last valid altitude for each grid cell and time step
# Dimensions: [lat, lon, time]
inversion_altitude <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1)) 
#surface_altitude <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))
#altitude_difference <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))
mean_cold_air_temperature <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))
inv_layer_temperature <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))
mean_cold_air_wind_speed <- array(NA, dim = c(nrow(temperature), ncol(temperature), length(time),1))

# plotting function
pp = function  (df){
  # data <- data.frame(
  #   altitude = altitude,  # Altitude (on the y-axis)
  #   temperature = temp_profile  # Temperature for the selected point (on the x-axis)
  # )
  # # Plot the temperature profile
  ggplot(df, aes(y = temperature, x = altitude)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    geom_text(aes(label = round(temperature,5)),  # Add temperature labels (rounded to 1 decimal)
              hjust = -0.5, vjust = -0.5, size = 3, color = "black") +
    
    ggtitle("Temperature Profile at Selected Point") +
    xlab("Altitude (m)") +
    ylab("Temperature (K)") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14)
    ) 
}



# Function to detect temperature inversion and the last NA altitude
detect_inversion <- function(temp_profile, alt_profile) {
  inversion_level <- NA
  surface_level <- NA
  
  # Iterate through the temperature profile to find inversion
  clean_data <- na.omit(data.frame(altitude, temp_profile))
    surface_level = 0
    deepest_temp_index <- which.min(clean_data$temp_profile)
    deepest_temp <- clean_data$temp_profile[deepest_temp_index]
    deepest_altitude <- clean_data$altitude[deepest_temp_index+1]    
    for (i in 2:length(temp_profile)) {    
    if (is.na(temp_profile[i-1]) && !is.na(temp_profile[i])) {
      surface_level <- alt_profile[i]  # Store the last valid altitude before NA
    } 
  }
  return(list(inversion_altitude = deepest_altitude, surface_altitude = surface_level))
}

t = 7
# Iterate over each latitude, longitude, and time step, and detect inversion
# Iterate over each latitude, longitude, and time step, and detect inversion
for (t in 7:7){
for (lat in 1:199) {
  for (lon in 1:145) {

      # Extract the temperature profile and altitude profile for the specific time, latitude, and longitude
      temp_profile <- temperature[lat, lon, , t]
      wind_speed_profile <- wind_speed[lat, lon, , t]
      alt_profile <- altitude  # Assuming altitude is consistent across latitudes for each time
      
      # Detect inversion and the last valid altitude
      result <- detect_inversion(temp_profile, alt_profile)
      
      # Check if both inversion_level and surface_level are not NA
      if (!is.na(result$inversion_altitude) && !is.na(result$surface_altitude) ) {
        # Find the indices of the inversion_level and surface_level in the altitude profile
        inversion_index <- which(alt_profile == result$inversion_altitude)
        last_na_index <- which(alt_profile == result$surface_altitude)
        
        # Check if valid indices exist
        if (length(inversion_index) > 0 && length(last_na_index) > 0) {
          # Ensure correct range of indices
          start_index <- min(inversion_index, last_na_index)
          end_index <- max(inversion_index, last_na_index)
          
          # Subset the temperature profile
          sub_data_temp <- na.omit(data.frame(
            altitude = alt_profile[start_index:end_index],
            temperature = temp_profile[start_index:end_index]
          ))
          
          # Calculate altitude differences and weighted mean temperature
          #if (nrow(sub_data_temp) > 1 & (sub_data_temp$temperature[1] <  sub_data_temp$temperature[2])){
            if (nrow(sub_data_temp) > 1){  
            sub_data_temp$altitude_diff <- c(diff(sub_data_temp$altitude), NA)
            valid_data_temp <- sub_data_temp[!is.na(sub_data_temp$altitude_diff), ]
            cold_air_mean_temp <- sum(valid_data_temp$temperature * valid_data_temp$altitude_diff) /  sum(valid_data_temp$altitude_diff)
          } 
          # Subset wind speed profile (if provided)
          sub_data_wind_speed <- na.omit(data.frame(
            altitude = alt_profile[start_index:end_index],
            WindSpd = wind_speed_profile[start_index:end_index]  # Assuming wind_speed is defined
          ))
          
          # Calculate altitude differences and weighted mean wind speed
          #if (nrow(sub_data_wind_speed) > 1 & (sub_data_temp$temperature[1] <  sub_data_temp$temperature[2])) {
          if (nrow(sub_data_wind_speed) > 1 ) {
            sub_data_wind_speed$altitude_diff <- c(diff(sub_data_wind_speed$altitude), NA)
            valid_data_wind_speed <- sub_data_wind_speed[!is.na(sub_data_wind_speed$altitude_diff), ]
            cold_air_mean_wind_speed <- sum(valid_data_wind_speed$WindSpd * valid_data_wind_speed$altitude_diff) / 
              sum(valid_data_wind_speed$altitude_diff)
          } 
          
          # Extract the temperature of the inversion level
          #if (nrow(sub_data_temp) > 1 && (sub_data_temp$temperature[1] < sub_data_temp$temperature[2]))
          if (nrow(sub_data_temp) > 1)
          inv_temp <- temp_profile[inversion_index]  
          
          # Store the computed values for this grid cell and time step
          mean_cold_air_temperature[lat, lon, t,] <- cold_air_mean_temp
          inv_layer_temperature[lat, lon, t,] <- inv_temp
          mean_cold_air_wind_speed[lat, lon, t,] <- cold_air_mean_wind_speed
          
          # Optional: Print results for debugging
          if (printout) {
            print(paste("Mean temp at lat", lat, "lon", lon, "time", t, ":", cold_air_mean_temp))
            print(paste("Inversion temp at lat", lat, "lon", lon, "time", t, ":", inv_temp))
            print(paste("Inversion altitude at lat", lat, "lon", lon, "time", t, ":", result$inversion_altitude))
            print(paste("Last valid altitude at lat", lat, "lon", lon, "time", t, ":", result$surface_altitude))
            
          }

        }
        
      }
    }

  # Convert the array for the first time step to a raster
  #plot(rast(altitude_difference[, , timestart,]))
  #plot(rast(inversion_altitude[, , timestart,])) 

  
}
  r = min(rast(mean_cold_air_temperature[, , t,]) - rast(inv_layer_temperature[, , t,]))
  plot(t(flip(r, direction = "horizontal",)),  main = paste(" Temp diff cold air/inv ts", "(", t, ")"), 
       col = terrain.colors(100))  

}
crplan = clamp(rplan, lower = -0.5, upper = 0.0, values = FALSE)  
crorig = clamp(rorig, lower = -0.5, upper = 0.0, values = FALSE)
dif=rplan-rorig
plot(clamp(t(flip(dif, direction = "horizontal",)), lower = -0.5, upper = 0.5, values = FALSE))
plot(t(flip(dif, direction = "horizontal",)))
plot(t(flip(crplan, direction = "horizontal",))- t(flip(crorig, direction = "horizontal",)))
# mask=t(flip(rast(mean_cold_air_temperature[, , t,]), direction = "horizontal",))
# r = rast(inv_layer_temperature[, , t,])
# plot(mask(t(flip(r, direction = "horizontal",)),mask),  main = paste(" Inversion Temp ts", "(", t, ")"), 
#      col = terrain.colors(100))
# r= rast(mean_cold_air_temperature[, , t,])
# plot(t(flip(r, direction = "horizontal",)),  main = paste(" Mean Temp cold air ts", "(", t, ")"), 
#      col = terrain.colors(100))
r = min(rast(mean_cold_air_temperature[, , t,]) - rast(inv_layer_temperature[, , t,]))
plot(t(flip(r, direction = "horizontal",)),,  main = paste(" Temp diff cold air/inv ts", "(", t, ")"), 
     col = terrain.colors(100))

# # Convert the array for the first time step to a raster
# #plot(rast(altitude_difference[, , timestart,]))
# #plot(rast(inversion_altitude[, , timestart,])) 
# mask=t(flip(rast(mean_cold_air_temperature[, , timestart,]), direction = "horizontal",))
# r = rast(inv_layer_temperature[, , timestart,])
# plot(mask(t(flip(r, direction = "horizontal",)),mask),  main = paste(" Inversion Temp ts", "(", t, ")"), 
#      col = terrain.colors(100))
# r= rast(mean_cold_air_temperature[, , timestart,])
# plot(t(flip(r, direction = "horizontal",)),  main = paste(" Mean Temp cold air ts", "(", t, ")"), 
#      col = terrain.colors(100))
# r = min(rast(inv_layer_temperature[, , timestart,])-rast(mean_cold_air_temperature[, , timestart,]))
# plot(t(flip(r, direction = "horizontal",)))
# r = min(rast(mean_cold_air_wind_speed[, , timestart,]))
# plot(t(flip(r, direction = "horizontal",)))

#### tipping point

# identify_tipping_point <- function(tempfile) {
#   # Reverse the tempfile to search from the deepest point toward the surface
#   reversed_temp <- rev(temp_profile)
#   
#   # Iterate through the reversed temperature profile
#   for (i in 2:length(reversed_temp)) {
#     # Ensure both values are non-NA
#     if (!is.na(reversed_temp[i]) && !is.na(reversed_temp[i - 1])) {
#       # Check if temperature starts increasing (from bottom to top in reversed_temp)
#       if (reversed_temp[i - 1] < reversed_temp[i]) {
#         # Return the original index of the tipping point (convert from reversed index)
#         return(length(tempfile) - i + 1)
#       }
#     }
#   }
#   
#   return(NA)  # Return NA if no tipping point is found
# }
# 
# tipping_point_index <- identify_tipping_point(temp_profile)
# 
# # Output the tipping point and its temperature value
# if (!is.na(tipping_point_index)) {
#   print(paste("Tipping point index:", tipping_point_index))
#   print(paste("Temperature at tipping point:", tempfile[tipping_point_index]))
# } else {
#   print("No tipping point found.")
# }