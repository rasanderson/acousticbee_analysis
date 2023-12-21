# Calculate night vs day
# Code from Steve
# rmax=rmax*1000000      # was in MJ/m2 now in j  for whole day        ## needs to be in J/min/cm^2
# rmax=rmax/10000  # in j/cm^2
# light_dark<-function(today, hour, latitude){
#   hours_of_light=matrix(1, nrow=24, ncol=1)
#   hours_light=round(daylength(JDay = today, latitude = latitude)$Sunset)
#   hours_around=round(hours_light/2)
#   before_noon=12-hours_around
#   post_noon=12+hours_around
#   # divide around midday
#   for(t in 1:12){
#     if(t<before_noon){
#       hours_of_light[t,1]=0
#     }
#   }
#   for(t in 12:24){
#     if(t>post_noon){
#       hours_of_light[t,1]=0
#     }
#   }
#   light=  hours_of_light[hour,1]
#   return(light)
# }
# 
# Steve's code stalling so will use alternative
# This is clunky but works
light_dark<-function(day_no, hour, latitude){
  # Returns night or day
  all_res <- NULL
  sunrise_hr <- Sunrise(day_no, latitude)
  sunset_hr  <- Sunset(day_no, latitude)
  for(i in 1:length(sunrise_hr)){
  if(hour[i] <= 12){
    # Before noon
    if(hour[i] > abs(sunrise_hr[i])){
      tmp <- "daylight"
    } else {
      tmp <- "nightime"
    }
  } else{
    # After noon
    if(hour[i]-12 < abs(sunset_hr[i])){
      tmp <- "daylight"
    } else {
      tmp <- "nightime"
    }
  }
    all_res <- c(all_res, tmp)
  }
  return(all_res)
}