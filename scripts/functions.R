#### functions for reading and screening data form the histological images 

read_profile_data_from_csv= function (listpaths=Sys.glob("Summary_starch_*.csv")){
  file <- lapply(listpaths, read.csv)
  
  file <-do.call("rbind", file)%>%
    separate(Slice, c("ID", "radius", "deep_class", "photo", "depth_mm", 
                      "mm", "phototype", "measurement_type"))
  
  file[is.na(file)] <- 0
  file$ID=paste(file$ID, file$radius)
  
  return(file)
}

profile_means=function (data){
  Profi_means=as.data.frame(group_by(data, ID, depth_mm, measurement_type)%>%
                              summarise(mean = mean(X.Area, na.rm = T),
                                        sd = sd(X.Area, na.rm = T),
                                        sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area))))
}

general_profile_average=function (data){
  profile_average=as.data.frame(group_by(data, depth_mm, measurement_type)%>%
                                  summarise(mean = mean(X.Area, na.rm = T),
                                            sd = sd(X.Area, na.rm = T),
                                            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area))))
}