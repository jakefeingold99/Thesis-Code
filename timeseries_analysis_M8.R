###turning shepherd data into a time series and merging it with existing time series to extend the time range

shepherd_timeseries <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/shep_timeseries_M8.csv") %>%
  filter(shepID != 0)


####setting up time series based on hierarchical cluster IDs
shepherd_timeseries_hCID <- shepherd_timeseries %>%  ##collapsing duplicate years
  arrange(hCID, SURVEY_YR) %>%
  group_by(hCID, SURVEY_YR) %>%
  summarise() %>%
  ungroup() %>%
  mutate(rad_dist = 1) ##note that this one doesnt mean anything, its just to signify that there is defoliation


##filling in the time series with zero years and empty cluster IDs
cluster_id_list <- c(1:8) 
nclusters <- length(cluster_id_list)
start <- min(shepherd_timeseries_hCID$SURVEY_YR)
finish <- 1990

shep_timeseries_complete <- tibble()

for (i in 1:nclusters) {
  tmp_defol <- shepherd_timeseries_hCID %>%
    filter(hCID == cluster_id_list[i])
  all <- tibble(hCID=cluster_id_list[i],
                SURVEY_YR= start:finish)
  
  tmp_defol<-full_join(tmp_defol,all,by=c("SURVEY_YR","hCID")) %>%
    dplyr::mutate(rad_dist=ifelse(is.na(rad_dist),0,rad_dist)
                  #, 
                  #scat_dist_base=ifelse(is.na(scat_dist_base),0,scat_dist_base), 
                  #scat_dist_mod=ifelse(is.na(scat_dist_mod),0,scat_dist_mod)
    )
  
  shep_timeseries_complete <- bind_rows(shep_timeseries_complete, tmp_defol)
  
}

shep_timeseries_complete <- shep_timeseries_complete %>%
  arrange(hCID, SURVEY_YR)



##assigning outbreak years
shep_timeseries_complete_outbrkyrs <- tibble()




for (i in 1:length(unique(shep_timeseries_complete$hCID))) {  ##this loop assigns each cell in shep_timeseries_complete a yr_of_outbreak value
  
  filtered <- shep_timeseries_complete %>%
    filter(hCID == unique(shep_timeseries_complete$hCID)[i]) %>%
    mutate(yr_of_outbreak = 0)
  
  for (j in 1:length(filtered$hCID)) {
    if (j==1) {
      filtered$yr_of_outbreak[j] <- ifelse(filtered$rad_dist[j] != 0, 1, 0)
    }
    
    if (j==2) {
      if (filtered$yr_of_outbreak[j-1] != 0) {
        if (filtered$rad_dist[j]==0) {
          if (filtered$rad_dist[j+1] ==0) {
            filtered$yr_of_outbreak[j] <- 0
          }
          if (filtered$rad_dist[j+1] != 0) {
            filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
          }
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
        }
      }
      if (filtered$yr_of_outbreak[j-1] == 0) {
        if (filtered$rad_dist[j] == 0) {
          filtered$yr_of_outbreak[j] <- 0
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- 1
        }
      }
    }
    
    if (j != 2 & j != 1 & j != 75) {
      if (filtered$rad_dist[j-2] == 0 & filtered$rad_dist[j-1] == 0) {
        if (filtered$rad_dist[j] == 0) {
          filtered$yr_of_outbreak[j] <- 0
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- 1
        }
      }
      if (filtered$yr_of_outbreak[j-1] != 0) {
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
        }
        if (filtered$rad_dist[j] == 0) {
          if (filtered$rad_dist[j+1] == 0) {
            filtered$yr_of_outbreak[j] <- 0
          }
          if (filtered$rad_dist[j+1] != 0) {
            filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
          }
        }
      }
    }
    if (j == 75) {
      if (filtered$yr_of_outbreak[j-1] != 0) {
        if (filtered$rad_dist[j] == 0) {
          filtered$yr_of_outbreak[j] <- 0
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
        }
      }
      if (filtered$yr_of_outbreak[j-1] == 0) {
        if (filtered$rad_dist[j] == 0) {
          filtered$yr_of_outbreak[j] <- 0
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- 1
        }
      }
    }
  }
  
  shep_timeseries_complete_outbrkyrs <-  bind_rows(shep_timeseries_complete_outbrkyrs, filtered) 
  
  print(i) 
}



#assigning outbreak IDs
shep_timeseries_complete_outbreakIDs <- tibble()
filtered <- shep_timeseries_complete_outbrkyrs %>%
  filter(yr_of_outbreak == 1) %>%
  mutate(outbreak_ID = 0)

for (j in 1:length(filtered$yr_of_outbreak)) {
  filtered$outbreak_ID[j] <- j
}

filtered <- filtered %>%
  select(hCID, SURVEY_YR, yr_of_outbreak, outbreak_ID)

shep_timeseries_complete_outbreakIDs <- shep_timeseries_complete_outbrkyrs %>%
  left_join(filtered, by = c("hCID", "yr_of_outbreak", "SURVEY_YR"))


require(zoo)

shep_timeseries_complete_outbreakIDs$outbreak_ID[1]<-0
shep_timeseries_complete_outbreakIDs$outbreak_ID<-na.locf(shep_timeseries_complete_outbreakIDs$outbreak_ID)

shep_timeseries_complete_outbreakIDs <- shep_timeseries_complete_outbreakIDs %>%
  dplyr::mutate(outbreak_ID_final = ifelse(yr_of_outbreak == 0, 0, outbreak_ID)) %>%
  select(-outbreak_ID)


###writing this defol summary as a csv
write_csv(shep_timeseries_complete_outbreakIDs, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/shep_timeseries_hCID_complete_M8.csv")











###REPLICATING ABOVE PIPELINE FOR NAL3CODE

####setting up time series based on hierarchical cluster IDs
shepherd_timeseries_NAL3CODE <- shepherd_timeseries %>%  ##collapsing duplicate years
  arrange(NA_L3CODE, SURVEY_YR) %>%
  group_by(NA_L3CODE, SURVEY_YR) %>%
  summarise() %>%
  ungroup() %>%
  mutate(rad_dist = 1) ##note that this one doesnt mean anything, its just to signify that there is defoliation


##filling in the time series with zero years and empty cluster IDs
cluster_id_list <- unique(shepherd_timeseries_NAL3CODE$NA_L3CODE) 
nclusters <- length(cluster_id_list)
start <- min(shepherd_timeseries_NAL3CODE$SURVEY_YR)
finish <- 1990

shep_timeseries_complete_NAL3CODE <- tibble()

for (i in 1:nclusters) {
  tmp_defol <- shepherd_timeseries_NAL3CODE %>%
    filter(NA_L3CODE == cluster_id_list[i])
  all <- tibble(NA_L3CODE=cluster_id_list[i],
                SURVEY_YR= start:finish)
  
  tmp_defol<-full_join(tmp_defol,all,by=c("SURVEY_YR","NA_L3CODE")) %>%
    dplyr::mutate(rad_dist=ifelse(is.na(rad_dist),0,rad_dist)
                  #, 
                  #scat_dist_base=ifelse(is.na(scat_dist_base),0,scat_dist_base), 
                  #scat_dist_mod=ifelse(is.na(scat_dist_mod),0,scat_dist_mod)
    )
  
  shep_timeseries_complete_NAL3CODE <- bind_rows(shep_timeseries_complete_NAL3CODE, tmp_defol)
  
}

shep_timeseries_complete_NAL3CODE <- shep_timeseries_complete_NAL3CODE %>%
  arrange(NA_L3CODE, SURVEY_YR)



##assigning outbreak years
shep_timeseries_complete_NAL3CODE_outbrkyrs <- tibble()




for (i in 1:length(unique(shep_timeseries_complete_NAL3CODE$NA_L3CODE))) {  ##this loop assigns each cell in shep_timeseries_complete a yr_of_outbreak value
  
  filtered <- shep_timeseries_complete_NAL3CODE %>%
    filter(NA_L3CODE == unique(shep_timeseries_complete_NAL3CODE$NA_L3CODE)[i]) %>%
    mutate(yr_of_outbreak = 0)
  
  for (j in 1:length(filtered$NA_L3CODE)) {
    if (j==1) {
      filtered$yr_of_outbreak[j] <- ifelse(filtered$rad_dist[j] != 0, 1, 0)
    }
    
    if (j==2) {
      if (filtered$yr_of_outbreak[j-1] != 0) {
        if (filtered$rad_dist[j]==0) {
          if (filtered$rad_dist[j+1] ==0) {
            filtered$yr_of_outbreak[j] <- 0
          }
          if (filtered$rad_dist[j+1] != 0) {
            filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
          }
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
        }
      }
      if (filtered$yr_of_outbreak[j-1] == 0) {
        if (filtered$rad_dist[j] == 0) {
          filtered$yr_of_outbreak[j] <- 0
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- 1
        }
      }
    }
    
    if (j != 2 & j != 1 & j != 75) {
      if (filtered$rad_dist[j-2] == 0 & filtered$rad_dist[j-1] == 0) {
        if (filtered$rad_dist[j] == 0) {
          filtered$yr_of_outbreak[j] <- 0
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- 1
        }
      }
      if (filtered$yr_of_outbreak[j-1] != 0) {
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
        }
        if (filtered$rad_dist[j] == 0) {
          if (filtered$rad_dist[j+1] == 0) {
            filtered$yr_of_outbreak[j] <- 0
          }
          if (filtered$rad_dist[j+1] != 0) {
            filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
          }
        }
      }
    }
    if (j == 75) {
      if (filtered$yr_of_outbreak[j-1] != 0) {
        if (filtered$rad_dist[j] == 0) {
          filtered$yr_of_outbreak[j] <- 0
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- filtered$yr_of_outbreak[j-1] + 1
        }
      }
      if (filtered$yr_of_outbreak[j-1] == 0) {
        if (filtered$rad_dist[j] == 0) {
          filtered$yr_of_outbreak[j] <- 0
        }
        if (filtered$rad_dist[j] != 0) {
          filtered$yr_of_outbreak[j] <- 1
        }
      }
    }
  }
  
  shep_timeseries_complete_NAL3CODE_outbrkyrs <-  bind_rows(shep_timeseries_complete_NAL3CODE_outbrkyrs, filtered) 
  
  print(i) 
}



#assigning outbreak IDs
shep_timeseries_complete_NAL3CODE_outbreakIDs <- tibble()
filtered <- shep_timeseries_complete_NAL3CODE_outbrkyrs %>%
  filter(yr_of_outbreak == 1) %>%
  mutate(outbreak_ID = 0)

for (j in 1:length(filtered$yr_of_outbreak)) {
  filtered$outbreak_ID[j] <- j
}

filtered <- filtered %>%
  select(NA_L3CODE, SURVEY_YR, yr_of_outbreak, outbreak_ID)

shep_timeseries_complete_NAL3CODE_outbreakIDs <- shep_timeseries_complete_NAL3CODE_outbrkyrs %>%
  left_join(filtered, by = c("NA_L3CODE", "yr_of_outbreak", "SURVEY_YR"))


require(zoo)

shep_timeseries_complete_NAL3CODE_outbreakIDs$outbreak_ID[1]<-ifelse(shep_timeseries_complete_NAL3CODE_outbreakIDs$yr_of_outbreak[1]==1, 1, 0)
shep_timeseries_complete_NAL3CODE_outbreakIDs$outbreak_ID<-na.locf(shep_timeseries_complete_NAL3CODE_outbreakIDs$outbreak_ID)

shep_timeseries_complete_NAL3CODE_outbreakIDs <- shep_timeseries_complete_NAL3CODE_outbreakIDs %>%
  dplyr::mutate(outbreak_ID_final = ifelse(yr_of_outbreak == 0, 0, outbreak_ID)) %>%
  select(-outbreak_ID)


###writing this defol summary as a csv
write_csv(shep_timeseries_complete_NAL3CODE_outbreakIDs, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/shep_timeseries_NAL3CODE_complete.csv")









##NOW THAT WE HAVE THE TIME SERIES FOR THE HIERARCHICAL CLUSTERS AND THE ECOLEVEL 3 REGIONS, WE HAVE TO JOIN THEM TO THE RECENT TIME SERIES TO COMPLETE THE FULL TIME SERIES

post1990_timeseries_hCID <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/defol_summary_hCID_M8_post1988.csv")
shep_timeseries_hCID <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/shep_timeseries_hCID_complete_M8.csv")
post1990_timeseries_NAL3CODE <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/defol_summary_ecolevel_complete.csv")
shep_timeseries_NAL3CODE <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/shep_timeseries_NAL3CODE_complete.csv")


hCID_timeseries_complete <- bind_rows(shep_timeseries_hCID, post1990_timeseries_hCID) %>%
  arrange(hCID, SURVEY_YR) %>%
  select(-outbreak_ID_final)

NAL3CODE_timeseries_complete <- bind_rows(shep_timeseries_NAL3CODE, post1990_timeseries_NAL3CODE) %>%
  arrange(NA_L3CODE, SURVEY_YR)%>%
  select(-outbreak_ID_final)


#reassigning outbreak IDs
hCID_timeseries_complete_outbreakIDs <- tibble()   ###for hCID
filtered <- hCID_timeseries_complete %>%
  filter(yr_of_outbreak == 1) %>%
  mutate(outbreak_ID = 0)

for (j in 1:length(filtered$yr_of_outbreak)) {
  filtered$outbreak_ID[j] <- j
}

filtered <- filtered %>%
  select(hCID, SURVEY_YR, yr_of_outbreak, outbreak_ID)

hCID_timeseries_complete_outbreakIDs <- hCID_timeseries_complete %>%
  left_join(filtered, by = c("hCID", "yr_of_outbreak", "SURVEY_YR"))


require(zoo)

hCID_timeseries_complete_outbreakIDs$outbreak_ID[1]<-ifelse(hCID_timeseries_complete_outbreakIDs$yr_of_outbreak[1]==1, 1, 0)
hCID_timeseries_complete_outbreakIDs$outbreak_ID<-na.locf(hCID_timeseries_complete_outbreakIDs$outbreak_ID)

hCID_timeseries_complete_outbreakIDs <- hCID_timeseries_complete_outbreakIDs %>%
  dplyr::mutate(outbreak_ID_final = ifelse(yr_of_outbreak == 0, 0, outbreak_ID)) %>%
  select(-outbreak_ID)

##writing the final time series as a csv
write_csv(hCID_timeseries_complete_outbreakIDs, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/hCID_timeseries_complete_outbreakIDs_M8.csv")




NAL3CODE_timeseries_complete_outbreakIDs <- tibble()   ###for ecoregion level 3
filtered <- NAL3CODE_timeseries_complete %>%
  filter(yr_of_outbreak == 1) %>%
  mutate(outbreak_ID = 0)

for (j in 1:length(filtered$yr_of_outbreak)) {
  filtered$outbreak_ID[j] <- j
}

filtered <- filtered %>%
  select(NA_L3CODE, SURVEY_YR, yr_of_outbreak, outbreak_ID)

NAL3CODE_timeseries_complete_outbreakIDs <- NAL3CODE_timeseries_complete %>%
  left_join(filtered, by = c("NA_L3CODE", "yr_of_outbreak", "SURVEY_YR"))


require(zoo)

NAL3CODE_timeseries_complete_outbreakIDs$outbreak_ID[1]<-ifelse(NAL3CODE_timeseries_complete_outbreakIDs$yr_of_outbreak[1]==1, 1, 0)
NAL3CODE_timeseries_complete_outbreakIDs$outbreak_ID<-na.locf(NAL3CODE_timeseries_complete_outbreakIDs$outbreak_ID)

NAL3CODE_timeseries_complete_outbreakIDs <- NAL3CODE_timeseries_complete_outbreakIDs %>%
  dplyr::mutate(outbreak_ID_final = ifelse(yr_of_outbreak == 0, 0, outbreak_ID)) %>%
  select(-outbreak_ID)

##writing the final time series as a csv
write_csv(NAL3CODE_timeseries_complete_outbreakIDs, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/NAL3CODE_timeseries_complete_outbreakIDs.csv")

