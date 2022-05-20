##running the main pipeline for the level 2 ecoregions

library(tidyverse)
library(factoextra)

setwd("~/Desktop/DwyerLab_Thesis/for office mac from personal mac_ 9-20-21/")

centroids_M8_Ecolevel2 <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_NA_L1-2CODE.csv") %>%
  filter(SURVEY_YR > 1989 & NA_L2CODE > -1)

####code from previous pipeline
##calculating "final defol data" for the entire polygons set at once
defol_summary <- centroids_M8_Ecolevel2 %>%
  group_by(NA_L2CODE, SURVEY_YR) %>%
  dplyr::summarise(n_pol = n(),
                   defol_total = sum(area4),
                   rad_dist= sqrt(sum(area4)/pi)
                   #,      the scattering stuff is commeted because we don't have hublines for the hierarchical clusters yet
                   #scat_dist_base = (pi*rad_dist^2)*mean(weighted_dist),
                   #scat_dist_mod = scat_dist_base/sum(inv_area)
  ) %>%
  arrange(NA_L2CODE, SURVEY_YR) %>%
  ungroup()


##adding in years with 0 defoliation to complete the time series

defol_summary_complete <- tibble()
start <- 1991
finish <- 2019


cluster_ids<-unique(defol_summary$NA_L2CODE)
nclusters<-length(cluster_ids)
all <- tibble(NA_L2CODE=sort(rep((cluster_ids),finish-start+1)),
              SURVEY_YR=rep((start:finish),nclusters))

defol_summary_complete<-full_join(defol_summary,all,by=c("SURVEY_YR","NA_L2CODE")) %>%
  dplyr::mutate(n_pol=ifelse(is.na(n_pol),0,n_pol), 
                defol_total=ifelse(is.na(defol_total),0,defol_total), 
                rad_dist=ifelse(is.na(rad_dist),0,rad_dist)
                #, 
                #scat_dist_base=ifelse(is.na(scat_dist_base),0,scat_dist_base), 
                #scat_dist_mod=ifelse(is.na(scat_dist_mod),0,scat_dist_mod)
  )




defol_summary_complete <- defol_summary_complete %>%
  arrange(NA_L2CODE, SURVEY_YR)



###adding delta R and delta S
defol_summary_complete <- defol_summary_complete %>%
  group_by(NA_L2CODE) %>%
  dplyr::mutate(delta_rad=ifelse(is.na(lag(rad_dist, default=NA)), 0, rad_dist-lag(rad_dist, default = 0))
                #, 
                #delta_scat_base=ifelse(is.na(lag(scat_dist_base, default=NA)),0, scat_dist_base-lag(scat_dist_base, default = 0)),
                #delta_scat_mod = ifelse(is.na(lag(scat_dist_mod, default=NA)),0, scat_dist_mod-lag(scat_dist_mod, default = 0))
  ) %>%
  ungroup()



###assigning year of outbreak and outbreak IDs

##looking at autocorrelation within each grouping level to determine cutoff for number of zeros in a row that determines when an outbreak has ended
library(timetk)
defol_summary_complete <- defol_summary_complete %>%
  dplyr::mutate(yr_as_date=as.Date(SURVEY_YR-1990, origin = "1991-01-01")) ##this date doesn't actually mean anything, it's just so that we can feed it into tk_acf_diagnostics



acf_table<-defol_summary_complete %>%
  group_by(NA_L2CODE) %>%
  arrange(yr_as_date) %>%
  tk_acf_diagnostics(yr_as_date, rad_dist) %>%
  dplyr::mutate(signif_onezero = ifelse(ACF > .white_noise_lower & ACF < .white_noise_upper, 0, 1), lvl_M = lvl_M_list[i]) %>%
  ungroup()



##plotting acfs as a function of lag, colored by lvl_M
pl<-acf_table %>%
  ggplot() +
  aes(x = lag, y = ACF) +
  ggtitle("Auto-Correlation vs Lag") +
  geom_point() +
  geom_hline(yintercept = unique(acf_table$.white_noise_upper)) +
  geom_hline(yintercept = unique(acf_table$.white_noise_lower)) +
  scale_color_viridis_c()
ggsave("./useful graphics for presentation/ACFvsLag_Ecolevel3", device = "pdf", plot = pl)



###conclusion: no more than 1 zero in a row is tolerated for defining a continuous outbreak

#assigning outbreak years

defol_summary_complete_outbrkyrs <- tibble()



for (i in 1:length(unique(defol_summary_complete$NA_L2CODE))) {  ##this loop assigns each cell in defol_summary_complete a yr_of_outbreak value
  
  filtered <- defol_summary_complete %>%
    filter(NA_L2CODE == unique(defol_summary_complete$NA_L2CODE)[i]) %>%
    mutate(yr_of_outbreak = 0)
  
  for (j in 1:length(filtered$NA_L2CODE)) {
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
    
    if (j != 2 & j != 1 & j != 29) {
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
    if (j == 29) {
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
  
  defol_summary_complete_outbrkyrs <-  bind_rows(defol_summary_complete_outbrkyrs, filtered) 
  
  print(i) 
}


###assigning outbreak IDs for individual outbreaks
defol_summary_complete_outbreakIDs <- tibble()



filtered <- defol_summary_complete_outbrkyrs %>%
  filter(yr_of_outbreak == 1) %>%
  mutate(outbreak_ID = 0)

for (j in 1:length(filtered$yr_of_outbreak)) {
  filtered$outbreak_ID[j] <- j
}

filtered <- filtered %>%
  select(NA_L2CODE, SURVEY_YR, yr_of_outbreak, outbreak_ID)

defol_summary_complete_outbreakIDs <- defol_summary_complete_outbrkyrs %>%
  left_join(filtered, by = c("NA_L2CODE", "yr_of_outbreak", "SURVEY_YR"))




require(zoo)

defol_summary_complete_outbreakIDs$outbreak_ID[1]<-0
defol_summary_complete_outbreakIDs$outbreak_ID<-na.locf(defol_summary_complete_outbreakIDs$outbreak_ID)

defol_summary_complete_outbreakIDs <- defol_summary_complete_outbreakIDs %>%
  dplyr::mutate(outbreak_ID_final = ifelse(yr_of_outbreak == 0, 0, outbreak_ID)) %>%
  select(-outbreak_ID)


###writing this defol summary as a csv
write_csv(defol_summary_complete_outbreakIDs, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/defol_summary_level2eco_post1988.csv")
defol_summary_complete_outbreakIDs_level2eco <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/defol_summary_level2eco_post1988.csv")%>%
  mutate(eco_level = 2)
colnames(defol_summary_complete_outbreakIDs_level2eco)[1]<-"Ecoregion_Code"
defol_summary_complete_outbreakIDs_level2eco$Ecoregion_Code <- as.character(defol_summary_complete_outbreakIDs_level2eco$Ecoregion_Code)


