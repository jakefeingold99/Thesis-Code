##attemping hierarchical clustering on the polygon data set

library(tidyverse)
library(factoextra)

setwd("~/Desktop/DwyerLab_Thesis/for office mac from personal mac_ 9-20-21/")


###the code below calculates the clusters

##attempting hierarchical clustering on the full set of centroids
centroids_XY <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/102039_centroids_allCIDs_dissolved_byyear_area4_XY.csv") %>%
  filter(SURVEY_YR > 1989)

centroids_XY_formatted_forclust<-centroids_XY %>%
  select(xcoord, ycoord)

lvl_M_list <- c(3, 4, 5, 7, 8, 10, 13, 16, 20, 25, 31, 38, 47, 59, 73, 90, 112, 139, 172, 213, 263, 326, 404, 499)


##creating a long data frame of all the polygons
polygons_hCIDs_long <- tibble()
for (i in 1:length(lvl_M_list)) {  ##this for loop clusters the polygons, adds their cluster ID and clustering level as a field, then joins that frame to an empty data frame
  h_clusters <- hcut(centroids_XY_formatted_forclust, k = lvl_M_list[i])
  tmp_polygons_hCIDs <- centroids_XY %>%
    mutate(hCID = 0, lvl_M = lvl_M_list[i])
  tmp_polygons_hCIDs$hCID <- h_clusters$cluster
  polygons_hCIDs_long <- bind_rows(polygons_hCIDs_long, tmp_polygons_hCIDs)
  
  
  pl<-fviz_cluster(h_clusters)
  ggsave(filename = paste0("./useful graphics for presentation/hierarchical clustering visualizations/hclust_M",lvl_M_list[i]), device = "pdf", plot = pl)
  
  
  
  print(i)
}

polygons_hCIDs_long <- polygons_hCIDs_long %>% #selecting the variables we want
  mutate(inv_area = 1/area4) %>%
  select(SURVEY_YR, ID, area4, xcoord, ycoord, hCID, lvl_M, inv_area)


#writing the M16 and M8 clusters and all the clusters as csvs to read in to QGIS
polygons_hCIDs_M16 <- polygons_hCIDs_long %>%
  filter(lvl_M == 16)

polygons_hCIDs_M8 <- polygons_hCIDs_long %>%
  filter(lvl_M == 8)

write_csv(polygons_hCIDs_M16, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_hCID_M16.csv")
write_csv(polygons_hCIDs_M8, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_hCID_M8.csv")
write_csv(polygons_hCIDs_long, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_hCIDs_allclusters.csv")





##THIS CODE CAN BE RUN ONCE YOU HAVE THE CLUSTERS


####code from previous pipeline
##calculating "final defol data" for the entire polygons set at once
polygons_hCIDs_long <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_hCIDs_allclusters.csv")




defol_summary <- polygons_hCIDs_long %>%
  group_by(lvl_M, hCID, SURVEY_YR) %>%
  dplyr::summarise(n_pol = n(),
                   defol_total = sum(area4),
                   rad_dist= sqrt(sum(area4)/pi)
                   #,      the scattering stuff is commeted because we don't have hublines for the hierarchical clusters yet
                   #scat_dist_base = (pi*rad_dist^2)*mean(weighted_dist),
                   #scat_dist_mod = scat_dist_base/sum(inv_area)
                   ) %>%
  arrange(lvl_M,  hCID, SURVEY_YR) %>%
  ungroup()


##adding in years with 0 defoliation to complete the time series

defol_summary_complete <- tibble()
start <- 1991
finish <- 2019

for (i in 1:length(lvl_M_list)) {
  tmp_defol <- defol_summary %>%
    filter(lvl_M == lvl_M_list[i])
  cluster_ids<-unique(tmp_defol$hCID)
  nclusters<-length(cluster_ids)
  all <- tibble(lvl_M=lvl_M_list[i], hCID=sort(rep((cluster_ids),finish-start+1)),
                SURVEY_YR=rep((start:finish),nclusters))
  
  tmp_defol<-full_join(tmp_defol,all,by=c("lvl_M","SURVEY_YR","hCID")) %>%
    dplyr::mutate(n_pol=ifelse(is.na(n_pol),0,n_pol), 
                  defol_total=ifelse(is.na(defol_total),0,defol_total), 
                  rad_dist=ifelse(is.na(rad_dist),0,rad_dist)
                  #, 
                  #scat_dist_base=ifelse(is.na(scat_dist_base),0,scat_dist_base), 
                  #scat_dist_mod=ifelse(is.na(scat_dist_mod),0,scat_dist_mod)
                  )
  
  defol_summary_complete <- bind_rows(defol_summary_complete, tmp_defol)
  
}

defol_summary_complete <- defol_summary_complete %>%
  arrange(lvl_M,  hCID, SURVEY_YR)



###adding delta R and delta S
defol_summary_complete <- defol_summary_complete %>%
  group_by(lvl_M, hCID) %>%
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


acf_table <- tibble()

for (i in 1:length(lvl_M_list)) {
  acf_tmp<-defol_summary_complete %>%
    filter(lvl_M == lvl_M_list[i]) %>%
    group_by(hCID) %>%
    arrange(yr_as_date) %>%
    tk_acf_diagnostics(yr_as_date, rad_dist) %>%
    dplyr::mutate(signif_onezero = ifelse(ACF > .white_noise_lower & ACF < .white_noise_upper, 0, 1), lvl_M = lvl_M_list[i]) %>%
    ungroup()
  
  acf_table <- bind_rows(acf_table, acf_tmp)
}

acf_summary <- acf_table %>%
  group_by(lvl_M, lag) %>%
  dplyr::summarise(pct_signif = mean(signif_onezero))

##plotting acfs as a function of lag, colored by lvl_M
pl<-acf_table %>%
  ggplot() +
  aes(x = lag, y = ACF, color = lvl_M) +
  ggtitle("Auto-Correlation vs Lag") +
  geom_point() +
  geom_hline(yintercept = unique(acf_table$.white_noise_upper)) +
  geom_hline(yintercept = unique(acf_table$.white_noise_lower)) +
  scale_color_viridis_c()
ggsave("./useful graphics for presentation/ACFvsLag_allMs_hclust", device = "pdf", plot = pl)

library(reshape2)

pl<-acf_summary %>%
  ggplot() +
  aes(x = lag, y = pct_signif, color = lvl_M, group = lvl_M) +
  ggtitle("% of significant ACF values vs Lag") +
  ylab("% of significant ACF values") +
  xlab("lag (years)") +
  labs(color = "Clustering Level") +
  geom_line() +
  geom_point()+
  scale_color_viridis_c() 
ggsave("./Figures for Final Draft/Figure 4.pdf", device = "pdf", plot = pl)

pl<-acf_summary %>%
  filter(lag != 0) %>%
  ggplot() +
  aes(x = lvl_M, y = pct_signif, color = lag, group = lag) +
  ggtitle("% of ACF values that are significant vs clustering level") +
  ylab("% of significant ACF values") + 
  xlab("Number of Clusters") +
  labs(color = "lag (years)") +
  geom_line() +
  geom_point() +
  scale_color_viridis_c() 
ggsave("./Figures for Final Draft/Figure 5.pdf", device = "pdf", plot = pl)

acf_summary_2 <- acf_summary %>%
  arrange(lag)
###conclusion: no more than 1 zero in a row is tolerated for defining a continuous outbreak


defol_summary_complete_outbrkyrs <- tibble()

for (i in 1:length(lvl_M_list)) { ##this loop assigns each cell in defol_summary_complete a yr_of_outbreak value
  
  tmp_defol <- defol_summary_complete %>%
    filter(lvl_M == lvl_M_list[i])
  
  for (i in 1:length(unique(tmp_defol$hCID))) {  
    
    filtered <- tmp_defol %>%
      filter(hCID == unique(tmp_defol$hCID)[i]) %>%
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
    
    
  }
  print(i)
}


###assigning outbreak IDs for individual outbreaks
defol_summary_complete_outbreakIDs <- tibble()

for (i in 1:length(lvl_M_list)) {
  tmp_defol <- defol_summary_complete_outbrkyrs %>%
    filter(lvl_M == lvl_M_list[i])
  
  filtered <- tmp_defol %>%
    filter(yr_of_outbreak == 1) %>%
    mutate(outbreak_ID = 0)
  
  for (j in 1:length(filtered$yr_of_outbreak)) {
    filtered$outbreak_ID[j] <- j
  }
  
  filtered <- filtered %>%
    select(lvl_M, hCID, SURVEY_YR, yr_of_outbreak, outbreak_ID)
  
  tmp_defol <- tmp_defol %>%
    left_join(filtered, by = c("lvl_M", "hCID", "yr_of_outbreak", "SURVEY_YR"))
  
  defol_summary_complete_outbreakIDs <- defol_summary_complete_outbreakIDs %>%
    bind_rows(tmp_defol)
  
}

require(zoo)

defol_summary_complete_outbreakIDs$outbreak_ID[1]<-0
defol_summary_complete_outbreakIDs$outbreak_ID<-na.locf(defol_summary_complete_outbreakIDs$outbreak_ID)

defol_summary_complete_outbreakIDs <- defol_summary_complete_outbreakIDs %>%
  dplyr::mutate(outbreak_ID_final = ifelse(yr_of_outbreak == 0, 0, outbreak_ID)) %>%
  select(-outbreak_ID)


###writing this defol summary as a csv
write_csv(defol_summary_complete_outbreakIDs, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/defol_summary_hCIDs_allMs.csv")

defol_summary_complete_outbreakIDs <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/defol_summary_hCIDs_allMs.csv")

#writing just the 16 and the 8 clusters level as a csv for other use
defol_summary_complete_outbreakIDs_M16 <- defol_summary_complete_outbreakIDs %>%
  filter(lvl_M == 16)
write_csv(defol_summary_complete_outbreakIDs_M16, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/defol_summary_hCID_M16_post1988.csv")

defol_summary_complete_outbreakIDs_M8 <- defol_summary_complete_outbreakIDs %>%
  filter(lvl_M == 8)
write_csv(defol_summary_complete_outbreakIDs_M8, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/defol_summary_hCID_M8_post1988.csv")







#### calculating the lengths of individual outbreaks, then calculating the average length of an outbreak within each clustering level, then plotting as a function of clustering level
outbreak_length_summary <- defol_summary_complete_outbreakIDs %>%
  group_by(lvl_M, outbreak_ID_final) %>%
  dplyr::summarise(length_of_outbreak = n()) %>%
  ungroup() %>%
  filter(outbreak_ID_final != 0) %>%
  group_by(lvl_M) %>%
  dplyr::mutate(mean_outbreak_length = mean(length_of_outbreak)) %>%
  dplyr::mutate(SD_outbreak_length = sd(length_of_outbreak)) %>%
  dplyr::mutate(CV_outbreak_length = SD_outbreak_length / mean_outbreak_length) %>%
  ungroup() %>%
  group_by(lvl_M, length_of_outbreak) %>% ###this is for coloring the graph below by the number of data points at that value (so we can tell when points are overlapping and by how much)
  dplyr::mutate(n_outbreaks_with_length_x = n()) %>%
  ungroup() %>%
  group_by(lvl_M) %>%
  dplyr::mutate(total_outbreaks_atlvlM = max(outbreak_ID_final)) %>%
  ungroup() %>%
  mutate(pct_outbreaks_w_lengthX_atlvlM = n_outbreaks_with_length_x / total_outbreaks_atlvlM)

pl<-outbreak_length_summary %>%
  ggplot() +
  geom_point(aes(x = lvl_M, y = length_of_outbreak, color = pct_outbreaks_w_lengthX_atlvlM)) +
  scale_x_log10() +
  scale_color_viridis_c() 
ggsave("./useful graphics for presentation/length_of_outbreak_vs_lvlM_hclust.pdf", device = "pdf", plot = pl)


pl<-outbreak_length_summary %>%
  ggplot() +
  geom_point(aes(x = lvl_M, y = mean_outbreak_length)) +
  geom_line(aes(x = lvl_M, y = mean_outbreak_length)) +
  scale_x_log10()
ggsave("./useful graphics for presentation/mean_length_of_outbreak_vs_lvlM_hclust.pdf", device = "pdf", plot = pl)


outbreak_length_summary %>%
  ggplot() +
  aes(x = lvl_M, y = SD_outbreak_length) +
  geom_point() +
  geom_line() +
  scale_x_log10()

outbreak_length_summary %>%
  ggplot() +
  aes(x = lvl_M, y = CV_outbreak_length) +
  geom_point() +
  geom_line() +
  scale_x_log10()













####summarizing and visualizing the frequencies of number of outbreaks observed in clusters across groups
outbreak_summary <- defol_summary_complete_outbrkyrs %>%
  group_by(lvl_M, hCID) %>%
  summarise(n_outbreaks = length(which(yr_of_outbreak==1)), nyrs_defol = length(which(defol_total != 0))) %>%
  ungroup() %>%
  group_by(lvl_M) %>%
  mutate(n_clust_revised = length(unique(hCID))) %>%
  ungroup()


###making histogram

outbreak_summary_2<-outbreak_summary %>%
  group_by(lvl_M, n_clust_revised, n_outbreaks) %>%
  summarise(nclust_w_noutbrks = n()) %>%
  ungroup() %>%
  mutate(pctclust_w_noutbrks = nclust_w_noutbrks/n_clust_revised)



pl<-outbreak_summary_2 %>%
  ggplot() +
  aes(x=lvl_M, y = pctclust_w_noutbrks, fill = n_outbreaks) +
  ggtitle("% of clusters with n outbreaks vs clustering level") +
  xlab("Clustering Level") +
  ylab("% of clusters") +
  geom_bar(stat = "identity") +
  scale_x_log10() +
  scale_fill_viridis_c()
ggsave("./useful graphics for presentation/pctclust_w_noutbrks_vs_lvlM_hclust", device = "pdf", plot = pl)



###making histogram
outbreak_summary_3<-outbreak_summary %>%
  group_by(lvl_M, n_clust_revised, nyrs_defol) %>%
  summarise(nclust_w_nyrs_defol = n()) %>%
  ungroup() %>%
  mutate(pctclust_w_nyrs_defol = nclust_w_nyrs_defol/n_clust_revised)



pl<-outbreak_summary_3 %>%
  ggplot() +
  aes(x=lvl_M, y = pctclust_w_nyrs_defol, fill = nyrs_defol) +
  ggtitle("% of clusters with n years of defoliation vs clustering level") +
  xlab("Clustering Level") +
  ylab("% of Clusters") +
  geom_bar(stat = "identity", color = "white") +
  scale_x_log10() +
  scale_fill_viridis_c()
ggsave("./useful graphics for presentation/pctclust_w_nyrs_defol_vs_lvlM_hclust", device = "pdf", plot = pl)




outbreaks_per_cluster_summary <- defol_summary_complete_outbreakIDs %>%
  group_by(lvl_M) %>%
  summarise(n_outbreaks = (length(unique(outbreak_ID_final))-1)/unique(lvl_M)) %>%
  ungroup() 

pl <- outbreaks_per_cluster_summary %>%
  ggplot() +
  aes(x=lvl_M, y = n_outbreaks) +
  geom_point() +
  geom_line() +
  xlab("Number of Clusters") +
  ylab("Mean Outbreaks per Cluster") +
  ggtitle("Number of Outbreaks per Cluster vs Number of Clusters") +
  scale_x_log10()
