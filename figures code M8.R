##making some figures for thesis
setwd("~/Desktop/DwyerLab_Thesis/for office mac from personal mac_ 9-20-21/")
library(tidyverse)

##heatmap
#hierarchical clustering
polygons_hCID_M8 <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_hCID_M8.csv") %>%
  arrange(ID)
#connection_matrix_M16 <- matrix(nrow = length(polygons_hCID_M16$ID), ncol = length(polygons_hCID_M16$ID))
#for (i in 1:length(polygons_hCID_M16$ID)) {
#for (j in 1:length(polygons_hCID_M16$ID)) {
#if (polygons_hCID_M16$hCID[i] == polygons_hCID_M16$hCID[j]) {
#connection_matrix_M16[i,j] <- 1
#}
#else {
#connection_matrix_M16[i,j] <- 0
#}
#}
#}

test <- expand.grid(X=polygons_hCID_M16$hCID, Y=polygons_hCID_M16$hCID)
test <- test %>%
  mutate(connection = ifelse(X == Y, 1, 0))
test %>%
  ggplot() +
  geom_tile(aes(X, Y, fill = connection))

#ecoregions
centroids_hCID_M16_kCIDs_Ecolevel3 <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_hCID-M16_kCIDs_Ecolevel3.csv") %>%
  arrange(ID)
#connection_matrix_ecoregions <- matrix(nrow = length(centroids_hCID_M16_kCIDs_Ecolevel3$ID), ncol = length(centroids_hCID_M16_kCIDs_Ecolevel3$ID))
#for (i in 1:length(centroids_hCID_M16_kCIDs_Ecolevel3$ID)) {
#for (j in 1:length(centroids_hCID_M16_kCIDs_Ecolevel3$ID)) {
#if (centroids_hCID_M16_kCIDs_Ecolevel3$ID[i] == centroids_hCID_M16_kCIDs_Ecolevel3$ID[j]) {
#connection_matrix_ecoregions[i,j] <- 1
#}
#else {
#connection_matrix_ecoregions[i,j] <- 0
#}
#}
#}









##no longer doing heatmap, doing this instead

centroids_hCID_M16_kCIDs_Ecolevel3 <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_hCID-M16_kCIDs_Ecolevel3.csv")
centroids_hCID_ecoregion_summary <- centroids_hCID_M16_kCIDs_Ecolevel3 %>%
  group_by(CID_M8) %>%
  mutate(cluster_size = n()) %>%
  ungroup() %>%
  group_by(CID_M8, NA_L3CODE, cluster_size) %>%
  summarise(overlap = n()) %>%
  ungroup() %>%
  mutate(overlap_pct = overlap/cluster_size)
hCIDs <- sort(unique(centroids_hCID_M16_kCIDs_Ecolevel3$CID_M8))
ecoregions <- sort(unique(centroids_hCID_M16_kCIDs_Ecolevel3$NA_L3CODE))

all <- tibble(CID_M8=sort(rep((hCIDs),length(ecoregions))),
              NA_L3CODE=rep((ecoregions),length(hCIDs)))
centroids_hCID_ecoregion_summary <- full_join(centroids_hCID_ecoregion_summary, all, by = c("CID_M8", "NA_L3CODE"))
centroids_hCID_ecoregion_summary <- centroids_hCID_ecoregion_summary %>%
  arrange(CID_M8, NA_L3CODE)
centroids_hCID_ecoregion_summary$overlap[is.na(centroids_hCID_ecoregion_summary$overlap)] <- 0

pl<-centroids_hCID_ecoregion_summary %>%
  ggplot() +
  geom_tile(aes(CID_M8, NA_L3CODE, fill = overlap_pct), color = "white", lwd = .5, linetype = 1) +
  coord_fixed() +
  ylab("Level III Ecoregion Code") +
  xlab("Cluster ID") +
  labs(fill = "Overlap %") + 
  ggtitle("Distribution of Clusters over Ecoregions")
ggsave("./Figures for Final Draft/Figure 8.pdf", device = "pdf", plot = pl)

centroids_hCID_ecoregion_summary %>%
  ggplot() +
  geom_tile(aes(hCID, NA_L3CODE, fill = log(overlap)), color = "white", lwd = .5, linetype = 1) +
  coord_fixed() +
  ylab("Level III Ecoregion Code") +
  xlab("Cluster ID") +
  ggtitle("Distribution of Clusters over Ecoregions")




##Figure 3- tracking radial distance through time
hCID_timeseries_post1991 <- hCID_timeseries %>%
  filter(SURVEY_YR>1990)

pl<-hCID_timeseries_post1991 %>%
  filter(hCID == 6) %>%
  ggplot() +
  aes(x = SURVEY_YR, y = rad_dist) +
  xlab("Year") +
  ylab("Radial Distance of Defoliation (m)") +
  geom_point() +
  geom_line()
ggsave("./Figures for Final Draft/Figure 3_new.png", device = "png", plot = pl)



##plotting a histogram of number of ecoregions inhabited by each clusters
overlap_data <- centroids_hCID_ecoregion_summary %>%
  filter(overlap_pct > 0) %>%
  group_by(CID_M8) %>%
  summarise(n_ecoregions = n())

pl <- overlap_data %>%
  ggplot() +
  aes(x = n_ecoregions) +
  geom_histogram(bins = 7, color = "black", fill = "grey") +
  xlab("Number of Ecoregions Occupied") +
  ylab("Frequency") +
  ggtitle("Histogram of Number of Ecoregions Occupied by Clusters")
ggsave("./Figures for Final Draft/Figure 8_new.png", device = "png", plot = pl)


##plotting a histogram of number of outbreaks in each cluster and ecoregion
summary_outbreaks_percluster<-hCID_timeseries %>%
  group_by(hCID) %>%
  summarise(n_outbreaks = length(unique(outbreak_ID_final))-1)
summary_outbreaks_percluster %>% 
  ggplot() +
  aes(x = n_outbreaks) %>%
  geom_histogram(bins = 10, color = "black", fill = "grey") +
  xlab("Number of Outbreaks in a Cluster") +
  ylab("Frequency") +
  ggtitle("Histogram of Number of Outbreaks per Cluster")
ggsave("./Figures for Final Draft/Figure 8.pdf", device = "pdf", plot = pl)



##creating figure 6 for ecoregions
defol_summary_complete_outbreakIDs_all_ecolevels <- bind_rows(defol_summary_complete_outbreakIDs_level1eco, defol_summary_complete_outbreakIDs_level2eco, defol_summary_complete_outbreakIDs_level3eco)

outbreak_length_summary <- defol_summary_complete_outbreakIDs_all_ecolevels %>%
  group_by(eco_level, outbreak_ID_final) %>%
  dplyr::summarise(length_of_outbreak = n()) %>%
  ungroup() %>%
  filter(outbreak_ID_final != 0) %>%
  group_by(eco_level) %>%
  dplyr::mutate(mean_outbreak_length = mean(length_of_outbreak)) %>%
  dplyr::mutate(SD_outbreak_length = sd(length_of_outbreak)) %>%
  dplyr::mutate(CV_outbreak_length = SD_outbreak_length / mean_outbreak_length) %>%
  ungroup() %>%
  group_by(eco_level, length_of_outbreak) %>% ###this is for coloring the graph below by the number of data points at that value (so we can tell when points are overlapping and by how much)
  dplyr::mutate(n_outbreaks_with_length_x = n()) %>%
  ungroup() %>%
  group_by(eco_level) %>%
  dplyr::mutate(total_outbreaks_atlvlM = max(outbreak_ID_final)) %>%
  ungroup() %>%
  mutate(pct_outbreaks_w_lengthX_atlvlM = n_outbreaks_with_length_x / total_outbreaks_atlvlM)

pl<-outbreak_length_summary %>%
  ggplot() +
  geom_point(aes(x = eco_level, y = length_of_outbreak, color = pct_outbreaks_w_lengthX_atlvlM)) +
  scale_x_log10() +
  scale_color_viridis_c() 
ggsave("./useful graphics for presentation/length_of_outbreak_vs_lvlM_hclust.pdf", device = "pdf", plot = pl)


pl<-outbreak_length_summary %>%
  ggplot() +
  geom_point(aes(x = eco_level, y = mean_outbreak_length)) +
  geom_line(aes(x = eco_level, y = mean_outbreak_length)) +
  scale_x_log10()
ggsave("./useful graphics for presentation/mean_length_of_outbreak_vs_lvlM_hclust.pdf", device = "pdf", plot = pl)


outbreak_length_summary %>%
  ggplot() +
  aes(x = eco_level, y = SD_outbreak_length) +
  geom_point() +
  geom_line() +
  scale_x_log10()

outbreak_length_summary %>%
  ggplot() +
  aes(x = eco_level, y = CV_outbreak_length) +
  geom_point() +
  geom_line() +
  scale_x_log10()













####summarizing and visualizing the frequencies of number of outbreaks observed in clusters across groups
outbreak_summary <- defol_summary_complete_outbreakIDs_all_ecolevels %>%
  group_by(eco_level, Ecoregion_Code) %>%
  summarise(n_outbreaks = length(which(yr_of_outbreak==1)), nyrs_defol = length(which(defol_total != 0))) %>%
  ungroup() %>%
  group_by(eco_level) %>%
  mutate(n_clust_revised = length(unique(Ecoregion_Code))) %>%
  ungroup()


###making histogram

outbreak_summary_2<-outbreak_summary %>%
  group_by(eco_level, n_clust_revised, n_outbreaks) %>%
  summarise(nclust_w_noutbrks = n()) %>%
  ungroup() %>%
  mutate(pctclust_w_noutbrks = nclust_w_noutbrks/n_clust_revised)



pl<-outbreak_summary_2 %>%
  ggplot() +
  aes(x=eco_level, y = pctclust_w_noutbrks, fill = n_outbreaks) +
  ggtitle("% of Ecoregions with n outbreaks vs Ecoregion Level") +
  xlab("Ecoregion Level") +
  ylab("% of Ecoregions") +
  geom_bar(stat = "identity") +
  scale_x_log10() +
  scale_fill_viridis_c()
ggsave("./Figures for Final Draft/Figure 6B.png", device = "png", plot = pl)




###Replacing figures 11-12 with different figures
NAL3CODE_timeseries_outbreaksummary <- NAL3CODE_timeseries_outbreaksummary %>%
  mutate(period = length_of_outbreak + time_to_next_outbreak)

NAL3CODE_timeseries_outbreaksummary %>%
  filter(time_to_next_outbreak != 0) %>%
  ggplot() +
  aes(x = length_of_outbreak, y = time_to_next_outbreak) +
  geom_point()


pl <- NAL3CODE_timeseries_outbreaksummary %>%
  filter(time_to_next_outbreak != 0) %>%
  ggplot() +
  aes(x = period) +
  geom_histogram(color = "black", fill = "grey") +
  xlim(0, 30) +
  xlab("Period (years)") +
  ylab("Frequency") +
  ggtitle("Histogram of Outbreak Periods for Level III Ecoregions") +
  geom_text(label = paste("Mean period:", mean(NAL3CODE_timeseries_outbreaksummary$period)), x = 20, y = 4) +
  geom_text(label = paste("Period St. Dev:", sd(NAL3CODE_timeseries_outbreaksummary$period)), x = 20, y = 3.7)
ggsave("./Figures for Final Draft/Figure 11_new.png", device = "png", plot = pl) 




hCID_timeseries_outbreaksummary <- hCID_timeseries %>%
  group_by(outbreak_ID_final) %>%
  summarise(startyear = SURVEY_YR[1],
            hCID = mean(hCID),
            length_of_outbreak = n()) %>%
  filter(outbreak_ID_final != 0) %>%
  mutate(endyear = startyear + (length_of_outbreak - 1), time_to_next_outbreak = 0)

##adding time til next outbreak
for (i in 1:length(hCID_timeseries_outbreaksummary$time_to_next_outbreak)) {
  if (i < length(hCID_timeseries_outbreaksummary$time_to_next_outbreak)) {
    if (hCID_timeseries_outbreaksummary$endyear[i] < hCID_timeseries_outbreaksummary$startyear[i+1]) {
      hCID_timeseries_outbreaksummary$time_to_next_outbreak[i] <- hCID_timeseries_outbreaksummary$startyear[i+1] - hCID_timeseries_outbreaksummary$endyear[i]
    }
  }
}

hCID_timeseries_outbreaksummary <-hCID_timeseries_outbreaksummary %>%
  mutate(period = length_of_outbreak + time_to_next_outbreak)

pl <- hCID_timeseries_outbreaksummary %>%
  filter(time_to_next_outbreak != 0) %>%
  ggplot() +
  aes(x = period) +
  geom_histogram(color = "black", fill = "grey") +
  xlim(0, 30) +
  xlab("Period (years)") +
  ylab("Frequency") +
  ggtitle("Histogram of Outbreak Periods for Clusters") +
  geom_text(label = paste("Mean Period:", mean(hCID_timeseries_outbreaksummary$period)), x = 23, y = 14) +
  geom_text(label = paste("Period St. Dev:", sd(hCID_timeseries_outbreaksummary$period)), x = 23, y = 13)
ggsave("./Figures for Final Draft/Figure 10_new.png", device = "png", plot = pl) 



##ALTERNATIVE TO FIGURE 6
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
ggsave("./Figures for Final Draft/Figure 6_new_V2.png", device = "png", plot = pl)



outbreaks_per_ecoregion_summary <- defol_summary_complete_outbreakIDs_all_ecolevels %>%
  group_by(eco_level) %>%
  summarise(n_outbreaks = (length(unique(outbreak_ID_final))-1)/unique(eco_level)) %>%
  ungroup() 

pl <- outbreaks_per_ecoregion_summary %>%
  ggplot() +
  aes(x=eco_level, y = n_outbreaks) +
  geom_point() +
  geom_line() +
  xlab("Ecoregion Level") +
  ylab("Mean Outbreaks per Ecoregion") +
  ylim(0,20) +
  ggtitle("Number of Outbreaks per Ecoregion vs Ecoregion Level") 
ggsave("./Figures for Final Draft/Figure 6_new_V2_2.png", device = "png", plot = pl)
