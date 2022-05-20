###using centroids with hierarchical clustering IDs to compute averages of different variables across clusters
###this is for semi-replicating the Shepherd 88 analysis of outbreak dynamics across groups.

setwd("~/Desktop/DwyerLab_Thesis/for office mac from personal mac_ 9-20-21/")
library(tidyverse)

centroids_hCID_M16_kCIDs_Ecolevel3 <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/centroids_hCID-M16_kCIDs_Ecolevel3.csv")

##reading in the data frame that has just the relevant host1 data
centroids_host1 <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/102039_alldefol_dissolved_hCID_kCID_ecolevel3_joined-to-undissolved.csv")
centroids_host1 <- centroids_host1 %>%    ###calculating the number of elements in each cluster for each type of clustering
  group_by(CID_M8) %>%
  mutate(n_hCID = n()) %>%
  ungroup() %>%
  group_by(NA_L3CODE) %>%
  mutate(n_NAL3Code = n()) %>%
  ungroup() %>%
  group_by(CID_M13) %>%
  mutate(n_kCID_M13 = n()) %>%
  ungroup()

centroids_host1$HOST1_2[is.na(centroids_host1$HOST1_2)] <- 0 #setting NAs in host1 to 0

##summarizing host1 data by these groups
host1list<-unique(centroids_host1$HOST1_2)

hCID_M16_host1_summary <- centroids_host1 %>% ##summarizing for hierarchical clusters
  group_by(CID_M8, n_hCID) %>%
  dplyr::summarise(n_host_neg1 = sum(HOST1_2 == -1),
                   n_host_202 = sum(HOST1_2 == 202),
                   n_host_15 = sum(HOST1_2 == 15),
                   n_host_999 = sum(HOST1_2 == 999),
                   n_host_299 = sum(HOST1_2 == 299),
                   n_host_9107 = sum(HOST1_2 == 9107),
                   n_host_10 = sum(HOST1_2 == 10),
                   n_host_19 = sum(HOST1_2 == 19),
                   n_host_9102 = sum(HOST1_2 == 9102),
                   n_host_122 = sum(HOST1_2 == 122),
                   n_host_113 = sum(HOST1_2 == 113),
                   n_host_93 = sum(HOST1_2 == 93),
                   n_host_101 = sum(HOST1_2 == 101),
                   n_host_2 = sum(HOST1_2 == 2),
                   n_host_96 = sum(HOST1_2 == 96),
                   n_host_NA = sum(HOST1_2 == 0)
  ) %>%
  ungroup() %>%
  mutate(pct_host_neg1 = n_host_neg1/n_hCID,
         pct_host_202 = n_host_202/n_hCID,
         pct_host_15 = n_host_15/n_hCID,
         pct_host_999 = n_host_999/n_hCID,
         pct_host_299 = n_host_299/n_hCID,
         pct_host_9107 = n_host_9107/n_hCID,
         pct_host_10 = n_host_10/n_hCID,
         pct_host_19 = n_host_19/n_hCID,
         pct_host_9102 = n_host_9102/n_hCID,
         pct_host_122 = n_host_122/n_hCID,
         pct_host_113 = n_host_113/n_hCID,
         pct_host_93 = n_host_93/n_hCID,
         pct_host_101 = n_host_101/n_hCID,
         pct_host_2 = n_host_2/n_hCID,
         pct_host_96 = n_host_96/n_hCID,
         pct_host_NA = n_host_NA/n_hCID
  )



ecolevel3_host1_summary <- centroids_host1 %>% ##summarizing for ecolevel 3 groups
  group_by(NA_L3CODE, n_NAL3Code) %>%
  dplyr::summarise(n_host_neg1 = sum(HOST1_2 == -1),
                   n_host_202 = sum(HOST1_2 == 202),
                   n_host_15 = sum(HOST1_2 == 15),
                   n_host_999 = sum(HOST1_2 == 999),
                   n_host_299 = sum(HOST1_2 == 299),
                   n_host_9107 = sum(HOST1_2 == 9107),
                   n_host_10 = sum(HOST1_2 == 10),
                   n_host_19 = sum(HOST1_2 == 19),
                   n_host_9102 = sum(HOST1_2 == 9102),
                   n_host_122 = sum(HOST1_2 == 122),
                   n_host_113 = sum(HOST1_2 == 113),
                   n_host_93 = sum(HOST1_2 == 93),
                   n_host_101 = sum(HOST1_2 == 101),
                   n_host_2 = sum(HOST1_2 == 2),
                   n_host_96 = sum(HOST1_2 == 96),
                   n_host_NA = sum(HOST1_2 == 0)
  ) %>%
  ungroup() %>%
  mutate(pct_host_neg1 = n_host_neg1/n_NAL3Code,
         pct_host_202 = n_host_202/n_NAL3Code,
         pct_host_15 = n_host_15/n_NAL3Code,
         pct_host_999 = n_host_999/n_NAL3Code,
         pct_host_299 = n_host_299/n_NAL3Code,
         pct_host_9107 = n_host_9107/n_NAL3Code,
         pct_host_10 = n_host_10/n_NAL3Code,
         pct_host_19 = n_host_19/n_NAL3Code,
         pct_host_9102 = n_host_9102/n_NAL3Code,
         pct_host_122 = n_host_122/n_NAL3Code,
         pct_host_113 = n_host_113/n_NAL3Code,
         pct_host_93 = n_host_93/n_NAL3Code,
         pct_host_101 = n_host_101/n_NAL3Code,
         pct_host_2 = n_host_2/n_NAL3Code,
         pct_host_96 = n_host_96/n_NAL3Code,
         pct_host_NA = n_host_NA/n_NAL3Code
  )






#looking at average latitude and longitude and altitude
centroids_4326_XY <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/4326_centroids_hCID-M16_kCIDs_Ecolevel3_XY.csv") 
altitude_data <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/lat-long-alt_data.csv")
colnames(altitude_data) <- c("ycoord4326","xcoord4326","alt_m")


centroids_4326_XYZ <- full_join(centroids_4326_XY, altitude_data, by = c("xcoord4326","ycoord4326"))


centroids_4326_XY_hCID_summary <- centroids_4326_XYZ %>% ##calculating average x and ys based on hierarchical clustering
  group_by(CID_M8) %>%
  summarise(avg_x = mean(xcoord4326),
            avg_y = mean(ycoord4326),
            avg_alt = mean(alt_m, na.rm = T)) %>%
  ungroup()

centroids_4326_XY_NAL3Code_summary <- centroids_4326_XYZ %>% #calculating average x and y based on ecolevel clustering
  group_by(NA_L3CODE) %>%
  summarise(avg_x = mean(xcoord4326),
            avg_y = mean(ycoord4326),
            avg_alt = mean(alt_m, na.rm = T)) %>%
  ungroup()


#for altitude: extracting just x and y coords from centroids to be fed into converter website
XYcoords<-centroids_4326_XY %>%
  select(ycoord4326,xcoord4326)
colnames(XYcoords) <- c("latitude", "longitude")
write_csv(XYcoords, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/XYcoords_forconverter.csv")


##joining average x and y frames with host 1 breakdown frames

hCID_M8_summary <- left_join(hCID_M16_host1_summary, centroids_4326_XY_hCID_summary, by = c("CID_M8"))
NAL3Code_summary <- left_join(ecolevel3_host1_summary, centroids_4326_XY_NAL3Code_summary, by = c("NA_L3CODE"))

#writing summaries as csvs
write_csv(hCID_M8_summary, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/hCID_M8_summary.csv")
write_csv(NAL3Code_summary, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/NAL3CODE_summary.csv")


hCID_M8_summary <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/hCID_M8_summary.csv")
NAL3Code_summary <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/NAL3CODE_summary.csv")





#extracting length of outbreak and time between outbreak data from time series summaries
hCID_timeseries <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/hCID_timeseries_complete_outbreakIDs_M8.csv")
NAL3CODE_timeseries <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/NAL3CODE_timeseries_complete_outbreakIDs.csv")



#length of outbreak and length of gap for hCID
#length of outbreak
hCID_timeseries <- hCID_timeseries %>%
  mutate(length_of_outbreak =0)

for (i in 1:length(hCID_timeseries$length_of_outbreak)) {
  if (hCID_timeseries$outbreak_ID_final[i] != 0) {
    hCID_timeseries$length_of_outbreak[i] <- length(which(hCID_timeseries$outbreak_ID_final == hCID_timeseries$outbreak_ID_final[i]))
  }
  if (hCID_timeseries$outbreak_ID_final[i] == 0) {
    hCID_timeseries$length_of_outbreak[i] <- 0
  }
}

#calculating outbreak summary
hCID_timeseries_outbreaksummary <- hCID_timeseries %>%
  group_by(outbreak_ID_final) %>%
  summarise(startyear = SURVEY_YR[1],
            hCID = mean(hCID),
            length_of_outbreak = mean(length_of_outbreak)) %>%
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



#length of outbreak and length of gap for NAL3CODE
#length of outbreak
defol_summary_complete_outbreakIDs_level3eco <- defol_summary_complete_outbreakIDs_level3eco %>%
  mutate(length_of_outbreak =0)

for (i in 1:length(defol_summary_complete_outbreakIDs_level3eco$length_of_outbreak)) {
  if (defol_summary_complete_outbreakIDs_level3eco$outbreak_ID_final[i] != 0) {
    defol_summary_complete_outbreakIDs_level3eco$length_of_outbreak[i] <- length(which(defol_summary_complete_outbreakIDs_level3eco$outbreak_ID_final == defol_summary_complete_outbreakIDs_level3eco$outbreak_ID_final[i]))
  }
  if (defol_summary_complete_outbreakIDs_level3eco$outbreak_ID_final[i] == 0) {
    defol_summary_complete_outbreakIDs_level3eco$length_of_outbreak[i] <- 0
  }
}

#calculating outbreak summary
NAL3CODE_timeseries_outbreaksummary <- defol_summary_complete_outbreakIDs_level3eco %>%
  group_by(outbreak_ID_final) %>%
  summarise(startyear = SURVEY_YR[1],
            Ecoregion_Code = unique(Ecoregion_Code),
            length_of_outbreak = mean(length_of_outbreak)) %>%
  filter(outbreak_ID_final != 0) %>%
  mutate(endyear = startyear + (length_of_outbreak - 1), time_to_next_outbreak = 0)

##adding time til next outbreak
for (i in 1:length(NAL3CODE_timeseries_outbreaksummary$time_to_next_outbreak)) {
  if (i < length(NAL3CODE_timeseries_outbreaksummary$time_to_next_outbreak)) {
    if (NAL3CODE_timeseries_outbreaksummary$endyear[i] < NAL3CODE_timeseries_outbreaksummary$startyear[i+1]) {
      NAL3CODE_timeseries_outbreaksummary$time_to_next_outbreak[i] <- NAL3CODE_timeseries_outbreaksummary$startyear[i+1] - NAL3CODE_timeseries_outbreaksummary$endyear[i]
    }
  }
}

#writing outbreak summaries as csv's
write_csv(hCID_timeseries_outbreaksummary, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/hCID_timeseries_outbreaksummary_M8.csv")
write_csv(NAL3CODE_timeseries_outbreaksummary, file = "./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/NAL3CODE_timeseries_outbreaksummary.csv")









##running ANOVA across clusters to see if outbreak dynamics are differing

#hCID
anova_hCID_outbreak_length <- aov(length_of_outbreak ~ as.factor(hCID), data = hCID_timeseries_outbreaksummary)
summary(anova_hCID_outbreak_length)

hCID_timeseries_outbreaksummary<-hCID_timeseries_outbreaksummary %>%
  group_by(hCID, length_of_outbreak) %>%
  mutate(n_length_hCID = n()) %>%
  ungroup() %>%
  group_by(hCID, time_to_next_outbreak) %>%
  mutate(n_gap_hCID = n())

pl<-hCID_timeseries_outbreaksummary %>%
  ggplot() +
  aes(x = hCID, y = length_of_outbreak, group = hCID, size = n_length_hCID) +
  xlab("Cluster ID") +
  ylab("Length of Outbreak (years)")+
  labs(size = "# of outbreaks") +
  geom_point()
ggsave("./Figures for Final Draft/Figure 11.pdf", device = "pdf", plot = pl)

anova_hCID_gap_length <- aov(time_to_next_outbreak ~ as.factor(hCID), data = hCID_timeseries_outbreaksummary)
summary(anova_hCID_gap_length)

pl<-hCID_timeseries_outbreaksummary %>%
  ggplot() +
  aes(x = hCID, y = time_to_next_outbreak, group = hCID, size = n_gap_hCID) +
  xlab("Cluster ID") +
  ylab("Time to Next Outbreak (years)")+
  labs(size = "# of outbreaks") +
  geom_point()
ggsave("./Figures for Final Draft/Figure 12.pdf", device = "pdf", plot = pl)



#NAL3CODE
anova_NAL3CODE_outbreak_length <- aov(length_of_outbreak ~ as.factor(NA_L3CODE), data = NAL3CODE_timeseries_outbreaksummary)
summary(anova_NAL3CODE_outbreak_length)

NAL3CODE_timeseries_outbreaksummary %>%
  ggplot() +
  aes(x = NA_L3CODE, y = length_of_outbreak, group = NA_L3CODE) +
  geom_boxplot()

anova_NAL3CODE_gap_length <- aov(time_to_next_outbreak ~ as.factor(NA_L3CODE), data = NAL3CODE_timeseries_outbreaksummary)
summary(anova_NAL3CODE_gap_length)

NAL3CODE_timeseries_outbreaksummary %>%
  ggplot() +
  aes(x = NA_L3CODE, y = time_to_next_outbreak, group = NA_L3CODE) +
  geom_boxplot()


hCID_timeseries_outbreaksummary %>%
  ggplot() +
  aes(x = startyear, y = endyear, color=hCID) +
  scale_color_viridis_c() +
  geom_point() +
  geom_abline(slope = 1, intercept = 4)

pl<-hCID_timeseries_outbreaksummary %>%
  ggplot() +
  geom_point(aes(x = startyear, y=hCID)) +
  xlab("Start Year of Outbreak") +
  ylab("Cluster ID")
#geom_point(aes(x = endyear, y= hCID, color = 2))
ggsave("./Figures for Final Draft/Figure 10.pdf", device = "pdf", plot = pl)

pl<-NAL3CODE_timeseries_outbreaksummary %>%
  ggplot() +
  geom_point(aes(x = startyear, y=Ecoregion_Code)) +
  xlab("Start Year of Outbreak") +
  ylab("Ecoregion Code")
#geom_point(aes(x = endyear, y= NA_L3CODE, color = 2))
ggsave("./Figures for Final Draft/Figure 9.pdf", device = "pdf", plot = pl)





hCID_timeseries_outbreaksummary <- read_csv("./Cluster analysis/dissolved/clusters_dispersal.gdb/defol_summaries/hCID_timeseries_outbreaksummary.csv")
hCID_timeseries_outbreaksummary_merged <- full_join(hCID_timeseries_outbreaksummary, hCID_M16_summary, by = c("hCID"))

pl<-hCID_M8_summary %>%
  ggplot() +
  aes(x = avg_y, y= avg_alt) +
  scale_color_viridis_c() +
  xlab("Average Latitude of Population (Decimal Degrees)") +
  ylab("Average Altitude of Defoliation (m)") +
  geom_point() +
  geom_text()
ggsave("./Figures for Final Draft/Figure 13_new.png", device = "png", plot = pl)

