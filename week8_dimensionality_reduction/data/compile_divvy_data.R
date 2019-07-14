library(jsonlite)
stations <- read_json("divvy_stations.json", simplifyVector = FALSE)

all_stations <- NULL
for (st in stations$stationBeanList){
  if (is.null(all_stations)){
    all_stations <- tibble(id = st$id, name = st$stationName, latitude = st$latitude, longitude = st$longitude)
  } else {
    all_stations <- all_stations %>% add_row(id = st$id, name = st$stationName, latitude = st$latitude, longitude = st$longitude)
  }
}

distance_matrix <- as.matrix(dist(all_stations %>% dplyr::select(latitude, longitude)))
colnames(distance_matrix) <- all_stations$id
rownames(distance_matrix) <- all_stations$id
save(distance_matrix, file = "divvy_stations_distances.RData")
write_csv(all_stations, path = "divvy_stations.csv")
