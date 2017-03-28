library(plyr)
library(dplyr, lib = "/home/bmayer/CMVGrowth/libraryFolder")
library(doParallel)

my_session_info <- devtools::session_info()

platform <- my_session_info[[1]]
packages <- my_session_info[[2]]

# TABLE 1
values = matrix(unlist(platform), nrow = length(platform))

my_session_info1 <- data.frame(
  name = names(platform),
  value = values[,1])

names(my_session_info1)[2] <- "Session info."

my_session_info2 <- as.data.frame(
  matrix(unlist(packages), ncol = length(packages))
) %>% select(-V2) # Only want attached packages
names(my_session_info2) = names(packages)[-2]


write.table(my_session_info1, "r-server-session-info.csv", sep = "\t", row.names = F)
write.table(my_session_info2, "r-server-package-info.csv", sep = "\t", row.names = F)

