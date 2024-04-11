# we have the addresses, and need to stem them into regions

library(stringr)
library(readxl)

data <- read_excel("data_all_00to23.xlsx")
data <- data[-which(is.na(data$Abstract)), ]
data$Time_slices <- data$Publication.Year-1999

tmp <- str_split(data$Address, "; ")
for (i in 1:length(tmp)) {
  tmp2 <- str_split_i(tmp[[i]], ", ", -1)
  # There are states in the address in USA, change them all into USA
  tmp2[which(str_detect(tmp2, "U.S.A"))] <- "USA"
  tmp2[which(str_detect(tmp2, "USA"))] <- "USA"
  tmp2[which(str_detect(tmp2, "\\."))] <- ""
  # Combine four countries in UK into one
  tmp2[which(str_detect(tmp2, "England"))] <- "UK"
  tmp2[which(str_detect(tmp2, "North Ireland"))] <- "UK"
  tmp2[which(str_detect(tmp2, "Wales"))] <- "UK"
  tmp2[which(str_detect(tmp2, "Scotland"))] <- "UK"
  tmp2 <- unique(tmp2)
  data$Regions[i] <- str_c(tmp2, collapse = "; ")
}

# Top six countries: USA 5854, UK 1076, PRC 922, France 479, Canada 472, Germany 458
Regionss <- summary(factor(c(str_split(data$Regions, "; ", simplify = TRUE))))

data$USA <- str_detect(data$Regions, "USA")
data$UK <- str_detect(data$Regions, "UK")
data$PRC <- str_detect(data$Regions, "Peoples R China")
data$France <- str_detect(data$Regions, "France")
data$Canada <- str_detect(data$Regions, "Canada")
data$Germany <- str_detect(data$Regions, "Germany")

write.csv(data[,-c(5,8,9)], file = "data_all.csv", row.names = FALSE)
