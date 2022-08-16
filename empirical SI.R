rm(list = ls())

allpairs <- read.csv("C:/Users/Dongxuan Chen/SPH Dropbox/Dongxuan Chen/Transmissibility_Methods/Data/China/allpairs0908.csv")

library(ggplot2)
library(ggpubr)

# format date

colexpo <- grep(c("expo"), colnames(allpairs)) 
for(i in colexpo){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}

colonset <- grep(c("Onset"), colnames(allpairs)) 
for(i in colonset){
  allpairs[, i] <- as.Date(allpairs[, i], "%m/%d/%Y")
}


## forward
startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-15"), to = as.Date("2020-01-30"), by = "1 day"))
endseq <- c(as.Date("2020-01-20"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-04"), by = "1 day"),
            as.Date("2020-02-29"))


n <- length(startseq)

samplesize <- numeric(n)
empSI <- vector("list", n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  
  # i <- 1
  subdata <- subset(allpairs, Onset_Infector >= startseq[i] & Onset_Infector <= endseq[i])
  samplesize[i] <- nrow(subdata)
  empSI[[i]] <- as.numeric(subdata$Onset_Infectee - subdata$Onset_Infector)
  setTxtProgressBar(progressbar, i)
}

timewindow <- c()
for(i in 1:n){
  timewindow <- c(timewindow, rep(paste0(format(startseq[i], format =  "%b %d"), " - ",
                                         format(endseq[i], format =  "%b %d")), samplesize[i]))
}

dfempSI <- data.frame(
  empSI = unlist(empSI),
  timewindow = timewindow
)
write.csv(dfempSI, "forwardempSI.csv")


## backward

startseq <- c(as.Date("2020-01-01"), seq.Date(from = as.Date("2020-01-21"), to = as.Date("2020-02-05"), by = "1 day"))
endseq <- c(as.Date("2020-01-26"), seq.Date(from = as.Date("2020-01-27"), to = as.Date("2020-02-10"), by = "1 day"),
            as.Date("2020-02-29"))


n <- length(startseq)

samplesize <- numeric(n)
empSI <- vector("list", n)

progressbar <- txtProgressBar(min = 0, max = n, style = 3)

for(i in 1:n){
  
  # i <- 1
  subdata <- subset(allpairs, Onset_Infectee >= startseq[i] & Onset_Infectee <= endseq[i])
  samplesize[i] <- nrow(subdata)
  empSI[[i]] <- as.numeric(subdata$Onset_Infectee - subdata$Onset_Infector)
  setTxtProgressBar(progressbar, i)
}

timewindow <- c()
for(i in 1:n){
  timewindow <- c(timewindow, rep(paste0(format(startseq[i], format =  "%b %d"), " - ",
                                         format(endseq[i], format =  "%b %d")), samplesize[i]))
}

dfempSI <- data.frame(
  empSI = unlist(empSI),
  timewindow = timewindow
)

summary(dfempSI$empSI)

write.csv(dfempSI, "backwardempSI.csv")






