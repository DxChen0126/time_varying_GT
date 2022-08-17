######### check with raw data
######## Xu's raw data's drawbacks in exposure information

library(readxl)

pair1 <- read_excel("allpossiblepairs.xlsx", sheet = 2)
pair2 <- read_excel("allpossiblepairs.xlsx", sheet = 3)


ID.infector1 <- unique(pair1$`ID(I)`)
ID.infectee1 <- unique(pair1$`ID(S)`)
ID.data1 <- unique(union(ID.infector1, ID.infectee1))

ID.infector2 <- unique(pair2$`ID(I)`)
ID.infectee2 <- unique(pair2$`ID(S)`)

ID.data2 <- unique(union(ID.infector2, ID.infectee2))

ID.all <- unique(union(ID.data1, ID.data2))
length(ID.all) # 2989


### ID(I) Age(I) Gender(I) ID(S) Age(S) Gender(S)

length(ID.infector1) # 810
length(ID.infector2) # 532
length(ID.infectee1) # 1274
length(ID.infectee2) # 717

age.infector1 <- sex.infector1 <- numeric(810)

for(i in 1:810){
  ind.tmp <- which(pair1$`ID(I)` == ID.infector1[i])[1]
  age.infector1[i] <- pair1$`Age(I)`[ind.tmp]
  sex.infector1[i] <- pair1$`Gender(I)`[ind.tmp]
}

age.infector2 <- sex.infector2 <- numeric(532)

for(i in 1:532){
  ind.tmp <- which(pair2$`ID(I)` == ID.infector2[i])[1]
  age.infector2[i] <- pair2$`Age(I)`[ind.tmp]
  sex.infector2[i] <- pair2$`Gender(I)`[ind.tmp]
}

age.infectee1 <- sex.infectee1 <- numeric(1274)

for(i in 1:1274){
  ind.tmp <- which(pair1$`ID(S)` == ID.infectee1[i])[1]
  age.infectee1[i] <- pair1$`Age(S)`[ind.tmp]
  sex.infectee1[i] <- pair1$`Gender(S)`[ind.tmp]
}


age.infectee2 <- sex.infectee2 <- numeric(717)

for(i in 1:717){
  ind.tmp <- which(pair2$`ID(S)` == ID.infectee2[i])[1]
  age.infectee2[i] <- pair2$`Age(S)`[ind.tmp]
  sex.infectee2[i] <- pair2$`Gender(S)`[ind.tmp]
}

allcase.combine <- data.frame(
  ID = c(ID.infector1, ID.infector2, ID.infectee1, ID.infectee2),
  Age = c(age.infector1, age.infector2, age.infectee1, age.infectee2),
  Sex = c(sex.infector1, sex.infector2, sex.infectee1, sex.infectee2)
)

# remove duplicates
ind.uniq <- numeric(2989)

for(i in 1:2989){
  ind.uniq[i] <- which(allcase.combine$ID == ID.all[i])[1]
}

allcase.combine.uniq <- allcase.combine[ind.uniq,]

summary(as.numeric(allcase.combine.uniq$Age)) # med 46 mean 45.5 1stq 33 3rdq 58 NA 133
table(allcase.combine.uniq$Sex) # male 120 + 1378 = 1498 female 62 + 1385 = 1447 NA 44

allpairs <- read_excel("allpossiblepairs.xlsx", sheet = 4)

ID.I <- unique(allpairs$Infector.ID)
ID.S <- unique(allpairs$Infectee.ID)

infector.ind <- numeric(428)

for(i in 1:428){
  infector.ind[i] <- which(allpairs$Infector.ID == ID.I[i])[1]
}

summary(as.numeric(allpairs$age.infectee)) # med 49 mean 47.3 1stq 34 3rdq 60.5 NA 14
table(allpairs$infecteesex) # female 13 + 320 = 333 male 18 + 278 = 296

summary(as.numeric(allpairs[infector.ind,]$age.infector)) # med 47 mean 47.3 1stq 37 3rdq 57 NA 7
table(allpairs[infector.ind,]$infectorsex) # female 170 + 3 = 173 male 20 + 234 = 254 NA 1





