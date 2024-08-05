localDir = "."
data.dir = file.path(localDir, "data")
model.dir = file.path(localDir, "models")
functions.dir <- file.path(localDir, "functions")
figures.dir <- file.path(localDir, "figures")

dat <- read.csv(file=file.path(data.dir, "wpfg_db_V2.csv"))
head(dat)

unique(dat$pfg)
#  [1] "ATe"     "Tdr"     "Tda"     "ATw"     "ATl"     "AT"      "Se"      "ARp"     "ARf"     "T"       "Sk"      "Sr"      "Ate"     "S"       "sk"      "AR"     
# [17] "Atw"     "Atl"     "A"       "T/ATw"   "Tda "    "Arf"     "Tda/ATe"

dat$pfg[which(dat$pfg=="sk")] <- "Sk"
dat$pfg[which(dat$pfg=="Atw")] <- "ATw"
dat$pfg[which(dat$pfg=="Arf")] <- "ARf"
dat$pfg[which(dat$pfg=="Atl")] <- "ATl"
dat$pfg[which(dat$pfg=="Ate")] <- "ATe"
dat$pfg[which(dat$pfg=="Tda ")] <- "Tda"
dat$pfg[which(dat$pfg=="Tda/ATe")] <- "Tda"
dat$pfg[which(dat$pfg=="T/ATw")] <- "T"


bad <- which(dat$pfg == "A" |dat$pfg == "S"|dat$pfg == "AR"| dat$pfg == "T"|dat$pfg == "AT")
dat.bad <- dat[bad,]
dat.gd <- dat[-bad,]

# how many spp records? ####
tab <- table(dat.gd$match_pfg_nom)
length(tab) # 3413

xclass <- c()
for(i in 1:8) xclass[i] <- (sum(tab == i))
xclass # [1] 1890  838  335  167  102   41   19   13

names(tab)[tab==8]

xtab <- table(dat.gd$match_pfg_nom, dat.gd$pfg)
cnts <- apply(xtab>0,1, sum)
xtab[which(cnts>1),]

unique(dat.gd$pfg)

# reshape unique records to wide
library(reshape2)
db1 <- dcast(dat.gd, match_pfg_nom ~ pfg, fun.aggregate=length, value.var="pfg")
raw_database <- db1
save(raw_database, file=file.path(data.dir,"raw_database_V2.RData"))

head(raw_database)
dim(raw_database)

gen.only <- grep(pattern = "sp.", x = db1[,1])

spp_database <- raw_database[-gen.only,]

rownames(spp_database) <- spp_database[,1]
spp_database <- spp_database[,-1]
head(spp_database)
save(spp_database, file=file.path(data.dir,"spp_database.RData"))

rownames(spp_database)

pfg.assign <- apply(spp_database, 1, sum)
pfg.count <- apply(spp_database>0, 1, sum); table(pfg.count)
# pfg.count
#    1    2    3 
# 2191  268   13 
2191+268 + 13

hist(pfg.assign)
hist(pfg.count)


disp2.pfg <- spp_database[which(pfg.count == 2),]
write.csv(disp2.pfg, file = file.path(data.dir,"disp2pfg.csv"))

library(plot.matrix)

d2.pa <- ifelse(disp2.pfg>0,1,0)
plot(as.matrix(d2.pa))
pfg <- colnames(d2.pa)
combs <- combn(pfg, m=2)

counts<- numeric()
for(i in 1:ncol(combs)){
  # i = 2
  x <- combs[,i]
  c1 <- which(colnames(disp2.pfg) == x[1])
  c2 <- which(colnames(disp2.pfg) == x[2])
  counts[i] <- length(which(disp2.pfg[,c1] > 0 & disp2.pfg[,c2] >0)) 
}

wpfgX2 <- data.frame(grp1 = combs[1,], grp2 = combs[2,], count = counts)
wpfgX2 <- wpfgX2[order(wpfgX2$count, decreasing =TRUE),]
write.csv(wpfgX2, file=file.path(data.dir, "wpfgX2.csv"))

# 3 different WPFG
pfgX3 <- rownames(spp_database)[which(pfg.count==3)]
disp3.pfg <- spp_database[which(rownames(spp_database) %in% pfgX3),]
write.csv(disp3.pfg, file = file.path(data.dir,"disp3pfg.csv"))

d3.pa <- ifelse(disp3.pfg>0,1,0)
plot(d3.pa, main ="", xlab="WPFG",ylab="Species", col = c("grey","red"), key= NULL)#as.matrix(disp3.pfg[disp3.pfg>0,1,0]))

round(apply(disp2.pfg>0, 2, sum)/sum(apply(disp2.pfg>0,2,sum)),2)

# ARf   ARp  ATe  ATl  ATw   Se   Sk   Sr  Tda  Tdr 
#  8     41   94   47   14   12    8    7   164  141 
# 0.01 0.08 0.18 0.09 0.03 0.02 0.01 0.01  0.31 0.26 

db.consensus <- spp_database[which(pfg.count == 1 & pfg.assign>1),]


db.conf <- apply(spp_database, 2, function(x,y) x/y, y = pfg.assign)

head(db.conf)
dim(db.conf)
save(db.conf, file=file.path(data.dir,"db_conf.RData"))
head(db.conf)

