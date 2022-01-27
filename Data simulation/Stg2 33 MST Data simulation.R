################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#                             DIF in CORE ONLY 
#                             Magnitude = .25
#                           for 11% Amount of DIF 
################################################################################
#Install and call package (this version was specifically created for Leslie Rutkowski et al.)
#install.packages("C:/Users/mbvaldiv/Box/DIF in MST/Simulation MST/mstRforLR/v1.1/mstRforLR_1.1.zip", repos = NULL, type = "win.binary")
install.packages("G:/DIF in MST Routing in ILSA/Simulation MST/mstRforLR/v1.2/mstRforLR_1.2.zip", repos = NULL, type = "win.binary")
library("mstRforLR") 
setwd("G:\\DIF in MST Routing in ILSA\\Simulation MST\\mstR first change")
source("randomMSTmv1.R")
install.packages("R.utils")
library("R.utils")
#install.packages("TAM")
#library("TAM")
#################################################################################
#                                 Data call                                     #
#################################################################################
#Set pathway 
path <- "G:\\DIF in MST Routing in ILSA\\Simulation MST\\Data"
setwd(path)


# Call true item parameters
itempar = read.csv("Test assemble/core.test.25.csv", header = T)
itembank = read.csv("item bank for modules.csv", header = T)
Slope = itembank[, 4]
Location = itempar[, 4] # this is for selecting the DIF introduced to only 11%
it1 = data.frame(Slope, Location)
it1$c = 0
it1$d = 1
colnames(it1)[colnames(it1)=="Slope"] <- "a"
colnames(it1)[colnames(it1)=="Location"] <- "b"
it = it1

#################################################################################
#                              Assignment of modules                            #
#################################################################################
modules <- matrix(0, 72, 6)
modules[1:12, 1] <- modules[13:24, 2] <- modules[25:36, 3] <- modules[37:48, 4] <- modules[49:60, 5] <- modules[61:72, 6]<- 1

colSums(modules)

### Attach the module membership to each item
it_mod <- cbind(it, mod=0)
it_mod[1:12, "mod"]<- 1
it_mod[13:24, "mod"]<- 2
it_mod[25:36, "mod"]<- 3
it_mod[37:48, "mod"]<- 4
it_mod[49:60, "mod"]<- 5
it_mod[61:72, "mod"]<- 6

it_mod <- data.frame(it_mod)


# Cration of the transition matrix to define a 123 MST
trans <- matrix(0, 6, 6)
trans[1, 2:3] <- trans[2, 4:5] <- trans[3, 5:6]<-1
trans

plot.mst(trans)


# List for estimation of  theta
start<- list(fixModule = 1)
test1 <- list(moduleSelect="MFI", method = "EAP",  prob=c(1))
final<- list(method = "EAP")

#res1 <- randomMSTmv1(trueTheta = .99, itemBank = it, modules = modules, transMatrix = trans,
#                  start = start, test = test1, final = final, allTheta = TRUE)

#res1$testItems




#data[ ii, paste0("i.", c(item_pool[prov$testItems, "id"])) ] <- prov$pattern
#    data[ ii, c("OP1", "OP2")] <- prov$best.module
#    data[ ii, c("MOD1", "MOD2", "MOD3")] <- c( prov$selected.modules )

#data[1, c("TH0")] <- res1$trueTheta
#data[1, paste0("i.", c(res1$testItems)) ] <- res1$pattern
#data[1, c("OP1", "OP2")] <- res1$best.module
#data[1, c("MOD1", "MOD2", "MOD3")] <- res1$selected.modules
#data[1, c("th1", "th2", "th3")] <- res1$thetaProv


#################################################################################
#                                 MST SIMULATION                                #
#################################################################################

#===============================================================================
#                    ROUTING PROBABILITY = 1
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = ".25", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = "1", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(1), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.25.33.p1", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm

#===============================================================================
#                    ROUTING PROBABILITY = .70
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = ".25", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = ".70", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(.7), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.25.33.p2", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm

#===============================================================================
#                    ROUTING PROBABILITY = .50
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = ".25", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = ".50", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(.5), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.25.33.p3", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm

################################################################################
################################################################################
################################################################################
################################################################################
#                             DIF in CORE ONLY 
#                             Magnitude = .50
#                           for 11% Amount of DIF 
################################################################################
#Install and call package (this version was specifically created for Leslie Rutkowski et al.)

#library("mstRforLR") 
#setwd("G:\\DIF in MST Routing in ILSA\\Simulation MST\\mstR first change")
#source("randomMSTmv1.R")

#install.packages("R.utils")
#library("R.utils")
#install.packages("TAM")
#library("TAM")
#################################################################################
#                                 Data call                                     #
#################################################################################
#Set pathway 
path <- "G:\\DIF in MST Routing in ILSA\\Simulation MST\\Data"
setwd(path)


# Call true item parameters
itempar = read.csv("Test assemble/core.test.50.csv", header = T)
itembank = read.csv("item bank for modules.csv", header = T)
Slope = itembank[, 4]
Location = itempar[, 4] # this is for selecting the DIF introduced to only 11%
it1 = data.frame(Slope, Location)
it1$c = 0
it1$d = 1
colnames(it1)[colnames(it1)=="Slope"] <- "a"
colnames(it1)[colnames(it1)=="Location"] <- "b"
it = it1

#################################################################################
#                              Assignment of modules                            #
#################################################################################
modules <- matrix(0, 72, 6)
modules[1:12, 1] <- modules[13:24, 2] <- modules[25:36, 3] <- modules[37:48, 4] <- modules[49:60, 5] <- modules[61:72, 6]<- 1

colSums(modules)

### Attach the module membership to each item
it_mod <- cbind(it, mod=0)
it_mod[1:12, "mod"]<- 1
it_mod[13:24, "mod"]<- 2
it_mod[25:36, "mod"]<- 3
it_mod[37:48, "mod"]<- 4
it_mod[49:60, "mod"]<- 5
it_mod[61:72, "mod"]<- 6

it_mod <- data.frame(it_mod)


# Cration of the transition matrix to define a 123 MST
trans <- matrix(0, 6, 6)
trans[1, 2:3] <- trans[2, 4:5] <- trans[3, 5:6]<-1
trans

plot.mst(trans)


# List for estimation of  theta
start<- list(fixModule = 1)
test1 <- list(moduleSelect="MFI", method = "EAP",  prob=c(1))
final<- list(method = "EAP")

#res1 <- randomMSTmv1(trueTheta = .99, itemBank = it, modules = modules, transMatrix = trans,
#                  start = start, test = test1, final = final, allTheta = TRUE)

#res1$testItems




#data[ ii, paste0("i.", c(item_pool[prov$testItems, "id"])) ] <- prov$pattern
#    data[ ii, c("OP1", "OP2")] <- prov$best.module
#    data[ ii, c("MOD1", "MOD2", "MOD3")] <- c( prov$selected.modules )

#data[1, c("TH0")] <- res1$trueTheta
#data[1, paste0("i.", c(res1$testItems)) ] <- res1$pattern
#data[1, c("OP1", "OP2")] <- res1$best.module
#data[1, c("MOD1", "MOD2", "MOD3")] <- res1$selected.modules
#data[1, c("th1", "th2", "th3")] <- res1$thetaProv


#################################################################################
#                                 MST SIMULATION                                #
#################################################################################

#===============================================================================
#                    ROUTING PROBABILITY = 1
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = ".50", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = "1", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(1), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.50.33.p1", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm

#===============================================================================
#                    ROUTING PROBABILITY = .70
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = ".50", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = ".70", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(.7), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.50.33.p2", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm

#===============================================================================
#                    ROUTING PROBABILITY = .50
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = ".50", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = ".50", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(.5), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.50.33.p3", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm

################################################################################
################################################################################
################################################################################
################################################################################
#                             DIF in CORE ONLY 
#                             Magnitude = 1
#                           for 11% Amount of DIF 
################################################################################
#Install and call package (this version was specifically created for Leslie Rutkowski et al.)
#install.packages("C:/Users/mbvaldiv/Box/DIF in MST/Simulation MST/mstRforLR/v1.1/mstRforLR_1.1.zip", repos = NULL, type = "win.binary")

#library("mstRforLR") 
#install.packages("R.utils")
#library("R.utils")
#install.packages("TAM")
#library("TAM")
#################################################################################
#                                 Data call                                     #
#################################################################################
#Set pathway 
path <- "G:\\DIF in MST Routing in ILSA\\Simulation MST\\Data"
setwd(path)


# Call true item parameters
itempar = read.csv("Test assemble/core.test.1.csv", header = T)
itembank = read.csv("item bank for modules.csv", header = T)
Slope = itembank[, 4]
Location = itempar[, 4] # this is for selecting the DIF introduced to only 11%
it1 = data.frame(Slope, Location)
it1$c = 0
it1$d = 1
colnames(it1)[colnames(it1)=="Slope"] <- "a"
colnames(it1)[colnames(it1)=="Location"] <- "b"
it = it1

#################################################################################
#                              Assignment of modules                            #
#################################################################################
modules <- matrix(0, 72, 6)
modules[1:12, 1] <- modules[13:24, 2] <- modules[25:36, 3] <- modules[37:48, 4] <- modules[49:60, 5] <- modules[61:72, 6]<- 1

colSums(modules)

### Attach the module membership to each item
it_mod <- cbind(it, mod=0)
it_mod[1:12, "mod"]<- 1
it_mod[13:24, "mod"]<- 2
it_mod[25:36, "mod"]<- 3
it_mod[37:48, "mod"]<- 4
it_mod[49:60, "mod"]<- 5
it_mod[61:72, "mod"]<- 6

it_mod <- data.frame(it_mod)


# Cration of the transition matrix to define a 123 MST
trans <- matrix(0, 6, 6)
trans[1, 2:3] <- trans[2, 4:5] <- trans[3, 5:6]<-1
trans

plot.mst(trans)


# List for estimation of  theta
start<- list(fixModule = 1)
test1 <- list(moduleSelect="MFI", method = "EAP",  prob=c(1))
final<- list(method = "EAP")

#res1 <- randomMSTmv1(trueTheta = .99, itemBank = it, modules = modules, transMatrix = trans,
#                  start = start, test = test1, final = final, allTheta = TRUE)

#res1$testItems




#data[ ii, paste0("i.", c(item_pool[prov$testItems, "id"])) ] <- prov$pattern
#    data[ ii, c("OP1", "OP2")] <- prov$best.module
#    data[ ii, c("MOD1", "MOD2", "MOD3")] <- c( prov$selected.modules )

#data[1, c("TH0")] <- res1$trueTheta
#data[1, paste0("i.", c(res1$testItems)) ] <- res1$pattern
#data[1, c("OP1", "OP2")] <- res1$best.module
#data[1, c("MOD1", "MOD2", "MOD3")] <- res1$selected.modules
#data[1, c("th1", "th2", "th3")] <- res1$thetaProv


#################################################################################
#                                 MST SIMULATION                                #
#################################################################################

#===============================================================================
#                    ROUTING PROBABILITY = 1
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = "1", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = "1", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(1), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.1.33.p1", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm

#===============================================================================
#                    ROUTING PROBABILITY = .70
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = "1", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = ".70", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(.7), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.1.33.p2", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm

#===============================================================================
#                    ROUTING PROBABILITY = .70
#===============================================================================
## Time
ptm <- proc.time()

path <- "G:/DIF in MST Results/"
Ntheta = 1000
NoReps <- c(3000)
nrSamples = 100
theta <- list(mode="vector",length=nrSamples)
data <- expand.grid(rep = NA,
                    Group = "Focal",
                    Country = c(rep("MAR", 1000), rep("ZAF", 1000), rep("SAU", 1000)), 
                    DIFmag = "1", 
                    DIFper = "33",
                    DIFloc = "Stg 2", 
                    Prob = ".50", 
                    TH0 = NA)

data <- data.frame(data,
                   th1 = NA, th2 = NA, th3 = NA, 
                   seT = NA,
                   set1 = NA, set2 = NA, set3 = NA,
                   OP1 = NA, OP2 = NA, 
                   MOD1 = NA, MOD2 = NA, MOD3 = NA, 
                   i = matrix(NA, ncol = 72))


# MST loop 
for (i in 1:100){ 
  set.seed(i*1992)
  th <- theta[[i]] <- c(rnorm(Ntheta, mean= -1.16, sd = .80 ), 
                        rnorm(Ntheta, mean= -1.28, sd = .87),
                        rnorm(Ntheta, mean= -1.32, sd = .86)) 
  ressim <- NULL
  for (k in 1:NoReps[1]){
    test2 <- list(method = "EAP", moduleSelect= "MFI", prob=c(.5), seed.prob = (i*k+47405))
    res2 <- randomMSTmv1(trueTheta = th[k], itemBank = it, modules = modules, transMatrix = trans,
                      start = start, 
                      test = test2, 
                      final = final)
    
    data[k, c("rep")] <- i
    data[k, c("TH0")] <- res2$trueTheta
    data[k, paste0("i.", c(res2$testItems)) ] <- res2$pattern
    data[k, c("OP1", "OP2")] <- res2$best.module
    data[k, c("MOD1", "MOD2", "MOD3")] <- res2$selected.modules
    data[k, c("th1", "th2", "th3")] <- res2$thetaProv
    data[k, c("seT")] <- res2$seFinal
    data[k, c("set1", "set2", "set3")] <- res2$seProv
    #ressim = cbind(res, data)
  }
  saveObject(data, file = paste0(path, "stg2.1.33.p3", "_i", i, ".Rbin"))
  
}

#Stop the clock
proc.time( ) - ptm






