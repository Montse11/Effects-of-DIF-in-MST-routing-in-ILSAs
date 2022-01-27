#############################################################################################
#                                     NoDIF
#############################################################################################


###############################################################################
#                             DIF in MST
#                         Item Parameter Recovery
###############################################################################


#install.packages("TAM")
library("TAM")
library("R.utils")
library("tidyverse")
library("reshape2")
library("parallel")
library("CDM")
library("doParallel")
library("dplyr")
#=========================================================================================
#                           NoDIF.1.11.p1
#=========================================================================================
rm(list=ls())

#Call data



#data1= loadObject("D:/DIF in MST Routing in ILSA/DIF in MST Results/NoDIF/Reference/Reference.p1/NoDIF.R.p1_i50.Rbin")

#acountry = loadObject("D:/DIF in MST Routing in ILSA/DIF in MST Analysis/Results Latent Regression/Person Parameter Results/Country/NoDIF.1_p3.Rbin")

# Call true item parameters
# itempar = read.csv("D:/DIF in MST Routing in ILSA/Simulation MST/Data/Test assemble/NoDIF.test.1.csv", header = T)
# itembank = read.csv("D:/DIF in MST Routing in ILSA/Simulation MST/Data/item bank for modules.csv", header = T)
# Slope = itembank[, 4]
# Location = itempar[, 3] # this is for selecting the DIF introduced to only 11%
# it = data.frame(a = Slope, b = Location)
#call background data
background_10 = read.csv("D:/DIF in MST Routing in ILSA/Latent regression/TIMSS 2015/Background variables 10.csv", header = TRUE)
#background_10$country = factor(background_10$country, levels = c("SGP", "KOR", "TWN", "SWE", "MLT", "MYS", "MAR", "ZAF", "SAU"))

#==============================================================================
##          Latent Item regression
#==============================================================================

latent_reg = function(data, j, k, l, m){
  data1 <- NULL
  data1 <-subset(data, rep==j)
  data1 <- data1[order(data1$Country),]
  data1$mstpath = ifelse(data1$MOD2 == 2 & data1$MOD3 == 4, "Core-Low-Low",
                         ifelse(data1$MOD2 == 2 & data1$MOD3 == 5, "Core-Low-Med",
                                ifelse(data1$MOD2 == 3 & data1$MOD3 == 6, "Core-High-High",
                                       ifelse(data1$MOD2 == 3 & data1$MOD3 == 5, "Core-High-Med",
                                              ifelse(data1$MOD2 == 2 & data1$MOD3 == 6, "Core-Low-High",
                                                     ifelse(data1$MOD2 == 3 & data1$MOD3 == 4, "Core-High-Low", NA))))))
  
  #estimate item parameters
  tam.mod <- TAM::tam.mml.2pl( data1[, paste0("i.", 1:72)],
                               irtmodel = "2PL", 
                               est.variance = TRUE,
                               group = data1$C_mstpath,
                               control = list(maxiter = 10000), 
                               verbose = F)
  
  # rescale estimated item parameters using mean/sigma method
  ## extract estimated parameters
  # extract (estimated) item parameters https://embracingheterogeneity.github.io/website/01a_LR.html#the-tam-version
  est_a <- tam.mod$item$B.Cat1.Dim1 #estimated descrimination/slope
  est_b <- tam.mod$item$AXsi_.Cat1 / tam.mod$item$B.Cat1.Dim1  #estimated difficulty/location
  
  u <- mean(est_a)/mean(it$a)
  v <- mean(it$b) - u*mean(est_b)
  
  rs_est_a <- est_a/u          #rescaled item discrimination
  rs_est_b <- u*est_b + v      #rescaled item difficulty
  
  
  
  ## Specify item parameters (1) from TAM; (2) from rescaling
  ## Save for scoring and for parameter recovery
  tam_ip <- tam.mod$item
  tam_ip$b <- rs_est_b
  tam_ip$a <- rs_est_a
  tam_ip$xsi.fixed.estimated <- tam_ip$a*tam_ip$b
  tam_ip$B.fixed.estimated <- tam_ip$a
  
  #
  tam.mod$xsi.fixed.estimated[, 2] <-  tam_ip$xsi.fixed.estimated
  tam.mod$B.fixed.estimated[seq(2, 2* nrow(it), by = 2), 4] <- tam_ip$B.fixed.estimated 
  
  # tam.mod$xsi.fixed.estimated[, 2] <- rs_est_a * rs_est_b
  # tam.mod$B.fixed.estimated[seq(2, 2* nrow(it), by = 2), 4] <- rs_est_a
  # 
  estimated_items = data.frame(rep = rep(j, 72),
                               item = tam_ip$item,
                               DIFmag = rep(m, 72),
                               DIFloc = rep("Stage 2", 72),
                               DIFper = rep(l, 72),
                               Prob = rep(k, 72),
                               it$a,
                               tam_a = est_a,
                               rs_est_a = tam_ip$a,
                               it$b,
                               tam_b = est_b,
                               rs_est_b= tam_ip$b)
  
  #Save both item parameters and the object
  path5 <- "D:/DIF in MST Routing in ILSA/DIF in MST Analysis/Results Latent Regression/Rescaled Item Parameters/"
  saveObject(estimated_items, file = paste0(path5,"Stg2", m, "_", l, "_", k, "_p", j,".Rbin"))
  path6 <- "D:/DIF in MST Routing in ILSA/DIF in MST Analysis/Results Latent Regression/Item Parameters/"
  saveObject(tam.mod, file = paste0(path6,"Stg2", m, "_", l, "_", k, "_p", j,".Rbin"))
  # # setting the fixed parameters in format
  
  tam.mod.mg <- TAM::tam.mml.2pl( data1[, paste0("i.", 1:72)], 
                                  irtmodel = "2PL",
                                  est.variance = TRUE, 
                                  group = data1$C_mstpath,
                                  xsi.inits = tam.mod$xsi.fixed.estimated,
                                  xsi.fixed = tam.mod$xsi.fixed.estimated,
                                  B.fixed = tam.mod$B.fixed.estimated,
                                  control = list(maxiter = 10000), verbose = F)
  
  se1 <- TAM::tam.se( tam.mod.mg ) 
  # #get the likelihood
  # like <- CDM::IRT.likelihood(tam.mod.mg)
  # X <- background_10[2:10]
  # group <- background_10[,11]
  # 
  # #build the tam model
  # latent_model <- TAM::tam.latreg(like, 
  #                                 Y=X, 
  #                                 group = group,
  #                                 pid=as.numeric(rownames(background_10)), 
  #                                 verbose = FALSE)
  # 
  #Save all the objects
  #path2 <- "D:/DIF in MST Routing in ILSA/DIF in MST Analysis/Results Latent Regression/Latent model object/"
  #saveObject(latent_model, file = paste0(path2,"Stg2", m, "_", l, "_", k, "_p", j,".Rbin"))
  
  #PLAUSIBLE VALUES
  pvmod <- TAM::tam.pv( tam.mod.mg, 
                        nplausible = 10, #10 plausible values
                        samp.regr = TRUE)
  
  data1 <- cbind(data1, pvmod$pv)
  #=============================================================================
  # For data Analysis
  #=============================================================================
  ## Summarizing by country and mstpath
  pv <- stats::aggregate(data1[, paste0("PV", 1:10, ".Dim1")], 
                         by = list(data1$Country, data1$mstpath), FUN = mean)
  
  pv$PV <- rowMeans( pv[, paste0("PV", 1:10, ".Dim1")])
  
  pv_n <- stats::aggregate(data1[, paste0("PV", 1, ".Dim1")], 
                           by = list(data1$Country, data1$mstpath), FUN = length)
  pv <- merge(pv_n, pv, by = c("Group.1", "Group.2"))
  
  th0 <- stats::aggregate(data1[, "TH0"], 
                          by = list(data1$Country, data1$mstpath), FUN = mean)
  th0andpv <- merge(th0, pv, by = c("Group.1", "Group.2"))
  colnames(th0andpv)[1] <- "Country"
  colnames(th0andpv)[2] <- "mstpath"
  colnames(th0andpv)[3] <- "TH0"
  colnames(th0andpv)[4] <- "N"
  
  
  pv1 = NULL
  pv1 <- data.frame(Country = pv$Group.1,
                    mstpath = pv$Group.2,
                    PV1sumS = (pv$PV1.Dim1 - pv$PV)^2,
                    PV2sumS = (pv$PV2.Dim1 - pv$PV)^2,
                    PV3sumS = (pv$PV3.Dim1 - pv$PV)^2,
                    PV4sumS = (pv$PV4.Dim1 - pv$PV)^2,
                    PV5sumS = (pv$PV5.Dim1 - pv$PV)^2,
                    PV6sumS = (pv$PV6.Dim1 - pv$PV)^2,
                    PV7sumS = (pv$PV7.Dim1 - pv$PV)^2,
                    PV8sumS = (pv$PV8.Dim1 - pv$PV)^2,
                    PV9sumS = (pv$PV9.Dim1 - pv$PV)^2,
                    PV10sumS = (pv$PV10.Dim1 - pv$PV)^2)
  
  pv1$Mvar_PV = rowSums(pv1[ , 3:12])/9
  
  th0pvsd <- merge(th0andpv[,c(1:4, 15)], pv1[, c(1,2,13)], by = c("Country", "mstpath") )
  #we need to order in the same order se1 is
  th0pvsd$interaction = interaction(th0pvsd$Country, th0pvsd$mstpath)
  th0pvsd = th0pvsd[order(th0pvsd$interaction),]
  
  # ## Summarizing by country
  # 
  # pv2 <- stats::aggregate(data1[, paste0("PV", 1:10, ".Dim1")], 
  #                        by = list(data1$Country), FUN = mean)
  # pv2$PV <- rowMeans( pv2[, paste0("PV", 1:10, ".Dim1")])
  # th0_2 <- stats::aggregate(data1[, "TH0"], 
  #                         by = list(data1$Country), FUN = mean)
  # th0andpv_2 <- merge(th0_2, pv2[, c("Group.1", "PV")], by = c("Group.1"))
  # colnames(th0andpv_2)[1] <- "Country"
  # colnames(th0andpv_2)[2] <- "TH0"
  # 
  # pv3 = NULL
  # pv3 <- data.frame(Country = pv2$Group.1,
  #                   PV1sumS = (pv2$PV1.Dim1 - pv2$PV)^2,
  #                   PV2sumS = (pv2$PV2.Dim1 - pv2$PV)^2,
  #                   PV3sumS = (pv2$PV3.Dim1 - pv2$PV)^2,
  #                   PV4sumS = (pv2$PV4.Dim1 - pv2$PV)^2,
  #                   PV5sumS = (pv2$PV5.Dim1 - pv2$PV)^2,
  #                   PV6sumS = (pv2$PV6.Dim1 - pv2$PV)^2,
  #                   PV7sumS = (pv2$PV7.Dim1 - pv2$PV)^2,
  #                   PV8sumS = (pv2$PV8.Dim1 - pv2$PV)^2,
  #                   PV9sumS = (pv2$PV9.Dim1 - pv2$PV)^2,
  #                   PV10sumS = (pv2$PV10.Dim1 - pv2$PV)^2)
  # 
  # pv3$Mvar_PV = rowSums(pv3[ , 2:11])/9
  # 
  # th0pvsd_2 <- merge(th0andpv_2, pv3[, c(1,12)], by = c("Country") )
  # 
  # th0pvsd_2 = th0pvsd_2[order(th0pvsd_2$Country),]
  
  
  #data1 <-  NoDIF1[which(NoDIF1$rep==j & NoDIF1$Prob==Prob[p] &  NoDIF1$DIFper == per[q]), 1:20]
  # 
  # theta <- aggregate( TH0 ~ Country, data = data1, mean)
  # est_theta <- stats::aggregate(latent_model[["person"]]$EAP, 
  #                               by = list(data2$Country, data2$mstpath), data = data1, mean ) 
  # 
  # colnames(est_theta)[2]  <- "Eth"
  # sd <- stats::aggregate(latent_model[["person"]]$SD.EAP ~ Country, data = data1, FUN = stats::sd )
  # colnames(sd)[2] <- "Esd"
  
  
  Results = NULL
  
  
  Results <- data.frame(rep = rep(j, nrow(th0pvsd)),
                        DIFmag = rep(m, nrow(th0pvsd)),
                        DIFloc = rep("Stage 2", nrow(th0pvsd)),
                        DIFper = rep(l, nrow(th0pvsd)),
                        Prob = rep(k, nrow(th0pvsd)),
                        th0pvsd, 
                        se1$beta)
  
  Results$bias_pv= Results$PV - Results$TH0  
  Results$abs.bias_pv = abs(Results$PV - Results$TH0)
  Results$RMSE = Results$Mvar_PV + Results$bias_pv^2
  
  # Results_bycountry = NULL
  # Results_bycountry <- data.frame(rep = rep(j, 9),
  #                       DIFmag = rep(m, 9),
  #                       DIFloc = rep("Stage 2", 9),
  #                       DIFper = rep(l, 9),
  #                       Prob = rep(k, 9),
  #                       th0pvsd_2,
  #                       se1$beta)
  # 
  # Results_bycountry$bias_pv = Results_bycountry$PV - Results_bycountry$TH0
  # Results_bycountry$abs.bias_pv = abs(Results_bycountry$PV - Results_bycountry$TH0)
  # Results_bycountry$RMSE = Results_bycountry$Mvar_PV + Results_bycountry$bias_pv^2
  
  
  
  
  #Save what you will use for your results
  path4 <- "D:/DIF in MST Routing in ILSA/DIF in MST Analysis/Results Latent Regression/Person Parameter Results/"
  saveObject(Results, file = paste0(path4, "Country and mstpath/" ,"Stg2", m, "_", l, "_", k, "_p", j,".Rbin"))
  # saveObject(Results_bycountry, file = paste0(path4, "Country/" ,"Stg2",  k, "_p", j,".Rbin"))
  
  
  
}



##=======================================================================##
## Run the function                                                 ##
##=======================================================================##

# Call true item parameters for CORE magnitude 1
itempar = read.csv("D:/DIF in MST Routing in ILSA/Simulation MST/Data/Test assemble/stg2.test.1.csv",
                   header = T)
itembank = read.csv("D:/DIF in MST Routing in ILSA/Simulation MST/Data/item bank for modules.csv", 
                    header = T)
Slope = itembank[, 4]
Location = itempar[, 3] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

# Call simulated results
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.11/stg2.1.11.p1/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)


# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-2)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())


parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = 1, l = 11, m = 1 ))

# stop the core division
stopCluster(cl)

#=====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.11/stg2.1.11.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())


parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = .7, l = 11, m = 1 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.11/stg2.1.11.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())


parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = .5, l = 11, m = 1 ))

# stop the core division
stopCluster(cl)

#=================
#=================
it = NULL
stg2 = NULL
Location = itempar[, 4] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.33/stg2.1.33.p1/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = 1, l = 33, m = 1 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.33/stg2.1.33.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = .7, l = 33, m = 1 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.33/stg2.1.33.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())



parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = .5, l = 33, m = 1 ))
# stop the core division
stopCluster(cl)

#=================
#=================
it = NULL
stg2 = NULL
Location = itempar[, 5] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.56/stg2.1.56.p1/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = 1, l = 56, m = 1 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.56/stg2.1.56.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)


# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())


parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = .7, l = 56, m = 1 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.1/stg2.1.56/stg2.1.56.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)
# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .5, l = 56, m = 1 ))

# stop the core division
stopCluster(cl)

#=======================================================================

# Call true item parameters for CORE magnitude 1
itempar = read.csv("D:/DIF in MST Routing in ILSA/Simulation MST/Data/Test assemble/stg2.test.50.csv",
                   header = T)
itembank = read.csv("D:/DIF in MST Routing in ILSA/Simulation MST/Data/item bank for modules.csv", 
                    header = T)
Slope = itembank[, 4]
Location = itempar[, 3] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.11/stg2.50.11.p1/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = 1, l = 11, m = .50 ))
# stop the core division
stopCluster(cl)

#====
stg2 = NULL

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.11/stg2.50.11.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .7, l = 11, m = .50 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.11/stg2.50.11.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .5, l = 11, m = .50 ))

# stop the core division
stopCluster(cl)

#=============================
stg2 = NULL
it = NULL
Location = itempar[, 4] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.33/stg2.50.33.p1/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = 1, l = 33, m = .50 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.33/stg2.50.33.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .7, l = 33, m = .50 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.33/stg2.50.33.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .5, l = 33, m = .50 ))

# stop the core division
stopCluster(cl)

#=============================
stg2 = NULL
it = NULL
Location = itempar[, 5] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.56/stg2.50.56.p1/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = 1, l = 56, m = .50 ))
# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.56/stg2.50.56.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)
# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, 1:100, function(j) latent_reg(stg2, j, k = .7, l = 56, m = .50 ))
# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.50/stg2.50.56/stg2.50.56.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .5, l = 56, m = .50 ))
# stop the core division
stopCluster(cl)


#=====================================================================================================
#=====================================================================================================

# Call true item parameters for CORE magnitude 1
itempar = read.csv("D:/DIF in MST Routing in ILSA/Simulation MST/Data/Test assemble/stg2.test.25.csv",
                   header = T)
itembank = read.csv("D:/DIF in MST Routing in ILSA/Simulation MST/Data/item bank for modules.csv", 
                    header = T)
Slope = itembank[, 4]
Location = itempar[, 3] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.11/stg2.25.11.p1/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = 1, l = 11, m = .25 ))
# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.11/stg2.25.11.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .7, l = 11, m = .25 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.11/stg2.25.11.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())


parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .5, l = 11, m = .25 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL
it = NULL
Location = itempar[, 4] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.33/stg2.25.33.p1/")


files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())


parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = 1, l = 33, m = .25 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.33/stg2.25.33.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())


parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .7, l = 33, m = .25 ))

# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.33/stg2.25.33.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .5, l = 33, m = .25 ))
# stop the core division
stopCluster(cl)

#====
stg2 = NULL
it = NULL
Location = itempar[, 5] # this is for selecting the NODIF item parameters
it = data.frame(a = Slope, b = Location)

setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.56/stg2.25.56.p1/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = 1, l = 56, m = .25 ))
# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.56/stg2.25.56.p2/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())

parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .7, l = 56, m = .25 ))
# stop the core division
stopCluster(cl)

#====
stg2 = NULL
setwd("D:/DIF in MST Routing in ILSA/DIF in MST Results/stg2/stg2.25/stg2.25.56/stg2.25.56.p3/")

files <- dir(pattern = "*.Rbin")
files

stg2 <- files %>%
  map(loadObject) %>%
  reduce(rbind)

# Set parallel 
no_cores <- detectCores()        #number of cores that will be used
cl <- makeCluster(no_cores-4)    #creates copies of R to run on cores (i.e., computational clusters)
registerDoParallel(cl)           #register parallel backend 

clusterExport(cl,list('latent_reg', 'stg2', 'it', 'background_10', 'R.utils',
                      'saveObject'), envir=environment())


parSapply(cl, seq_along(1:100), function(j) latent_reg(stg2, j, k = .5, l = 56, m = .25 ))

# stop the core division
stopCluster(cl)



