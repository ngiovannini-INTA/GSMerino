library(AlphaSimR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(gdata)

out.i <- data.frame()
cor.i <- vector()

system.time(
  for (rep in 1:50) {
    cat('Repeticion:',rep,'\n')

    # Simulate founder genomes ----------------------------------------------------------------------
    #founderGenomes = quickHaplo(nInd=3000, nChr=26, segSites=3000)
    system.time(
      founderGenomes <- runMacs2(nInd = 5000,                                  #Estimate that approximately half are female
                                 nChr = 26,                                    #Sciencie 26
                                 segSites = 3000,                              #(QTLs + SNPs) per Chr
                                 Ne = 120,                                     #vozzi et al, 120
                                 bp = 100000000,                               #genome size divided by the number of chromosomes, https://www.sheephapmap.org/
                                 genLen=1,                                     #genome size in morgans
                                 mutRate = 2.5e-08,                            #Check
                                 #histNe = c(500, 1500, 6000, 12000, 1e+05),   #default
                                 #histGen = c(100, 1000, 10000, 1e+05, 1e+06), #default number of generations ago for effective population sizes given in histNe
                                 histNe = c(140, 240, 500, 1000, 5000),
                                 histGen = c(10, 100, 500, 1000, 2000),
                                 inbred = FALSE,                               #should founder individuals be inbred?
                                 split = NULL,                                 #an optional historic population split in terms of generations ago
                                 ploidy = 2L, 
                                 returnCommand = FALSE, 
                                 nThreads = NULL)
    )
    
    #Global simulation parameters
    SP = SimParam$new(founderGenomes)
    
    SP$setSexes('yes_sys')
    
    #TRAITS ORDER: DES PCD PVL PDF PCA
    #DES: Survival to weaning
    #PCD: Weaning weight
    #PVL: Clean wool weight
    #PDF: Average fiber diameter
    #PCA: Adult weight
    
    phenoVar = c(0.14,11.90,0.09,2.07,20.25)
    
    heritability = c(0.04,0.2,0.3,0.45,0.35) 
    
    genVar = phenoVar*heritability
    
    genCor = matrix(c( 1.00, 0.05, 0.00, 0.01, 0.09,
                       0.05, 1.00, 0.15, 0.00, 0.65,
                       0.00, 0.15, 1.00, 0.25, 0.25,
                       0.01, 0.00, 0.25, 1.00, 0.10,
                       0.09, 0.65, 0.25, 0.10, 1.00), ncol=5)
    
    SP$addTraitA(nQtlPerChr = 1000, mean=c(0.7,23,2,18,45), var=genVar, corA = genCor)
    
    SP$addSnpChip(nSnpPerChr = 2000) #52k (2000)    
    
    #Base Population of founders
    founders = newPop(founderGenomes)
    
    
    #Phenotype the founders
    founders = setPheno(pop=founders, h2=heritability)
    
    umb <- mean(founders@pheno[,1]) - sd(founders@pheno[,1])*1.25 #threshold value 1ds = 16% of area
    
    founders@pheno[,1] <- ifelse(founders@pheno[,1] < umb,0,1)
    founders@pheno[,2] <- ifelse(founders@pheno[,1]==0,NA,founders@pheno[,2])
    founders@pheno[,3] <- ifelse(founders@pheno[,1]==0,NA,founders@pheno[,3])
    founders@pheno[,4] <- ifelse(founders@pheno[,1]==0,NA,founders@pheno[,4])
    founders@pheno[,5] <- ifelse(founders@pheno[,1]==0,NA,founders@pheno[,5])
    
    #Random reproduction cycles
    #Assign year of birth for the founders
    year = 0
    
    founders@misc$YearOfBirth <- year
    
    # Generate Initial parent populations ----------------------------------------------------------------
    #Size of the breeding nucleus (without lambs)
    nFemales = 1000
    nMales = 50
    
    #Males
    sobM = 0.97 #male survival
    
    repM = round(nMales*.35) #male replacement
    
    nMales0 = nFemales*0.5 #BT, male lambs. I consider mortality up to weaning within the dynamics
    nMales1 = round(repM) #2T
    nMales2 = round(repM * sobM) #4T
    nMales3 = round(repM * sobM^2) #6T
    
    sum(nMales1,nMales2,nMales3)
    c(nMales1,nMales2,nMales3)
    nMales0
    
    males = selectInd(pop = founders, nInd = sum(nMales0,nMales1,nMales2,nMales3), use='rand', sex='M')
    
    #Sires3
    (start = 1)
    (end = nMales3)
    
    sires3 = males[start:end] #6T
    sires3@misc$YearOfBirth <- rep(-3,nMales3)
    
    #sires3 = setMisc(x=sires3,
    
    #Sires2
    (start = end + 1)
    (end = start - 1 + nMales2)
    
    sires2 = males[start:end] #4T
    sires2@misc$YearOfBirth <- rep(-2,nMales2)
    
    #Sires1
    (start = end + 1)
    (end = start - 1 + nMales1)
    
    sires1 = males[start:end] #2T
    sires1@misc$YearOfBirth <- rep(-1,nMales1)
    
    #Sires0
    (start = end + 1)
    (end = start - 1 + nMales0)
    
    sires0 = males[start:end] #DL
    sires0@misc$YearOfBirth <- rep(0,nMales0)

    
    #sires0 only have DES data, they are phenotyped in the loop
    sires0@pheno[,2] <- NA
    sires0@pheno[,3] <- NA
    sires0@pheno[,4] <- NA
    sires0@pheno[,5] <- NA
    
    sires = c(sires1,sires2,sires3)
    #table(unlist(getMisc(x=sires, node='YearOfBirth')))
    table(sires@misc$YearOfBirth)
    
    #Females
    sobH = 0.97 #ewe survival
    repH = round(nFemales*.21) #ewe replacement
    
    nDams0 = nFemales*0.5 #BT, femele lambs. I consider mortality up to weaning within the dynamics
    nDams1 = round(repH) #2T
    nDams2 = round(repH * sobH) #4T
    nDams3 = round(repH * sobH^2) #6T
    nDams4 = round(repH * sobH^3) #8T
    nDams5 = round(repH * sobH^4) #HT
    
    sum(nDams1,nDams2,nDams3,nDams4,nDams5)
    c(nDams1,nDams2,nDams3,nDams4,nDams5)
    nDams0
    
    females = selectInd(pop = founders, nInd = sum(nDams0,nDams1,nDams2,nDams3,nDams4,nDams5), use='rand', sex='F')
    
    #Dams5
    (start = 1)
    (end = nDams5)
    
    dams5 = females[start:end]
    dams5@misc$YearOfBirth <- rep(-5,nDams5)
    nInd(dams5)
    
    #Dams4
    (start = end + 1)
    (end = start - 1 + nDams4)
    
    dams4 = females[start:end]
    dams4@misc$YearOfBirth <- rep(-4,nDams4)

    nInd(dams4)
    
    #Dams3
    (start = end + 1)
    (end = start - 1 + nDams3)
    
    dams3 = females[start:end]
    dams3@misc$YearOfBirth <- rep(-3,nDams3)

    nInd(dams3)
    
    #Dams2
    (start = end + 1)
    (end = start - 1 + nDams2)
    
    dams2 = females[start:end]
    dams2@misc$YearOfBirth <- rep(-2,nDams2)

    nInd(dams2)
    
    #Dams1
    (start = end + 1)
    (end = start - 1 + nDams1)
    
    dams1 = females[start:end]
    dams1@misc$YearOfBirth <- rep(-1,nDams1)

    nInd(dams1)
    
    #Dams0
    (start = end + 1)
    (end = start - 1 + nDams0)
    
    dams0 = females[start:end]
    dams0@misc$YearOfBirth <- rep(0,nDams0)
    
    #Dams0 only have DES data, they are phenotyped in the loop
    dams0@pheno[,2] <- NA
    dams0@pheno[,3] <- NA
    dams0@pheno[,4] <- NA
    dams0@pheno[,5] <- NA
    
    dams = c(dams5,dams4,dams3,dams2,dams1)
    nInd(dams)
    
    #table(unlist(getMisc(x = dams, node = 'YearOfBirth')))
    table(dams@misc$YearOfBirth)
    
    # Data recording ----------------------------------------------------------
    #function to record and colect data
    recordData <- function(data=NULL, pop, YearOfUse=NA){
      popData = data.frame(id = pop@id,
                           father = pop@father,
                           mother = pop@mother,
                           sex = pop@sex,
                           gvDES = pop@gv[,'Trait1'],
                           gvPCD = pop@gv[,'Trait2'],
                           gvPVL = pop@gv[,'Trait3'],
                           gvPDF = pop@gv[,'Trait4'],
                           gvPCA = pop@gv[,'Trait5'],
                           pDES = pop@pheno[,'Trait1'],
                           pPCD = pop@pheno[,'Trait2'],
                           pPVL = pop@pheno[,'Trait3'],
                           pPDF = pop@pheno[,'Trait4'],
                           pPCA = pop@pheno[,'Trait5'],
                           #YearOfBirth = unlist(getMisc(x = pop, node = 'YearOfBirth')),
                           YearOfBirth = pop@misc$YearOfBirth,
                           YearOfUse = as.numeric(YearOfUse),
                           geno = pullSnpGeno(pop)
      )
      
      #manage first instance of calling this function, when data is NULL
      if (is.null(data)){
        ret = popData
      } else { ret = rbind(data,popData)
      }
      return(ret)
      
    }
    
    data4AllAnimals = recordData(pop=founders)
    data4NewParents = recordData(pop=c(sires1,dams1))
    
    # Simulation of years of dynamics prior to the scenarios ------------------
    
    #I define economic values of selection index
    econWt = c(927.4,29.3,531.1,-28.3,2.9) #(a) de "Actualización de objetivos... .pdf"
    
    varP = matrix(c(0.140,	0.232,	-0.013,	0.000,	0.556,
                    0.232,	11.903,	0.311,	0.000,	6.986,
                    -0.013,	0.311,	0.090,	0.108,	0.405,
                    0.000,	0.000,	0.108,	2.074,	0.648,
                    0.556,	6.986,	0.405,	0.648,	20.250), ncol=5)
    
    
    varG = matrix(c(0.006,	0.006,	0.000,	0.001,	0.018,
                    0.006,	2.381,	0.040,	0.000,	0.131,
                    0.000,	0.040,	0.031,	0.042,	0.116,
                    0.001,	0.000,	0.042,	0.933,	0.257,
                    0.018,	0.131,	0.116,	0.257,	7.088), ncol=5)
    
    b <- as.vector(smithHazel(econWt,varG=varG ,varP = varP))
    
    #Dynamics
    for (year in 1:10) {
      cat('Working on the year:',year,'\n')
      
      #generate progeny from current sires and dams
      candidates = randCross2(males = sires, females = dams, nCrosses = nInd(dams))
      candidates@misc$YearOfBirth = rep(year, nInd(dams)) 
      #candidates = setMisc(x = candidates, node = 'YearOfBirth', value = year)
      #candidates = attrition(candidates, p=0.25) #This is where the animals that do not reach phenotyping are filtered out
      candidates = setPheno(candidates, h2 = 0.04, traits = 1) #Phenotype DES only lambs
      candidates@pheno[,1] <- ifelse(candidates@pheno[,1] < umb,0,1)
      
      #The animals are phenotyped at 1 year of age
      sires0 = setPheno(sires0, h2 = c(0.20, 0.30, 0.45, 0.35), traits = c(2,3,4,5))
      dams0 = setPheno(dams0, h2 = c(0.20, 0.30, 0.45, 0.35), traits = c(2,3,4,5))
      
      #record data for all newborn animals
      data4AllAnimals = recordData(data = data4AllAnimals,
                                   pop = c(dams0,sires0))
      
      #Update and select sires
      sires3 = selectInd(pop = sires2, nInd = nMales3, use = 'rand')
      sires2 = selectInd(pop = sires1, nInd = nMales2, use = 'rand')
      sires1 = selectInd(pop = sires0, nInd = nMales1, use = 'pheno', trait=selIndex, b=b)
      sires0  = selectInd(pop = candidates, nInd = sum(isMale(candidates) & candidates@pheno[,1]==1), trait = 1, use= "pheno", sex = 'M')
      sires0b = selectInd(pop = candidates, nInd = sum(isMale(candidates) & candidates@pheno[,1]==0), trait = 1, use= "pheno", selectTop = F, sex = 'M') #I'm only including these to add to the database later for BLUP (they're for informational purposes)
      sires   = c(sires3,sires2,sires1)
      
      #Update and select dams
      dams5 = selectInd(pop = dams4, nInd = nDams5, use = 'rand')
      dams4 = selectInd(pop = dams3, nInd = nDams4, use = 'rand')
      dams3 = selectInd(pop = dams2, nInd = nDams3, use = 'rand')
      dams2 = selectInd(pop = dams1, nInd = nDams2, use = 'rand')
      dams1 = selectInd(pop = dams0, nInd = nDams1, use = 'rand')
      dams0  = selectInd(pop = candidates, nInd = sum(isFemale(candidates) & candidates@pheno[,1]==1), trait = 1, use = 'pheno', sex = 'F')
      dams0b = selectInd(pop = candidates, nInd = sum(isFemale(candidates) & candidates@pheno[,1]==0), trait = 1, use = 'pheno', selectTop = F, sex = 'F')#I'm only including these to add to the database later for BLUP (they're for informational purposes)
      dams  = c(dams5,dams4,dams3,dams2,dams1)
      
      #record data for dead ones
      data4AllAnimals = recordData(data = data4AllAnimals,
                                   pop = c(sires0b,dams0b))
      
      #record data for de newly selected sires and dams (just the new ones)
      data4NewParents = recordData(data = data4NewParents,
                                   pop = c(sires1,dams1))
      
    }
    
# Scenario 5 GSp -------------------------------
#Presel + GS
    
    #Reference Population (animals from last 4 generations of the last stage n=2029)
    SNPref <- data4NewParents %>% dplyr::filter(YearOfBirth >= 3)
    pre.snp.ref <- apply(SNPref[,17:ncol(SNPref)], MARGIN = 1, FUN = paste, sep = '', collapse= '')
    id.ref <- SNPref$id
    
   pre.snp.i <- cbind(id.ref, pre.snp.ref)

          for (year in 11:21) {
          cat('Working on the year:',year,'\n')
          
          #generate progeny from current sires and dams
          candidates = randCross2(males = sires, females = dams, nCrosses = nInd(dams))
          candidates@misc$YearOfBirth = rep(year, nInd(dams)) 
          candidates = setPheno(candidates, h2 = 0.04, traits = 1) #DES phenotype only lambs
          candidates@pheno[,1] <- ifelse(candidates@pheno[,1] < umb,0,1)
          
          #The animals are phenotyped at 1 year of age
          sires0 = attrition(sires0, p=0.40) #phenotyped males 60% of DES=1
          sires0 = setPheno(sires0, h2 = c(0.20, 0.30, 0.45, 0.35), traits = c(2,3,4,5))
          
          dams0 = attrition(dams0, p=0.05) #phenotyped females 95% of DES=1 
          dams0 = setPheno(dams0, h2 = c(0.20, 0.30, 0.45, 0.35), traits = c(2,3,4,5))
          
          
          #record data for all newborn animals
          data4AllAnimals = recordData(data = data4AllAnimals,
                                       pop = c(dams0,sires0))
          
          ped <- data.frame(id = data4AllAnimals$id, 
                            father = data4AllAnimals$father, 
                            mother = data4AllAnimals$mother)
          
          dat <- data.frame(id = data4AllAnimals$id, 
                            sex = data4AllAnimals$sex, 
                            AN = data4AllAnimals$YearOfBirth, 
                            DES = data4AllAnimals$pDES, 
                            PCD = data4AllAnimals$pPCD, 
                            PVL = data4AllAnimals$pPVL,
                            PDF = data4AllAnimals$pPDF,
                            PCA = data4AllAnimals$pPCA)
          
          dat <- subset(dat, AN >= 3 ) #define a partir de que anio hay registros 
          
                  
          #Linux server
          #write.table(ped, file = "/home/corva/NicoG/doc/ped.dat", sep=" ", row.names = F, col.names = F, quote=F, na='-999')
          #write.table(dat, file = "/home/corva/NicoG/doc/dat.dat", sep=" ", row.names = F, col.names = F, quote=F, na='-999')
          #write.fwf(pre.snp.i, file = "/home/corva/NicoG/doc/snp.dat", sep="", rownames = F, colnames = F, width = c(7,52000))
          
          #Local windows
          write.table(ped, file = "d:/PROVINO/temp/e5/ped.dat", sep=" ", row.names = F, col.names = F, quote=F, na='-999')
          write.table(dat, file = "d:/PROVINO/temp/e5/dat.dat", sep=" ", row.names = F, col.names = F, quote=F, na='-999')
          write.fwf(pre.snp.i, file = "d:/PROVINO/temp/e5/snp.dat", sep="", rownames = F, colnames = F, width = c(7,52000))
          
          #run BLUP1 in local windows (preselection)
          setwd('d:/PROVINO/temp/e5/')
          shell("renumf90 renum.par")
          shell("blupf90 renf90.par")
          
          #run BLUP1 in Linux server (preselection)
          #setwd('/home/corva/NicoG/doc/')
          #system(command = "renumf90 renum.par")
          #system(command = "blupf90 renf90.par")
          
          #Imports EBV's
          VC <-read.delim("solutions", skip=1, header=F, sep='')
          names(VC)=c("trait","effect","ID.r","solution","se")
          VC <- subset(VC,VC$effect==3)
          
          ped.r <-read.delim("renadd03.ped", header=F,sep='')
          names(ped.r)=c("ID.r","IDUP.r","IDUM.r","v4","V5","V6",'V7',"n.pad","n.mad","id")
          
          VC.j <- merge(VC,ped.r, by='ID.r')
          VC.j$id <- as.character(VC.j$id)
          
          #I identify sires0 to select superior 25%
          id.Mblup1 <- data4AllAnimals %>% dplyr::filter(sex == 'M' & YearOfBirth == year-1 & pDES==1) %>% select(id)
          
          EBV.Mblup1 <- dplyr::filter(VC.j, id %in% id.Mblup1$id) %>%
            select(id,trait, solution)
          
          EBV.Mblup1 <- spread(EBV.Mblup1,"trait","solution")
          names(EBV.Mblup1)=c("id","Trait1","Trait2","Trait3","Trait4","Trait5")
          
          EBV.Mblup1 <- EBV.Mblup1 %>% dplyr::mutate(IND = Trait1*econWt[1] + 
                                                Trait2*econWt[2] +
                                                Trait3*econWt[3] +
                                                Trait4*econWt[4] +
                                                Trait5*econWt[5]) %>% slice_max(IND, n = 50)
          
          #these are 20% superior males for genotyping
          SNPcand <- dplyr::filter(data4AllAnimals, id %in% EBV.Mblup1$id)
          pre.snp.cand <- apply(SNPcand[,17:ncol(SNPcand)], MARGIN = 1, FUN = paste, sep = '', collapse= '')
          id.cand <- SNPcand$id
          snp1 <- cbind(id.cand, pre.snp.cand)
          
          pre.snp.i <- rbind(pre.snp.i,snp1)
          
          #run BLUP2 windows 
          setwd('d:/PROVINO/temp/e5/')
          shell("renumf90 renum.par")
          shell("blupf90 renf90.par")
          
          #run BLUP2 Linux
          #setwd('/home/corva/NicoG/doc/')
          #system(command = "renumf90 renum.par")
          #system(command = "blupf90 renf90.par")
          
          #Imports EBV's
          VC <-read.delim("solutions", skip=1, header=F, sep='')
          names(VC)=c("trait","effect","ID.r","solution","se")
          VC <- subset(VC,VC$effect==3)
          
          ped.r <-read.delim("renadd03.ped", header=F,sep='')
          names(ped.r)=c("ID.r","IDUP.r","IDUM.r","v4","V5","V6",'V7',"n.pad","n.mad","id")
          
          VC.j <- merge(VC,ped.r, by='ID.r')
          VC.j$id <- as.character(VC.j$id)
          
          EBV.pob <- dplyr::filter(VC.j, id %in% data4AllAnimals$id) %>%
            select(id,trait, solution)
          
          EBV.pob <- spread(EBV.pob,"trait","solution")
          names(EBV.pob)=c("id","Trait1","Trait2","Trait3","Trait4","Trait5")
          
          ev.pob <- merge(data4AllAnimals,EBV.pob,by = 'id')
          
          cor <- c(rep,
                   year,
                   cor(ev.pob$gvDES,ev.pob$Trait1),
                   cor(ev.pob$gvPCD,ev.pob$Trait2),
                   cor(ev.pob$gvPVL,ev.pob$Trait3),
                   cor(ev.pob$gvPDF,ev.pob$Trait4),
                   cor(ev.pob$gvPCA,ev.pob$Trait5))
          
          cor.i <- rbind(cor.i,cor)
          
          
          #Extract EBVs from PedBLUP output (sires0)
          EBV <-dplyr::filter(VC.j, id %in% sires0@id) %>%
            select(id,trait, solution)
          
          EBV <- spread(EBV,"trait","solution")
          
          EBV <- arrange(EBV,as.numeric(id))
          names(EBV)=c("id","Trait1","Trait2","Trait3","Trait4","Trait5")
          
          # Insert externally created EBV into the population
          sires0 <- sires0[order(sires0@id),]
          sires0@ebv <- as.matrix(EBV[,c(2:6)]) #Assign EBVs to population 
          
          #Extract EBVs from PedBLUP output (dams0)
          EBV <- dplyr::filter(VC.j, id %in% dams0@id) %>%
            select(id,trait, solution)
          
          EBV <- spread(EBV,"trait","solution")
          
          EBV <- arrange(EBV,as.numeric(id))
          names(EBV)=c("id","Trait1","Trait2","Trait3","Trait4","Trait5")
          
          # Insert externally created EBV into the population
          dams0 <- dams0[order(dams0@id),]
          dams0@ebv <- as.matrix(EBV[,c(2:6)]) #Assign EBVs to population  
          
          ##Aging and selection of males
          sires3 = selectInd(pop = sires2, nInd = nMales3, use = 'rand') #4T to 6T
          sires2 = selectInd(pop = sires1, nInd = nMales2, use = 'rand') #2T to 4T
          sires1 = selectInd(pop = sires0, nInd = nMales1, use= "ebv", sex = 'M', trait = selIndex, b=econWt) #BT to 2D
          sires0  = selectInd(pop = candidates, nInd = sum(isMale(candidates)  & candidates@pheno[,1]==1), trait=1, sex = 'M', use = 'pheno')
          sires0b = selectInd(pop = candidates, nInd = sum(isMale(candidates)  & candidates@pheno[,1]==0), trait=1, sex = 'M', selectTop = F, use = 'pheno') 
          sires = c(sires3,sires2,sires1)
          
          #Aging and selection of females
          dams5 = selectInd(pop = dams4, nInd = nDams5, use = 'rand')
          dams4 = selectInd(pop = dams3, nInd = nDams4, use = 'rand')
          dams3 = selectInd(pop = dams2, nInd = nDams3, use = 'rand')
          dams2 = selectInd(pop = dams1, nInd = nDams2, use = 'rand')
          dams1 = selectInd(pop = dams0, nInd = nDams1, use = 'ebv', sex = 'F', trait=selIndex, b=econWt) #BT to 2T, best are selected and pass to be breeding females on next generation
          dams0  = selectInd(pop = candidates, nInd = sum(isFemale(candidates) & candidates@pheno[,1]==1), trait=1, sex ='F', use = 'pheno')
          dams0b = selectInd(pop = candidates, nInd = sum(isFemale(candidates) & candidates@pheno[,1]==0), trait=1, sex ='F', selectTop = F, use = 'pheno')
          
          dams = c(dams5,dams4,dams3,dams2,dams1)
          
          #record data for dead ones
          data4AllAnimals = recordData(data = data4AllAnimals,
                                       pop = c(sires0b,dams0b))
          
          #record data for de newly selected sires and dams (just the new ones)
          data4NewParents = recordData(data = data4NewParents,
                                       pop = c(sires1,dams1))
          
        }
    
    out <- data4AllAnimals %>%
      dplyr::select(1:16) %>%
      #group_by(YearOfBirth) %>%
      #summarise(mean = mean(gvPCD)) %>%
      dplyr::mutate(rep=rep)
    
    #save inbreeding file - linux server
    #file.rename('/home/corva/NicoG/doc/renf90.inb',paste('/home/corva/NicoG/doc/renf90',rep,'.inb',sep=''))
    
    #save inbreeding file - local windows
    file.rename('d:/PROVINO/temp/e5/renf90.inb',paste('d:/PROVINO/temp/e5/renf90',rep,'.inb',sep=''))
    
    out.i <- rbind(out.i,out)
    
    #temp files
    #Linux server
    #write.table(out.i,paste('/home/corva/NicoG/doc/','out_',rep,'.csv',sep=''), sep=';', row.names = F)
    #write.table(cor.i,paste('/home/corva/NicoG/doc/','cor_',rep,'.csv',sep=''), sep=';', row.names = F)
    
    #local windows
    write.table(out.i,paste('d:/PROVINO/temp/e5/','out_',rep,'.csv', sep=''),sep=';',row.names = F)
    write.table(cor.i,paste('d:/PROVINO/temp/e5/','cor_',rep,'.csv', sep=''),sep=';',row.names = F)
    
    rm(list=setdiff(ls(), c("out.i","cor.i")))
  }
)  

#Linux server
#write.table(out.i,"/home/corva/NicoG/doc/out.csv", sep=';', row.names = F)
#write.table(cor.i,"/home/corva/NicoG/doc/cor.csv", sep=';', row.names = F)

#local windows
write.table(out.i,"d:/PROVINO/temp/e5/out.csv", sep=';')
write.table(cor.i,"d:/PROVINO/temp/e5/cor.csv", sep=';')


