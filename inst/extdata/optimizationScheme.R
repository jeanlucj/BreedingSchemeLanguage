initializePopulation(sEnv, nInd=parmList$nFounder) # popID 0
phenotype(sEnv, nRep=parmList$nRepC0)
if (parmList$yesGeno) genotype(sEnv)
predictValue(sEnv)
select(sEnv) # popID 1
cross(sEnv, nProgeny=parmList$nProg) # popID 2
phenotype(sEnv, nRep=parmList$nRepC1)
if (parmList$yesGeno) genotype(sEnv)
predictValue(sEnv)
select(sEnv) # popID 3
