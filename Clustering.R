require(rJava)
options(java.parameters = "-Xmx8000m")
require(xlsx)
require(chemometrics)
require(rcdk)
require(cluster)
require(ChemmineR)
require(rgl)
require(vegan)
require(factoextra)
require(fingerprint)
require(NbClust)
require(ggplot2)
require(gridExtra)
require(fmcsR)
sdfset <- read.SDFset("ZINC_NP_Dataset.sdf")
view(sdfset)
view(sdfset[1:4])
as(sdfset[1:4], "list")
sdfid(sdfset)[1:4]
blockmatrix <- datablock2ma(datablocklist=datablock(sdfset))
numchar <- splitNumChar(blockmatrix=blockmatrix)
propma <- data.frame(MF=MF(sdfset), MW=MW(sdfset), atomcountMA(sdfset))
datablock(sdfset) <- propma

#AtomPair calculation for the dataset
apset <- sdf2ap(sdfset)
cmp.search(apset, apset[1], type=3, cutoff = 0.3, quiet=TRUE)
cmp.search(apset[1], apset[6:35858], type=3, cutoff = 0.3, quiet=TRUE)

#Fingerprint calculation for the dataset
showClass("FPset")
fpset <- desc2fp(apset)
view(fpset[1:2])
fpset[1:4]
length(fpset)
cid(fpset)
fpma <- as.matrix(fpset)
as(fpma, "FPset")
fpchar <- as.character(fpset)
as(fpchar, "FPset")

#Fingerprint-driven Tanimoto coefficient based calculation for the dataset
fpSim_comp1 <- fpSim(fpset[1], fpset, method="Tanimoto", cutoff=0.5, top=100)
fpSim_comp2 <- fpSim(fpset[2], fpset, method="Tanimoto", cutoff=0.5, top=100)
fpSim_comp3 <- fpSim(fpset[3], fpset, method="Tanimoto", cutoff=0.5, top=100)
fpSim_comp4 <- fpSim(fpset[4], fpset, method="Tanimoto", cutoff=0.5, top=100)
fpSim_comp5 <- fpSim(fpset[5], fpset, method="Tanimoto", cutoff=0.5, top=100)
params <- genParameters(fpset)
params
fpSim(fpset[[1]], fpset, top=50, parameters=params)
fpSim_comp1_param <- fpSim(fpset[[1]], fpset, cutoff=0.5, scoreType="similarity", parameters=params)
fpSim_comp2_param <- fpSim(fpset[[2]], fpset, cutoff=0.5, scoreType="similarity", parameters=params)
fpSim_comp3_param <- fpSim(fpset[[3]], fpset, cutoff=0.5, scoreType="similarity", parameters=params)
fpSim_comp4_param <- fpSim(fpset[[4]], fpset, cutoff=0.5, scoreType="similarity", parameters=params)
fpSim_comp5_param <- fpSim(fpset[[5]], fpset, cutoff=0.5, scoreType="similarity", parameters=params)

