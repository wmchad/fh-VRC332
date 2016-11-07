library(VRC332); library(ggplot2)

data("FcCoreAssay")
data("ICSStatsPreChallenge")
data("PreChallengeCOMPASS8Color")
data("PreChallengeCOMPASS18Color")
data("PostChallengeCOMPASS8Color")

Animals = FcCoreAssay[!duplicated(FcCoreAssay$AnimalID),
                      c("AnimalID", "GroupNm", "PeakVL", "SetpointVL",
                        "Challenges to Infection")]
names(Animals)[5] = "NoChallenges"
Animals$LogPeakVL = log10(as.numeric(paste(Animals$PeakVL)))
Animals$LogSetpointVL = log10(Animals$SetpointVL)

Animals = Animals[order(Animals$NoChallenges, decreasing=TRUE),]

anBasePlot = ggplot(Animals, aes(LogPeakVL, LogSetpointVL))

anBasePlot + geom_point(aes(color=GroupNm, size=NoChallenges))
anBasePlot + geom_point() + stat_smooth(method="lm", se=FALSE)
anBasePlot + geom_point(aes(color=GroupNm, size=NoChallenges)) +
  stat_smooth(mapping=aes(color=GroupNm), method="lm", se=FALSE)

## Change AnimalIDs to match other tables

Animals$AnimalID2 = gsub("-", "", Animals$AnimalID)
FcCoreAssay$AnimalID2 = gsub("-", "", FcCoreAssay$AnimalID)
ICSStatsPreChallenge$AnimalID2 = gsub("-", "", ICSStatsPreChallenge$AnimalID)
PostChallengeCOMPASS8Color$AnimalID2 = gsub("-", "", PostChallengeCOMPASS8Color$AnimalID)
PreChallengeCOMPASS8Color$AnimalID2 = gsub("-", "", PreChallengeCOMPASS8Color$AnimalID)
PreChallengeCOMPASS18Color$AnimalID2 = gsub("-", "", PreChallengeCOMPASS18Color$AnimalID)

animalIds = unique(Animals$AnimalID2)
## 87
animalIds.ics = unique(ICSStatsPreChallenge$AnimalID2)
## 100
animalIds.post = unique(PostChallengeCOMPASS8Color$AnimalID2)
## 92
animalIds.pre8 = unique(PreChallengeCOMPASS8Color$AnimalID2)
## 106
animalIds.pre18 = unique(PreChallengeCOMPASS18Color$AnimalID2)
## 59

sum(animalIds %in% animalIds.post)
## 82
animalIds[!(animalIds %in% animalIds.post)]
## "08D0130" "08D220e" "08D140"  "08D227"  "08D235" 
animalIds.post[!(animalIds.post %in% animalIds)]
## "08D130" "08D220" "8116"   "AV45X"  "DBMK"   "DBNF"   "DBNL"   "DBPJ"   "DBPP"   "ZE24"
animalIds[!(animalIds %in% animalIds.pre8)]
## None
animalIds.pre8[!(animalIds.pre8 %in% animalIds)]
## "DBMK"       "DBNF"       "DBNL"       "DBPJ"       "DBPP"       "T4276_blip" "08D225"     "8116"      
## "82"         "820"        "870"        "AS42"       "AS70"       "AV45X"      "T476"       "ZE24"      
## "ZG05"       "ZG12"       "ZG68"   
length(animalIds[animalIds %in% animalIds.pre18])
## 53

## Look at FcCore readout values

FcCore.readout = FcCoreAssay[!is.na(FcCoreAssay$readout),]
ggplot(FcCore.readout) + geom_density(aes(readout))
ggplot(FcCore.readout) + geom_density(aes(log(readout)))
max(FcCore.readout$readout)
## 217786
ggplot(FcCore.readout) + geom_density(aes(readout, color=factor(Timepoint)))
sum(FcCore.readout$readout > 50000)
## 9144
sum(FcCore.readout$readout > 100000)
## 3400
sum(FcCore.readout$readout > 200000)
## 27

## Box-Cox Transformation
bc = boxcox(lm(FcCore.readout$readout ~ 1))
lambda = bc$x[which.max(bc$y)]
## -0.2626263
FcCore.readout$bcReadout = (FcCore.readout$readout^lambda-1)/lambda
ggplot(FcCore.readout) + geom_density(aes(bcReadout))

FcCore.readout2 = FcCore.readout[FcCore.readout$readout < 25000,]
ggplot(FcCore.readout2) + geom_density(aes(readout, color=factor(Timepoint)))
ggplot(FcCore.readout) + geom_density(aes(log(readout), color=factor(Timepoint)))
table(FcCore.readout$Timepoint)

FcCore.readout.tp8 = FcCore.readout[FcCore.readout$Timepoint == 8,]
FcCore.readout.tp8$NoChallenges = FcCore.readout.tp8$`Challenges to Infection`
FcCore.readout.tp8$TF = FcCore.readout.tp8$`T/F number`
FcCore.readout.tp8$AK = FcCore.readout.tp8$`A|K Infection Type`
## 14964 observations

ggplot(FcCore.readout.tp8) + geom_density(aes(log(readout), color=GroupNm))
ggplot(FcCore.readout.tp8) + geom_density(aes(log(readout), color=factor(NoChallenges)))
ggplot(FcCore.readout.tp8) + geom_density(aes(log(readout), color=factor(Antigen)))
ggplot(FcCore.readout.tp8) + geom_density(aes(log(readout), color=factor(DetectionReagent)))
ggplot(FcCore.readout.tp8) + geom_density(aes(log(readout), color=factor(TF)))
ggplot(FcCore.readout.tp8) + geom_density(aes(log(readout), color=factor(AK)))

FcCore.readout.tp8.siv = FcCore.readout.tp8[substr(FcCore.readout.tp8$Antigen, 0, 3) == "SIV",]
ggplot(FcCore.readout.tp8.siv) + geom_density(aes(log(readout), color=factor(Antigen)))
ggplot(FcCore.readout.tp8.siv) + geom_density(aes(log(readout), color=factor(DetectionReagent)))

table(FcCore.readout.tp8.siv$Antigen[FcCore.readout.tp8.siv$DetectionReagent=="MBL"])
## SIVmac239.gp140.AVI.His 
## 87 
table(FcCore.readout.tp8.siv$Antigen[FcCore.readout.tp8.siv$DetectionReagent=="R2A.2"])
## SIVmac239.gp140 
## 87

