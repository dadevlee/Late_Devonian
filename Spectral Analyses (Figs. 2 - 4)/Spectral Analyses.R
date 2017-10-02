rm(list=ls())
library(astrochron)

setwd("/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Monte Carlo R1/")

# Copy the "optimal" pointers obtained in the Monte Carlo script (Monte_Carlo_Devonian.R, look for variable "Pointers_sum_best") ----
Pointers_optimal=c(2719.81061, 2455.09067, 2313.04357, 2122.42466, 2060.64356, 1830.25260, 1570.54610, 1301.58676,  913.05088,  773.97566,  571.46070, 339.81020,  242.32673,  107.65829,   60.37882,  -20.00000, -488.48499, -653.01997)

# Read the "rough" tie-points on Figures 5 and 6
Pointers_H32=read.csv("Pointers_H32.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
colnames(Pointers_H32)=c('Height',"Relative Age older than F-F boundary (kyr)")
Pointers_H32[,2]=Pointers_optimal[c(1,2,4,6:15,17:18)] # Replace by "the best astronomical fit" Pointers. Frasnian are positive ages, Famennian negative ages.

Pointers_CG1=read.csv("Pointers_CG1.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
colnames(Pointers_H32)=c('Height',"Relative Age older than F-F boundary (kyr)")
Pointers_CG1[,2]=Pointers_optimal[c(1:13)] # Replace by "the best astronomical fit" Pointers. Frasnian are positive ages, Famennian negative ages.

Pointers_SectionC=read.csv("Pointers_SectionC.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
colnames(Pointers_H32)=c('Height',"Relative Age older than F-F boundary (kyr)")
Pointers_SectionC[,2]=Pointers_optimal[c(1:18)] # Replace by "the best astronomical fit" Pointers. Frasnian are positive ages, Famennian negative ages.

Pointers_Sinsin=read.csv("Pointers_Sinsin.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
colnames(Pointers_H32)=c('Height',"Relative Age older than F-F boundary (kyr)")
Pointers_Sinsin[,2]=Pointers_optimal[c(10:18)] # Replace by "the best astronomical fit" Pointers. Frasnian are positive ages, Famennian negative ages.

Pointers_Fuhe=read.csv("Pointers_Fuhe.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
colnames(Pointers_H32)=c('Height',"Relative Age older than F-F boundary (kyr)")
Pointers_Fuhe[,2]=Pointers_optimal[c(10:17)] # Replace by "the best astronomical fit" Pointers. Frasnian are positive ages, Famennian negative ages.

Pointers_Kowala=read.csv("Pointers_Kowala.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
colnames(Pointers_H32)=c('Height',"Relative Age older than F-F boundary (kyr)")
Pointers_Kowala[,2]=Pointers_optimal[c(10,12:14,16:17)] # Replace by "the best astronomical fit" Pointers. Frasnian are positive ages, Famennian negative ages.

# Calibrate the proxy records using the "tune" function in astrochron ----
setwd("/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Data_depth_domain/")

# H-32 ----
H32_d13C_carb=read.csv("H32_d13C_run2.csv",header=T)
H32_d13C_carb=H32_d13C_carb[,c(1,2)]
H32_d13C_carb[,1]=H32_d13C_carb[,1]/100
H32_d13C_carb_T1=tune(H32_d13C_carb,Pointers_H32,extrapolate=T) # Convert depth to time according to the "best astronomical fit" Pointers
write.csv(H32_d13C_carb_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/H32_d13C_carb_T1.csv",row.names = F)

H32_d13C_org=read.delim("d13Corg_H32.txt",header=T)
H32_d13C_org=H32_d13C_org[,c(3,2)]
H32_d13C_org_T1=tune(H32_d13C_org,Pointers_H32,extrapolate=T)
write.csv(H32_d13C_org_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/H32_d13C_org_T1.csv",row.names = F)

H32_MS=read.csv("H32_MS.csv",header=T)
H32_MS=trim(H32_MS[,c(1,2)])
H32_MS[,2]=H32_MS[,2]*10^8
H32_MS[,1]=H32_MS[,1]/100
H32_MS_T1=tune(H32_MS,Pointers_H32,extrapolate=F)
write.csv(H32_MS_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/H32_MS_T1.csv",row.names = F)

# CG-1 ----
CG1_d13C_carb=read.delim("isotopes_CG1.txt",header=T)
CG1_d13C_carb=sortNave(CG1_d13C_carb[,c(1,2)])
CG1_d13C_carb_T1=tune(CG1_d13C_carb,Pointers_CG1,extrapolate = T)
write.csv(CG1_d13C_carb_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/CG1_d13C_carb_T1.csv",row.names = F)

CG1_MS=read.delim("MS_CG1.txt",header=T)
CG1_MS=sortNave(CG1_MS[,c(1,2)])
CG1_MS_T1=tune(CG1_MS,Pointers_CG1,extrapolate = T)
write.csv(CG1_MS_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/CG1_MS_T1.csv",row.names = F)

# Section C ----
SectionC_d13C_org=read.csv("AlbertaC_d13Corg.csv",header=F)
SectionC_d13C_org=sortNave(SectionC_d13C_org[,c(1,2)])
SectionC_d13C_org_T1=tune(SectionC_d13C_org,Pointers_SectionC,extrapolate = T)
write.csv(SectionC_d13C_org_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/SectionC_d13C_org_T1.csv",row.names = F)

SectionC_MS=read.csv("AlbertaC_MS.csv",header=T)
SectionC_MS=sortNave(SectionC_MS[,c(1,2)])
SectionC_MS_T1=tune(SectionC_MS,Pointers_SectionC,extrapolate = F)
write.csv(SectionC_MS_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/SectionC_MS_T1.csv",row.names = F)

# Sinsin ----
Sinsin_d13C_carb=read.delim("Sinsin_Isotopes.txt",header=T)
Sinsin_d13C_carb=sortNave(Sinsin_d13C_carb[,c(1,2)])
Sinsin_d13C_carb[,1]=Sinsin_d13C_carb[,1]/100
Sinsin_d13C_carb_T1=tune(Sinsin_d13C_carb,Pointers_Sinsin,extrapolate = T)
write.csv(Sinsin_d13C_carb_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Sinsin_d13C_carb_T1.csv",row.names = F)

Sinsin_MS=read.delim("Sinsin_MS.txt",header=T)
Sinsin_MS=sortNave(Sinsin_MS[,c(1,2)])
Sinsin_MS[,2]=Sinsin_MS[,2]*10^8
Sinsin_MS[,1]=Sinsin_MS[,1]/100
Sinsin_MS_T1=tune(Sinsin_MS,Pointers_Sinsin,extrapolate = T)
write.csv(Sinsin_MS_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Sinsin_MS_T1.csv",row.names = F)

# Fuhe ----
Fuhe_d13C_carb_highres=read.csv("Fuhe_d13C_carb_highres.csv",header=F)
Fuhe_d13C_carb_highres=sortNave(Fuhe_d13C_carb_highres[,c(1,2)])
Fuhe_d13C_carb_highres_T1=tune(Fuhe_d13C_carb_highres,Pointers_Fuhe,extrapolate = T)
write.csv(Fuhe_d13C_carb_highres_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Fuhe_d13C_carb_highres_T1.csv",row.names = F)

Fuhe_d13C_carb_lowres=read.csv("Fuhe_d13C_carb_lowres.csv",header=F)
Fuhe_d13C_carb_lowres=sortNave(Fuhe_d13C_carb_lowres[,c(1,2)])
Fuhe_d13C_carb_lowres_T1=tune(Fuhe_d13C_carb_lowres,Pointers_Fuhe,extrapolate = T)
write.csv(Fuhe_d13C_carb_lowres_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Fuhe_d13C_carb_lowres_T1.csv",row.names = F)

Fuhe_d13C_org_highres=read.csv("Fuhe_d13C_org_highres.csv",header=F)
Fuhe_d13C_org_highres=sortNave(Fuhe_d13C_org_highres[,c(1,2)])
Fuhe_d13C_org_highres_T1=tune(Fuhe_d13C_org_highres,Pointers_Fuhe,extrapolate = T)
write.csv(Fuhe_d13C_org_highres_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Fuhe_d13C_org_highres_T1.csv",row.names = F)

Fuhe_d13C_org_lowres=read.csv("Fuhe_d13C_org_lowres.csv",header=F)
Fuhe_d13C_org_lowres=sortNave(Fuhe_d13C_org_lowres[,c(1,2)])
Fuhe_d13C_org_lowres_T1=tune(Fuhe_d13C_org_lowres,Pointers_Fuhe,extrapolate = T)
write.csv(Fuhe_d13C_org_lowres_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Fuhe_d13C_org_lowres_T1.csv",row.names = F)

Fuhe_MS_lower=read.csv("Data_Zenghui.csv",header = F)
Fuhe_MS_lower=sortNave(Fuhe_MS_lower[,c(1,2)])
Fuhe_MS_lower_T1=tune(Fuhe_MS_lower,Pointers_Fuhe,extrapolate = T)
write.csv(Fuhe_MS_lower_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Fuhe_MS_lower_T1.csv",row.names = F)

Fuhe_MS_upper=read.csv("MagneticSusceptibility_upper.csv",header = F)
Fuhe_MS_upper=sortNave(Fuhe_MS_upper[,c(1,2)])
Fuhe_MS_upper[,2]=Fuhe_MS_upper[,2]*10^8
Fuhe_MS_upper_T1=tune(Fuhe_MS_upper,Pointers_Fuhe,extrapolate = T)
write.csv(Fuhe_MS_upper_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Fuhe_MS_upper_T1.csv",row.names = F)

# Kowala ----
Kowala_d13C_carb_highres=read.csv("Kowala_d13C_carb_highres.csv",header=F)
Kowala_d13C_carb_highres=sortNave(Kowala_d13C_carb_highres[,c(1,2)])
Kowala_d13C_carb_highres_T1=tune(Kowala_d13C_carb_highres,Pointers_Kowala,extrapolate = T)
write.csv(Kowala_d13C_carb_highres_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Kowala_d13C_carb_highres_T1.csv",row.names = F)

Kowala_d13C_carb_lowres=read.csv("Kowala_d13C_carb_lowres.csv",header=F)
Kowala_d13C_carb_lowres=sortNave(Kowala_d13C_carb_lowres[,c(1,2)])
Kowala_d13C_carb_lowres_T1=tune(Kowala_d13C_carb_lowres,Pointers_Kowala,extrapolate = T)
write.csv(Kowala_d13C_carb_lowres_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Kowala_d13C_carb_lowres_T1.csv",row.names = F)

Kowala_d13C_org=read.csv("Kowala_d13C_org_lowres.csv",header=F)
Kowala_d13C_org=sortNave(Kowala_d13C_org[,c(1,2)])
Kowala_d13C_org_T1=tune(Kowala_d13C_org,Pointers_Kowala,extrapolate = T)
write.csv(Kowala_d13C_org_T1,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Kowala_d13C_org_T1.csv",row.names = F)

############ MTM ######### ----
# 
setwd('/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/')
H32_MS_T1_i=linterp(H32_MS_T1,5)
H32_MS_T1_i=detrend(H32_MS_T1_i)
H32_MS_T1_MTM=mtm(H32_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1)
H32_MS_T1_eha=eha(H32_MS_T1_i,tbw=2,demean = F, detrend = T,fmax = 0.05, genplot=3)
write.csv(H32_MS_T1_MTM,"H32_MS_T1_MTM.csv",row.names = F)

CG1_MS_T1_i=linterp(CG1_MS_T1,5)
CG1_MS_T1_i=detrend(CG1_MS_T1_i)
CG1_MS_T1_MTM=mtm(CG1_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1)
CG1_MS_T1_eha=eha(CG1_MS_T1_i,tbw=2,demean = F, detrend = T,fmax = 0.05, genplot=3)
write.csv(CG1_MS_T1_MTM,"CG1_MS_T1_MTM.csv",row.names = F)

Sinsin_MS_T1_i=linterp(Sinsin_MS_T1,5)
Sinsin_MS_T1_i=detrend(Sinsin_MS_T1_i)
Sinsin_MS_T1_MTM=mtm(Sinsin_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1)
Sinsin_MS_T1_eha=eha(Sinsin_MS_T1_i,tbw=2,demean = F, detrend = T,fmax = 0.05, genplot=3)
write.csv(Sinsin_MS_T1_MTM,"Sinsin_MS_T1_MTM.csv",row.names = F)

Sinsin_d13C_carb_T1_i=linterp(Sinsin_d13C_carb_T1,7)
Sinsin_d13C_carb_T1_i=detrend(Sinsin_d13C_carb_T1_i)
Sinsin_d13C_carb_T1_i_MTM=mtm(Sinsin_d13C_carb_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1)
write.csv(Sinsin_d13C_carb_T1_i_MTM,"Sinsin_d13C_carb_T1_MTM.csv",row.names = F)
setwd('/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Figure_4/')
Sinsin_obliquity_UKE=bandpass(Sinsin_d13C_carb_T1_i,flow = 0.02, fhigh = 0.04,xmax = 0.07)
write.csv(Sinsin_obliquity_UKE,"Sinsin_obliquity_UKE.csv",row.names = F)

SectionC_MS_T1_i=linterp(SectionC_MS_T1,7)
SectionC_MS_T1_i=detrend(SectionC_MS_T1_i)
SectionC_MS_T1_MTM=mtm(SectionC_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1)
write.csv(SectionC_MS_T1_MTM,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/SectionC_MS_T1_MTM.csv",row.names = F)

Fuhe_MS_T1=rbind(Fuhe_MS_lower_T1,Fuhe_MS_upper_T1)
Fuhe_MS_T1=sortNave(Fuhe_MS_T1)
Fuhe_MS_T1_i=linterp(Fuhe_MS_T1,5)
Fuhe_MS_T1_i=detrend(Fuhe_MS_T1_i)
Fuhe_MS_T1_MTM=mtm(Fuhe_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1)
write.csv(Fuhe_MS_T1_MTM,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Fuhe_MS_T1_MTM.csv",row.names = F)

Fuhe_d13C_carb_T1_i=linterp(Fuhe_d13C_carb_highres_T1,3)
Fuhe_d13C_carb_T1_i=detrend(Fuhe_d13C_carb_T1_i)
Fuhe_d13C_carb_T1_MTM=mtm(Fuhe_d13C_carb_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1)
write.csv(Fuhe_d13C_carb_T1_MTM,"/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/Fuhe_d13C_carb_T1_MTM.csv",row.names = F)
setwd('/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Figure_4/')
Fuhe_obliquity_UKE=bandpass(Fuhe_d13C_carb_T1_i,flow = 0.02, fhigh = 0.04,xmax = 0.07)
write.csv(Fuhe_obliquity_UKE,"Fuhe_obliquity_UKE.csv",row.names = F)

# Data displayed in Figure 4 ---- 
dev.off()
Fuhe_d13C_carb_T1_eha=eha(Fuhe_d13C_carb_T1_i,tbw=2,win=150,output=2,step=3)
Fuhe_obliquity_power=c()
Fuhe_total_power=c()
Fuhe_obliquity_total=c()
Position=seq(-158.45, 96.55, by = 3)
for (i in 2:length(Fuhe_d13C_carb_T1_eha[1,])){
  Fuhe_obliquity_power[i-1]=sum(Fuhe_d13C_carb_T1_eha[c(7:16),i])
  Fuhe_total_power[i-1]=sum(Fuhe_d13C_carb_T1_eha[c(1:39),i])
  Fuhe_obliquity_total[i-1]=Fuhe_obliquity_power[i-1]/Fuhe_total_power[i-1]
}
Fuhe_obliquity_power=cbind(Position,Fuhe_obliquity_power)
Fuhe_obliquity_total=cbind(Position,Fuhe_obliquity_total)
plot(Fuhe_obliquity_total)
write.csv(Fuhe_obliquity_power,"Fuhe_obliquity_power.csv",row.names = F) # Black line in Figure 4 for the Fuhe section
write.csv(Fuhe_obliquity_total,"Fuhe_obliquity_total.csv",row.names = F) # Red line in Figure 4 for the Fuhe section

setwd('/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Data_time_domain/')
Kowala_d13C_carb_T1_i=linterp(Kowala_d13C_carb_highres_T1,5)
Kowala_d13C_carb_T1_i=detrend(Kowala_d13C_carb_T1_i)
Kowala_d13C_carb_T1_MTM=mtm(Kowala_d13C_carb_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07)
write.csv(Kowala_d13C_carb_T1_MTM,"Kowala_d13C_carb_T1_MTM.csv",row.names = F)
setwd('/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Figure_4/')
Kowala_obliquity_UKE=bandpass(Kowala_d13C_carb_T1_i,flow = 0.02, fhigh = 0.04,xmax = 0.07)
write.csv(Kowala_obliquity_UKE,"Kowala_obliquity_UKE.csv",row.names = F) # Obliquity filter shown on Figure 4

dev.off()
Kowala_d13C_carb_T1_eha=eha(Kowala_d13C_carb_T1_i,tbw=2,win=300,output=2,step=5)
Kowala_obliquity_power=c()
Kowala_total_power=c()
Kowala_obliquity_total=c()
Position=seq(-168.92, 616.08, by = 5)
for (i in 2:length(Kowala_d13C_carb_T1_eha[1,])){
  Kowala_obliquity_power[i-1]=sum(Kowala_d13C_carb_T1_eha[c(11:26),i])
  Kowala_total_power[i-1]=sum(Kowala_d13C_carb_T1_eha[,i])
  Kowala_obliquity_total[i-1]=Kowala_obliquity_power[i-1]/Kowala_total_power[i-1]
}
Kowala_obliquity_power=cbind(Position,Kowala_obliquity_power)
Kowala_obliquity_total=cbind(Position,Kowala_obliquity_total)
plot(Kowala_obliquity_total)
write.csv(Kowala_obliquity_power,"Kowala_obliquity_power.csv",row.names = F) # Black line in Figure 4 for the Kowala section
write.csv(Kowala_obliquity_total,"Kowala_obliquity_total.csv",row.names = F) # Red line in Figure 4 for the Kowala section.

