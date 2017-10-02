rm(list=ls())
library(astrochron)
setwd('/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Monte Carlo R1/') #Set Working Directory

# Start from "rough" tie-point ages. These are the ages listed in BLACK on Figures 5 and 6 ----
Pointers_original=c(2750,2530,2430,2280,2130,1930,1700,1315,860,700,500,400,300,150,100,-20,-340,450)
offset=c(200,100,150,150,200,230,385,455,160,200,100,100,150,50,120,320,110) #The time difference between the tie-point ages in the line above

# Initiate some variables ----
astropower_sum=c()
Pointers_sum_best=c()
astropower_sum_best=1
offset_matrix=matrix(0,10000,17)
offset_thisloop=c()
start_i=1
sd_disturbance=0.02


# Monte Carlo loop ----
for (i in 1:5000) {
  if (i %in% c(401,801,1201,1601,2001,2401,2801,3201,3601,4001,4401,4801,5201,5601,6001,6401,6801,7201,7601,8001,8401,8801,9201,9601)==TRUE) {
    start_i=i
    offset=c(220,110,200,100,270,250,300,450,150,200,200,100,130,60,80,350,130)
  }
  Pointers=c()
  Pointers[1]=Pointers_original[1]
  for (j in 1:17) {
    Pointers[j+1]=Pointers[j]-offset[j]*rnorm(1,1,sd_disturbance) # Disturb the time-differences between consecutive tie-points
  }
  
  Pointers_H32=read.csv("Pointers_H32.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
  colnames(Pointers_H32)=c('Height',"Relative Age (kyr)")
  Pointers_H32[,2]=Pointers[c(1,2,4,6:15,17:18)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  
  Pointers_CG1=read.csv("Pointers_CG1.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
  colnames(Pointers_CG1)=c('Height',"Relative Age (kyr)")
  Pointers_CG1[,2]=Pointers[c(1:13)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  
  Pointers_SectionC=read.csv("Pointers_SectionC.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
  colnames(Pointers_SectionC)=c('Height',"Relative Age (kyr)")
  Pointers_SectionC[,2]=Pointers[c(1:18)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  
  Pointers_Sinsin=read.csv("Pointers_Sinsin.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
  colnames(Pointers_Sinsin)=c('Height',"Relative Age (kyr)")
  Pointers_Sinsin[,2]=Pointers[c(10:18)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  
  Pointers_Fuhe=read.csv("Pointers_Fuhe.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
  colnames(Pointers_Fuhe)=c('Height',"Relative Age (kyr)")
  Pointers_Fuhe[,2]=Pointers[c(10:17)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  
  Pointers_Kowala=read.csv("Pointers_Kowala.csv",header=T) # Read the "rough" tie-point ages for H-32 stratigraphic levels
  colnames(Pointers_Kowala)=c('Height',"Relative Age (kyr)")
  Pointers_Kowala[,2]=Pointers[c(10,12:14,16:17)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  
  # Calibrate the proxy records using the "tune" function in astrochron ----
  setwd("/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Data_depth_domain/")
  # H-32 ----
  H32_d13C_carb=read.csv("H32_d13C_run2.csv",header=T)
  H32_d13C_carb=H32_d13C_carb[,c(1,2)]
  H32_d13C_carb[,1]=H32_d13C_carb[,1]/100
  H32_d13C_carb_T1=tune(H32_d13C_carb,Pointers_H32,extrapolate=T) # Convert depth to time according to the "disturbed" tie-point ages
  
  H32_d13C_org=read.delim("d13Corg_H32.txt",header=T)
  H32_d13C_org=H32_d13C_org[,c(3,2)]
  H32_d13C_org_T1=tune(H32_d13C_org,Pointers_H32,extrapolate=T)
  
  H32_MS=read.csv("H32_MS.csv",header=T)
  H32_MS=trim(H32_MS[,c(1,2)])
  H32_MS[,2]=H32_MS[,2]*10^8
  H32_MS[,1]=H32_MS[,1]/100
  H32_MS_T1=tune(H32_MS,Pointers_H32,extrapolate=F)
  
  # CG-1 ----
  CG1_d13C_carb=read.delim("isotopes_CG1.txt",header=T)
  CG1_d13C_carb=sortNave(CG1_d13C_carb[,c(1,2)])
  CG1_d13C_carb_T1=tune(CG1_d13C_carb,Pointers_CG1,extrapolate = T)
  
  CG1_MS=read.delim("MS_CG1.txt",header=T)
  CG1_MS=sortNave(CG1_MS[,c(1,2)])
  CG1_MS_T1=tune(CG1_MS,Pointers_CG1,extrapolate = T)
  
  # Section C ----
  SectionC_d13C_org=read.csv("AlbertaC_d13Corg.csv",header=F)
  SectionC_d13C_org=sortNave(SectionC_d13C_org[,c(1,2)])
  SectionC_d13C_org_T1=tune(SectionC_d13C_org,Pointers_SectionC,extrapolate = T)
  
  SectionC_MS=read.csv("AlbertaC_MS.csv",header=T)
  SectionC_MS=sortNave(SectionC_MS[,c(1,2)])
  SectionC_MS_T1=tune(SectionC_MS,Pointers_SectionC,extrapolate = F)
  
  # Sinsin ----
  Sinsin_d13C_carb=read.delim("Sinsin_Isotopes.txt",header=T)
  Sinsin_d13C_carb=sortNave(Sinsin_d13C_carb[,c(1,2)])
  Sinsin_d13C_carb[,1]=Sinsin_d13C_carb[,1]/100
  Sinsin_d13C_carb_T1=tune(Sinsin_d13C_carb,Pointers_Sinsin,extrapolate = T)
  
  Sinsin_MS=read.delim("Sinsin_MS.txt",header=T)
  Sinsin_MS=sortNave(Sinsin_MS[,c(1,2)])
  Sinsin_MS[,2]=Sinsin_MS[,2]*10^8
  Sinsin_MS[,1]=Sinsin_MS[,1]/100
  Sinsin_MS_T1=tune(Sinsin_MS,Pointers_Sinsin,extrapolate = T)
  
  # Fuhe ----
  Fuhe_d13C_carb_highres=read.csv("Fuhe_d13C_carb_highres.csv",header=F)
  Fuhe_d13C_carb_highres=sortNave(Fuhe_d13C_carb_highres[,c(1,2)])
  Fuhe_d13C_carb_highres_T1=tune(Fuhe_d13C_carb_highres,Pointers_Fuhe,extrapolate = T)
  
  Fuhe_d13C_carb_lowres=read.csv("Fuhe_d13C_carb_lowres.csv",header=F)
  Fuhe_d13C_carb_lowres=sortNave(Fuhe_d13C_carb_lowres[,c(1,2)])
  Fuhe_d13C_carb_lowres_T1=tune(Fuhe_d13C_carb_lowres,Pointers_Fuhe,extrapolate = T)
  
  Fuhe_d13C_org_highres=read.csv("Fuhe_d13C_org_highres.csv",header=F)
  Fuhe_d13C_org_highres=sortNave(Fuhe_d13C_org_highres[,c(1,2)])
  Fuhe_d13C_org_highres_T1=tune(Fuhe_d13C_org_highres,Pointers_Fuhe,extrapolate = T)
  
  Fuhe_d13C_org_lowres=read.csv("Fuhe_d13C_org_lowres.csv",header=F)
  Fuhe_d13C_org_lowres=sortNave(Fuhe_d13C_org_lowres[,c(1,2)])
  Fuhe_d13C_org_lowres_T1=tune(Fuhe_d13C_org_lowres,Pointers_Fuhe,extrapolate = T)
  
  Fuhe_MS_lower=read.csv("Data_Zenghui.csv",header = F)
  Fuhe_MS_lower=sortNave(Fuhe_MS_lower[,c(1,2)])
  Fuhe_MS_lower_T1=tune(Fuhe_MS_lower,Pointers_Fuhe,extrapolate = T)
  
  Fuhe_MS_upper=read.csv("MagneticSusceptibility_upper.csv",header = F)
  Fuhe_MS_upper=sortNave(Fuhe_MS_upper[,c(1,2)])
  Fuhe_MS_upper[,2]=Fuhe_MS_upper[,2]*10^8
  Fuhe_MS_upper_T1=tune(Fuhe_MS_upper,Pointers_Fuhe,extrapolate = T)
  
  # Kowala ----
  Kowala_d13C_carb_highres=read.csv("Kowala_d13C_carb_highres.csv",header=F)
  Kowala_d13C_carb_highres=sortNave(Kowala_d13C_carb_highres[,c(1,2)])
  Kowala_d13C_carb_highres_T1=tune(Kowala_d13C_carb_highres,Pointers_Kowala,extrapolate = T)
  
  Kowala_d13C_carb_lowres=read.csv("Kowala_d13C_carb_lowres.csv",header=F)
  Kowala_d13C_carb_lowres=sortNave(Kowala_d13C_carb_lowres[,c(1,2)])
  Kowala_d13C_carb_lowres_T1=tune(Kowala_d13C_carb_lowres,Pointers_Kowala,extrapolate = T)
  
  Kowala_d13C_org=read.csv("Kowala_d13C_org_lowres.csv",header=F)
  Kowala_d13C_org=sortNave(Kowala_d13C_org[,c(1,2)])
  Kowala_d13C_org_T1=tune(Kowala_d13C_org,Pointers_Kowala,extrapolate = T)
  
  ############ MTM ######### ----
  setwd('/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Monte Carlo R1/')
  H32_MS_T1_i=linterp(H32_MS_T1,7)
  H32_MS_T1_i=detrend(H32_MS_T1_i)
  H32_MS_T1_MTM=mtm(H32_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07) # Calculate MTM power spectrum for H-32, according to the "disturbed" tie-point ages
  
  idx1=which(H32_MS_T1_MTM$Frequency>0.008 & H32_MS_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(H32_MS_T1_MTM$Frequency>0.025 & H32_MS_T1_MTM$Frequency<0.035) # Obliquity band
  idx3=which(H32_MS_T1_MTM$Frequency>0.045 & H32_MS_T1_MTM$Frequency<0.065) # Precession band

  H32_MS_100peak_idx=which(H32_MS_T1_MTM$Power==max(H32_MS_T1_MTM$Power[idx1])) 
  H32_MS_100_misfit=abs(H32_MS_T1_MTM$Frequency[H32_MS_100peak_idx]-0.0105)/0.0105
  H32_MS_100_height=(100-H32_MS_T1_MTM$AR1_CL[H32_MS_100peak_idx])/100
  H32_MS_100peak_Fidx=which(H32_MS_T1_MTM$Harmonic_CL==max(H32_MS_T1_MTM$Harmonic_CL[idx1])) 
  H32_MS_100_Fmisfit=abs(H32_MS_T1_MTM$Frequency[H32_MS_100peak_Fidx]-0.0105)/0.0105
  H32_MS_100_Fheight=(100-H32_MS_T1_MTM$Harmonic_CL[H32_MS_100peak_Fidx])/100
  
  H32_MS_obpeak_idx=which(H32_MS_T1_MTM$Power==max(H32_MS_T1_MTM$Power[idx2])) 
  H32_MS_ob_misfit=abs(H32_MS_T1_MTM$Frequency[H32_MS_obpeak_idx]-0.031)/0.031
  H32_MS_ob_height=(100-H32_MS_T1_MTM$AR1_CL[H32_MS_obpeak_idx])/100
  H32_MS_obpeak_Fidx=which(H32_MS_T1_MTM$Harmonic_CL==max(H32_MS_T1_MTM$Harmonic_CL[idx2])) 
  H32_MS_ob_Fmisfit=abs(H32_MS_T1_MTM$Frequency[H32_MS_obpeak_Fidx]-0.031)/0.031
  H32_MS_ob_Fheight=(100-H32_MS_T1_MTM$Harmonic_CL[H32_MS_obpeak_Fidx])/100
  
  H32_MS_P1peak_idx=which(H32_MS_T1_MTM$Power==max(H32_MS_T1_MTM$Power[idx3]))
  H32_MS_P1_misfit=abs(H32_MS_T1_MTM$Frequency[H32_MS_P1peak_idx]-0.055)/0.055
  H32_MS_P1_height=(100-H32_MS_T1_MTM$AR1_CL[H32_MS_P1peak_idx])/100
  H32_MS_P1peak_Fidx=which(H32_MS_T1_MTM$Harmonic_CL==max(H32_MS_T1_MTM$Harmonic_CL[idx3])) 
  H32_MS_P1_Fmisfit=abs(H32_MS_T1_MTM$Frequency[H32_MS_P1peak_Fidx]-0.055)/0.055
  H32_MS_P1_Fheight=(100-H32_MS_T1_MTM$Harmonic_CL[H32_MS_P1peak_Fidx])/100
  
  H32_MS_misfit=mean(c(H32_MS_100_misfit,H32_MS_ob_misfit,H32_MS_P1_misfit,H32_MS_100_height,H32_MS_ob_height,H32_MS_P1_height))
  
  CG1_MS_T1_i=linterp(CG1_MS_T1,5)
  CG1_MS_T1_i=detrend(CG1_MS_T1_i)
  CG1_MS_T1_MTM=mtm(CG1_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07) # Calculate MTM power spectrum for CG1, according to the "disturbed" tie-point ages
  
  idx1=which(CG1_MS_T1_MTM$Frequency>0.008 & CG1_MS_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(CG1_MS_T1_MTM$Frequency>0.025 & CG1_MS_T1_MTM$Frequency<0.035) # Obliquity band
  idx3=which(CG1_MS_T1_MTM$Frequency>0.045 & CG1_MS_T1_MTM$Frequency<0.065) # Precession 1 band

  CG1_MS_100peak_idx=which(CG1_MS_T1_MTM$Power==max(CG1_MS_T1_MTM$Power[idx1])) 
  CG1_MS_100_misfit=abs(CG1_MS_T1_MTM$Frequency[CG1_MS_100peak_idx]-0.0105)/0.0105
  CG1_MS_100_height=(100-CG1_MS_T1_MTM$AR1_CL[CG1_MS_100peak_idx])/100
  CG1_MS_100peak_Fidx=which(CG1_MS_T1_MTM$Harmonic_CL==max(CG1_MS_T1_MTM$Harmonic_CL[idx1])) 
  CG1_MS_100_Fmisfit=abs(CG1_MS_T1_MTM$Frequency[CG1_MS_100peak_Fidx]-0.0105)/0.0105
  CG1_MS_100_Fheight=(100-CG1_MS_T1_MTM$Harmonic_CL[CG1_MS_100peak_Fidx])/100
  
  CG1_MS_obpeak_idx=which(CG1_MS_T1_MTM$Power==max(CG1_MS_T1_MTM$Power[idx2])) 
  CG1_MS_ob_misfit=abs(CG1_MS_T1_MTM$Frequency[CG1_MS_obpeak_idx]-0.031)/0.031
  CG1_MS_ob_height=(100-CG1_MS_T1_MTM$AR1_CL[CG1_MS_obpeak_idx])/100
  CG1_MS_obpeak_Fidx=which(CG1_MS_T1_MTM$Harmonic_CL==max(CG1_MS_T1_MTM$Harmonic_CL[idx2])) 
  CG1_MS_ob_Fmisfit=abs(CG1_MS_T1_MTM$Frequency[CG1_MS_obpeak_Fidx]-0.031)/0.031
  CG1_MS_ob_Fheight=(100-CG1_MS_T1_MTM$Harmonic_CL[CG1_MS_obpeak_Fidx])/100
  
  CG1_MS_P1peak_idx=which(CG1_MS_T1_MTM$Power==max(CG1_MS_T1_MTM$Power[idx3]))
  CG1_MS_P1_misfit=abs(CG1_MS_T1_MTM$Frequency[CG1_MS_P1peak_idx]-0.055)/0.055
  CG1_MS_P1_height=(100-CG1_MS_T1_MTM$AR1_CL[CG1_MS_P1peak_idx])/100
  CG1_MS_P1peak_Fidx=which(CG1_MS_T1_MTM$Harmonic_CL==max(CG1_MS_T1_MTM$Harmonic_CL[idx3])) 
  CG1_MS_P1_Fmisfit=abs(CG1_MS_T1_MTM$Frequency[CG1_MS_P1peak_Fidx]-0.055)/0.055
  CG1_MS_P1_Fheight=(100-CG1_MS_T1_MTM$Harmonic_CL[CG1_MS_P1peak_Fidx])/100

  CG1_MS_misfit=mean(c(CG1_MS_100_misfit,CG1_MS_ob_misfit,CG1_MS_P1_misfit,CG1_MS_100_height,CG1_MS_ob_height,CG1_MS_P1_height))
  
  SectionC_MS_T1_i=linterp(SectionC_MS_T1,7)
  SectionC_MS_T1_i=detrend(SectionC_MS_T1_i)
  SectionC_MS_T1_MTM=mtm(SectionC_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07) # Calculate MTM power spectrum for Section C, according to the "disturbed" tie-point ages
  
  idx1=which(SectionC_MS_T1_MTM$Frequency>0.008 & SectionC_MS_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(SectionC_MS_T1_MTM$Frequency>0.025 & SectionC_MS_T1_MTM$Frequency<0.035) # Obliquity band
  idx3=which(SectionC_MS_T1_MTM$Frequency>0.045 & SectionC_MS_T1_MTM$Frequency<0.065) # Precession band

  SectionC_MS_100peak_idx=which(SectionC_MS_T1_MTM$Power==max(SectionC_MS_T1_MTM$Power[idx1])) 
  SectionC_MS_100_misfit=abs(SectionC_MS_T1_MTM$Frequency[SectionC_MS_100peak_idx]-0.0105)/0.0105
  SectionC_MS_100_height=(100-SectionC_MS_T1_MTM$AR1_CL[SectionC_MS_100peak_idx])/100
  SectionC_MS_100peak_Fidx=which(SectionC_MS_T1_MTM$Harmonic_CL==max(SectionC_MS_T1_MTM$Harmonic_CL[idx1])) 
  SectionC_MS_100_Fmisfit=abs(SectionC_MS_T1_MTM$Frequency[SectionC_MS_100peak_Fidx]-0.0105)/0.0105
  SectionC_MS_100_Fheight=(100-SectionC_MS_T1_MTM$Harmonic_CL[SectionC_MS_100peak_Fidx])/100
  
  SectionC_MS_obpeak_idx=which(SectionC_MS_T1_MTM$Power==max(SectionC_MS_T1_MTM$Power[idx2])) 
  SectionC_MS_ob_misfit=abs(SectionC_MS_T1_MTM$Frequency[SectionC_MS_obpeak_idx]-0.031)/0.031
  SectionC_MS_ob_height=(100-SectionC_MS_T1_MTM$AR1_CL[SectionC_MS_obpeak_idx])/100
  SectionC_MS_obpeak_Fidx=which(SectionC_MS_T1_MTM$Harmonic_CL==max(SectionC_MS_T1_MTM$Harmonic_CL[idx2])) 
  SectionC_MS_ob_Fmisfit=abs(SectionC_MS_T1_MTM$Frequency[SectionC_MS_obpeak_Fidx]-0.031)/0.031
  SectionC_MS_ob_Fheight=(100-SectionC_MS_T1_MTM$Harmonic_CL[SectionC_MS_obpeak_Fidx])/100
  
   SectionC_MS_misfit=mean(c(SectionC_MS_100_misfit,SectionC_MS_ob_misfit,SectionC_MS_100_height,SectionC_MS_ob_height))
  
  Sinsin_MS_T1_i=linterp(Sinsin_MS_T1,4)
  Sinsin_MS_T1_i=detrend(Sinsin_MS_T1_i)
  Sinsin_MS_T1_MTM=mtm(Sinsin_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07) # Calculate MTM power spectrum for Sinsin, according to the "disturbed" tie-point ages
  
  idx1=which(Sinsin_MS_T1_MTM$Frequency>0.008 & Sinsin_MS_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(Sinsin_MS_T1_MTM$Frequency>0.025 & Sinsin_MS_T1_MTM$Frequency<0.035) # Obliquity band
  idx3=which(Sinsin_MS_T1_MTM$Frequency>0.045 & Sinsin_MS_T1_MTM$Frequency<0.065) # Precession band

  Sinsin_MS_100peak_idx=which(Sinsin_MS_T1_MTM$Power==max(Sinsin_MS_T1_MTM$Power[idx1])) 
  Sinsin_MS_100_misfit=abs(Sinsin_MS_T1_MTM$Frequency[Sinsin_MS_100peak_idx]-0.0105)/0.0105
  Sinsin_MS_100_height=(100-Sinsin_MS_T1_MTM$AR1_CL[Sinsin_MS_100peak_idx])/100
  Sinsin_MS_100peak_Fidx=which(Sinsin_MS_T1_MTM$Harmonic_CL==max(Sinsin_MS_T1_MTM$Harmonic_CL[idx1])) 
  Sinsin_MS_100_Fmisfit=abs(Sinsin_MS_T1_MTM$Frequency[Sinsin_MS_100peak_Fidx]-0.0105)/0.0105
  Sinsin_MS_100_Fheight=(100-Sinsin_MS_T1_MTM$Harmonic_CL[Sinsin_MS_100peak_Fidx])/100
  
  Sinsin_MS_obpeak_idx=which(Sinsin_MS_T1_MTM$Power==max(Sinsin_MS_T1_MTM$Power[idx2])) 
  Sinsin_MS_ob_misfit=abs(Sinsin_MS_T1_MTM$Frequency[Sinsin_MS_obpeak_idx]-0.031)/0.031
  Sinsin_MS_ob_height=(100-Sinsin_MS_T1_MTM$AR1_CL[Sinsin_MS_obpeak_idx])/100
  Sinsin_MS_obpeak_Fidx=which(Sinsin_MS_T1_MTM$Harmonic_CL==max(Sinsin_MS_T1_MTM$Harmonic_CL[idx2])) 
  Sinsin_MS_ob_Fmisfit=abs(Sinsin_MS_T1_MTM$Frequency[Sinsin_MS_obpeak_Fidx]-0.031)/0.031
  Sinsin_MS_ob_Fheight=(100-Sinsin_MS_T1_MTM$Harmonic_CL[Sinsin_MS_obpeak_Fidx])/100
  
  Sinsin_MS_P1peak_idx=which(Sinsin_MS_T1_MTM$Power==max(Sinsin_MS_T1_MTM$Power[idx3]))
  Sinsin_MS_P1_misfit=abs(Sinsin_MS_T1_MTM$Frequency[Sinsin_MS_P1peak_idx]-0.055)/0.055
  Sinsin_MS_P1_height=(100-Sinsin_MS_T1_MTM$AR1_CL[Sinsin_MS_P1peak_idx])/100
  Sinsin_MS_P1peak_Fidx=which(Sinsin_MS_T1_MTM$Harmonic_CL==max(Sinsin_MS_T1_MTM$Harmonic_CL[idx3])) 
  Sinsin_MS_P1_Fmisfit=abs(Sinsin_MS_T1_MTM$Frequency[Sinsin_MS_P1peak_Fidx]-0.055)/0.055
  Sinsin_MS_P1_Fheight=(100-Sinsin_MS_T1_MTM$Harmonic_CL[Sinsin_MS_P1peak_Fidx])/100
  
  Sinsin_MS_misfit=mean(c(Sinsin_MS_100_misfit,Sinsin_MS_ob_misfit,Sinsin_MS_P1_misfit,Sinsin_MS_100_height,Sinsin_MS_ob_height,Sinsin_MS_P1_height))
  
  Fuhe_MS_T1=rbind(Fuhe_MS_lower_T1,Fuhe_MS_upper_T1)
  Fuhe_MS_T1=sortNave(Fuhe_MS_T1)
  Fuhe_MS_T1_i=linterp(Fuhe_MS_T1,5)
  Fuhe_MS_T1_i=detrend(Fuhe_MS_T1_i)
  Fuhe_MS_T1_MTM=mtm(Fuhe_MS_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07) # Calculate MTM power spectrum for H-32, according to the "disturbed" tie-point ages
  
  idx1=which(Fuhe_MS_T1_MTM$Frequency>0.008 & Fuhe_MS_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(Fuhe_MS_T1_MTM$Frequency>0.025 & Fuhe_MS_T1_MTM$Frequency<0.035) # Obliquity band
  idx3=which(Fuhe_MS_T1_MTM$Frequency>0.045 & Fuhe_MS_T1_MTM$Frequency<0.065) # Precession band

  Fuhe_MS_100peak_idx=which(Fuhe_MS_T1_MTM$Power==max(Fuhe_MS_T1_MTM$Power[idx1])) 
  Fuhe_MS_100_misfit=abs(Fuhe_MS_T1_MTM$Frequency[Fuhe_MS_100peak_idx]-0.0105)/0.0105
  Fuhe_MS_100_height=(100-Fuhe_MS_T1_MTM$AR1_CL[Fuhe_MS_100peak_idx])/100
  Fuhe_MS_100peak_Fidx=which(Fuhe_MS_T1_MTM$Harmonic_CL==max(Fuhe_MS_T1_MTM$Harmonic_CL[idx1])) 
  Fuhe_MS_100_Fmisfit=abs(Fuhe_MS_T1_MTM$Frequency[Fuhe_MS_100peak_Fidx]-0.0105)/0.0105
  Fuhe_MS_100_Fheight=(100-Fuhe_MS_T1_MTM$Harmonic_CL[Fuhe_MS_100peak_Fidx])/100
  
  Fuhe_MS_obpeak_idx=which(Fuhe_MS_T1_MTM$Power==max(Fuhe_MS_T1_MTM$Power[idx2])) 
  Fuhe_MS_ob_misfit=abs(Fuhe_MS_T1_MTM$Frequency[Fuhe_MS_obpeak_idx]-0.031)/0.031
  Fuhe_MS_ob_height=(100-Fuhe_MS_T1_MTM$AR1_CL[Fuhe_MS_obpeak_idx])/100
  Fuhe_MS_obpeak_Fidx=which(Fuhe_MS_T1_MTM$Harmonic_CL==max(Fuhe_MS_T1_MTM$Harmonic_CL[idx2])) 
  Fuhe_MS_ob_Fmisfit=abs(Fuhe_MS_T1_MTM$Frequency[Fuhe_MS_obpeak_Fidx]-0.031)/0.031
  Fuhe_MS_ob_Fheight=(100-Fuhe_MS_T1_MTM$Harmonic_CL[Fuhe_MS_obpeak_Fidx])/100
  
  Fuhe_MS_P1peak_idx=which(Fuhe_MS_T1_MTM$Power==max(Fuhe_MS_T1_MTM$Power[idx3]))
  Fuhe_MS_P1_misfit=abs(Fuhe_MS_T1_MTM$Frequency[Fuhe_MS_P1peak_idx]-0.055)/0.055
  Fuhe_MS_P1_height=(100-Fuhe_MS_T1_MTM$AR1_CL[Fuhe_MS_P1peak_idx])/100
  Fuhe_MS_P1peak_Fidx=which(Fuhe_MS_T1_MTM$Harmonic_CL==max(Fuhe_MS_T1_MTM$Harmonic_CL[idx3])) 
  Fuhe_MS_P1_Fmisfit=abs(Fuhe_MS_T1_MTM$Frequency[Fuhe_MS_P1peak_Fidx]-0.055)/0.055
  Fuhe_MS_P1_Fheight=(100-Fuhe_MS_T1_MTM$Harmonic_CL[Fuhe_MS_P1peak_Fidx])/100
  
  Fuhe_MS_misfit=mean(c(Fuhe_MS_100_misfit,Fuhe_MS_ob_misfit,Fuhe_MS_P1_misfit,Fuhe_MS_100_height,Fuhe_MS_ob_height,Fuhe_MS_P1_height))
  
  
  # End average of "tied-in" sections (H32, CG1, Sinsin and Fuhe)
  astropower_sum[i]=mean(c(H32_MS_misfit,CG1_MS_misfit,Sinsin_MS_misfit,Fuhe_MS_misfit))

  for (k in 1:17) {offset_thisloop[k]=Pointers[k]-Pointers[k+1]}
  offset_matrix[i,c(1:17)]=offset_thisloop
  sd_disturbance=astropower_sum[i]^1.5

  if (astropower_sum[i]<astropower_sum_best){
    astropower_sum_best=astropower_sum[i]
    Pointers_sum_best=Pointers
    save.image(file = "Sum_best.Rdata")  }

  if (astropower_sum[i]<=min(astropower_sum[c(start_i:i)])){  
    for (k in 1:17) {offset[k]=Pointers[k]-Pointers[k+1]}
    }

}
plot(astropower_sum,type="l",log = "y")

