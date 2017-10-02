rm(list=ls())
library(astrochron)
setwd('/Users/pdoc3/Documents/Devonian/F-F global astrochronology paper/Nature Communications/R1/Monte Carlo demonstration/')

#### How the disturbed ETP series were created. ----
# etp=etp(tmin = 0, tmax = 3000, dt = 5, eWt = 1,oWt = 1,pWt = 1,esinw = T,genplot = F,standardize = T)
# colnames(etp)<-c("Time","ETP")
# 
# sedrate_S1=etp
# sedrate_S1[,2]=3+sedrate_S1$Time*0.0005
# colnames(sedrate_S1)<-c("Time","Sedrate_S1 (cm/kyr)")
# plot(sedrate_S1,type="l", ylim = c(3,4.5), xlab = "Time", ylab = "Sedimentation Rate (cm/kyr)")
# 
# sedrate_S2=etp
# sedrate_S2[,2]=3.75+0.75*sin(2*pi*sedrate_S2$Time/400)
# colnames(sedrate_S2)<-c("Time","Sedrate_S2 (cm/kyr)")
# lines(sedrate_S2,type="l",col="red")
# 
# sedrate_S3=etp
# sedrate_S3[1:201,2]=4.5
# sedrate_S3[202:401,2]=3.75
# sedrate_S3[402:601,2]=3
# colnames(sedrate_S3)<-c("Time","Sedrate_S3 (cm/kyr)")
# lines(sedrate_S3,type="l",col="blue")
# 
# 
# 
# signal_1=ar1etp(etpdat=etp,nsim=1,rho=1,wtAR=sqrt(1))
# for (i in 2:length(signal_1[,1])) {
#   signal_1[i,1]=signal_1[i-1,1]+((sedrate_S1[i,1]-sedrate_S1[i-1,1])*sedrate_S1[i,2]/100)
# }
# 
# signal_2=ar1etp(etpdat=etp,nsim=1,rho=1,wtAR=sqrt(1))
# for (i in 2:length(signal_2[,1])) {
#   signal_2[i,1]=signal_2[i-1,1]+((sedrate_S2[i,1]-sedrate_S2[i-1,1])*sedrate_S2[i,2]/100)
# }
# 
# signal_3=ar1etp(etpdat=etp,nsim=1,rho=1,wtAR=sqrt(1))
# for (i in 2:length(signal_3[,1])) {
#   signal_3[i,1]=signal_3[i-1,1]+((sedrate_S3[i,1]-sedrate_S3[i-1,1])*sedrate_S3[i,2]/100)
# }
#write.csv(signal_1,"signal_1.csv",row.names = F)        
#write.csv(signal_2,"signal_2.csv",row.names = F)        
#write.csv(signal_3,"signal_3.csv",row.names = F)        
#write.csv(etp,"etp.csv",row.names = F)

#### Load disturbed ETP signals ---- 
signal_1=read.csv("signal_1.csv",header=T)
signal_2=read.csv("signal_2.csv",header=T)
signal_3=read.csv("signal_3.csv",header=T)


# Create Supplementary Figure S3 ----
pl(r=1,c=4,mar=c(3,3,1,3),w = 100)

# plot (A) ETP
plot(etp[,2], etp[,1], col = "black", xlab = "ETP", ylab = "Time (ka)", main = "",type="l",lwd=1.2,ylim=c(3000,0), yaxs="i")
text(5,50,"(a)",cex=1.3,font=2)          
abline(h=210,col="blue")
abline(h=555,col="blue")
abline(h=855,col="blue")
abline(h=1015,col="blue")
abline(h=1300,col="blue")
abline(h=1730,col="blue")
abline(h=2020,col="blue")
abline(h=2222,col="blue")
abline(h=2445,col="blue")
abline(h=2955,col="blue")

# plot (B) Signal 1
plot(signal_1[,2], signal_1[,1], col = "black", xlab = "ETP", ylab = "Height (m)", main = "",type="l",lwd=1.2,ylim=c(113,0), yaxs="i")
text(3,2,"(b)",cex=1.3,font=2)          
abline(h=6.4,col="blue")
abline(h=17.3,col="blue")
abline(h=27.5,col="blue")
abline(h=33,col="blue")
abline(h=43.25,col="blue")
abline(h=59.25,col="blue")
abline(h=71,col="blue")
abline(h=79,col="blue")
abline(h=88.4,col="blue")
abline(h=110.5,col="blue")

# plot (C) Signal 2
plot(signal_2[,2], signal_2[,1], col = "black", xlab = "ETP", ylab = "Height (m)", main = "",type="l",lwd=1.2,ylim=c(113,0), yaxs="i")
text(3,2,"(c)",cex=1.3,font=2)          
abline(h=8.8,col="blue")
abline(h=21.75,col="blue")
abline(h=32.25,col="blue")
abline(h=39,col="blue")
abline(h=49.25,col="blue")
abline(h=65.5,col="blue")
abline(h=76,col="blue")
abline(h=84.25,col="blue")
abline(h=92,col="blue")
abline(h=111.75,col="blue")

# plot (D) Signal 3
plot(signal_3[,2], signal_3[,1], col = "black", xlab = "ETP", ylab = "Height (m)", main = "",type="l",lwd=1.2,ylim=c(113,0), yaxs="i")
text(3,2,"(d)",cex=1.3,font=2)          
abline(h=9.4,col="blue")
abline(h=25,col="blue")
abline(h=38.5,col="blue")
abline(h=45.75,col="blue")
abline(h=56.25,col="blue")
abline(h=72.25,col="blue")
abline(h=83.25,col="blue")
abline(h=89.1,col="blue")
abline(h=96,col="blue")
abline(h=111.2,col="blue")

# Start from "rough" tie-point ages. These are the ages listed in black on Suppl. Figure 3 ----
Pointers_original=c(3000,2500,2200,2000,1750,1300,1000,900,600,250)
offset=c(500,300,200,250,450,300,100,300,350) #The time difference between the tie-point ages in the line above
# Initiate some variables ----
astropower_sum=c()
astropower_product=c()
Pointers_product_best=c()
astropower_product_best=1
Pointers_sum_best=c()
astropower_sum_best=1
offset_matrix=matrix(0,10000,9)
offset_thisloop=c()
offset_fit=c()
start_i=1
sd_disturbance=0.02


# Monte Carlo loop ----
for (i in 1:5000) {
  if (i %in% c(401,801,1201,1601,2001,2401,2801,3201,3601,4001,4401,4801,5201,5601,6001,6401,6801,7201,7601,8001,8401,8801,9201,9601)==TRUE) {
    start_i=i
    offset=c(500,300,200,250,450,300,100,300,350)
  }
  Pointers=c()
  Pointers[1]=Pointers_original[1]
  for (j in 1:9) {
    Pointers[j+1]=Pointers[j]-offset[j]*rnorm(1,1,sd_disturbance) # Disturb the time-differences between consecutive tie-points
  }
  
  Pointers_S1=read.csv("Pointers_S1.csv",header=F) # Read the "rough" tie-point ages for H-32 stratigraphic levels
  colnames(Pointers_S1)=c('Height',"Relative Age (kyr)")
  Pointers_S1[,2]=Pointers[c(10:1)] # Replace the "rough" tie-point ages by "disturbed" tie-point ages
  
  Pointers_S2=read.csv("Pointers_S2.csv",header=F)
  colnames(Pointers_S2)=c('Height',"Relative Age (kyr)")
  Pointers_S2[,2]=Pointers[c(10:1)] 
  
  Pointers_S3=read.csv("Pointers_S3.csv",header=F)
  colnames(Pointers_S3)=c('Height',"Relative Age (kyr))")
  Pointers_S3[,2]=Pointers[c(10:1)] 
  
  # Calibrate the proxy records using the "tune" function in astrochron ----

  # Signal 1 ----
  signal_1=read.csv("signal_1.csv",header=T)
  signal_1_T1=tune(signal_1,Pointers_S1,extrapolate=T)
  
  # Signal 2 ----
  signal_2=read.csv("signal_2.csv",header=T)
  signal_2_T1=tune(signal_2,Pointers_S2,extrapolate=T)
  
  # Signal 3 ----
  signal_3=read.csv("signal_3.csv",header=T)
  signal_3_T1=tune(signal_3,Pointers_S3,extrapolate=T)
  
  ############ MTM ######### ----
  signal_1_T1_i=linterp(signal_1_T1,4)
  signal_1_T1_i=detrend(signal_1_T1_i)
  signal_1_T1_MTM=mtm(signal_1_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07) # Calculate MTM power spectrum for signal 1, according to the "disturbed" tie-point ages
  #write.csv(signal_1_T1_MTM,"Figure S3/signal_1_MTM_optimized.csv",row.names = F)
  
  idx1=which(signal_1_T1_MTM$Frequency>0.008 & signal_1_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(signal_1_T1_MTM$Frequency>0.02 & signal_1_T1_MTM$Frequency<0.03) # Obliquity band
  idx3=which(signal_1_T1_MTM$Frequency>0.04 & signal_1_T1_MTM$Frequency<0.048) # Precession 1 band
  idx4=which(signal_1_T1_MTM$Frequency>0.048 & signal_1_T1_MTM$Frequency<0.056) # Precession 1 band
  
  signal_1_100peak_idx=which(signal_1_T1_MTM$Power==max(signal_1_T1_MTM$Power[idx1])) 
  signal_1_100_misfit=abs(signal_1_T1_MTM$Frequency[signal_1_100peak_idx]-0.0105)/0.0105
  signal_1_100_height=(100-signal_1_T1_MTM$AR1_CL[signal_1_100peak_idx])/100
  signal_1_100peak_Fidx=which(signal_1_T1_MTM$Harmonic_CL==max(signal_1_T1_MTM$Harmonic_CL[idx1])) 
  signal_1_100_Fmisfit=abs(signal_1_T1_MTM$Frequency[signal_1_100peak_Fidx]-0.0105)/0.0105
  signal_1_100_Fheight=(100-signal_1_T1_MTM$Harmonic_CL[signal_1_100peak_Fidx])/100
  
  signal_1_obpeak_idx=which(signal_1_T1_MTM$Power==max(signal_1_T1_MTM$Power[idx2])) 
  signal_1_ob_misfit=abs(signal_1_T1_MTM$Frequency[signal_1_obpeak_idx]-0.0244)/0.0244
  signal_1_ob_height=(100-signal_1_T1_MTM$AR1_CL[signal_1_obpeak_idx])/100
  signal_1_obpeak_Fidx=which(signal_1_T1_MTM$Harmonic_CL==max(signal_1_T1_MTM$Harmonic_CL[idx2])) 
  signal_1_ob_Fmisfit=abs(signal_1_T1_MTM$Frequency[signal_1_obpeak_Fidx]-0.0244)/0.0244
  signal_1_ob_Fheight=(100-signal_1_T1_MTM$Harmonic_CL[signal_1_obpeak_Fidx])/100
  
  signal_1_P1peak_idx=which(signal_1_T1_MTM$Power==max(signal_1_T1_MTM$Power[idx3]))
  signal_1_P1_misfit=abs(signal_1_T1_MTM$Frequency[signal_1_P1peak_idx]-0.0435)/0.0435
  signal_1_P1_height=(100-signal_1_T1_MTM$AR1_CL[signal_1_P1peak_idx])/100
  signal_1_P1peak_Fidx=which(signal_1_T1_MTM$Harmonic_CL==max(signal_1_T1_MTM$Harmonic_CL[idx3])) 
  signal_1_P1_Fmisfit=abs(signal_1_T1_MTM$Frequency[signal_1_P1peak_Fidx]-0.0435)/0.0435
  signal_1_P1_Fheight=(100-signal_1_T1_MTM$Harmonic_CL[signal_1_P1peak_Fidx])/100
  
  signal_1_P2peak_idx=which(signal_1_T1_MTM$Power==max(signal_1_T1_MTM$Power[idx4]))
  signal_1_P2_misfit=abs(signal_1_T1_MTM$Frequency[signal_1_P2peak_idx]-0.0526)/0.0526
  signal_1_P2_height=(100-signal_1_T1_MTM$AR1_CL[signal_1_P2peak_idx])/100
  signal_1_P2peak_Fidx=which(signal_1_T1_MTM$Harmonic_CL==max(signal_1_T1_MTM$Harmonic_CL[idx4])) 
  signal_1_P2_Fmisfit=abs(signal_1_T1_MTM$Frequency[signal_1_P2peak_Fidx]-0.0526)/0.0526
  signal_1_P2_Fheight=(100-signal_1_T1_MTM$Harmonic_CL[signal_1_P2peak_Fidx])/100
  
  signal_1_misfit=mean(c(signal_1_100_misfit,signal_1_ob_misfit,signal_1_P1_misfit,signal_1_P2_misfit,signal_1_100_height,signal_1_ob_height,signal_1_P1_height,signal_1_P2_height,signal_1_100_Fmisfit,signal_1_ob_Fmisfit,signal_1_P1_Fmisfit,signal_1_P2_Fmisfit,signal_1_100_Fheight,signal_1_ob_Fheight,signal_1_P1_Fheight,signal_1_P2_Fheight))
  

  signal_2_T1_i=linterp(signal_2_T1,4)
  signal_2_T1_i=detrend(signal_2_T1_i)
  signal_2_T1_MTM=mtm(signal_2_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax=0.07) # Calculate MTM power spectrum for signal 2, according to the "disturbed" tie-point ages
  #write.csv(signal_2_T1_MTM,"Figure S3/signal_2_MTM_optimized.csv",row.names = F)
  
  idx1=which(signal_2_T1_MTM$Frequency>0.008 & signal_2_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(signal_2_T1_MTM$Frequency>0.02 & signal_2_T1_MTM$Frequency<0.03) # Obliquity band
  idx3=which(signal_2_T1_MTM$Frequency>0.04 & signal_2_T1_MTM$Frequency<0.048) # Precession 1 band
  idx4=which(signal_2_T1_MTM$Frequency>0.048 & signal_2_T1_MTM$Frequency<0.056) # Precession 2 band
  
  signal_2_100peak_idx=which(signal_2_T1_MTM$Power==max(signal_2_T1_MTM$Power[idx1])) 
  signal_2_100_misfit=abs(signal_2_T1_MTM$Frequency[signal_2_100peak_idx]-0.0105)/0.0105
  signal_2_100_height=(100-signal_2_T1_MTM$AR1_CL[signal_2_100peak_idx])/100
  signal_2_100peak_Fidx=which(signal_2_T1_MTM$Harmonic_CL==max(signal_2_T1_MTM$Harmonic_CL[idx1])) 
  signal_2_100_Fmisfit=abs(signal_2_T1_MTM$Frequency[signal_2_100peak_Fidx]-0.0105)/0.0105
  signal_2_100_Fheight=(100-signal_2_T1_MTM$Harmonic_CL[signal_2_100peak_Fidx])/100
  
  signal_2_obpeak_idx=which(signal_2_T1_MTM$Power==max(signal_2_T1_MTM$Power[idx2])) 
  signal_2_ob_misfit=abs(signal_2_T1_MTM$Frequency[signal_2_obpeak_idx]-0.0244)/0.0244
  signal_2_ob_height=(100-signal_2_T1_MTM$AR1_CL[signal_2_obpeak_idx])/100
  signal_2_obpeak_Fidx=which(signal_2_T1_MTM$Harmonic_CL==max(signal_2_T1_MTM$Harmonic_CL[idx2])) 
  signal_2_ob_Fmisfit=abs(signal_2_T1_MTM$Frequency[signal_2_obpeak_Fidx]-0.0244)/0.0244
  signal_2_ob_Fheight=(100-signal_2_T1_MTM$Harmonic_CL[signal_2_obpeak_Fidx])/100
  
  signal_2_P1peak_idx=which(signal_2_T1_MTM$Power==max(signal_2_T1_MTM$Power[idx3]))
  signal_2_P1_misfit=abs(signal_2_T1_MTM$Frequency[signal_2_P1peak_idx]-0.0435)/0.0435
  signal_2_P1_height=(100-signal_2_T1_MTM$AR1_CL[signal_2_P1peak_idx])/100
  signal_2_P1peak_Fidx=which(signal_2_T1_MTM$Harmonic_CL==max(signal_2_T1_MTM$Harmonic_CL[idx3])) 
  signal_2_P1_Fmisfit=abs(signal_2_T1_MTM$Frequency[signal_2_P1peak_Fidx]-0.0435)/0.0435
  signal_2_P1_Fheight=(100-signal_2_T1_MTM$Harmonic_CL[signal_2_P1peak_Fidx])/100
  
  signal_2_P2peak_idx=which(signal_2_T1_MTM$Power==max(signal_2_T1_MTM$Power[idx4]))
  signal_2_P2_misfit=abs(signal_2_T1_MTM$Frequency[signal_2_P2peak_idx]-0.0526)/0.0526
  signal_2_P2_height=(100-signal_2_T1_MTM$AR1_CL[signal_2_P2peak_idx])/100
  signal_2_P2peak_Fidx=which(signal_2_T1_MTM$Harmonic_CL==max(signal_2_T1_MTM$Harmonic_CL[idx4])) 
  signal_2_P2_Fmisfit=abs(signal_2_T1_MTM$Frequency[signal_2_P2peak_Fidx]-0.0526)/0.0526
  signal_2_P2_Fheight=(100-signal_2_T1_MTM$Harmonic_CL[signal_2_P2peak_Fidx])/100
  
  signal_2_misfit=mean(c(signal_2_100_misfit,signal_2_ob_misfit,signal_2_P1_misfit,signal_2_P2_misfit,signal_2_100_height,signal_2_ob_height,signal_2_P1_height,signal_2_P2_height,signal_2_100_Fmisfit,signal_2_ob_Fmisfit,signal_2_P1_Fmisfit,signal_2_P2_Fmisfit,signal_2_100_Fheight,signal_2_ob_Fheight,signal_2_P1_Fheight,signal_2_P2_Fheight))
  
  
  signal_3_T1_i=linterp(signal_3_T1,4)
  signal_3_T1_i=detrend(signal_3_T1_i)
  signal_3_T1_MTM=mtm(signal_3_T1_i,tbw=2,padfac = 5,demean = F, detrend = T, ar1 = T, output = 1,xmax = 0.07) # Calculate MTM power spectrum for signal 3, according to the "disturbed" tie-point ages
 #write.csv(signal_3_T1_MTM,"Figure S3/signal_3_MTM_optimized.csv",row.names = F)
  
  idx1=which(signal_3_T1_MTM$Frequency>0.008 & signal_3_T1_MTM$Frequency<0.013) # 100-kyr eccentricity band
  idx2=which(signal_3_T1_MTM$Frequency>0.02 & signal_3_T1_MTM$Frequency<0.03) # Obliquity band
  idx3=which(signal_3_T1_MTM$Frequency>0.04 & signal_3_T1_MTM$Frequency<0.048) # Precession 1 band
  idx4=which(signal_3_T1_MTM$Frequency>0.048 & signal_3_T1_MTM$Frequency<0.056) # Precession 2 band
  
  signal_3_100peak_idx=which(signal_3_T1_MTM$Power==max(signal_3_T1_MTM$Power[idx1])) 
  signal_3_100_misfit=abs(signal_3_T1_MTM$Frequency[signal_3_100peak_idx]-0.0105)/0.0105
  signal_3_100_height=(100-signal_3_T1_MTM$AR1_CL[signal_3_100peak_idx])/100
  signal_3_100peak_Fidx=which(signal_3_T1_MTM$Harmonic_CL==max(signal_3_T1_MTM$Harmonic_CL[idx1])) 
  signal_3_100_Fmisfit=abs(signal_3_T1_MTM$Frequency[signal_3_100peak_Fidx]-0.0105)/0.0105
  signal_3_100_Fheight=(100-signal_3_T1_MTM$Harmonic_CL[signal_3_100peak_Fidx])/100
  
  signal_3_obpeak_idx=which(signal_3_T1_MTM$Power==max(signal_3_T1_MTM$Power[idx2])) 
  signal_3_ob_misfit=abs(signal_3_T1_MTM$Frequency[signal_3_obpeak_idx]-0.0244)/0.0244
  signal_3_ob_height=(100-signal_3_T1_MTM$AR1_CL[signal_3_obpeak_idx])/100
  signal_3_obpeak_Fidx=which(signal_3_T1_MTM$Harmonic_CL==max(signal_3_T1_MTM$Harmonic_CL[idx2])) 
  signal_3_ob_Fmisfit=abs(signal_3_T1_MTM$Frequency[signal_3_obpeak_Fidx]-0.0244)/0.0244
  signal_3_ob_Fheight=(100-signal_3_T1_MTM$Harmonic_CL[signal_3_obpeak_Fidx])/100
  
  signal_3_P1peak_idx=which(signal_3_T1_MTM$Power==max(signal_3_T1_MTM$Power[idx3]))
  signal_3_P1_misfit=abs(signal_3_T1_MTM$Frequency[signal_3_P1peak_idx]-0.0435)/0.0435
  signal_3_P1_height=(100-signal_3_T1_MTM$AR1_CL[signal_3_P1peak_idx])/100
  signal_3_P1peak_Fidx=which(signal_3_T1_MTM$Harmonic_CL==max(signal_3_T1_MTM$Harmonic_CL[idx3])) 
  signal_3_P1_Fmisfit=abs(signal_3_T1_MTM$Frequency[signal_3_P1peak_Fidx]-0.0435)/0.0435
  signal_3_P1_Fheight=(100-signal_3_T1_MTM$Harmonic_CL[signal_3_P1peak_Fidx])/100
  
  signal_3_P2peak_idx=which(signal_3_T1_MTM$Power==max(signal_3_T1_MTM$Power[idx4]))
  signal_3_P2_misfit=abs(signal_3_T1_MTM$Frequency[signal_3_P2peak_idx]-0.0526)/0.0526
  signal_3_P2_height=(100-signal_3_T1_MTM$AR1_CL[signal_3_P2peak_idx])/100
  signal_3_P2peak_Fidx=which(signal_3_T1_MTM$Harmonic_CL==max(signal_3_T1_MTM$Harmonic_CL[idx4])) 
  signal_3_P2_Fmisfit=abs(signal_3_T1_MTM$Frequency[signal_3_P2peak_Fidx]-0.0526)/0.0526
  signal_3_P2_Fheight=(100-signal_3_T1_MTM$Harmonic_CL[signal_3_P2peak_Fidx])/100
  
  signal_3_misfit=mean(c(signal_3_100_misfit,signal_3_ob_misfit,signal_3_P1_misfit,signal_3_P2_misfit,signal_3_100_height,signal_3_ob_height,signal_3_P1_height,signal_3_P2_height,signal_3_100_Fmisfit,signal_3_ob_Fmisfit,signal_3_P1_Fmisfit,signal_3_P2_Fmisfit,signal_3_100_Fheight,signal_3_ob_Fheight,signal_3_P1_Fheight,signal_3_P2_Fheight))
  
  # End sum
  astropower_sum[i]=sum(signal_3_misfit,signal_2_misfit,signal_1_misfit)
  #astropower_product[i]=signal_3_misfit*signal_2_misfit*signal_1_misfit
  for (k in 1:9) {offset_thisloop[k]=Pointers[k]-Pointers[k+1]}
  offset_matrix[i,c(1:9)]=offset_thisloop
  offset_fit[i]=sum(abs(offset_thisloop-c(510, 223, 202, 290, 430, 285, 160, 300, 345)))
  sd_disturbance=astropower_sum[i]^2

  if (astropower_sum[i]<astropower_sum_best){
    astropower_sum_best=astropower_sum[i]
    Pointers_sum_best=Pointers
    save.image(file = "Sum_best.Rdata")  }

  if (astropower_sum[i]<=min(astropower_sum[c(start_i:i)])){  
    for (k in 1:9) {offset[k]=Pointers[k]-Pointers[k+1]}
    }

  if (astropower_product[i]<astropower_product_best){
    astropower_product_best=astropower_product[i]
    Pointers_product_best=Pointers
    #save.image(file = "Product_best.Rdata")
  }
  
}

# Create Supplementary Figure S4 ----
plot(astropower_sum,type="l",log = "y")
fit=lm(astropower_sum~offset_fit)
par(xaxs="r",yaxs="i")
plot(offset_fit,astropower_sum,xlim=c(0,500),ylim=c(0.05,0.25),log = "y")
points(0, 0.073,col='red')
abline(h=0.073,col="red")
abline(fit, col="grey")

