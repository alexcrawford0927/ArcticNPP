library(leaps)
library(car)
library(lmtest)

################ READ IN DATA ##################
v = 5
ct = '75'

dfa <- read.csv(paste0("/Volumes/Prospero/Arrigo/Regional",ct,"/ArrigoData_Regional_Annual_v",v,".csv"))
dfm <- read.csv(paste0("/Volumes/Prospero/Arrigo/Regional",ct,"/ArrigoData_Regional_Monthly_v",v,".csv"))

###################################
########## FOR LOOP IT ############
###################################

# Set up
ms = c(0,5,6,7,8,9)
regions = c(7,8,9,10,11,13,16,17,18)
reglist <- sapply(regions,rep,times=length(ms))
dim(reglist) <- NULL

##### Standard Model with HWF, uWper, SST, and Open Days for NPPTOT ####
#####

dfr1 <- data.frame(month=rep(ms,length(regions)),reg=reglist,
                  r2=0,BIC=0,AIC=0,bp_p=NA,sw_p=NA,max_vif=NA,
                  int=0,icefreedays=0,oisst=0,hwf=0,uWper=0,
                  int_p=0,icefreedays_p=0,oisst_p=0,hwf_p=0,uWper_p=0,
                  int_e=0,icefreedays_e=0,oisst_e=0,hwf_e=0,uWper_e=0)

# The For Loop to Run the Models
for(m in ms){
  if (m == 0) {df <- dfa}
  if (m > 0) {df <- dfm[dfm$Month == m,]}
  for(reg in regions){
    # Create Model
    lm1 <- lm(npptot~icefreedays+oisst+hwf+uWper, data.frame(scale(df[df$Reg == reg,])) )
    
    # Assign Coefficients
    dfr1[dfr1$month == m & dfr1$reg == reg,c("int","icefreedays","oisst","hwf","uWper")] <- lm1$coefficients
    
    # Assign p-values & errors
    pvals <- summary(lm1)$coefficients[,4]
    pvals <- c(pvals,rep(NA,length(lm1$coefficients)-length(pvals)))
    dfr1[dfr1$month == m & dfr1$reg == reg,c("int_p","icefreedays_p","oisst_p","hwf_p","uWper_p")] <- pvals
    
    errors <- summary(lm1)$coefficients[,2]
    errors <- c(errors,rep(NA,length(lm1$coefficients)-length(errors)))
    dfr1[dfr1$month == m & dfr1$reg == reg,c("int_e","icefreedays_e","oisst_e","hwf_e","uWper_e")] <- errors
    
    # Assign Summary Statistics
    dfr1[dfr1$month == m & dfr1$reg == reg,"r2"] <- summary(lm1)$r.squared
    dfr1[dfr1$month == m & dfr1$reg == reg,"AIC"] <- AIC(lm1)
    dfr1[dfr1$month == m & dfr1$reg == reg,"BIC"] <- BIC(lm1)
    dfr1[dfr1$month == m & dfr1$reg == reg,"bp_p"] <- bptest(lm1)$p.value
    if (length(lm1$residuals) > 2 & sum(abs(lm1$residuals) > 0)){
      dfr1[dfr1$month == m & dfr1$reg == reg,"sw_p"] <- shapiro.test(lm1$residuals)$p.value }
    if (sum(is.na(pvals)) == 0){
      dfr1[dfr1$month == m & dfr1$reg == reg,"max_vif"] <- max(vif(lm1)) }
  }}
write.csv(dfr1,paste0("/Volumes/Prospero/Arrigo/Regional",ct,"/ArrigoData_Regional",v,"_4VarModels_icefreedays_sst_hwf_uWper.csv"),row.names=F)


##### Best BIC Model with HWF, SST, Open Days, uWper, vSper for npptot ####
#####
dfr2 <- data.frame(month=rep(ms,length(regions)),reg=reglist,
                   r2=0,BIC=0,AIC=0,bp_p=NA,sw_p=NA,max_vif=NA,
                   int=NA,icefreedays=NA,oisst=NA,hwf=NA,uWper=NA,
                   int_p=NA,icefreedays_p=NA,oisst_p=NA,hwf_p=NA,uWper_p=NA,
                   int_e=NA,icefreedays_e=NA,oisst_e=NA,hwf_e=NA,uWper_e=NA)

indices = c(12,9,10,14,16)

# The For Loop to Run the Models
for(m in ms){
  if (m == 0) {df <- dfa}
  if (m > 0) {df <- dfm[dfm$Month == m,c(-1)]}
  for(reg in regions){
    if (nrow(df[(df$Reg == reg),]) >= 5){
      # Regsubsets Extraction
      df <- df[is.finite(df$npprate) == 1,]
      dff <- data.frame(scale(df[(df$Reg == reg),]))
      regs <- regsubsets(npptot~., dff[,indices])
      
      xvarinds = indices[2:length(indices)][summary(regs)$which[(summary(regs)$bic == min(summary(regs)$bic,na.rm=T)) & (is.finite(summary(regs)$bic == 1))][2:length(indices)]]
      
      # Create Model
      lm1 <- lm(npptot~.,dff[,c(indices[1],xvarinds)])
      coef1 <- summary(lm1)$coefficients
      
      # Assign Coefficients & p-values
      dfr2[dfr2$month == m & dfr2$reg == reg,"int"] <- coef1[1,1]
      dfr2[dfr2$month == m & dfr2$reg == reg,paste0("int_e")] <- coef1[1,2]
      dfr2[dfr2$month == m & dfr2$reg == reg,paste0("int_p")] <- coef1[1,4]
      
      for (xvar in rownames(coef1)[2:length(rownames(coef1))]){
        dfr2[dfr2$month == m & dfr2$reg == reg,xvar] <- coef1[xvar,1]
        dfr2[dfr2$month == m & dfr2$reg == reg,paste0(xvar,"_e")] <- coef1[xvar,2]
        dfr2[dfr2$month == m & dfr2$reg == reg,paste0(xvar,"_p")] <- coef1[xvar,4]
      }
      
      # Assign Summary Statistics
      dfr2[dfr2$month == m & dfr2$reg == reg,"r2"] <- summary(lm1)$r.squared
      dfr2[dfr2$month == m & dfr2$reg == reg,"AIC"] <- AIC(lm1)
      dfr2[dfr2$month == m & dfr2$reg == reg,"BIC"] <- BIC(lm1)
      dfr2[dfr2$month == m & dfr2$reg == reg,"bp_p"] <- bptest(lm1)$p.value
      if (length(lm1$residuals) > 2 & sum(abs(lm1$residuals) > 0)){
        dfr2[dfr2$month == m & dfr2$reg == reg,"sw_p"] <- shapiro.test(lm1$residuals)$p.value }
      if (sum(is.finite(coef1[,4])) > 2){
        dfr2[dfr2$month == m & dfr2$reg == reg,"max_vif"] <- max(vif(lm1)) }
    }}}

write.csv(dfr2,paste0("/Volumes/Prospero/Arrigo/Regional",ct,"/ArrigoData_Regional",v,"_BestBICModels_4Var_icefreedays_sst_hwf_uWper.csv"),row.names=F)

######################################################################################
# EXPERIMENTATION #

# VIFS
vif(lm(npptot~oisst+hwf+uWper,dfa[dfa$Reg == 8,]))
vif(lm(npptot~icefreedays+oisst+uWper,dfa[dfa$Reg == 9,]))
vif(lm(npptot~icefreedays+oisst+hwf,dfa[dfa$Reg == 10,]))
vif(lm(npptot~icefreedays+oisst+hwf,dfa[dfa$Reg == 11,]))
vif(lm(npptot~icefreedays+oisst+hwf,dfa[dfa$Reg == 17,]))
vif(lm(npptot~icefreedays+oisst,dfa[dfa$Reg == 18,]))
vif(lm(npptot~icefreedays+uWper,dfa[dfa$Reg == 13,]))

# June
vif(lm(npptot~icefreedays+hwf,dfm[dfm$Reg == 10 & dfm$Month == 6,]))
vif(lm(npptot~icefreedays+oisst,dfm[dfm$Reg == 9 & dfm$Month == 6,]))
vif(lm(npptot~icefreedays+oisst+hwf+uWper,dfm[dfm$Reg == 18 & dfm$Month == 6,]))
vif(lm(npptot~oisst+uWper,dfm[dfm$Reg == 13 & dfm$Month == 6,]))

# July
vif(lm(npptot~icefreedays+oisst,dfm[dfm$Reg == 9 & dfm$Month == 7,]))
vif(lm(npptot~icefreedays+hwf,dfm[dfm$Reg == 10 & dfm$Month == 7,]))
vif(lm(npptot~icefreedays+uWper,dfm[dfm$Reg == 18 & dfm$Month == 7,]))
vif(lm(npptot~icefreedays+oisst+hwf+uWper,dfm[dfm$Reg == 13 & dfm$Month == 7,]))

# August
vif(lm(npptot~icefreedays+oisst,dfm[dfm$Reg == 10 & dfm$Month == 8,]))
vif(lm(npptot~icefreedays+oisst,dfm[dfm$Reg == 11 & dfm$Month == 8,]))
vif(lm(npptot~icefreedays+oisst+hwf,dfm[dfm$Reg == 17 & dfm$Month == 8,]))
vif(lm(npptot~icefreedays+oisst,dfm[dfm$Reg == 18 & dfm$Month == 8,]))
vif(lm(npptot~icefreedays+oisst+hwf+uWper,dfm[dfm$Reg == 13 & dfm$Month == 8,]))

# September
vif(lm(npptot~oisst+uWper,dfm[dfm$Reg == 10 & dfm$Month == 9,]))

