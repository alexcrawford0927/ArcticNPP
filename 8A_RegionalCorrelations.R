
# Data Prep
vars = c("npptot","npprate","chl","oisst","icefreedays","opendays","hwf","uWper","vSper")
cors <- data.frame(Reg=NA,Month=NA,v1=NA,v2=NA,r=NA,p=NA)
row = 0
method = 'pearson'

# ANNUAL
dfa <- read.csv("/Volumes/Prospero/Arrigo/Regional75/ArrigoData_Regional_Annual_v4.csv")

for (r in unique(dfa$Reg)){
  for (v1 in vars){
    for (v2 in vars){
      row = row+1 # Advance Row
      
      # Assign values to new data frame
      cors[row,"Reg"] <- r
      cors[row,"v1"] <- v1
      cors[row,"v2"] <- v2
      cors[row,"Month"] <- 0
      
      # Extract vectors for correlating
      x1 <- dfa[dfa$Reg == r,v1]
      x2 <- dfa[dfa$Reg == r,v2]
      
      # If there are enough observations, correlate
      if (sum(is.finite(x1)) > 5 & sum(is.finite(x2)) > 5){
        cor <- cor.test(x1,x2, method=method) # Correlation Test
        
        # Assign to new data frame
        cors[row,"r"] <- cor$estimate
        cors[row,"p"] <- cor$p.value
      }
    }}}

# MONTHLY
dfm <- read.csv("/Volumes/Prospero/Arrigo/Regional75/ArrigoData_Regional_Monthly_v4.csv")

for (m in 5:9){
  for (r in unique(dfa$Reg)){
    for (v1 in vars){
      for (v2 in vars){
        row = row+1 # Advance Row
        
        # Assign values to new data frame
        cors[row,"Reg"] <- r
        cors[row,"v1"] <- v1
        cors[row,"v2"] <- v2
        cors[row,"Month"] <- m
        
        # Extract vectors for correlating
        x1 <- dfm[dfm$Reg == r & dfm$Month == m,v1]
        x2 <- dfm[dfm$Reg == r & dfm$Month == m,v2]
        
        # If there are enough observations, correlate
        if (sum(is.finite(x1)) > 5 & sum(is.finite(x2)) > 5){
          cor <- cor.test(x1,x2,method=method) # Correlation Test
          
          # Assign to new data frame
          cors[row,"r"] <- cor$estimate
          cors[row,"p"] <- cor$p.value
        }
}}}}

write.csv(cors,paste0("/Volumes/Prospero/Arrigo/Regional75/ArrigoData_Regional4_Correlations_",method,".csv"),row.names = F)
