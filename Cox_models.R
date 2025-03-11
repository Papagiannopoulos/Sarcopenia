
#Libraries
library(tidyverse)
library(survival); library(survminer)
library(biostat3)#For linear combinations
library(plotRCS); library(ggplot2)
library(xlsx); library(readxl)

#Import your data including your covariates, exposures and outcomes.
#The user needs to include the time to event feature named as 'time' and the stratification factor (in our research we used the 'sex' feature (binary))
data <- read_excel("")
#str(data)
#data.frame:	~ 400,000 obs. of  N variables:
# Covariate1      : num  25.1 25.8 36.9 24.8 30.3 ...
# Covariate2      : int  59 63 57 60 68 41 49 59 51 63 ...
# Covariate3      : Factor  0 0 1 0 0 0 1 0 2 0 ...
# Exposure1       : int  0 0 0 1 0 0 0 1 0 1 ...
# Exposure2       : int  0 1 0 0 0 1 0 0 0 0 ...
# Outcome1        : int  0 0 0 0 1 0 0 1 0 0 ...
# Outcome2        : int  0 1 0 0 0 1 1 0 0 1 ...
# Time to event   : num  12.2 12.6 11.2 11.9 10.7 ...
#.
#.
#.

#Create a vector of strings with elements the variable names of every Exposure in your data
exposure <- c("Exposure1","Exposure2")

#Create the corresponding abbreviations of your Exposure variable (with the same order as in exposure variable)
abr <- c("ABR - Exposure1","ABR - Exposure2")

#Create a vector of strings with elements the variable names of every Exposure in your data
covariates = c("Covariate1","Covariate2","Covariate3")

#Create a vector of strings with elements the variable names of every Outcome in your data
outcome <- c("Outcome1","Outcome2")
outcome <- c("Bladder","Brain","Breast","Colorectal","Endometrial","Gastric","Hematologic",
             "Kidney","Liver","Lung","Melanoma","Oesophageal","Oral",
             "Ovarian","Overall","Pancreatic","Prostate")

data <- data[, c(covariates, exposure, outcome, "time")]

crt <- coxph.control(1e-11/2, iter.max = 50)
#### Cox model - Model 1 ####
results <- list()
for(j in 1:length(exposure)){
  results[[j]] <- matrix(NA, nrow = length(outcome),ncol = 9)
  colnames(results[[j]]) <- c("Cancer","Exposure","PVal","Female","PVal","Male","PVal","Interaction", "PH")
  results[[j]][,1] <- outcome
  
  for(i in 1:length(outcome)){
    print(paste0(exposure[j], " --> ", outcome[i]))
    
    if(outcome[i] %in% c("Prostate")){
      data1 <- data %>% filter(sex == 1)#Keep only males
      form <- as.formula(paste0("Surv(time,",outcome[i],")~",exposure[j], "+", paste0(covariates, collapse = "+")))
      cox <- coxph(form, data = data1,ties = "efron",control = crt)
      #No interactions by sex in sex-specific level
      coxInter <- cox
    }else if(outcome[i] %in% c("Breast","Cervical","Endometrial","Ovarian")){
      data1 <- data %>% filter(sex == 0)#keep only females
      #if(sum(outcome[i] %in% c("Cervical","Endometrial")) > 0){
      #  data1 <- data1 %>% filter(exclude_cervical == 0)#Binary feature (0: Keep, 1:exclude)
      #}else if(outcome[i] %in% c("Ovarian")){
      #  data1 <- data1 %>% filter(exclude_ovarian == 0)#Binary feature (0: Keep, 1:exclude)
      #}
      form <- as.formula(paste0("Surv(time,",outcome[i],")~",exposure[j], "+", paste0(covariates, collapse = "+")))
      cox <- coxph(form,data = data1,ties = "efron",control = crt)
      #No interactions by sex in sex-specific level
      coxInter <- cox
    }else{
      form <- as.formula(paste0("Surv(time,",outcome[i],")~",exposure[j], "+", paste0(covariates, collapse = "+")))
      cox <- coxph(form,data = data,ties = "efron",control = crt)
      #Interaction Model by Sex
      formInter <- as.formula(paste0("Surv(time,",outcome[i],")~",exposure[j],"*sex+", paste0(covariates, collapse = "+")))
      coxInter <- coxph(formInter,data = data,ties = "efron",control = crt)
    }
    
    #Results - Without interaction by sex
    row <- grepl(exposure[j],rownames(coef(summary(cox))),fixed = T)
    coefs <- apply(coef(summary(cox)), 2, round, 2)
    coefs_intervals <- apply(exp(confint(cox)), 2, round,2)
    #Overall
    results[[j]][i,2] <- paste0(coefs[row,"exp(coef)"],"(", coefs_intervals[row,"2.5 %"],",", coefs_intervals[row,"97.5 %"],")")
    results[[j]][i,3] <- round(coefs[row,"Pr(>|z|)"],3)
    
    #Results - With interaction by sex
    row_inter <- grepl(exposure[j],rownames(coef(summary(coxInter))),fixed = T)
    coefs_inter <- apply(coef(summary(coxInter)), 2, round, 2)
    coefs_intervals_inter <- apply(exp(confint(coxInter)), 2, round,2)
    #Females
    results[[j]][i,4] <- ifelse(outcome[i] %in% "Prostate",NA, paste0(coefs_inter[,"exp(coef)"][1],"(", coefs_intervals_inter[row_inter,"2.5 %"][1],",", coefs_intervals_inter[row_inter,"97.5 %"][1],")"))
    results[[j]][i,5] <- ifelse(outcome[i] %in% "Prostate",NA, round(coefs_inter[row_inter,"Pr(>|z|)"][1],3))
    #Males
    results[[j]][i,6] <- ifelse(outcome[i] %in% c("Breast","Cervical","Endometrial","Ovarian") | sum(is.na(coef(coxInter)[row_inter])) > 0,NA,#Check if is NA the effect (due to zero cases in those categories)
                                paste0(round(lincom(coxInter,paste(rownames(coefs_inter)[row_inter],collapse = '+'), eform = T,singular.ok = T)["Estimate"][[1]],2),"(",
                                      round(lincom(coxInter,paste(rownames(coefs_inter)[row_inter],collapse = '+'), eform = T,singular.ok = T)["2.5 %"],2),",",
                                      round(lincom(coxInter,paste(rownames(coefs_inter)[row_inter],collapse = '+'), eform = T,singular.ok = T)["97.5 %"],2),")"))
    results[[j]][i,7] <- ifelse(outcome[i] %in% c("Breast","Cervical","Endometrial","Ovarian") |
                                  sum(is.na(coef(coxInter)[row_inter])) > 0,NA,#Check if is NA the effect (due to zero cases in those categories)
                                round(lincom(coxInter,paste(rownames(coefs_inter)[row_inter],collapse = '+'), eform = T,singular.ok = T)["Pr(>Chisq)"],3)[[1]])
    #P-Value - Interaction term
    results[[j]][i,8] <- ifelse(outcome[i] %in% c("Breast","Cervical","Endometrial","Ovarian","Prostate"),NA,
                                coefs_inter[grepl(":sex",rownames(coefs_inter),fixed = T),"Pr(>|z|)"])
    
    #Proportional hazard
    zph <- cox.zph(coxInter, transform = "identity")$table
    results[[j]][i,9] <- round(zph[which(rownames(zph) %in% exposure[j]),"p"],3)

  }
  results[[j]] <- data.frame(cbind(results[[j]],exposure = ""))
  colnames(results[[j]])[length(results[[j]])] <- abr[j]
}

unajusted_results <- results
for(i in 1:length(results)){
  k <- c(which(grepl("PVal", names(results[[i]]))),length(results[[i]]))
  for(j in k){
    results[[i]][,j] <- ifelse(p.adjust(as.numeric(results[[i]][,j]),method = "fdr") < 0.05,paste0(results[[i]][,j],"*"),results[[i]][,j])
  }
  
  if(!dir.exists("../Results")){ dir.create("../Results")}
  write.xlsx(results[[i]], paste0("../Results/", exposure[i],"_HR.xlsx"), row.names = F, )
}
