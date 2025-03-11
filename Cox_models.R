
#Libraries
library(tidyverse)
library(survival); library(survminer)
library(biostat3)#For linear combinations
library(plotRCS); library(ggplot2)
 
data <- read.csv("C:/Users/user/Desktop/Papers/Sarcopenia/Git/for_Git.csv")
outcome <- c("pre_sarcopenic","sarcopenic","muscle_mass_index","cat_muscle_mass_index","grip_strength","Cat_grip_strength","sarcopenic_obesity","walking_pace")
nnn <- c("Probable Sarcopenia","Sarcopenia","Muscle Mass Index", "Muscle Mass Index Categories", "Grip Strength", "Grip Strength Categories", "Sarcopenic Obesity","Walking Pace")
check_cancer <- c("Bladder","Brain","Breast","Colorectal","Endometrial","Gastric","Hematologic",
                  "Kidney","Liver","Lung","Melanoma","Oesophageal","Oral",
                  "Ovarian","Overall","Pancreatic","Prostate")
crt <- coxph.control(1e-11/2,iter.max = 50)

covariates = c("bmi","age","smoking","smoke_intensity","alcohol","physical_activity","family_history","TDI")

#### Cox model - Model 1 ####
results <- list()
for(j in 1:length(outcome)){
  print(outcome[j])
  
  results[[j]] <- matrix(NA, nrow = length(check_cancer),ncol = 9)
  colnames(results[[j]]) <- c("Cancer","Exposure","PVal","Female","PVal","Male","PVal","Interaction", "PH")
  results[[j]][,1] <- check_cancer
  
  for(i in 1:length(check_cancer)){
    print(check_cancer[i])
    
    if(check_cancer[i] %in% c("Prostate")){
      data1 <- data %>% filter(gender == 1)#Keep only males
      form <- as.formula(paste0("Surv(time,",check_cancer[i],")~",outcome[j], "+", paste0(covariates, collapse = "+")))
      cox <- coxph(form, data = data1,ties = "efron",control = crt)
      #No interactions by sex in sex-specific level
      coxInter <- cox
    }else if(check_cancer[i] %in% c("Breast","Cervical","Endometrial","Ovarian")){
      data1 <- data %>% filter(gender == 0)#keep only females
      if(sum(check_cancer[i] %in% c("Cervical","Endometrial")) > 0){
        data1 <- data1 %>% filter(exclude_cervical == 0)#Binary feature (0: Keep, 1:exclude)
      }else if(check_cancer[i] %in% c("Ovarian")){
        data1 <- data1 %>% filter(exclude_ovarian == 0)#Binary feature (0: Keep, 1:exclude)
      }
      form <- as.formula(paste0("Surv(time,",check_cancer[i],")~",outcome[j], "+", paste0(covariates, collapse = "+")))
      cox <- coxph(form,data = data1,ties = "efron",control = crt)
      #No interactions by sex in sex-specific level
      coxInter <- cox
    }else{
      form <- as.formula(paste0("Surv(time,",check_cancer[i],")~",outcome[j], "+", paste0(covariates, collapse = "+")))
      cox <- coxph(form,data = data,ties = "efron",control = crt)
      #Interaction Model by Sex
      formInter <- as.formula(paste0("Surv(time,",check_cancer[i],")~",outcome[j],"*sex+", paste0(covariates, collapse = "+")))
      coxInter <- coxph(formInter,data = data,ties = "efron",control = crt)
    }
    
    #Results - Without interaction by sex
    row <- grepl(outcome[j],rownames(coef(summary(cox))),fixed = T)
    coefs <- apply(coef(summary(cox)), 2, round, 2)
    coefs_intervals <- apply(exp(confint(cox)), 2, round,2)
    #Overall
    results[[j]][i,2] <- paste0(coefs[row,"exp(coef)"],"(", coefs_intervals[row,"2.5 %"],",", coefs_intervals[row,"97.5 %"],")")
    results[[j]][i,3] <- round(coefs[row,"Pr(>|z|)"],3)
    
    #Results - With interaction by sex
    row_inter <- grepl(outcome[j],rownames(coef(summary(coxInter))),fixed = T)
    coefs_inter <- apply(coef(summary(coxInter)), 2, round, 2)
    coefs_intervals_inter <- apply(exp(confint(coxInter)), 2, round,2)
    #Females
    results[[j]][i,4] <- ifelse(check_cancer[i] %in% "Prostate",NA, paste0(coefs_inter[,"exp(coef)"][1],"(", coefs_intervals_inter[row_inter,"2.5 %"][1],",", coefs_intervals_inter[row_inter,"97.5 %"][1],")"))
    results[[j]][i,5] <- ifelse(check_cancer[i] %in% "Prostate",NA, round(coefs_inter[row_inter,"Pr(>|z|)"][1],3))
    #Males
    results[[j]][i,6] <- ifelse(check_cancer[i] %in% c("Breast","Cervical","Endometrial","Ovarian") | sum(is.na(coef(coxInter)[row_inter])) > 0,NA,#Check if is NA the effect (due to zero cases in those categories)
                                paste0(round(lincom(coxInter,paste(rownames(coefs_inter)[row_inter],collapse = '+'), eform = T,singular.ok = T)["Estimate"][[1]],2),"(",
                                      round(lincom(coxInter,paste(rownames(coefs_inter)[row_inter],collapse = '+'), eform = T,singular.ok = T)["2.5 %"],2),",",
                                      round(lincom(coxInter,paste(rownames(coefs_inter)[row_inter],collapse = '+'), eform = T,singular.ok = T)["97.5 %"],2),")"))
    results[[j]][i,7] <- ifelse(check_cancer[i] %in% c("Breast","Cervical","Endometrial","Ovarian") |
                                  sum(is.na(coef(coxInter)[row_inter])) > 0,NA,#Check if is NA the effect (due to zero cases in those categories)
                                round(lincom(coxInter,paste(rownames(coefs_inter)[row_inter],collapse = '+'), eform = T,singular.ok = T)["Pr(>Chisq)"],3)[[1]])
    #P-Value - Interaction term
    results[[j]][i,8] <- ifelse(check_cancer[i] %in% c("Breast","Cervical","Endometrial","Ovarian","Prostate"),NA, coefs_inter[grepl(paste0(":sex1"),rownames(coefs_inter),fixed = T),"Pr(>|z|)"])

  }
  results[[j]] <- data.frame(cbind(results[[j]],outcome = ""))
  colnames(results[[j]])[length(results[[j]])] <- nnn[j]
}

unajusted_results <- results
for(i in 1:length(results)){
  k <- c(which(grepl("PVal", names(results[[i]]))),length(results[[i]]))
  for(j in k){
    results[[i]][,j] <- ifelse(p.adjust(as.numeric(results[[i]][,j]),method = "fdr") < 0.05,paste0(results[[i]][,j],"*"),results[[i]][,j])
  }
  write.csv(results[[i]], paste0(outcome[i],"_HR.csv"), row.names = F, na = "")
}
