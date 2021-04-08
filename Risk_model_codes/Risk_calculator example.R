source(paste0(getwd(), "//PLCOm2012.R")) # load PLCOm2012 model code
source(paste0(getwd(), "//PLCOm2012noRace.R")) # load PLCOm2012noRace model code
source(paste0(getwd(), "//LLP2.R")) # load LLPv2 model code
source(paste0(getwd(), "//LLP3.R")) # load LLPv3 model code
source(paste0(getwd(), "//BachLCDiag_adapted.R")) # load Bach model code for lung cancer diagnosis
source(paste0(getwd(), "//BachLCDeath_adapted.R")) # load Bach model code for other-cause death
source(paste0(getwd(), "//Bachmultiyear.R")) # load Bach model code for iterative predictions
library(lcrisks) #Load LCrisks package (NB: have to load this manually; problems with R versions 4+)



#Example persons: 55-year old high school graduated White male, current smoker, who smoked 15 cigarettes per day for 40 years,
#has a BMI of 23, no COPD, no asbestos exposure, no personal history of cancer, no personal history of pneumonia,
#and no family history of lung cancer 
Age=55 #Continuous value for age at the start of hte prediction
Gender = 0 #Gender; 0 represents men, 1 represents women
Education = 2 #Education: This variable has 6 levels: (1 = less than high school grad, 2 = high school grad, 3 = post high school training, 4 = some college, 5 = college grad, 6 = postgraduate/professional)
Race = 1 #Race/Ethnicity indicator. This variable has 6 levels: (1 = White, 2 = Black (non-Hispanic), 3 = Hispanic, 4 = Asian, 5 = American Indian/Alaskan Native, 6 = Native Hawaiian/Pacific Islander)
Smokingstatus = 1 #Binary indicator for smoking status (1= current smoker, 0= former)
CPD = 15 #Value for the average number of cigarettes smoked over the persons' smoking history
Yearssmoked = 40 #Value for number of years the person has smoked
Yearsquit = 0 #Value for the number of years since smoking cessation. This value should be set to 0 if the person is a current smoker
BMI = 23 #Value for BMI
COPD = 0 #Binary indicator for COPD/emphysema/chronic bronchitis
Asbestos = 0 #Binary indicator for asbestos exposure
PHcancer = 0 #Binary indicator for personal history of cancer
PHpneumonia = 0 #Binary indicator for personal history of cancer
FHlungcancer = 0 #Binary indicator for family history of lung cancer 


#Rebind variables for use in the LC(D)RAT models
Race2 = Race 
Race2 = if(Race ==1)
{Race2 =0 }
age = c(Age)
bmi = c(BMI)
cpd <- c(CPD)
emp <- c(COPD)
fam.lung.trend <- c(FHlungcancer )
female <- c(Gender )
smkyears <- c(Yearssmoked )
 
qtyears <- c(Yearsquit )
race <- c(Race2 )
edu6 <- c(Education )
years <- 6


persons=cbind(age,
female,
smkyears,
qtyears,
cpd,
race,
emp,
fam.lung.trend,
bmi,
edu6)

#Bach model: 6 years risk for lung cancer incidence 
Bach=Bachmultiyear(6,Smokingstatus,Age ,Gender,Yearssmoked,CPD ,Yearsquit,Asbestos)

#LCRAT and LCDRAT model risks (5-year risks)
personspredictions <- lcrisk(persons,years)
LCRAT=(personspredictions[1,5]/1000)
LCDRAT=(personspredictions[1,3]/1000)

#LLP model risks (5-year risks)
LLP2=model.LLP2(Age,Yearssmoked,FHlungcancer,PHcancer,PHpneumonia ,Asbestos ,Gender )
LLP3=model.LLP3(Age,Yearssmoked,FHlungcancer,PHcancer,PHpneumonia ,Asbestos ,Gender )

#PLCOm2012 model risk (6-year risk)
PLCOm2012=model.PLCOm2012(Age,Education,BMI,COPD,PHcancer,FHlungcancer,Race,Smokingstatus,CPD,Yearssmoked,Yearsquit)
PLCOm2012noRace=model.PLCOm2012NoRace(Age,Education,BMI,COPD,PHcancer,FHlungcancer,Smokingstatus,CPD,Yearssmoked,Yearsquit)  








