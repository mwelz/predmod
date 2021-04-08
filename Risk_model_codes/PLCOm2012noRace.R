# PLCOm2012NoRace model 
# Author          Kevin ten Haaf
# Organization    Erasmus Medical Center Rotterdam
#Last adjusted: September 01 2020
model.PLCOm2012NoRace = function(Age,education,BMI,COPD,personalhistory,familyhistory,Currentsmoker,Smokingintensity,Smokingduration,yearsquitsmoking)  
{
#Input types: 
#Age: number 
#Education: number corresponding to education level: 1 = Less than high school grad, 2= High school grad, 3= Post high school training, 4= Some college, 5= College grad, 6= Postgraduate/professional
#BMI: number
#COPD: boolean  0 = no, 1 = yes
#personalhistory: boolean 0 = no, 1 = yes
#familyhistory: boolean  0 = no, 1 = yes
#Currentsmoker: boolean  0 = former smoker, 1 = current smoker
#Smokingintensity: number
#Smokingduration: number
#yearsquitsmoking: number 


#Coefficients corresponding to inputs: Age,education,BMI,COPD,personalhistory,familyhistory,race,Currentsmoker,Smokingintensity,Smokingduration,yearsquitsmoking
Coeffs=c(
0.0778895,
-0.0811569,
-0.0251066,
0.3606082,
0.4683545,
0.584541,
0.2675539,
-1.767578,
0.031949,
-0.0312719)


#First the center values for the variables are defined
agecentervalue = 62
educationcentervalue = 4
bmicentervalue = 27
CPDcentervalue =0.4021541613
Smokingdurationcentervalue = 27
Smokingcessationcentervalue = 10

#Then each model parameter's contribution is calculated
Modelconstant=-4.536696
Agecontribution = (Age-agecentervalue)*Coeffs[1]
Educationcontribution = (education-educationcentervalue)*Coeffs[2]
Bmicontribution = (BMI-bmicentervalue)*Coeffs[3]
Copdcontribution = COPD*Coeffs[4]
Personalhistorycontribution = personalhistory*Coeffs[5]
Familyhistorycontribution= familyhistory*Coeffs[6]
Smokingstatuscontribution= Currentsmoker*Coeffs[7]
CPDcontribution = (((Smokingintensity/10)^-1)-CPDcentervalue)*Coeffs[8]
Smokingdurationcontribution = (Smokingduration-Smokingdurationcentervalue )*Coeffs[9]
Smokingcessationcontribution = (yearsquitsmoking-Smokingcessationcentervalue)*Coeffs[10]


#The individual contributions are summed and the 6-year probability is returned
Sumvalues = Modelconstant+Agecontribution+Educationcontribution+Bmicontribution+Copdcontribution+Personalhistorycontribution+Familyhistorycontribution+Smokingstatuscontribution+CPDcontribution+Smokingdurationcontribution+Smokingcessationcontribution
 
Sixyearprobability = exp(Sumvalues)/(1+exp(Sumvalues))
return(Sixyearprobability)
}




model.PLCOm2012NoRace(63,1,36.01108,1,1,1,1,10,12,0)
#risk: 0.01827511


model.PLCOm2012NoRace(68,5,30.5241,0,1,1,1,20,33,0)
#0.06997445

model.PLCOm2012NoRace(61,5,34.25606,1,1,1,0,20,21,9)
#0.02192871

model.PLCOm2012NoRace(70,4,37.87401,0,0,0,1,20,30,0)
#0.0245205
