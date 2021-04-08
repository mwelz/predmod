# PLCOm2012 model (Tammemagi,NEJM,2013)
# Author          Kevin ten Haaf
# Organization    Erasmus Medical Center Rotterdam
#Last adjusted: December 3 2014
model.PLCOm2012 = function(currentage,education,bmi,copd,personalhistory,familyhistory,race,currentsmokingstatus,averageCPD,currentsmokingduration,currentyearsquit)  
{
#coeffs: age, education, bmi, copd, personal history, family history, smoking status, cpd, smoking duration, years since cessation
Coeffs=c(0.0778868,-0.0812744,-0.0274194,0.3553063,0.4589971,0.587185,0.2597431,-1.822606,0.0317321,-0.0308572)

#coeffs: white, black, hispanic, asian, native pacific, 
Racecoeffs=c(0,0.3944778,-0.734744,-0.466585,0,1.027152)


#First the center values for the variables are defined
agecentervalue = 62
educationcentervalue = 4
bmicentervalue = 27
CPDcentervalue =0.4021541613
Smokingdurationcentervalue = 27
Smokingcessationcentervalue = 10

#Then each model parameter's contribution is calculated
Modelconstant=-4.532506
Agecontribution = (currentage-agecentervalue)*Coeffs[1]
Educationcontribution = (education-educationcentervalue)*Coeffs[2]
Bmicontribution = (bmi-bmicentervalue)*Coeffs[3]
Copdcontribution = copd*Coeffs[4]
Personalhistorycontribution = personalhistory*Coeffs[5]
Familyhistorycontribution= familyhistory*Coeffs[6]
Smokingstatuscontribution= currentsmokingstatus*Coeffs[7]
CPDcontribution = (((averageCPD/10)^-1)-CPDcentervalue)*Coeffs[8]
Smokingdurationcontribution = (currentsmokingduration-Smokingdurationcentervalue )*Coeffs[9]
Smokingcessationcontribution = (currentyearsquit-Smokingcessationcentervalue)*Coeffs[10]
Racecontribution = Racecoeffs[race]

#The individual contributions are summed and the 6-year probability is returned
Sumvalues = Modelconstant+Agecontribution+Educationcontribution+Bmicontribution+Copdcontribution+Personalhistorycontribution+Familyhistorycontribution+Smokingstatuscontribution+CPDcontribution+Smokingdurationcontribution+Smokingcessationcontribution+Racecontribution    
 
Sixyearprobability = exp(Sumvalues)/(1+exp(Sumvalues))
return(Sixyearprobability)
}

#Example 1 70 year old, college grad, bmi of 28, no copd, no personal history, family history, white, former smoker smoked 15 cpd for 30 years, quit 20 years ago
#model.PLCOm2012(70,5,28,0,0,1,1,0,15,30,20)




