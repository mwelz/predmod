#Bach multiyear probability calculator
# Author          Kevin ten Haaf
# Organization    Erasmus Medical Center Rotterdam
# Last updated    April 8th, 2014


#This function calculates the risk of developing lung cancer over a pre-specified number of years using the 
#Bach incidence and mortality (in the absence of lung cancer) models iteratively



Bachmultiyear = function(Numberofyears,Currentsmokingstatus, AGE,SEX,SMK,CPD,QUIT,ASB)
{

   # Numberofyears : number of years for which the prospective probability should be calculated 
   #Currentsmokingstatus: The person's current smokingstatus, 1 if the person is currently smoking, 0 otherwise
   # age: age in years
   # sex: gender M = 0, F = 1
   # smk: smoking duration in years
   # cpd : number of cigarettes per day
   # quit: years since quitted smoking
   # asb: asbestos exposure  Y = 1, N =0
 

#A number of initial persons is generated
Startingnrpersons = 1000000

#The number of non-dead/non-lung cancer diagnosed persons at each year
Healthypersons=0

#The total number of persons diagnosed with lung cancer 
LCpersons=0

#The total number of persons who died in the absence of a lung cancer diagnosis
nonLCDeads=0

#The number of persons newly diagnosed with lung cancer 
newLCpersons=0

#The number of persons newly died in the absence of a lung cancer diagnosis
newnonLCDeads=0

#The transitions for the first year are generated, first the number of persons diagnosed with lung cancer:
LCpersons=LCpersons+( BachLCDiag(AGE,SEX,SMK,CPD,QUIT,ASB) * Startingnrpersons )
Healthypersons=Startingnrpersons -LCpersons

#Then the number of persons who die in the absence of lung cancer are generated
nonLCDeads=nonLCDeads+(Healthypersons*BachLCDeath(AGE,SEX,SMK,CPD,QUIT,ASB))
Healthypersons=Healthypersons-nonLCDeads

#This process is repeated until the number of years for which the probability should be calculated has been reached,
#updating the age, smoking duration and smoking cessation variables where necessary.
for (x in 2:Numberofyears)
{
AGE = AGE+1
if (Currentsmokingstatus == 0)
{
QUIT =QUIT+1
}

if (Currentsmokingstatus == 1)
{
SMK =SMK +1
}
newLCpersons=BachLCDiag(AGE,SEX,SMK,CPD,QUIT,ASB) * Healthypersons
LCpersons=LCpersons+newLCpersons
Healthypersons=Healthypersons-newLCpersons

newnonLCDeads=Healthypersons*BachLCDeath(AGE,SEX,SMK,CPD,QUIT,ASB)
nonLCDeads=nonLCDeads+newnonLCDeads
Healthypersons=Healthypersons-newnonLCDeads

}

#The probability to be diagnosed in Numberofyears is the number of persons who have been diagnosed with lung cancer
#divided by the number of starting persons.

Probdiagnosedinxyears = LCpersons/Startingnrpersons 
return(Probdiagnosedinxyears)


}







