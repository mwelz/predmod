# LLP model V2 (based on parameters provided by John Field)
# Author          Kevin ten Haaf
# Organization    Erasmus Medical Center Rotterdam
#Last adjusted: June 19 2017


#Note: Family history: 0 is no family history, 1 is early onset (age <60), 2 is late onset (age>=60)
model.LLP2 = function(Age,Smokingduration,Familyhistory,Personalhistory,Pneumonia,Asbestos,Gender)  
{
  #Assume Age is Age+0.5
  Agevalue= Age+0.5
  
  
  
  #Counter for number of years
  Agecounter = 5
  
  #Age coefficients for 40-44-, 45-49, <50-54, 55-59, 60-64, 65-69, 70-74, 75-79 and 80-84+
  Maleagecoeffs = c(-9.06,-8.16,-7.31,-6.63,-5.97,-5.56,-5.31,-4.83,-4.68)
  Femaleagecoeffs= c(-9.90,-8.06,-7.46, -6.50, -6.22, -5.99, -5.49,-5.23,-5.42)
 
  Time44min = 0
  Time4549 = 0
  Time5054 = 0
  Time5559 = 0
  Time6064 = 0
  Time6569 = 0
  Time7074 = 0
  Time7579 = 0
  Time80plus = 0
  
  if (Gender == 0)
  {
    if ((Agevalue<45) & (Agecounter>0)) 
      {
      Time44min = min(max((45-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time44min 
      Agecounter = Agecounter-Time44min 
    }
    if ((Agevalue>=45) &   (Agevalue<50)      &   (Agecounter>0)) 
    {
      Time4549 = min(max((50-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time4549 
      Agecounter = Agecounter-Time4549 
    }


    if ((Agevalue>=50) &   (Agevalue<55)      &   (Agecounter>0)) 
    {
      Time5054 = min(max((55-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time5054 
      Agecounter = Agecounter-Time5054 
    }

    if ((Agevalue>=55) &   (Agevalue<60)      &   (Agecounter>0)) 
    {
      Time5559 = min(max((60-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time5559
      Agecounter = Agecounter-Time5559
    }
    if ((Agevalue>=60) &   (Agevalue<65)      &   (Agecounter>0)) 
    {
      Time6064 = min(max((65-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time6064
      Agecounter = Agecounter-Time6064
    }
    if ((Agevalue>=65) &   (Agevalue<70)      &   (Agecounter>0)) 
    {
      Time6569 = min(max((70-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time6569
      Agecounter = Agecounter-Time6569
    }
    if ((Agevalue>=70) &   (Agevalue<75)      &   (Agecounter>0)) 
    {
      Time7074 = min(max((75-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time7074
      Agecounter = Agecounter-Time7074
    }
    if ((Agevalue>=75) &   (Agevalue<80)      &   (Agecounter>0)) 
    {
      Time7579 = min(max((80-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time7579
      Agecounter = Agecounter-Time7579
    }    
        if ((Agevalue>=80) & (Agecounter>0)) 
    {
      Time80plus= Agecounter
      Agevalue = Agevalue+Time80plus
      Agecounter = Agecounter-Agecounter
    }    
    
    
    
    Agecontribution = ((Time44min*Maleagecoeffs[1])+(Time4549*Maleagecoeffs[2])+(Time5054*Maleagecoeffs[3])+(Time5559*Maleagecoeffs[4])+(Time6064*Maleagecoeffs[5])+(Time6569*Maleagecoeffs[6])+(Time7074*Maleagecoeffs[7])+(Time7579*Maleagecoeffs[8])+(Time80plus*Maleagecoeffs[9]))/5
    
  }
  

  
  
  if (Gender == 1)
  {
  
       if ((Agevalue<45) & (Agecounter>0)) 
      {
      Time44min = min(max((45-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time44min 
      Agecounter = Agecounter-Time44min 
    }
    if ((Agevalue>=45) &   (Agevalue<50)      &   (Agecounter>0)) 
    {
      Time4549 = min(max((50-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time4549 
      Agecounter = Agecounter-Time4549 
    }


    if ((Agevalue>=50) &   (Agevalue<55)      &   (Agecounter>0)) 
    {
      Time5054 = min(max((55-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time5054 
      Agecounter = Agecounter-Time5054 
    }

    if ((Agevalue>=55) &   (Agevalue<60)      &   (Agecounter>0)) 
    {
      Time5559 = min(max((60-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time5559
      Agecounter = Agecounter-Time5559
    }
    if ((Agevalue>=60) &   (Agevalue<65)      &   (Agecounter>0)) 
    {
      Time6064 = min(max((65-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time6064
      Agecounter = Agecounter-Time6064
    }
    if ((Agevalue>=65) &   (Agevalue<70)      &   (Agecounter>0)) 
    {
      Time6569 = min(max((70-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time6569
      Agecounter = Agecounter-Time6569
    }
    if ((Agevalue>=70) &   (Agevalue<75)      &   (Agecounter>0)) 
    {
      Time7074 = min(max((75-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time7074
      Agecounter = Agecounter-Time7074
    }
    if ((Agevalue>=75) &   (Agevalue<80)      &   (Agecounter>0)) 
    {
      Time7579 = min(max((80-Agevalue),0),Agecounter)
      Agevalue = Agevalue+Time7579
      Agecounter = Agecounter-Time7579
    }    
        if ((Agevalue>=80) & (Agecounter>0)) 
    {
      Time80plus= Agecounter
      Agevalue = Agevalue+Time80plus
      Agecounter = Agecounter-Agecounter
    }    
    
    
    
    Agecontribution = ((Time44min*Femaleagecoeffs[1])+(Time4549*Femaleagecoeffs[2])+(Time5054*Femaleagecoeffs[3])+(Time5559*Femaleagecoeffs[4])+(Time6064*Femaleagecoeffs[5])+(Time6569*Femaleagecoeffs[6])+(Time7074*Femaleagecoeffs[7])+(Time7579*Femaleagecoeffs[8])+(Time80plus*Femaleagecoeffs[9]))/5
     
    
  }
  
  
  
  #Duration coefficients #never, 1-19 years, 20-39, 40-59, 60+
  Smokingdurationcoeffs = c(0,0.7692,1.4516,2.5072,2.7243)
  if (Smokingduration==0) {Smokingdurationindicator=1}
  if ((Smokingduration<20) & (Smokingduration>0)) {Smokingdurationindicator=2}
  if ((Smokingduration<40) & (Smokingduration>19)) {Smokingdurationindicator=3}
  if ((Smokingduration<60) & (Smokingduration>39)) {Smokingdurationindicator=4}
  if (Smokingduration>59) {Smokingdurationindicator=5} 
  
  
  #Family history contributions: never, early onset (<60), late onset (60+)
  Familyhistorycoeffs = c(0,0.7034,0.1677)
  

  Smokingdurationcontribution = Smokingdurationcoeffs[Smokingdurationindicator]
  Pneumoniacontribution = Pneumonia*0.6025
  Asbestoscontribution = Asbestos*0.6343
  Personalhistoryconstribution = Personalhistory*0.6754
  Familyhistoryconstribution = Familyhistorycoeffs[Familyhistory+1]
  
  
  Total = Agecontribution+Smokingdurationcontribution+Pneumoniacontribution+Asbestoscontribution+Personalhistoryconstribution+Familyhistoryconstribution
  
  probability=1/(1+exp(-Total))
  
  return(probability)
}


