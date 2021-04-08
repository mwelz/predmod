# Author          Rafael Meza
# Organization    University of Michigan
# Last updated    August 15, 2013

#Minor adaptations to original code by Kevin ten Haaf
# Organization    Erasmus Medical Center Rotterdam
# Last updated    March 28, 2014

 BachLCDeath<-function(AGE,SEX,SMK,CPD,QUIT,ASB)
{  
   # function that computes the probability of death in the absence lung cancer in the next year
   # given age and several risk factors.
   # Based on the model by Bach el at, JNCI 2003 (Variations in Lung Cancer Risk Among Smokers) 
   # age: continuous, positive between 50-75
   # sex: gender ('M', 'F') M = 0, F = 1
   # smk: smoking duration. continuous, positive between 25-55
   # cpd : continuous, positive between 10-60
   # quit: continuous, years since quitted smoking
   # asb: asbestos exposure ('Y','N') Y = 1, N =0
   # Example: BachLCDeath(70,0,40,38.86577601,6,0)
   
  S0=0.9917663
  model = -7.2036219 + (0.015490665 * CPD)
  
  if( CPD>15) { model = model -(0.00001737645 * (CPD -15)^3)} #for all values CPD>15
  if( CPD>20) { model = model + (0.000021924149 * (CPD - 20.185718)^3)} #for all values CPD>20 between 20 and 40 or for all values greator than 20
  if( CPD>40 ) { model = model- (0.0000045476985 * (CPD - 40)^3)} #for all values CDP>40
  
  model=model+ (0.020041889 * SMK) 
  if( SMK>27) {model = model + (0.0000065443781 * (SMK - 27.6577)^3)}#for all values SMK>27
  if( SMK>40) { model = model - (0.000013947696 * (SMK - 40)^3)} #for all values SMK>40
  if( SMK>50) { model = model + (0.0000074033175 * (SMK-50.910335)^3)} #for all values SMK>5
  
  model = model - (0.023358962 * QUIT) + (0.0019208669 * QUIT^3) #for all values
  
  if( QUIT>0) { model = model - (0.0020031611 * (QUIT - 0.50513347)^3)} #for all values QUIT>0
  if( QUIT>12) { model=model + (0.000082294194 * (QUIT - 12.295688)^3)} #for all values QUIT>12
  
  model= model +  (0.099168033 * AGE) 
  if( AGE>53) { model = model  + (0.0000062174577 * (AGE - 53.459001)^3)} #for all values AGE>53
  if( AGE>61) { model =model- (0.000012115774 * (AGE - 61.954825) ^3)} #for all values AGE>61
  if( AGE>70) { model=model+(0.0000058983164 * (AGE - 70.910335) ^3)} #for all values AGE>70
  
  
  if( ASB==1) {model = model + (0.06084611)}
#Sex == 1 is female, sex == 0 is male
  if( SEX==1) {model=model- (0.49042298)}
  
  pDeath=1-S0^exp(model)
  return(pDeath)
}