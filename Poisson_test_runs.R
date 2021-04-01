source(paste0(getwd(), "//Poissontests.R")) 

#test runs. Startinglifeyears is the base number of lifeyears, diseasereduction is the 1-% reduction in lifeyears due to the disease, relative_effect is the relative reduction due to treatment effect

Basecase = Poisson_runs(startinglifeyears=5, diseasereduction = 0.9,relative_effect=0.7)
Basecase$generated.ate  #-0.1621716
Basecase$generated.rel  #0.7
Basecase$risk.ate       #0.1637023
Basecase$risk.rel       #0.6981297
Basecase$effect.ate     #0.1765962
Basecase$effect.rel     #0.6789596
Basecase$generated.reduction.rate.1000 #36.24553
Basecase$estimated.reduction.effect.rate.1000 #  37.14827
Basecase$estimated.reduction.risk.rate.1000 # 34.43595

Higher_start_years = Poisson_runs(startinglifeyears=20, diseasereduction = 0.9,relative_effect=0.7) 
Higher_start_years$generated.ate  # -0.1621716
Higher_start_years$generated.rel  #0.7
Higher_start_years$risk.ate       #0.1666739
Higher_start_years$risk.rel       #0.6936545
Higher_start_years$effect.ate     #0.1764739
Higher_start_years$effect.rel     #0.6791385
Higher_start_years$generated.reduction.rate.1000 #8.850901
Higher_start_years$estimated.reduction.effect.rate.1000 # 9.256348
Higher_start_years$estimated.reduction.risk.rate.1000 # 8.742323

Greater_reduction = Poisson_runs(startinglifeyears=5, diseasereduction = 0.7,relative_effect=0.7)
Greater_reduction$generated.ate  #-0.1621716
Greater_reduction$generated.rel  #0.7
Greater_reduction$risk.ate       #0.1645277
Greater_reduction$risk.rel       #0.6968831
Greater_reduction$effect.ate     #0.1985196
Greater_reduction$effect.rel     #0.647824
Greater_reduction$generated.reduction.rate.1000 #45.00641
Greater_reduction$estimated.reduction.effect.rate.1000 # 46.57756
Greater_reduction$estimated.reduction.risk.rate.1000 # 38.60224


Greater_relative = Poisson_runs(startinglifeyears=5, diseasereduction = 0.9,relative_effect=0.4)
Greater_relative$generated.ate  #-0.3243431
Greater_relative$generated.rel  #0.4
Greater_relative$risk.ate       #0.3260974
Greater_relative$risk.rel       # 0.3984008
Greater_relative$effect.ate     #0.3444715
Greater_relative$effect.rel     #0.3798701
Greater_relative$generated.reduction.rate.1000 #70.97633
Greater_relative$estimated.reduction.effect.rate.1000 # 71.78776
Greater_relative$estimated.reduction.risk.rate.1000 # 67.95858



Greater_relative_and_disease = Poisson_runs(startinglifeyears=5, diseasereduction = 0.7,relative_effect=0.4)
Greater_relative_and_disease$generated.ate  #-0.3243431
Greater_relative_and_disease$generated.rel  #0.4
Greater_relative_and_disease$risk.ate       #0.3265975
Greater_relative_and_disease$risk.rel       #0.3978791
Greater_relative_and_disease$effect.ate     #0.380554
Greater_relative_and_disease$effect.rel     #0.3470335
Greater_relative_and_disease$generated.reduction.rate.1000 #84.87121
Greater_relative_and_disease$estimated.reduction.effect.rate.1000 #86.56662
Greater_relative_and_disease$estimated.reduction.risk.rate.1000 #74.29285



