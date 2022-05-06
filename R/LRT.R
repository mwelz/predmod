#' Likelihood ratio test for multiple imputed results
#' Based on the D2 method of pooled Chi2 statistics (Li, Meng, et al. (1991)
#' If models without imputations are provided the general likelihood ratio test is applied
#' @param noninteract list of (imputed) models without an interaction effect 
#' @param interact list of (imputed) models with an interaction effect  
#' @param params_testeddegrees of freedom (generally 1) 
#' @return 
#' 
#'


LRT <- function(noninteract, interact, params_tested){
  nr_of_imputations = length(noninteract)
  
  if(nr_of_imputations = 1)
  {
    
    classtester = class(interact[[1]])
    classtest=classtester[1]
    if(classtest!='coxph')
    {  

        temp=lmtest::lrtest(interact[[1]], noninteract[[1]])
        p.value=temp$`Pr(>Chisq)`[2]

    }
    
    if(classtest=='coxph')
    {  

        temp=survival:::anova.coxph(interact[[1]], noninteract[[1]])
        p.value=temp$`P(>|Chi|)`[2]

    }
    
  
  return(list(p.value=p_value))
  }
  
  d=vector(length=nr_of_imputations)	
	classtester = class(interact[[1]])
  classtest=classtester[1]
  if(classtest!='coxph')
  {  
     for (i in 1:nr_of_imputations)
	    {
	    temp=lmtest::lrtest(interact[[i]], noninteract[[i]])
	    d[i]=temp$Chisq[2]
	    }
  }
  
  
  if(classtest=='coxph')
  {  
    for (i in 1:nr_of_imputations)
    {
      temp=survival:::anova.coxph(interact[[i]], noninteract[[i]])
      d[i]=temp$Chisq[2]
    }
  }

average_d = sum(d)/nr_of_imputations
root_d = sum(sqrt(d))/nr_of_imputations
sum_dstat = 0
	
	for(i in 1:nr_of_imputations)
	{
	sum_dstat=sum_dstat+((sqrt(d[i])-root_d)^2)
	}

R_two = (1+(1/nr_of_imputations))*((1/(nr_of_imputations-1))*sum_dstat)
D_two = (average_d - ((params_tested*(nr_of_imputations-1))/(nr_of_imputations+1)*R_two))  / (params_tested*(1+R_two))
v_two = (params_tested^(-3/nr_of_imputations)) * (nr_of_imputations-1) * (1+((R_two^(-1))^2))


	p_value=1-pf(D_two, params_tested, v_two)

  # return
  return(list(p.value=p_value))

} # FUN

