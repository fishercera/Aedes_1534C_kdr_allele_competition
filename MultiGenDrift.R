MultiGenDrift<-function(SIMS, gens, size_of_pop, R_start_freq, R_end_freq)
{
  result<-NULL;#declare results vector
  library(plyr)
  for (i in 1:SIMS)#Simulations
  {	
    #declare the results vector to track allele frequencies each gen
    drift_result<-NULL; 
    gen<-gens; #no. of generations
    pop_size<-size_of_pop; #size_of_pop
    
    #starting with the input R allele frequency
    R_alleles<-R_start_freq;#R_start_freq;
    #cat("R allele frequency at start of gen 1 is", R_alleles,'\n');
    pop_alleles<-c(rep('R',2*pop_size*R_alleles),rep('S',2*pop_size*(1-R_alleles)));
    #take a sample from the whole population to make next generation
    pop_alleles<-sample(pop_alleles,2*pop_size,replace=TRUE);
    R_alleles<-(length(pop_alleles[pop_alleles=='R']))/(2*pop_size);#freq of R allele
    drift_result<-c(drift_result, R_alleles); #vector with R allele freq added 
    
    #each loop, R freq ending from previous loop used as the start R freq
    for(i in 2:gen)
    {
      R_alleles<-R_alleles; #uses new R freq
      # cat("R allele frequency at start of gen",i, "is", R_alleles,'\n');
      pop_alleles<-c(rep('R',2*pop_size*R_alleles),rep('S',2*pop_size*(1-R_alleles)));
      #take a sample from the whole population to make next generation
      pop_alleles<-sample(pop_alleles,2*pop_size,replace=TRUE);
      R_alleles<-(length(pop_alleles[pop_alleles=='R']))/(2*pop_size);#freq of R allele
      drift_result<-c(drift_result, R_alleles);
      #cat(drift_result, '\n');
    }
    
    drift_result<-as.vector(drift_result);
    #cat("R allele frequency at end of each generation is", drift_result, '\n');	
    
    #grabs the frequency of R at generation of interest from vector after loops
    #are completed
    res<-drift_result[gen];
    #cat(res, '\ ');
    result<-c(result, res);
    #cat(result, '\ ');
    
  }	
  
  #makes a new vector of the values of R at generation of interest for all simulations
  result<-as.vector(result);
 # cat("Generation", gen, "R allele frequency of simulations are", '\n', result, '\n');
  #hist(result)
  
  
  #counting number of times R allele frequency at generation of interest is
  #equal to or more extreme than observed frequency. need to designate which
  #direction of the extreme values is wanted
  empirical_p = NA
  if(R_end_freq < R_start_freq)
    ## Cera edits to make this more legible
    ## Also, we want to keep result vector to work with it later, so 
    ## we don't want to keep saving new values over it. 
  {
    MatchingResults <- result[result <= R_end_freq]
    # That is, the drift simulation results that are smaller than the end freq
    # (which is smaller than the start freq, allele freq decreased)
    NumberMatchingResults <- length(MatchingResults)
    # How many of such simulation results were there?
    #result<-result/SIMS;
    ProbOfSmallerResults <- NumberMatchingResults / SIMS
    empirical_p <- ProbOfSmallerResults
    # The frequency of smaller simulated results = the probability of the results
    # sort of, it would need an additional test
    cat("probability of getting end R freq or smaller is", ProbOfSmallerResults, '\n');
  }else{ # r_end_freq > r_start_freq
    MatchingResults <- result[result >= R_end_freq]
    # That is, the drift simulation results that are larger than the end freq
    # (which is smaller than the start freq, allele freq decreased)
    NumberMatchingResults <- length(MatchingResults)
    # How many of such simulation results were there?
    #result<-result/SIMS;
    ProbOfLargerResults <- NumberMatchingResults / SIMS
    empirical_p <- ProbOfLargerResults
    # The frequency of larger simulated results = the probability of the results
    # sort of, it would need an additional test
    cat("probability of getting end R freq or larger is", ProbOfLargerResults, '\n');
  }
  
  return(list(result, empirical_p))  
}
# 
# ### trouble-shooting ###
# nsims = 100
# popsize = 800
# gens = 2
# Rstart=0.5
# Rend=0.47
# 
# list <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
# SIMS = nsims
# size_of_pop=popsize
# gens=gens
# R_start_freq=Rstart
# R_end_freq=Rend
