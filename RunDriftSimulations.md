# Run Drift Simulations

#### Cera Fisher

#### 9/13/2022

What follows is a boring R session of running the same function  multiple times. Repeating anything more than twice is absurd, but here  we are! The results from this simulation run are the results we report  in our paper. 

Results are given as in-line comments of the “empirical p”, which is  the frequency of simulations as extreme or more than the real results.  P-value significance threshold is 0.05; a value higher than means that  the observed change in R allele frequency is not significant and may be  attributed to drift.

(Please be aware that running these simulations takes a very long  time, and that the results will be slightly different each time, because this is an exercise in random sampling, just like drift; **do not attempt to rerun this code!**)

```R
genos <- read.table(file="../data/AlleleCompetition_GenotypeCounts.txt", header=TRUE, sep="\t")
```

Get r allele frequency for each row.

```R
genos <- genos %>% mutate(
  r.freq = ((RR*2) + SR)/(n*2)
)
```

Check structure of data table.

```R
str(genos)
## 'data.frame':    26 obs. of  9 variables:
##  $ Generation: chr  "F1" "F1" "F2" "F2" ...
##  $ Cross     : chr  "A" "B" "A" "A" ...
##  $ Replicate : chr  "A" "B" "A1" "A2" ...
##  $ Full.Name : chr  "F1A" "F1B" "F2A1" "F2A2" ...
##  $ n         : int  88 87 87 85 84 87 82 84 84 85 ...
##  $ SS        : int  0 0 19 25 22 26 24 25 25 21 ...
##  $ SR        : int  88 87 45 44 51 48 35 38 44 49 ...
##  $ RR        : int  0 0 23 16 11 13 23 21 15 15 ...
##  $ r.freq    : num  0.5 0.5 0.523 0.447 0.435 ...
source("MultiGenDrift.R")
```



### F1 to F2

```R
## Population A1 
nsims = 50000
gens = 1 ## From F1 to F2, only one generation happens. 
popsize = 800
Rstart = 0.5
```

#### F2A1

```R
Rend = genos$r.freq[genos$Full.Name=="F2A1"]

ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or larger is 0.0348
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.0348
Rend > Rstart
## [1] TRUE
table(ResultsDF$Result >= Rend)
## 
## FALSE  TRUE 
## 48260  1740
```

#### F2A2

```R
Rend = genos$r.freq[genos$Full.Name=="F2A2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 2e-05
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 2e-05
Rend < Rstart
## [1] TRUE
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49999     1
```

#### F2A3

```R
Rend = genos$r.freq[genos$Full.Name=="F2A3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F2B1

```R
Rend = genos$r.freq[genos$Full.Name=="F2B1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F2B2

```R
Rend = genos$r.freq[genos$Full.Name=="F2B2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.31594
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.31594
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 34203 15797
```

#### F2B3

```R
Rend = genos$r.freq[genos$Full.Name=="F2B3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.02796
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.02796
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 48602  1398
```

### F1 to F3

```R
## F1 to F3 is two generations
gens = 2
```

#### F3A1

```R
Rend = genos$r.freq[genos$Full.Name=="F3A1"]

ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.00024
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.00024
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49988    12
```

#### F3A2

```R
Rend = genos$r.freq[genos$Full.Name=="F3A2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.023
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.023
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 48850  1150
```

#### F3A3

```R
Rend = genos$r.freq[genos$Full.Name=="F3A3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.00056
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.00056
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49972    28
```

#### F3B1

```R
Rend = genos$r.freq[genos$Full.Name=="F3B1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F3B2

```R
Rend = genos$r.freq[genos$Full.Name=="F3B2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F3B3

```R
Rend = genos$r.freq[genos$Full.Name=="F3B3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.00038
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.00038
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49981    19
```

### F1 to F5

```
gens = 4
### F1 to F5 is four generations
```

#### F5A1

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F5A1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F5A2

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F5A2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.00036
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.00036
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49982    18
```

#### F5A3

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F5A3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F5B1

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F5B1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F5B2

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F5B2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F5B3

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F5B3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

### F1 to F7

```
gens = 6
## F1s to F7 is 6 generations
```

#### F7A1

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F7A1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F7A2

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F7A2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.00054
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.00054
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49973    27
```

#### F7A3

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F7A3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F7B1

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F7B1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F7B2

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F7B2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result >= Rend)
## 
##  TRUE 
## 50000
```

#### F7B3

```R
Rstart = 0.5
Rend = genos$r.freq[genos$Full.Name=="F7B3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

### F2 to F3

```R
gens = 1
#### F2 to F3 is 1 generation
```

#### F3A1

```R
Rstart = genos$r.freq[genos$Full.Name=="F2A1"]
Rend = genos$r.freq[genos$Full.Name=="F3A1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F3A2

```R
Rstart = genos$r.freq[genos$Full.Name=="F2A2"]
Rend = genos$r.freq[genos$Full.Name=="F3A2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or larger is 0.07964
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.07964
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
##  3982 46018
```

#### F3A3

```R
Rstart = genos$r.freq[genos$Full.Name=="F2A3"]
Rend = genos$r.freq[genos$Full.Name=="F3A3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or larger is 0.25732
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.25732
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 12866 37134
```

#### F3B1

```R
Rstart = genos$r.freq[genos$Full.Name=="F2B1"]
Rend = genos$r.freq[genos$Full.Name=="F3B1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F3B2

```R
Rstart = genos$r.freq[genos$Full.Name=="F2B2"]
Rend = genos$r.freq[genos$Full.Name=="F3B2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F3B3

```R
Rstart = genos$r.freq[genos$Full.Name=="F2B3"]
Rend = genos$r.freq[genos$Full.Name=="F3B3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.00272
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.00272
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49864   136
```

### F3s to F5s

```R
gens = 2
### F3s to F5s is 2 generations
```

#### F5A1

```R
Rstart = genos$r.freq[genos$Full.Name=="F3A1"]
Rend = genos$r.freq[genos$Full.Name=="F5A1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F5A2

```R
Rstart = genos$r.freq[genos$Full.Name=="F3A2"]
Rend = genos$r.freq[genos$Full.Name=="F5A2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.00284
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.00284
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49858   142
```

#### F5A3

```R
Rstart = genos$r.freq[genos$Full.Name=="F3A3"]
Rend = genos$r.freq[genos$Full.Name=="F5A3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F5B1

One of the few times that the R allele goes up between two generations

```R
Rstart = genos$r.freq[genos$Full.Name=="F3B1"]
Rend = genos$r.freq[genos$Full.Name=="F5B1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or larger is 0.08274
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.08274
Rend > Rstart
## [1] TRUE
table(ResultsDF$Result >= Rend)
## 
## FALSE  TRUE 
## 45863  4137
```

#### F5B2

```R
Rstart = genos$r.freq[genos$Full.Name=="F3B2"]
Rend = genos$r.freq[genos$Full.Name=="F5B2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result >= Rend)
## 
##  TRUE 
## 50000
```

#### F5B3

```R
Rstart = genos$r.freq[genos$Full.Name=="F3B3"]
Rend = genos$r.freq[genos$Full.Name=="F5B3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 4e-05
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 4e-05
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49998     2
```

### F5 to F7

```R
gens = 2
## F5 to F7 is two generations
```

#### F7A1

Here, again, R allele goes up instead of down

```R
Rstart = genos$r.freq[genos$Full.Name=="F5A1"]
Rend = genos$r.freq[genos$Full.Name=="F7A1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or larger is 0.14512
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.14512
Rend > Rstart
## [1] TRUE
table(ResultsDF$Result >= Rend)
## 
## FALSE  TRUE 
## 42744  7256
```

#### F7A2

```R
Rstart = genos$r.freq[genos$Full.Name=="F5A2"]
Rend = genos$r.freq[genos$Full.Name=="F7A2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.18576
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.18576
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 40712  9288
```

#### F7A3

```R
Rstart = genos$r.freq[genos$Full.Name=="F5A3"]
Rend = genos$r.freq[genos$Full.Name=="F7A3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F7B1

```R
Rstart = genos$r.freq[genos$Full.Name=="F5B1"]
Rend = genos$r.freq[genos$Full.Name=="F7B1"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0.43002
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0.43002
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 28499 21501
```

#### F7B2

```R
Rstart = genos$r.freq[genos$Full.Name=="F5B2"]
Rend = genos$r.freq[genos$Full.Name=="F7B2"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 0
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 0
table(ResultsDF$Result <= Rend)
## 
## FALSE 
## 50000
```

#### F7B3

```R
Rstart = genos$r.freq[genos$Full.Name=="F5B3"]
Rend = genos$r.freq[genos$Full.Name=="F7B3"]
ResultsList <- MultiGenDrift(SIMS = nsims, gens = gens, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)
## probability of getting end R freq or smaller is 2e-04
empirical_p <- ResultsList[[2]]
ResultsDF <- data.frame(ResultsList[[1]])
colnames(ResultsDF) <- 'Result'
empirical_p
## [1] 2e-04
table(ResultsDF$Result <= Rend)
## 
## FALSE  TRUE 
## 49990    10
```