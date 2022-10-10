Statistical Methods
================

Supplementary methods for <b>A globally distributed insecticide
resistance allele confers a fitness cost in the absence of insecticide
in Aedes aegypti, the yellow fever mosquito.</b>

Authors: Cera R. Fisher, Anastacia E. Dressel, Juan J. Silva, Jeffrey G.
Scott

## Table of contents:

1.  [R allele frequency analysis](#Rfreq1)

-   [1.1 Horizontal regression](#1_1)

-   [1.2 First-order regression](#1_2)

-   [1.3 Exploratory plots](#1_3)

-   [1.4 Determine if Cross A is different from Cross B](#1_4)

-   [1.5 Determine if males are different from females](#1_5)

-   [1.6 ANOVA of pooled data and post-hoc Tukey HSD](#1_6)

2.  [Genetic drift analysis](#Drift_2)

-   [2.1 Explanation of simulation experiment](#2_1)

-   [2.2 Overall simulation](#2_2)

-   [2.3 Adjustment of p-values to correct for multiple
    comparisons](#2_3)

3.  [Genotype frequency analysis](#Geno_3)

-   [3.1 ANOVA and Tukey HSD of genotypes by generations](#3_1)

-   [3.2 Hardy-Weinberg equilibrium tests within generation](#3_2)

-   [3.3 Hardy-Weinberg equilibrium based on previous generation](#3_3)

-   [3.4 Fisher’s combined probability for Hardy-Weinberg
    equilibrium](#3_4)

<a id="Rfreq1"></a>

## 1. R allele frequency analysis

Read in data for genotype counts across all samples.

``` r
genos <- read.table(file="../data/AlleleCompetition_GenotypeCounts.txt", header=TRUE, sep="\t")
```

A bit of data manipulation for dummy coding the categorical data.

``` r
genos <- genos %>% mutate(gen = 
  case_when(
    Generation == "F1" ~ 1, 
    Generation == "F2" ~ 2,
    Generation == "F3" ~ 3, 
    Generation == "F5" ~ 5, 
    Generation == "F7" ~ 7
  )
)
```

Get r allele frequency for each row.

``` r
genos <- genos %>% mutate(
  r.freq = ((RR*2) + SR)/(n*2)
)
```

Check structure of data table.

``` r
str(genos)
```

    ## 'data.frame':    26 obs. of  10 variables:
    ##  $ Generation: chr  "F1" "F1" "F2" "F2" ...
    ##  $ Cross     : chr  "A" "B" "A" "A" ...
    ##  $ Replicate : chr  "A" "B" "A1" "A2" ...
    ##  $ Full.Name : chr  "F1A" "F1B" "F2A1" "F2A2" ...
    ##  $ n         : int  88 87 87 85 84 87 82 84 84 85 ...
    ##  $ SS        : int  0 0 19 25 22 26 24 25 25 21 ...
    ##  $ SR        : int  88 87 45 44 51 48 35 38 44 49 ...
    ##  $ RR        : int  0 0 23 16 11 13 23 21 15 15 ...
    ##  $ gen       : num  1 1 2 2 2 2 2 2 3 3 ...
    ##  $ r.freq    : num  0.5 0.5 0.523 0.447 0.435 ...

Now the input data is in the right shape and format, and we’ll be able
to do our allele frequency analyses.

<a id="1_1"></a>

### 1.1 Horizontal regression

Are the data different from a horizontal line? (We also are creating
these linear models as comparisons for the first order regressions.)

``` r
ln.mod.a <- lm(data = genos, formula = r.freq ~ 1)
summary(ln.mod.a)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ 1, data = genos)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.18147 -0.04718  0.01196  0.04606  0.12313 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.39986    0.01553   25.75   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.07919 on 25 degrees of freedom

The intercept is significant, so the data are not a horizontal line.

<a id="1_2"></a>

### 1.2 First-order regression

Fit a linear model with R allele frequency as response variable and
generation as predictor.

``` r
ln.mod2.a <- lm(data=genos, formula = r.freq ~ gen)
summary(ln.mod2.a)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ gen, data = genos)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.088200 -0.027159  0.006866  0.027431  0.094644 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.524225   0.020179  25.978  < 2e-16 ***
    ## gen         -0.031091   0.004495  -6.917 3.75e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04671 on 24 degrees of freedom
    ## Multiple R-squared:  0.6659, Adjusted R-squared:  0.652 
    ## F-statistic: 47.84 on 1 and 24 DF,  p-value: 3.746e-07

Generation is a significant predictor variable for both cross A and
cross B data.

<a id="1_3"></a>

### 1.3 Exploratory plots

Let’s see what our data and the linear model look like.

``` r
plot_a <- ggplot(genos, aes(x=gen, y=r.freq)) + 
  geom_point(aes(color=Replicate)) + 
  geom_smooth(method=lm)

plot_a
```

![](AlleleCompetitionDataAnalysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

We observe some data points that fall outside the confidence interval,
and in some cases, some populations seem to have an increase in R allele
frequency between measured generations. It would be useful to make a
plot where we can see each populations’ allele frequency over
generations.

#### 1.3.a Raw line plot of R frequency change over generations per replicate.

This plot uses some very fancy `geoms()` and labeling to make things
clear, because we like to be clear.

``` r
## Set up  labels
genos <- genos %>% mutate(
  label = if_else(gen == 7, as.character(Replicate), NA_character_)
)

## Set up some colors 

library(colorspace)
replicates = c("A1", "A2", "A3", "B1", "B2", "B3")
colsA = sequential_hcl(h1=15, c1=150, n=6)
colsB = sequential_hcl(h1=200, c1=150, n=6)
cols = c(colsA[1:3], colsB[1:3])
replicate_colors=cols
names(replicate_colors) = replicates

## Make a plot

grand_plot <- ggplot(genos[,], aes(x=gen, y=r.freq)) +
  geom_point(aes(color=Replicate)) + 
  geom_line(aes(group=Replicate, color=Replicate)) + 
  geom_label_repel(aes(label=label, color=Replicate), na.rm=TRUE) + 
  geom_line(data=genos[c(1,3),], aes(x=gen, y=r.freq), color=replicate_colors[1]) + 
  geom_line(data=genos[c(1,4),], aes(x=gen, y=r.freq), color=replicate_colors[2]) + 
  geom_line(data=genos[c(1,5),], aes(x=gen, y=r.freq), color=replicate_colors[3]) +
  geom_line(data=genos[c(2,6),], aes(x=gen, y=r.freq), color=replicate_colors[4]) +
  geom_line(data=genos[c(2,7),], aes(x=gen, y=r.freq), color=replicate_colors[5]) +
  geom_line(data=genos[c(2,8),], aes(x=gen, y=r.freq), color=replicate_colors[6]) +
  scale_color_manual(values = replicate_colors, guide="none") +
  theme_minimal() + 
    labs(title="Change in R allele frequency across generations by replicate", 
       x="Generation", 
       y="R allele frequency")


grand_plot
```

![](AlleleCompetitionDataAnalysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

That reveals the true shape of the data. Already, there’s an impression
that, while R allele frequency *does* decrease over all, the effect size
(related to selection strength) might be small.

We want to know if we can pool this data without violating the
assumption that the populations are not different from each other for
any important biological reasons, such as the sex of the initial
resistant parent population, or the sex of the individual mosquitos.

<a id="1_4"></a>

### 1.4 Determine if Cross A is different from Cross B

Does it matter if the population was founded from resistant males or
resistant females? (i.e., is there a mitochondrial DNA effect?)

Cross A: Female resistant and male susceptible Cross B: Female
susceptible and male resistant

We’re using both ANOVA and paired t.test to show that Cross A and Cross
B are not different.

``` r
## Pooled data - rfreq response variable, generation as predictor. 

lm.mod <- lm(formula=r.freq ~ gen, data=genos)
summary(lm.mod)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ gen, data = genos)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.088200 -0.027159  0.006866  0.027431  0.094644 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.524225   0.020179  25.978  < 2e-16 ***
    ## gen         -0.031091   0.004495  -6.917 3.75e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04671 on 24 degrees of freedom
    ## Multiple R-squared:  0.6659, Adjusted R-squared:  0.652 
    ## F-statistic: 47.84 on 1 and 24 DF,  p-value: 3.746e-07

``` r
## Pooled data - rfreq response variable, generation and cross as predictors. 

lm.mod.cross1 <- lm(formula=r.freq ~ Cross, data=genos)
summary(lm.mod.cross1)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ Cross, data = genos)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.19052 -0.05171  0.01196  0.04731  0.11408 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.40891    0.02226  18.367 1.23e-15 ***
    ## CrossB      -0.01809    0.03148  -0.575    0.571    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08027 on 24 degrees of freedom
    ## Multiple R-squared:  0.01357,    Adjusted R-squared:  -0.02753 
    ## F-statistic: 0.3301 on 1 and 24 DF,  p-value: 0.5709

``` r
## Pooled data - rfreq response variable, generation and a Cross*generation interaction term as predictor 
### This is more relevant, because the interaction term tracks Cross across generations. 

lm.mod.cross <- lm(formula=r.freq ~ gen + Cross:gen, data = genos)
summary(lm.mod.cross)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ gen + Cross:gen, data = genos)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.101492 -0.031598  0.006634  0.026393  0.084479 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.524225   0.020238  25.903  < 2e-16 ***
    ## gen         -0.029192   0.004951  -5.896  5.2e-06 ***
    ## gen:CrossB  -0.003798   0.004093  -0.928    0.363    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04685 on 23 degrees of freedom
    ## Multiple R-squared:  0.678,  Adjusted R-squared:   0.65 
    ## F-statistic: 24.21 on 2 and 23 DF,  p-value: 2.191e-06

These three linear models indicate that while generation remains a
significant predictor, Cross is not, nor is there an interaction between
Cross and generation.

#### 1.4.a ANOVA of linear models

``` r
#Fit analysis of variance model with a call to lm
aov_cross <- aov(data=genos, r.freq ~ Cross)
aov_cross
```

    ## Call:
    ##    aov(formula = r.freq ~ Cross, data = genos)
    ## 
    ## Terms:
    ##                      Cross  Residuals
    ## Sum of Squares  0.00212702 0.15463830
    ## Deg. of Freedom          1         24
    ## 
    ## Residual standard error: 0.08026994
    ## Estimated effects may be unbalanced

``` r
summary(aov_cross)
```

    ##             Df  Sum Sq  Mean Sq F value Pr(>F)
    ## Cross        1 0.00213 0.002127    0.33  0.571
    ## Residuals   24 0.15464 0.006443

``` r
#Compute anova table for variance model and print in 'pretty' way
print(anova(aov_cross))
```

    ## Analysis of Variance Table
    ## 
    ## Response: r.freq
    ##           Df   Sum Sq   Mean Sq F value Pr(>F)
    ## Cross      1 0.002127 0.0021270  0.3301 0.5709
    ## Residuals 24 0.154638 0.0064433

ANOVA shows that the means of the two crosses are not different from
each other.

#### 1.4.b T-tests of Cross A versus Cross B

Let’s also make this comparison with paired t-tests across generations.

``` r
library("plyr")
detach("package:plyr", unload=TRUE)
### This won't work if the plyr library is loaded over dplyr, so we load it and then detach it to make sure it's not masking dplyr.
means <- genos %>% 
  filter(!(Generation=="F1")) %>% 
  group_by(Cross, Generation) 

 
means <- means %>% 
  summarize(avg=mean(r.freq))

means
```

    ## # A tibble: 8 × 3
    ## # Groups:   Cross [2]
    ##   Cross Generation   avg
    ##   <chr> <chr>      <dbl>
    ## 1 A     F2         0.468
    ## 2 A     F3         0.449
    ## 3 A     F5         0.364
    ## 4 A     F7         0.324
    ## 5 B     F2         0.465
    ## 6 B     F3         0.402
    ## 7 B     F5         0.358
    ## 8 B     F7         0.302

``` r
means_a <- means$avg[means$Cross=="A"]
names(means_a) <- means$Generation[means$Cross=="A"]


means_b <- means$avg[means$Cross=="B"]
names(means_b) <- means$Generation[means$Cross=="B"]

means_a
```

    ##        F2        F3        F5        F7 
    ## 0.4681904 0.4492369 0.3642322 0.3236053

``` r
means_b
```

    ##        F2        F3        F5        F7 
    ## 0.4651268 0.4017654 0.3583116 0.3016728

``` r
t.test(means_a, means_b, paired=T, var.equal = T)
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  means_a and means_b
    ## t = 1.9256, df = 3, p-value = 0.1498
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.01279049  0.05198465
    ## sample estimates:
    ## mean difference 
    ##      0.01959708

``` r
t.test(means_a, means_b, paired=F, var.equal = T)
```

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  means_a and means_b
    ## t = 0.40215, df = 6, p-value = 0.7015
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.0996439  0.1388381
    ## sample estimates:
    ## mean of x mean of y 
    ## 0.4013162 0.3817191

Paired and unpaired t-tests indicate that the generation means within
each cross are not different from each other.

This is also clear because of the fact that their confidence intervals
overlap:

#### 1.4.c Smoothed linear model plot for two crosses combined

``` r
cross_cols = replicate_colors[c(2,5)]
names(cross_cols) = c("A", "B")

grand_plot <- ggplot(genos[,], aes(x=gen, y=r.freq, groups=Cross)) +
  geom_point(aes(color=Cross)) + 
  geom_smooth(data=subset(genos, Cross=="A"), method=lm, alpha=0.2,  aes(color=Cross, fill=Cross)) +
  geom_smooth(data=subset(genos, Cross=="B"), method=lm, alpha=0.2,  aes(color=Cross, fill=Cross)) +
  theme_minimal() + 
  scale_color_manual(values=cross_cols) + 
  labs(title="Decrease in R allele frequency across 7 generations", 
       x = "Generation", 
       y = "R allele frequency", 
       caption = "Shaded area represents 95% confidence interval based on linear regression")

grand_plot
```

![](AlleleCompetitionDataAnalysis_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

<a id="1_5"></a>

### 1.5 Determine if males are different from females

We also want to be sure that there is no sex-biased effect on the R
frequency or genotypes. The data segregated by sex was collated
separately in Excel, so we need to read that file in.

``` r
sex.genos <- read.table("../data/AlleleCompetition_GenotypeCounts_bySex.txt", header=T, sep="\t")

sex.genos <- sex.genos %>% mutate(gen = 
  case_when(
    Generation == "F1" ~ 1, 
    Generation == "F2" ~ 2,
    Generation == "F3" ~ 3, 
    Generation == "F5" ~ 5, 
    Generation == "F7" ~ 7
  )
)

sex.genos <- sex.genos %>% mutate(
  r.freq = ((RR*2) + SR)/(n*2)
)

str(sex.genos)
```

    ## 'data.frame':    52 obs. of  11 variables:
    ##  $ Generation: chr  "F1" "F1" "F1" "F1" ...
    ##  $ Cross     : chr  "A" "A" "B" "B" ...
    ##  $ Replicate : chr  "A" "A" "B" "B" ...
    ##  $ Full.Name : chr  "F1A" "F1A" "F1B" "F1B" ...
    ##  $ Sex       : chr  "M" "F" "M" "F" ...
    ##  $ n         : int  44 44 44 44 43 44 42 43 42 42 ...
    ##  $ SS        : int  0 0 0 0 8 11 14 11 12 10 ...
    ##  $ SR        : int  44 44 44 44 20 25 17 27 24 27 ...
    ##  $ RR        : int  0 0 0 0 15 8 11 5 6 5 ...
    ##  $ gen       : num  1 1 1 1 2 2 2 2 2 2 ...
    ##  $ r.freq    : num  0.5 0.5 0.5 0.5 0.581 ...

Now we have the data read in and we can run the same kinds of tests.

#### 1.5.a Fit linear models on sex segregated data

``` r
## Linear model with rfreq as response variable and generation as predictor. 
lm.fit.gen <- lm(formula=r.freq ~ gen, data=sex.genos)
summary(lm.fit.gen)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ gen, data = sex.genos)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.112709 -0.033688  0.002838  0.032429  0.119136 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.524807   0.016816  31.209  < 2e-16 ***
    ## gen         -0.031274   0.003746  -8.349 4.82e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05505 on 50 degrees of freedom
    ## Multiple R-squared:  0.5823, Adjusted R-squared:  0.574 
    ## F-statistic: 69.71 on 1 and 50 DF,  p-value: 4.82e-11

``` r
## Linear model with rfreq as response variable and sex as predictor. 
lm.fit.sex <- lm(formula=r.freq ~ Sex, data=sex.genos)
summary(lm.fit.sex)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ Sex, data = sex.genos)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.203776 -0.042779  0.008322  0.061166  0.178929 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 0.396957   0.016696  23.775   <2e-16 ***
    ## SexM        0.005509   0.023612   0.233    0.816    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.08513 on 50 degrees of freedom
    ## Multiple R-squared:  0.001087,   Adjusted R-squared:  -0.01889 
    ## F-statistic: 0.05443 on 1 and 50 DF,  p-value: 0.8165

``` r
## Linear model with rfreq as response variable, gen as predictor, and an interaction term of sex*gen as predictor. 
lm.fit.sexgen <- lm(formula=r.freq ~ gen + Sex:gen, data=sex.genos)
summary(lm.fit.sexgen)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ gen + Sex:gen, data = sex.genos)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.111971 -0.033582  0.002838  0.032271  0.118925 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  0.5248069  0.0169859  30.897  < 2e-16 ***
    ## gen         -0.0313792  0.0041553  -7.552 9.32e-10 ***
    ## gen:SexM     0.0002109  0.0034355   0.061    0.951    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05561 on 49 degrees of freedom
    ## Multiple R-squared:  0.5824, Adjusted R-squared:  0.5653 
    ## F-statistic: 34.16 on 2 and 49 DF,  p-value: 5.128e-10

#### 1.5.b ANOVA of sex-segregated data

There doesn’t appear to be any effect of sex. Let’s model with an ANOVA.

``` r
# Fit variance model with generation, sex, and a sex*gen term. 
aov_sex <- aov(data=sex.genos, formula=r.freq ~ Sex)
aov_sex
```

    ## Call:
    ##    aov(formula = r.freq ~ Sex, data = sex.genos)
    ## 
    ## Terms:
    ##                       Sex Residuals
    ## Sum of Squares  0.0003945 0.3623946
    ## Deg. of Freedom         1        50
    ## 
    ## Residual standard error: 0.08513455
    ## Estimated effects may be unbalanced

``` r
summary(aov_sex)
```

    ##             Df Sum Sq  Mean Sq F value Pr(>F)
    ## Sex          1 0.0004 0.000395   0.054  0.816
    ## Residuals   50 0.3624 0.007248

``` r
# Compute anova and print table in a 'pretty' way
print(anova_sex <- anova(aov_sex))
```

    ## Analysis of Variance Table
    ## 
    ## Response: r.freq
    ##           Df  Sum Sq   Mean Sq F value Pr(>F)
    ## Sex        1 0.00039 0.0003945  0.0544 0.8165
    ## Residuals 50 0.36239 0.0072479

ANOVA indicates that there is no difference between the group means when
the data is partitioned by sex.

#### 1.5.c T-tests of sex segregated data

Now let’s try the t-tests. First, filter and separate the data so we can
compare them.

``` r
# This won't work is plyr is loaded over dplyr
library("plyr")
detach("package:plyr", unload=TRUE)
means <- sex.genos %>%
  filter(!(Generation=="F1")) %>% 
  group_by(Sex, Generation) %>% 
  summarize(avg=mean(r.freq))


means_M <- means$avg[means$Sex=="M"]
names(means_M) <- means$Generation[means$Sex=="M"]
means_M
```

    ##        F2        F3        F5        F7 
    ## 0.4863530 0.4259178 0.3490433 0.3160393

``` r
means_F <- means$avg[means$Sex=="F"]
names(means_F) <- means$Generation[means$Sex=="F"]
means_F
```

    ##        F2        F3        F5        F7 
    ## 0.4469760 0.4258737 0.3732323 0.3074000

``` r
sex.diff <- t.test(means_F, means_M, var.equal = T, paired = T)
sex.diff
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  means_F and means_M
    ## t = -0.45462, df = 3, p-value = 0.6803
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.04774400  0.03580832
    ## sample estimates:
    ## mean difference 
    ##    -0.005967839

``` r
sex.diff <- t.test(means_F, means_M, var.equal = T, paired = F)
sex.diff
```

    ## 
    ##  Two Sample t-test
    ## 
    ## data:  means_F and means_M
    ## t = -0.12083, df = 6, p-value = 0.9078
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.1268215  0.1148859
    ## sample estimates:
    ## mean of x mean of y 
    ## 0.3883705 0.3943383

According to ANOVA and to t-tests, paired and unpaired, the R
frequencies are not different between male and female. We can safely
combine all this data.

<a id="1_6"></a>

### 1.6 ANOVA of pooled data and post-hoc Tukey HSD

We want to take the ANOVA over the whole set, which we know indicates
that generation is a significant term, and perform a Tukey Honestly
Significant Difference (post-hoc) test to determine which groups are
different from each other.

``` r
gens_lm <- lm(data=genos, r.freq ~ Generation)
summary(gens_lm)
```

    ## 
    ## Call:
    ## lm(formula = r.freq ~ Generation, data = genos)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.094248 -0.025854  0.002697  0.025285  0.088596 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.50000    0.03499  14.291 2.73e-12 ***
    ## GenerationF2 -0.03334    0.04040  -0.825 0.418487    
    ## GenerationF3 -0.07450    0.04040  -1.844 0.079335 .  
    ## GenerationF5 -0.13873    0.04040  -3.434 0.002492 ** 
    ## GenerationF7 -0.18736    0.04040  -4.638 0.000142 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04948 on 21 degrees of freedom
    ## Multiple R-squared:  0.672,  Adjusted R-squared:  0.6096 
    ## F-statistic: 10.76 on 4 and 21 DF,  p-value: 6.642e-05

``` r
gens_anova <- aov(data=genos,r.freq ~ Generation)
summary(gens_anova)
```

    ##             Df  Sum Sq  Mean Sq F value   Pr(>F)    
    ## Generation   4 0.10535 0.026338   10.76 6.64e-05 ***
    ## Residuals   21 0.05141 0.002448                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
print(anova(gens_anova))
```

    ## Analysis of Variance Table
    ## 
    ## Response: r.freq
    ##            Df   Sum Sq   Mean Sq F value    Pr(>F)    
    ## Generation  4 0.105352 0.0263381  10.758 6.642e-05 ***
    ## Residuals  21 0.051413 0.0024482                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## Please note here that we are using Generation as data type factor, and not generation as data type integer. This is because we're not trying to understand whether change over a continuous variable (the arrow of time!) is happening, we're trying to understand how groups are different from each other. The tukey test will need grouped data and that's what factors are.
## The way that the grouped data are different from each other will show us the "arrow of time"

tukey_rfreq <- TukeyHSD(x=gens_anova, 'Generation', conf.level=0.95)
tukey_rfreq
```

    ##   Tukey multiple comparisons of means
    ##     95% family-wise confidence level
    ## 
    ## Fit: aov(formula = r.freq ~ Generation, data = genos)
    ## 
    ## $Generation
    ##              diff        lwr         upr     p adj
    ## F2-F1 -0.03334143 -0.1536944  0.08701158 0.9197964
    ## F3-F1 -0.07449885 -0.1948519  0.04585416 0.3762995
    ## F5-F1 -0.13872811 -0.2590811 -0.01837510 0.0188432
    ## F7-F1 -0.18736095 -0.3077140 -0.06700794 0.0012018
    ## F3-F2 -0.04115742 -0.1262599  0.04394501 0.6096456
    ## F5-F2 -0.10538668 -0.1904891 -0.02028425 0.0106591
    ## F7-F2 -0.15401952 -0.2391220 -0.06891709 0.0002106
    ## F5-F3 -0.06422926 -0.1493317  0.02087317 0.2011663
    ## F7-F3 -0.11286210 -0.1979645 -0.02775967 0.0058795
    ## F7-F5 -0.04863284 -0.1337353  0.03646959 0.4537630

``` r
plot(tukey_rfreq)
```

![The family-wise confident interval plot is hard to read in native plot
size](Z:/Shared%20Documents/Cera%20Fisher/Research/Allele_Competition_1534C-RK/AC_1534C-RK_DataAnalysis/output/Tukey_RFreq_Generations_familywiseConfLevel.png)

``` r
# Generation-wise means 
library("plyr")
detach("package:plyr", unload=TRUE)
means_gen <- genos %>%
  group_by(Generation) %>% 
  summarize(avg=mean(r.freq))
means_gen 
```

    ## # A tibble: 5 × 2
    ##   Generation   avg
    ##   <chr>      <dbl>
    ## 1 F1         0.5  
    ## 2 F2         0.467
    ## 3 F3         0.426
    ## 4 F5         0.361
    ## 5 F7         0.313

``` r
### we're actually going to use the package agricolae for the tukey HSD, because
### it more transparently computes groupings
library(agricolae)

tukey_rfreq_2 <- HSD.test(gens_anova, trt="Generation")
tukey_rfreq_2
```

    ## $statistics
    ##       MSerror Df      Mean       CV
    ##   0.002448236 21 0.3998625 12.37417
    ## 
    ## $parameters
    ##    test     name.t ntr StudentizedRange alpha
    ##   Tukey Generation   5         4.212995  0.05
    ## 
    ## $means
    ##       r.freq        std r       Min       Max       Q25       Q50       Q75
    ## F1 0.5000000 0.00000000 2 0.5000000 0.5000000 0.5000000 0.5000000 0.5000000
    ## F2 0.4666586 0.03777438 6 0.4252874 0.5229885 0.4376576 0.4616246 0.4894744
    ## F3 0.4255011 0.03824376 6 0.3571429 0.4647059 0.4153516 0.4408263 0.4421907
    ## F5 0.3612719 0.03386668 6 0.3275862 0.4166667 0.3356742 0.3546816 0.3771780
    ## F7 0.3126390 0.07903259 6 0.2183908 0.4012346 0.2413607 0.3296650 0.3711310
    ## 
    ## $comparison
    ## NULL
    ## 
    ## $groups
    ##       r.freq groups
    ## F1 0.5000000      a
    ## F2 0.4666586      a
    ## F3 0.4255011     ab
    ## F5 0.3612719     bc
    ## F7 0.3126390      c
    ## 
    ## attr(,"class")
    ## [1] "group"

``` r
plot(tukey_rfreq_2)
```

![](AlleleCompetitionDataAnalysis_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

Groups a and c don’t overlap, though a and b do, and b and c do. Let’s
add that to our linear regression plot.

``` r
grp.df <- tukey_rfreq_2$groups
groups <- grp.df$groups
names(groups) <- row.names(grp.df)

cross_cols = replicate_colors[c(2,5)]
names(cross_cols) = c("A", "B")

## it's a shame to not be able to do this more programmatically, but we'll just have to hard code it for the publication plot. 
genos <- genos %>% mutate(
  tukey_labels = case_when(
    Generation == "F1" ~ "a", 
    Generation == "F2" ~ "a",
    Generation == "F3" ~ "ab", 
    Generation == "F5" ~ "bc", 
    Generation == "F7" ~ "c"
  )
) %>%
  mutate(
    tukey_labels = if_else(Replicate == "A1", as.character(tukey_labels), NA_character_)
  )

genos$tukey_labels[1] <- "a"

grand_plot <- ggplot(genos[,], aes(x=gen, y=r.freq)) +
  geom_point(aes(color=Cross)) + 
  geom_smooth(method=lm, alpha=0.2, fill="grey", color="grey") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_x_continuous(breaks=c(1,2,3,5,7), labels=c("F1","F2", "F3", "F5", "F7")) + 
  scale_color_manual(values=cross_cols) + 
  annotate(geom="text", x=genos$gen[1], y=genos$r.freq[1]+0.02, 
           label=genos$tukey_labels[1]) +
  annotate(geom="text", x=genos$gen[3], y=max(genos$r.freq[genos$Generation=="F2"])+0.02, 
           label=genos$tukey_labels[3]) +
  annotate(geom="text", x=genos$gen[9], y=max(genos$r.freq[genos$Generation=="F3"])+0.02, 
           label=genos$tukey_labels[9]) +
  annotate(geom="text", x=genos$gen[15], y=max(genos$r.freq[genos$Generation=="F5"])+0.02, 
           label=genos$tukey_labels[15]) +
  annotate(geom="text", x=genos$gen[21], y=max(genos$r.freq[genos$Generation=="F7"])+0.02, 
           label=genos$tukey_labels[21]) +
  labs(title="Decrease in R allele frequency across 7 generations", 
       x = "Generation", 
       y = "R allele frequency") 
  

grand_plot
```

![](AlleleCompetitionDataAnalysis_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

<a id="Drift_2"></a>

## 2. Genetic drift analysis by simulation

Are the changes in allele frequency just due to genetic drift? The best
way to determine this is by simulation.

<a id="2_1"></a>

### 2.1 Explanation of simulation experiment

The code for the function is provided separately in file
`MultiGenDrift.R` because it is pretty long.

``` r
source("MultiGenDrift.R")
```

<a id="2_2"></a>

### 2.2 Overall simulation

``` r
popsize = 800
Rstart = 0.5
Rend = mean(genos$r.freq[genos$Generation=="F7"])
nsims=100000
histo_vector <- MultiGenDrift(SIMS = nsims, gens = 6, size_of_pop = popsize, R_start_freq = Rstart, R_end_freq = Rend)

empirical_p <- histo_vector[[2]]
ResultsDF <- data.frame(histo_vector[[1]])
colnames(ResultsDF) <- 'Result'

hist_plot <- ggplot(data=ResultsDF, aes(x=Result)) + 
  geom_histogram(stat="bin", position="dodge", 
                 binwidth=0.005, fill="grey", alpha=0.5, color="black") +
  geom_vline(xintercept = Rstart, color="red") +
  geom_vline(xintercept = Rend, color="red") + 
  theme_minimal() +
  labs(title="Drift simulation results", 
       x=paste("Range of R allele frequencies at Generation", 7, sep=" "), 
       y=paste("Count of frequencies in 0.005 width bins", sep=""), 
       caption="Over 100,000 simulations of drift, the number of times an R allele
       frequency occurred less than or equal to our real result (0.315) was 0.")
hist_plot

# write.table(ResultsDF, file="../output/DriftSimulation_F1toF7_100k_Results.tab", 
#             row.names = F, quote=F, eol="\n", sep="\t")
```

That was 100,000 simulations, so we’re not going to re-run it in the
R-markdown generation, but here is the resulting plot at the time of
analysis:

![histogram plot of 100,000 genetic drift
simulations](../output/DriftSimulations_F1toF7.png)

That’s how we do this. But we actually need to compute this separately
for each population at each generation, so this really needs to be
crunched through via script.

This script is provided separately in the file `DriftSimulations.Rmd`.

<a id="2_3"></a>

### 2.3 Adjustment of p-values to correct for multiple comparisons

Alpha thresholds control for Type I error, but when we have a set of
repeated tests, we run the risk of hitting that alpha erroneously.
P-value correction is necessary when multiply testing. Bonferroni is the
classical correction, but the stat vignette says “There seems to be no
reason to use the unmodified Bonferroni correction because it is
dominated by Holm’s method, which is also valid under arbitrary
assumptions”. Therefore we are using `method="holm"`.

I collected the initial p-values generated at run time and put them into
a tab-delimited text file to read into R.

``` r
pvals <- read.table("../data/DriftSimulationPvalues.txt", sep="\t", header=T)
str(pvals)
```

    ## 'data.frame':    42 obs. of  2 variables:
    ##  $ Comparison: chr  "F1:F2A1" "F1:F2A2" "F1:F2A3" "F1:F2B1" ...
    ##  $ Pvalue    : num  0.03494 0.00002 0 0 0.31594 ...

``` r
pvals$padj <- p.adjust(p=pvals$Pvalue, method="holm")

write.table(pvals, "../output/DriftSimulationAdjustedPvalues.txt", row.names = F, 
            quote=F, eol="\n", sep="\t")

nonsig <- pvals %>% filter(padj > 0.05)

nonsig
```

    ##    Comparison  Pvalue    padj
    ## 1     F1:F2A1 0.03494 0.27952
    ## 2     F1:F2B2 0.31594 0.77196
    ## 3     F1:F2B3 0.02796 0.25164
    ## 4     F1:F3A2 0.02300 0.23000
    ## 5   F2A2:F3A2 0.07964 0.55748
    ## 6   F2A3:F3A3 0.25732 0.77196
    ## 7   F3B1:F5B1 0.08274 0.55748
    ## 8   F5A1:F7A1 0.14512 0.72560
    ## 9   F5A2:F7A2 0.18576 0.74304
    ## 10  F5B1:F7B1 0.43002 0.77196

After pval correction, there are 10 (of 42) comparisons that no longer
fall below the alpha threshold of 0.05.

<a id="Geno_3"></a>

## 3. Genotype frequency analysis

Does genotype frequency change over generations? If they do, are they
out of HW equilibrium? Being out of HW equilibrium indicates that the
assumptions of Hardy-Weinberg’s law are being violated, and one of the
four forces of evolution is active. Because we’re controlling for
everything but selection (and we tested for drift), if the populations
are out of equilibrium, it must be due to selection.

First, convert data to long format.

``` r
str(genos)
```

    ## 'data.frame':    26 obs. of  12 variables:
    ##  $ Generation  : chr  "F1" "F1" "F2" "F2" ...
    ##  $ Cross       : chr  "A" "B" "A" "A" ...
    ##  $ Replicate   : chr  "A" "B" "A1" "A2" ...
    ##  $ Full.Name   : chr  "F1A" "F1B" "F2A1" "F2A2" ...
    ##  $ n           : int  88 87 87 85 84 87 82 84 84 85 ...
    ##  $ SS          : int  0 0 19 25 22 26 24 25 25 21 ...
    ##  $ SR          : int  88 87 45 44 51 48 35 38 44 49 ...
    ##  $ RR          : int  0 0 23 16 11 13 23 21 15 15 ...
    ##  $ gen         : num  1 1 2 2 2 2 2 2 3 3 ...
    ##  $ r.freq      : num  0.5 0.5 0.523 0.447 0.435 ...
    ##  $ label       : chr  NA NA NA NA ...
    ##  $ tukey_labels: chr  "a" NA "a" NA ...

``` r
genos_long <- pivot_longer(genos, 
                           cols = 6:8, # operating on the genotype columns
                           names_to = "genotype", 
                           values_to = "genotype.count")
head(genos_long)
```

    ## # A tibble: 6 × 11
    ##   Gener…¹ Cross Repli…² Full.…³     n   gen r.freq label tukey…⁴ genot…⁵ genot…⁶
    ##   <chr>   <chr> <chr>   <chr>   <int> <dbl>  <dbl> <chr> <chr>   <chr>     <int>
    ## 1 F1      A     A       F1A        88     1    0.5 <NA>  a       SS            0
    ## 2 F1      A     A       F1A        88     1    0.5 <NA>  a       SR           88
    ## 3 F1      A     A       F1A        88     1    0.5 <NA>  a       RR            0
    ## 4 F1      B     B       F1B        87     1    0.5 <NA>  <NA>    SS            0
    ## 5 F1      B     B       F1B        87     1    0.5 <NA>  <NA>    SR           87
    ## 6 F1      B     B       F1B        87     1    0.5 <NA>  <NA>    RR            0
    ## # … with abbreviated variable names ¹​Generation, ²​Replicate, ³​Full.Name,
    ## #   ⁴​tukey_labels, ⁵​genotype, ⁶​genotype.count

``` r
genos_long <- subset(genos_long, Generation!="F1") # Exclude F1, because it is out of whack to start with!

genos_long$Generation <- factor(genos_long$Generation)
genos_long$genotype <- factor(genos_long$genotype)
genos_long <- mutate(genos_long, genot.freq.obs = genotype.count/n)
```

<a id="3_1"></a>

### 3.1 ANOVA of genotypes by generations

``` r
# Get means of genotype frequencies over generations 

library("plyr")
detach("package:plyr", unload=TRUE)
genos_means <- genos_long %>%
  filter(!(Generation=="F1")) %>% 
  group_by(Generation, genotype) %>% 
  summarize(avg=mean(genot.freq.obs))
genos_means
```

    ## # A tibble: 12 × 3
    ## # Groups:   Generation [4]
    ##    Generation genotype    avg
    ##    <fct>      <fct>     <dbl>
    ##  1 F2         RR       0.211 
    ##  2 F2         SR       0.512 
    ##  3 F2         SS       0.277 
    ##  4 F3         RR       0.170 
    ##  5 F3         SR       0.511 
    ##  6 F3         SS       0.319 
    ##  7 F5         RR       0.130 
    ##  8 F5         SR       0.463 
    ##  9 F5         SS       0.407 
    ## 10 F7         RR       0.0935
    ## 11 F7         SR       0.438 
    ## 12 F7         SS       0.468

``` r
lm.genos.2 <- lm(data=genos_long, formula=genot.freq.obs ~ genotype + Generation + Generation:genotype)
summary(lm.genos.2)
```

    ## 
    ## Call:
    ## lm(formula = genot.freq.obs ~ genotype + Generation + Generation:genotype, 
    ##     data = genos_long)
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.134853 -0.037446  0.002816  0.031364  0.129515 
    ## 
    ## Coefficients:
    ##                         Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)              0.21058    0.02323   9.065 7.68e-13 ***
    ## genotypeSR               0.30158    0.03285   9.180 4.92e-13 ***
    ## genotypeSS               0.06668    0.03285   2.030 0.046824 *  
    ## GenerationF3            -0.04047    0.03285  -1.232 0.222839    
    ## GenerationF5            -0.08081    0.03285  -2.460 0.016802 *  
    ## GenerationF7            -0.11711    0.03285  -3.565 0.000722 ***
    ## genotypeSR:GenerationF3  0.03909    0.04646   0.841 0.403525    
    ## genotypeSS:GenerationF3  0.08231    0.04646   1.772 0.081523 .  
    ## genotypeSR:GenerationF5  0.03166    0.04646   0.681 0.498282    
    ## genotypeSS:GenerationF5  0.21077    0.04646   4.537 2.80e-05 ***
    ## genotypeSR:GenerationF7  0.04330    0.04646   0.932 0.355061    
    ## genotypeSS:GenerationF7  0.30804    0.04646   6.630 1.07e-08 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0569 on 60 degrees of freedom
    ## Multiple R-squared:  0.8889, Adjusted R-squared:  0.8685 
    ## F-statistic: 43.63 on 11 and 60 DF,  p-value: < 2.2e-16

``` r
anova_genos <- aov(data = subset(genos_long, Generation!="F1"), formula=genot.freq.obs ~ genotype + Generation + Generation:genotype)
summary(anova_genos)
```

    ##                     Df Sum Sq Mean Sq F value   Pr(>F)    
    ## genotype             2 1.3507  0.6753  208.57  < 2e-16 ***
    ## Generation           3 0.0000  0.0000    0.00        1    
    ## genotype:Generation  6 0.2033  0.0339   10.46 6.47e-08 ***
    ## Residuals           60 0.1943  0.0032                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
print(anova(anova_genos))
```

    ## Analysis of Variance Table
    ## 
    ## Response: genot.freq.obs
    ##                     Df  Sum Sq Mean Sq F value    Pr(>F)    
    ## genotype             2 1.35067 0.67534 208.570 < 2.2e-16 ***
    ## Generation           3 0.00000 0.00000   0.000         1    
    ## genotype:Generation  6 0.20326 0.03388  10.462  6.47e-08 ***
    ## Residuals           60 0.19428 0.00324                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

There is a statistically significant difference between genotype
frequency means.

``` r
### Separate anovas for each genotype

lm_RR <- lm(data=subset(genos_long, genotype=="RR"), formula=genot.freq.obs~Generation)
lm_RR
```

    ## 
    ## Call:
    ## lm(formula = genot.freq.obs ~ Generation, data = subset(genos_long, 
    ##     genotype == "RR"))
    ## 
    ## Coefficients:
    ##  (Intercept)  GenerationF3  GenerationF5  GenerationF7  
    ##      0.21058      -0.04047      -0.08081      -0.11711

``` r
summary(lm_RR)
```

    ## 
    ## Call:
    ## lm(formula = genot.freq.obs ~ Generation, data = subset(genos_long, 
    ##     genotype == "RR"))
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.079626 -0.026554  0.000117  0.037972  0.069910 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.21058    0.01836  11.471    3e-10 ***
    ## GenerationF3 -0.04047    0.02596  -1.559 0.134726    
    ## GenerationF5 -0.08081    0.02596  -3.113 0.005483 ** 
    ## GenerationF7 -0.11711    0.02596  -4.511 0.000213 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.04496 on 20 degrees of freedom
    ## Multiple R-squared:  0.5325, Adjusted R-squared:  0.4624 
    ## F-statistic: 7.593 on 3 and 20 DF,  p-value: 0.001399

``` r
anova_RR <- aov(data=subset(genos_long, genotype=="RR"), formula=genot.freq.obs ~ Generation)
summary(anova_RR)
```

    ##             Df  Sum Sq  Mean Sq F value Pr(>F)   
    ## Generation   3 0.04606 0.015352   7.593 0.0014 **
    ## Residuals   20 0.04044 0.002022                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
print(anova(anova_RR))
```

    ## Analysis of Variance Table
    ## 
    ## Response: genot.freq.obs
    ##            Df   Sum Sq   Mean Sq F value   Pr(>F)   
    ## Generation  3 0.046055 0.0153518   7.593 0.001399 **
    ## Residuals  20 0.040437 0.0020218                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generation is significant for RR.

``` r
lm_SR <- lm(data=subset(genos_long, genotype=="SR"), formula = genot.freq.obs ~ Generation)
summary(lm_SR)
```

    ## 
    ## Call:
    ## lm(formula = genot.freq.obs ~ Generation, data = subset(genos_long, 
    ##     genotype == "SR"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.08533 -0.04458  0.00437  0.03113  0.09498 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.51216    0.02254  22.722 9.33e-16 ***
    ## GenerationF3 -0.00138    0.03188  -0.043   0.9659    
    ## GenerationF5 -0.04915    0.03188  -1.542   0.1387    
    ## GenerationF7 -0.07381    0.03188  -2.316   0.0313 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.05521 on 20 degrees of freedom
    ## Multiple R-squared:  0.2825, Adjusted R-squared:  0.1749 
    ## F-statistic: 2.625 on 3 and 20 DF,  p-value: 0.07864

``` r
anova_SR <- aov(data=subset(genos_long, genotype=="SR"), formula=genot.freq.obs ~ Generation)
summary(anova_SR)
```

    ##             Df  Sum Sq  Mean Sq F value Pr(>F)  
    ## Generation   3 0.02400 0.008001   2.625 0.0786 .
    ## Residuals   20 0.06097 0.003048                 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
print(anova(anova_SR))
```

    ## Analysis of Variance Table
    ## 
    ## Response: genot.freq.obs
    ##            Df   Sum Sq   Mean Sq F value  Pr(>F)  
    ## Generation  3 0.024004 0.0080014  2.6248 0.07864 .
    ## Residuals  20 0.060967 0.0030483                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generation is not signification for SR.

``` r
lm_SS <- lm(data=subset(genos_long, genotype=="SS"), formula=genot.freq.obs ~ Generation)
summary(lm_SS)
```

    ## 
    ## Call:
    ## lm(formula = genot.freq.obs ~ Generation, data = subset(genos_long, 
    ##     genotype == "SS"))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.13485 -0.03747 -0.00290  0.02067  0.12951 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   0.27726    0.02782   9.966 3.35e-09 ***
    ## GenerationF3  0.04185    0.03934   1.064  0.30016    
    ## GenerationF5  0.12996    0.03934   3.303  0.00355 ** 
    ## GenerationF7  0.19093    0.03934   4.853 9.65e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.06814 on 20 degrees of freedom
    ## Multiple R-squared:  0.5892, Adjusted R-squared:  0.5276 
    ## F-statistic: 9.561 on 3 and 20 DF,  p-value: 0.0004008

``` r
anova_SS <- aov(data=subset(genos_long, genotype=="SS"), formula=genot.freq.obs ~ Generation)
summary(anova_SS)
```

    ##             Df  Sum Sq Mean Sq F value   Pr(>F)    
    ## Generation   3 0.13320 0.04440   9.561 0.000401 ***
    ## Residuals   20 0.09287 0.00464                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
print(anova(anova_SS))
```

    ## Analysis of Variance Table
    ## 
    ## Response: genot.freq.obs
    ##            Df   Sum Sq  Mean Sq F value    Pr(>F)    
    ## Generation  3 0.133199 0.044400  9.5614 0.0004008 ***
    ## Residuals  20 0.092873 0.004644                      
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Generation is significant for SS.

#### 3.1.a Tukey HSD post-hoc test

``` r
library(agricolae)
tukey_genos_2 <- HSD.test(anova_genos, trt=c("Generation", "genotype"))
tukey_genos_2
```

    ## $statistics
    ##       MSerror Df      Mean       CV       MSD
    ##   0.003237935 60 0.3333333 17.07086 0.1117033
    ## 
    ## $parameters
    ##    test              name.t ntr StudentizedRange alpha
    ##   Tukey Generation:genotype  12         4.808477  0.05
    ## 
    ## $means
    ##       genot.freq.obs        std r        Min       Max        Q25       Q50
    ## F2:RR     0.21057810 0.06310928 6 0.13095238 0.2804878 0.15912779 0.2191176
    ## F2:SR     0.51216094 0.06557112 6 0.42682927 0.6071429 0.46859606 0.5174442
    ## F2:SS     0.27726096 0.03195014 6 0.21839080 0.2988506 0.26959930 0.2934003
    ## F3:RR     0.17011073 0.02278694 6 0.14285714 0.2068966 0.15454856 0.1705882
    ## F3:SR     0.51078083 0.05400563 6 0.42857143 0.5764706 0.48135525 0.5177187
    ## F3:SS     0.31910844 0.06216293 6 0.24705882 0.4285714 0.28616947 0.3097291
    ## F5:RR     0.12976851 0.03141482 6 0.08988764 0.1777778 0.10919540 0.1298851
    ## F5:SR     0.46300675 0.03391650 6 0.40229885 0.5056180 0.46130486 0.4662879
    ## F5:SS     0.40722474 0.04338406 6 0.34444444 0.4712644 0.38977273 0.4022472
    ## F7:RR     0.09346421 0.05097483 6 0.03370787 0.1444444 0.04621091 0.1061739
    ## F7:SR     0.43834967 0.06186155 6 0.36781609 0.5308642 0.38820894 0.4469823
    ## F7:SS     0.46818612 0.10866137 6 0.33333333 0.5977011 0.39880952 0.4468439
    ##             Q75
    ## F2:RR 0.2607759
    ## F2:SR 0.5432049
    ## F2:SS 0.2967437
    ## F3:RR 0.1780462
    ## F3:SR 0.5456583
    ## F3:SS 0.3333667
    ## F5:RR 0.1441288
    ## F5:SR 0.4750000
    ## F5:SS 0.4287098
    ## F7:RR 0.1345899
    ## F7:SR 0.4633721
    ## F7:SS 0.5629083
    ## 
    ## $comparison
    ## NULL
    ## 
    ## $groups
    ##       genot.freq.obs groups
    ## F2:SR     0.51216094      a
    ## F3:SR     0.51078083      a
    ## F7:SS     0.46818612      a
    ## F5:SR     0.46300675      a
    ## F7:SR     0.43834967      a
    ## F5:SS     0.40722474     ab
    ## F3:SS     0.31910844     bc
    ## F2:SS     0.27726096     cd
    ## F2:RR     0.21057810    cde
    ## F3:RR     0.17011073    def
    ## F5:RR     0.12976851     ef
    ## F7:RR     0.09346421      f
    ## 
    ## attr(,"class")
    ## [1] "group"

``` r
plot(tukey_genos_2)
```

This plot is hard to read as natively generated, so I’ll embed the
results with an expanded x axis.

![Genotype by generation Tukey groups and
ranges](../output/Tukey_Groups_and_Range_GenotypeFreqsInteractionTermWithGenerationAsFactor.png)
This pattern echoes what we saw with the adjusted p-values for the Drift
simulations. Differences between genotypes over generations are
difficult to detect when the generations are only one or two apart, but
are very obvious over 7 generations.

#### 3.1.b Build plot with significance groups

Get the group assignments from the Tukey HSD test so we can add it to
our genotype frequency plot.

``` r
data_summary <- function(data, varname, groupnames){
  require(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      N = length2(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  data_sum$se <- data_sum$sd / sqrt(data_sum$N)
 return(data_sum)
}

genos_long <- subset(genos_long, Generation!="F1")
errorbars <- data_summary(genos_long, varname="genot.freq.obs", groupnames=c("Generation", "genotype"))

tukey_genos_groups <- tukey_genos_2$groups
tukey_genos_groups <- tukey_genos_groups %>% 
  mutate(group_name = row.names(tukey_genos_groups)) %>% 
  mutate(Generation = gsub(x=group_name, ":.*", "")) %>% 
  mutate(genotype = gsub(x=group_name, ".*:", ""))

df <- errorbars %>% 
  mutate(gen=case_when(
    Generation=="F2" ~ 2, 
    Generation=="F3" ~ 3, 
    Generation=="F5" ~ 5, 
    Generation=="F7" ~ 7
  )) %>%
  mutate(group_name=paste(gen, genotype, sep=":")) %>%
  mutate(tukey_group=case_when(
    group_name=="3:SR" ~ "a", 
    group_name=="2:SR" ~ "a", 
    group_name=="5:SR" ~ "a",
    group_name=="7:SS" ~ "a", 
    group_name=="7:SR" ~ "a",
    group_name=="5:SS" ~ "ab",
    group_name=="3:SS" ~ "bc", 
    group_name=="2:SS" ~ "cd", 
    group_name=="2:RR" ~ "cde", 
    group_name=="3:RR" ~ "def", 
    group_name=="5:RR" ~ "ef", 
    group_name=="7:RR" ~ "f"
  ))


plot <- ggplot(data=df, aes(x=Generation, y=genot.freq.obs, groups=genotype)) + 
  # Make columns for the mean genotype frequency for each generation + genotype
  geom_col(aes(x=Generation, y=genot.freq.obs, fill=genotype), color="#000000", position=position_dodge()) + 
  # Add errorbars (currently standard deviation)
  geom_errorbar(aes(ymin=genot.freq.obs-se, ymax=genot.freq.obs+se), width=.2, position=position_dodge(.9))  + 
  # Add Tukey HSD groups
  geom_text(data=df, aes(x=Generation, y=genot.freq.obs+se+0.02, label=tukey_group), position=position_dodge(.9)) +
  theme_minimal() + 
  labs(title="Genotype frequency change over generations", 
       y="Observed genotype frequency")
  

plot
```

![](AlleleCompetitionDataAnalysis_files/figure-gfm/get-genos-tukey-groups-1.png)<!-- -->

This is the raw plot – we’ll, again, have to prettify it in Illustrator.

``` r
write.table(df, file="../data/GenotypeFrequenciesWithStandardErrorAndTukeyGroups.txt", sep="\t", eol="\n", quote=FALSE, row.names = FALSE)
```

<a id="3_2"></a>

### 3.2 HW Equilibrium test within generations

We will test for Hardy-Weinberg equilibrium within each generation using
a Chi-square goodness of fit test. If the p-value is below the alpha
threshold of 0.05, that would indicate that the genotypes observed are
not in equilibrium with the allele frequencies observed.

``` r
# start with original raw data and reshape it. 
HW_genos <- read.table("../data/AlleleCompetition_GenotypeCounts.txt", 
                       sep="\t", header=T)


str(HW_genos)
```

    ## 'data.frame':    26 obs. of  8 variables:
    ##  $ Generation: chr  "F1" "F1" "F2" "F2" ...
    ##  $ Cross     : chr  "A" "B" "A" "A" ...
    ##  $ Replicate : chr  "A" "B" "A1" "A2" ...
    ##  $ Full.Name : chr  "F1A" "F1B" "F2A1" "F2A2" ...
    ##  $ n         : int  88 87 87 85 84 87 82 84 84 85 ...
    ##  $ SS        : int  0 0 19 25 22 26 24 25 25 21 ...
    ##  $ SR        : int  88 87 45 44 51 48 35 38 44 49 ...
    ##  $ RR        : int  0 0 23 16 11 13 23 21 15 15 ...

Data manipulations done, let’s do the HW chiSq test for each generation
with its own R allele frequency (the classic method)

#### 3.4.a Within generation HWE

Read in and reshape data.

``` r
str(HW_genos)
```

    ## 'data.frame':    26 obs. of  8 variables:
    ##  $ Generation: chr  "F1" "F1" "F2" "F2" ...
    ##  $ Cross     : chr  "A" "B" "A" "A" ...
    ##  $ Replicate : chr  "A" "B" "A1" "A2" ...
    ##  $ Full.Name : chr  "F1A" "F1B" "F2A1" "F2A2" ...
    ##  $ n         : int  88 87 87 85 84 87 82 84 84 85 ...
    ##  $ SS        : int  0 0 19 25 22 26 24 25 25 21 ...
    ##  $ SR        : int  88 87 45 44 51 48 35 38 44 49 ...
    ##  $ RR        : int  0 0 23 16 11 13 23 21 15 15 ...

``` r
# Current data format is "long" data - each genotype gets its own row

HW_genos <- HW_genos %>%
  mutate(rfreq=(SR+(RR*2))/(2*n))

# Functions to get expected genotype counts
getSSexp <- function(rfreq, n) {
  sfreq = 1 - rfreq
  SSexp = (sfreq**2) * n
  return(SSexp)
}

getRRexp <- function(rfreq, n) {
  RRexp = (rfreq**2) * n
  return(RRexp)
}

getSRexp <- function(rfreq, n) {
  SRexp = (2*rfreq*(1-rfreq)) * n
  return(SRexp)
}

#Check to make sure we get expected results 
getSSexp(0.5, 100)
```

    ## [1] 25

``` r
getSRexp(0.5, 100)
```

    ## [1] 50

``` r
getRRexp(0.5, 100)
```

    ## [1] 25

``` r
## Add expected genotype counts to the DF 

HW_genos <- HW_genos %>%
  mutate(SSexp = round(getSSexp(rfreq, n), 2)) %>%
  mutate(SRexp = round(getSRexp(rfreq, n), 2)) %>%
  mutate(RRexp = round(getRRexp(rfreq, n), 2))
```

We are using the chi-square goodness-of-fit test.

The formula for $X^{2}$ is:

$\begin{equation} X^{2} = \frac{(n_{SS} - e_{SS})^{2}}{e_{SS}} + \frac{(n_{SR} - e_{SR})^{2}}{e_{SR}} + \frac{(n_{RR} - e_{RR})^{2}}{e_{RR}} \end{equation}$

``` r
chiSq = c(0)
for (i in seq(1, length(HW_genos$Full.Name), 1)) {
    x2 = 
           sum(
             ((HW_genos$SS[i] - HW_genos$SSexp[i])**2 / 
                HW_genos$SSexp[i]), 
              ((HW_genos$SR[i] - HW_genos$SRexp[i])**2 / 
                HW_genos$SRexp[i]), 
              ((HW_genos$RR[i] - HW_genos$RRexp[i])**2 / 
                HW_genos$RRexp[i]) 
           )
  chiSq = c(chiSq, x2) ## add each value to the chiSq vector
  
}

## Remove that pesky first value that we made when we initiated the variable
chiSq <- chiSq[-1]

HW_genos$x2 <- chiSq

#Calculate pvalues from the Chi-Squared distribution
pvaluesthisgen <- pchisq(q=as.numeric(HW_genos$x2), 2, lower.tail = FALSE)

HW_genos$x2.pvals <- pvaluesthisgen

HW_genos$Signif <- HW_genos$x2.pvals < 0.05

HW_genos[HW_genos$Signif,]
```

    ##   Generation Cross Replicate Full.Name  n SS SR RR rfreq SSexp SRexp RRexp x2
    ## 1         F1     A         A       F1A 88  0 88  0   0.5 22.00  44.0 22.00 88
    ## 2         F1     B         B       F1B 87  0 87  0   0.5 21.75  43.5 21.75 87
    ##       x2.pvals Signif
    ## 1 7.781132e-20   TRUE
    ## 2 1.282892e-19   TRUE

These results are the same as my Excel calculations. Other than the F1,
there aren’t any populations out of equilibrium. We certainly expect
that the F1 is out of equilibrium, because that population was not
panmictic - we forced cross-breeding of homozygous resistant and
homozygous susceptible.

However, we want to know if genotype frequencies are changing between
generations in a way that violated HWE, so we want to perform the
Chi-Square test using the expected genotype frequencies for this
generation if the allele frequency were the same as the previous
generation.

<a id="3_3"></a>

### 3.3 HW Equilibrium based on previous generation

Thus, we need to build a table that has observed HW genotype counts and
expected HW genotype counts for each population. I’ve calculated the
expected genotype frequencies for the current and last generations in
Excel, and put that into a tab delimited sheet.

``` r
genosHW_long <- read.table("../data/HWTableForR.txt", sep="\t", header=T)
str(genosHW_long)
```

    ## 'data.frame':    78 obs. of  9 variables:
    ##  $ genotype               : chr  "SS" "SR" "RR" "SS" ...
    ##  $ fullname               : chr  "F1A" "F1A" "F1A" "F1B" ...
    ##  $ n                      : int  88 88 88 88 88 88 87 87 87 85 ...
    ##  $ genot.count.obs        : int  0 88 0 0 88 0 19 45 23 25 ...
    ##  $ genot.freq.obs         : num  0 1 0 0 1 ...
    ##  $ genot.freq.exp.samegen : num  0.25 0.5 0.25 0.25 0.5 ...
    ##  $ genot.count.exp.samegen: num  22 44 22 22 44 ...
    ##  $ genot.freq.exp.lastgen : num  NA NA NA NA NA NA 0.25 0.5 0.25 0.25 ...
    ##  $ genot.count.exp.lastgen: num  NA NA NA NA NA ...

``` r
# Reshape data 

genosHW_wide <- genosHW_long %>% 
  pivot_wider(names_from = genotype, values_from = c(genot.freq.obs, genot.freq.exp.samegen, genot.freq.exp.lastgen, genot.count.obs, genot.count.exp.samegen, genot.count.exp.lastgen))

## The genot.count.exp fields (genotype count expected) are based on the previous generation's R allele frequency and the current generation's tested population size

## subset data 

genosHW_wide_lastgen <- genosHW_wide %>% select(
  fullname, n, genot.count.obs_SS, genot.count.obs_SR, genot.count.obs_RR, 
  genot.count.exp.lastgen_SS, genot.count.exp.lastgen_SR, genot.count.exp.lastgen_RR
) 
str(genosHW_wide_lastgen)
```

    ## tibble [26 × 8] (S3: tbl_df/tbl/data.frame)
    ##  $ fullname                  : chr [1:26] "F1A" "F1B" "F2A1" "F2A2" ...
    ##  $ n                         : int [1:26] 88 88 87 85 84 87 82 84 84 85 ...
    ##  $ genot.count.obs_SS        : int [1:26] 0 0 19 25 22 26 24 25 25 21 ...
    ##  $ genot.count.obs_SR        : int [1:26] 88 88 45 44 51 48 35 38 44 49 ...
    ##  $ genot.count.obs_RR        : int [1:26] 0 0 23 16 11 13 23 21 15 15 ...
    ##  $ genot.count.exp.lastgen_SS: num [1:26] NA NA 21.8 21.2 21 ...
    ##  $ genot.count.exp.lastgen_SR: num [1:26] NA NA 43.5 42.5 42 ...
    ##  $ genot.count.exp.lastgen_RR: num [1:26] NA NA 21.8 21.2 21 ...

``` r
chiSq = c(0)

for (i in seq(1, length(genosHW_wide_lastgen$fullname), 1)) {
    x2 = 
           sum(
             ((genosHW_wide_lastgen$genot.count.obs_SS[i] - genosHW_wide_lastgen$genot.count.exp.lastgen_SS[i])**2 / 
                genosHW_wide_lastgen$genot.count.exp.lastgen_SS[i]), 
              ((genosHW_wide_lastgen$genot.count.obs_SR[i] - genosHW_wide_lastgen$genot.count.exp.lastgen_SR[i])**2 / 
                genosHW_wide_lastgen$genot.count.exp.lastgen_SR[i]), 
              ((genosHW_wide_lastgen$genot.count.obs_RR[i] - genosHW_wide_lastgen$genot.count.exp.lastgen_RR[i])**2 / 
                genosHW_wide_lastgen$genot.count.exp.lastgen_RR[i]) 
           )
  chiSq = c(chiSq, x2) ## add each value to the chiSq vector
  
}

chiSq = chiSq[-1] ## Remove the blank that we used to initiate the vector
genosHW_wide_lastgen$x2 <- chiSq

#Calculate pvalues from the Chi-Squared distribution
pvalueslastgen <- pchisq(q=as.numeric(genosHW_wide_lastgen$x2), 2, lower.tail = FALSE)

genosHW_wide_lastgen$x2.pval <- pvalueslastgen

Signif <- genosHW_wide_lastgen$x2.pval < 0.05
genosHW_wide_lastgen$significant <- Signif

na.omit(genosHW_wide_lastgen[genosHW_wide_lastgen$significant == TRUE,])
```

    ## # A tibble: 5 × 11
    ##   fullname     n genot.c…¹ genot…² genot…³ genot…⁴ genot…⁵ genot…⁶    x2 x2.pval
    ##   <chr>    <int>     <int>   <int>   <int>   <dbl>   <dbl>   <dbl> <dbl>   <dbl>
    ## 1 F2A3        84        22      51      11    21      42     21     6.74 0.0344 
    ## 2 F5A1        87        38      40       9    27.2    42.8   17.0   8.18 0.0167 
    ## 3 F5A3        89        36      45       8    27.7    43.7   17.6   7.72 0.0211 
    ## 4 F7A3        87        52      32       3    37.7    38.9   10.4  11.9  0.00257
    ## 5 F7B2        89        53      33       3    40.4    39.0    9.57  9.34 0.00939
    ## # … with 1 more variable: significant <lgl>, and abbreviated variable names
    ## #   ¹​genot.count.obs_SS, ²​genot.count.obs_SR, ³​genot.count.obs_RR,
    ## #   ⁴​genot.count.exp.lastgen_SS, ⁵​genot.count.exp.lastgen_SR,
    ## #   ⁶​genot.count.exp.lastgen_RR

Happily, these values agree with what I calculated separately in Excel.
(Should make sure to include that Excel sheet in Supplementary Data).
(NAs are because the F1s don’t have a previous generation to compare
to.)

``` r
write.table(genosHW_wide_lastgen, file="../output/HW_ChiSqGoodnessOfFit_basedOnLastGen.txt", 
            row.names = F, quote = F, eol="\n", sep="\t")
```

<a id="3_4"></a>

### 3.4 Fisher’s combined probability

Finally, we want to perform a meta-analysis of these tests to get a
global p-value, assessing whether or not the evolving populations are in
Hardy-Weinberg Equilibrium. We’ll use Fisher’s Combined Probability test
for this. In effect, this is testing the null meta-hypothesis that the
null hypothesis was actually true for each and every one of the tests.
The alternative meta-hypothesis is that at least one of the separate
alternative hypotheses was true. In this case, a p-value below 0.005
means that we’re correct in saying that in at least one of these cases,
selection is happening.

<a id="3_5"></a>

#### 3.5.a Meta-analysis of previous generation HWE

``` r
P.fishers <- fisher(na.omit(genosHW_wide_lastgen$x2.pval))

P.fishers
```

    ## combined p-values with:      Fisher's method
    ## number of p-values combined: 24 
    ## test statistic:              88.13 ~ chi-square(df = 48) 
    ## adjustment:                  none 
    ## combined p-value:            0.0003683506

The pooled p-value is less than our alpha threshold of 0.05. This
indicates that we would be correct in rejecting at least one alternative
hypothesis.

#### 3.5.b Meta-analysis of within-generation HWE

This is more or less a positive control. Except for the F1s, all of the
populations were found to be in equilibrium based on their own R allele
frequencies. Fisher’s combined test here is asking the question “Are we
right to have rejected the alternative hypothesis (that the populations
are not in equilibrium) in every case?”

It’s worth noting that Fisher’s method amplifies statistical
significance while p-value adjustment diminishes statistical
significance. With Fisher’s method, if we have a group of p-values that
are marginal but not below our alpha threshold, it’s very possible to
get a global p-value that is below the alpha threshold. This would
probably reflect low power in our tests.

``` r
str(HW_genos)
```

    ## 'data.frame':    26 obs. of  15 variables:
    ##  $ Generation: chr  "F1" "F1" "F2" "F2" ...
    ##  $ Cross     : chr  "A" "B" "A" "A" ...
    ##  $ Replicate : chr  "A" "B" "A1" "A2" ...
    ##  $ Full.Name : chr  "F1A" "F1B" "F2A1" "F2A2" ...
    ##  $ n         : int  88 87 87 85 84 87 82 84 84 85 ...
    ##  $ SS        : int  0 0 19 25 22 26 24 25 25 21 ...
    ##  $ SR        : int  88 87 45 44 51 48 35 38 44 49 ...
    ##  $ RR        : int  0 0 23 16 11 13 23 21 15 15 ...
    ##  $ rfreq     : num  0.5 0.5 0.523 0.447 0.435 ...
    ##  $ SSexp     : num  22 21.8 19.8 26 26.9 ...
    ##  $ SRexp     : num  44 43.5 43.4 42 41.3 ...
    ##  $ RRexp     : num  22 21.8 23.8 17 15.9 ...
    ##  $ x2        : num  88 87 0.117 0.189 4.657 ...
    ##  $ x2.pvals  : num  7.78e-20 1.28e-19 9.43e-01 9.10e-01 9.74e-02 ...
    ##  $ Signif    : logi  TRUE TRUE FALSE FALSE FALSE FALSE ...

``` r
## We do not want to include F1 because they violate HW assumptions

HW_genos2 <- subset(HW_genos, Generation!="F1")
dim(HW_genos2)
```

    ## [1] 24 15

``` r
P.fishers <- fisher(na.omit(HW_genos2$x2.pvals))
P.fishers
```

    ## combined p-values with:      Fisher's method
    ## number of p-values combined: 24 
    ## test statistic:              18.077 ~ chi-square(df = 48) 
    ## adjustment:                  none 
    ## combined p-value:            0.9999738

Our global p value for the HW equilibrium test within generations is not
significant. This indicates that we would be incorrect if we had
rejected any of the null hypotheses.

### Save RData file for making publication plots

``` r
save.image(file="../ACData.RData")
```

### R Session Info

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23 ucrt)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19044)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] plyr_1.8.7       agricolae_1.3-5  colorspace_2.0-3 poolr_1.1-1     
    ## [5] ggrepel_0.9.1    ggplot2_3.3.6    tidyr_1.2.0      dplyr_1.0.9     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.9       lattice_0.20-45  assertthat_0.2.1 digest_0.6.29   
    ##  [5] utf8_1.2.2       mime_0.12        R6_2.5.1         AlgDesign_1.2.1 
    ##  [9] labelled_2.9.1   evaluate_0.16    highr_0.9        pillar_1.8.1    
    ## [13] rlang_1.0.5      rstudioapi_0.14  miniUI_0.1.1.1   Matrix_1.4-1    
    ## [17] combinat_0.0-8   rmarkdown_2.16   mathjaxr_1.6-0   labeling_0.4.2  
    ## [21] splines_4.2.1    stringr_1.4.1    questionr_0.7.7  munsell_0.5.0   
    ## [25] shiny_1.7.2      compiler_4.2.1   httpuv_1.6.5     xfun_0.32       
    ## [29] pkgconfig_2.0.3  mgcv_1.8-40      htmltools_0.5.3  tidyselect_1.1.2
    ## [33] tibble_3.1.8     fansi_1.0.3      withr_2.5.0      later_1.3.0     
    ## [37] MASS_7.3-57      grid_4.2.1       nlme_3.1-157     xtable_1.8-4    
    ## [41] gtable_0.3.1     lifecycle_1.0.2  DBI_1.1.3        magrittr_2.0.3  
    ## [45] scales_1.2.1     cli_3.3.0        stringi_1.7.8    farver_2.1.1    
    ## [49] promises_1.2.0.1 ellipsis_0.3.2   generics_0.1.3   vctrs_0.4.1     
    ## [53] klaR_1.7-1       tools_4.2.1      forcats_0.5.2    glue_1.6.2      
    ## [57] purrr_0.3.4      hms_1.1.2        fastmap_1.1.0    yaml_2.3.5      
    ## [61] cluster_2.1.3    knitr_1.40       haven_2.5.0
