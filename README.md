# MR-Tool Docs

### About:
This tool should be able to run MR regressions and also output them to an interactive plot that has both the regressions and the forest plot available for the user to see (and highlight specific ones)
The aim is to allow some interactivity to help explore the data and compare the given MR approaches with each other (ones can be clicked on in legend to be removed or added). While this can be used as a replacement to creating the forest plots and MR plots yourself, it may not go into sufficient detail for the analysis you have requested

---

## Quick Start:
<i>Please see example/Example_case.R for working example code.</i>
Ensure you have plotly installed:
```R
install.packages("plotly")
```
Make a dataframe of your exposure and outcome SNPs and read it into R
```R
df <- read.table("Your_File_location.tsv",header=T)
```
Update the strands to make sure that they're all the right way around (your exposure SNP effects will all be positive, and your outcome SNPs will be flipped to be inline with them)
Use strings for the required column headers ( you will need your exposure alleles, A1 and A2, and exposure effect sizes, as well as outcome alleles and effect sizes)
```R
df <- df %>% strand.update('x_ea','x_nea','x_beta','y_ea','y_nea','y_beta')
```
You can then get the MR regressions using (use $ to get columns here):
```R
my.regressions <- all.regressions(df$x_beta,df$ybeta,
                                  df$x_se, df$y_se)
```
Then run the plot maker function
```R
my.fig <- make.plot(df$x_beta,df$ybeta,
                    df$x_se, df$y_se,
                    my.regressions)
```
You can view the interactive plot in Rstudio using
```R
my.fig
```
Or you can save it as an html file (this can be opened in a browser and interacted with accordingly)
```R
htmlwidgets::saveWidget(as_widget(my.fig), 'MR-Tool_Interactive_Plot.html')
```
(the htmlwidgets package should have been installed with your plotly install)
If a static image is required, this can be done by clicking the camera icon <image src="https://cdn1.iconfinder.com/data/icons/ios-11-glyphs/30/camera-512.png" width="14" height="14"> within the interactive

---

### Further Information

#### Included functions(required arguments):
- [strand.update(merged_df, expo_A1_col, expo_A2_col, expo_effect_col, out_A1_col, out_A2_col, out_effect_col)](#strandupdate)
- [all.regressions(x, y, x_se, y_se)](#allregressions)
- [make.plot(x, y, x_se, y_se, regression_summary)](#makeplot)
- [IVW(x, y, x_se, y_se)](#ivw)
- [Egger(x, y, x_se, y_se)](#egger)
- [median.regression(x, y, x_se, y_se)](#medianregression)
- [weighted.median(x, weights)](#weightedmedian)
- [weights.penalised(x, y, x_se, y_se, betaIVW)](#weightspenalised)
- [t.val(beta, se)](#tval)
- [p.val(t, df)](#pval)
- [boot.strap(x, y, x_se, y_se, weighting, n)](#bootstrap)

#### strand.update
|Arguments||
|---|---|
|merged_df |dataframe: including the known alleles, effect sizes and standard errors|
|expo_A1_col |string: name of the exposure allele column|
|expo_A2_col |string: name of the exposure non-effect allele column|
|expo_effect_col |string: name of the exposure effect (beta) column|
|out_A1_col |string: name of the outcome allele column|
|out_A2_col |string: name of the outcome allele column|
|out_effect_col |string: name of the outcome effect (beta) column|

This function will return a modified version fo the merged dataframe that is input with an updated SNP list to make sure that all the exposure alleles have positive effect sizes and the outcome alleles are correctly aligned to the exposures.
This may also be considered "harmonization".

#### all.regressions
|Arguments||
|---|---|
|x|vector: exposure effect size (ensure same length as y)|
|y|vector: outcome effect size (ensure same length as x)|
|x_se|vector: exposure standard error for each SNP|
|y_se|vector: outcome standard error for each SNP|

This runs all the different regressions at once: IVW, MR-Egger, simple median, weighted median and penalised weighted median.
Returns a dataframe of all the regressions.

#### make.plot
|Arguments||
|---|---|
|x|vector: exposure effect size (ensure same length as y)|
|y|vector: outcome effect size (ensure same length as x)|
|x_se|vector: exposure standard error for each SNP|
|y_se|vector: outcome standard error for each SNP|
|regression_summary| dataframe: a dataframe of regressions|
|title| string: The title of the final interactive, <i>default</i>: "Mendelian Randomisation Plot"|

regression_summary has columns: type, beta, se, tval, p and intercept
<i>all.regressions</i> will output the appropriate dataframe format to be used with make.plot
If desired, the dataframe can be modified (or created) to plot different regression lines in the interactive output (type, beta, se and intercept required).
Returns a plotly interactive figure.

#### IVW
|Arguments||
|---|---|
|x|vector: exposure effect size (ensure same length as y)|
|y|vector: outcome effect size (ensure same length as x)|
|x_se|vector: exposure standard error for each SNP|
|y_se|vector: outcome standard error for each SNP|
|weighting| function(xw,yw): weighting formula for regression, <i>default</i>: yw<sup>-2</sup>|

Produces a linear regression model (can be viewed using <i>summary()</i>) in line with the IVW method
<i>weighting</i> can be altered if desired (would no longer be inverse-variance weighted)

#### Egger
|Arguments||
|---|---|
|x|vector: exposure effect size (ensure same length as y)|
|y|vector: outcome effect size (ensure same length as x)|
|x_se|vector: exposure standard error for each SNP|
|y_se|vector: outcome standard error for each SNP|
|weighting|function(xw,yw): weighting formula for regression, <i>default</i>: yw<sup>-2</sup>|

Produces a linear regression model in line with MR-Egger method
<i>weighting</i> can be altered if desired (would no longer be inverse-variance weighted)

#### median.regression
Arguments: x, y, x_se, y_se, weighting
|Arguments||
|---|---|
|x|vector: exposure effect size (ensure same length as y)|
|y|vector: outcome effect size (ensure same length as x)|
|x_se|vector: exposure standard error for each SNP|
|y_se|vector: outcome standard error for each SNP|
|weighting|function(xw,yw): weighting formula for regression, <i>default</i>: yw<sup>-2</sup>|

Produces a median regression model based on <i>weighted.median</i>
<i>weighting</i> can be altered (which it is within <i>all.regresions</i>)
Uses <i>boot.strap</i> to calculate the standard error of the betas and <i>t.val</i> and <i>p.val</i> to calculate further statistics based on the bootstrapped errors.
Returns a list of beta, se, tval and p.

#### weighted.median
Arguments: x, weights
|Arguments||
|---|---|
|x|vector: exposure effect size (ensure same length as y)|
|weights|vector: weights for median (ensure same length as x)|

Returns a weighted median

#### weights.penalised
Arguments: x, y, x_se, y_se, betaIVW
|Arguments||
|---|---|
|x|vector: exposure effect size (ensure same length as y)|
|y|vector: outcome effect size (ensure same length as x)|
|x_se|vector: exposure standard error for each SNP|
|y_se|vector: outcome standard error for each SNP|
|betaIVW|float: IVW coefficient|

Returns a vector of penalised weights

#### t.val
Arguments: beta, se
|Arguments||
|---|---|
|beta|float: coefficient of regression|
|se|float: standard error of regression|

Applies beta/se to calculate a t-value (used for bootstrapped functions, ie median regressions)

#### p.val
Arguments: t, df, dist
|Arguments||
|---|---|
|t|float: t-value of rgeression|
|df|float: degrees of freedom for the chi-squared distribution|
|dist|string: distribution to use, <i>default</i>; "t"|

Calculates the p-value of a given distribution given a certain t-value
set <i>dist</i> to "t" for the t-distribution and "norm" to the normal distribution (N(0,1))

#### boot.strap
Arguments: x, y, x_se, y_se, weighting, n
|Arguments||
|---|---|
|x|vector: exposure effect size (ensure same length as y)|
|y|vector: outcome effect size (ensure same length as x)|
|x_se|vector: exposure standard error for each SNP|
|y_se|vector: outcome standard error for each SNP|
|weighting|vector: weights for the <i>weighted.median</i> calculation when bootstrapping|
|n|integer: number of times to repeat for bootstrapping standard errors|

Estimates the standard error of the given weightings for the <i>median.regression</i> function.