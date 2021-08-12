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
If a static image is required, this can be done by clicking the camera icon <image src="https://cdn1.iconfinder.com/data/icons/ios-11-glyphs/30/camera-512.png" style="width:1em;height:1em"> within the interactive

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
Arguments: merged_df, expo_A1_col, expo_A2_col, expo_effect_col, out_A1_col, out_A2_col, out_effect_col
|||
|---|---|
|merged_df |dataframe: including the known alleles, effect sizes and standard errors|
|expo_A1_col |string: name of the exposure allele column|
|expo_A2_col |string: name of the exposure non-effect allele column|
|expo_effect_col |string: name of the effect (beta) column|
|out_A1_col |string: name of the outcome allele column|
|out_A2_col |string: name of the outcome allele column|
|out_effect_col |string: name of the outcome effect (beta) column|

#### all.regressions
Arguments: x, y, x_se, y_se

#### make.plot
Arguments: x, y, x_se, y_se, regression_summary, title

#### IVW
Arguments: x, y, x_se, y_se, weighting

#### Egger
Arguments: x, y, x_se, y_se, weighting

#### median.regression
Arguments: x, y, x_se, y_se, weighting

#### weighted.median
Arguments: x, weights

#### weights.penalised
Arguments: x, y, x_se, y_se, betaIVW

#### t.val
Arguments: beta, se
Applies beta/se to calculate a t-value (used for bootstrapped functions, ie median regressions)

#### p.val
Arguments: t, df, dist
Calculates the p-value of a given distribution given a certain t-value
set 'dist' to 't' for the t-distribution and 'norm' to the normal distribution (N(0,1))

#### boot.strap
Arguments: x, y, x_se, y_se, weighting, n

