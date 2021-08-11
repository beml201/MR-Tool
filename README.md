# MR-Tool Docs

### About:
An MR tool that will plot the regressions with a first plot and numerical output in an interactive html

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

