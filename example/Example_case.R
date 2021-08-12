## Source the MRVis.R file,, which contains all the functions
## (this will also require the plotly package to be installed)
source('../src/MRVis.R')

## Read your data in
dataframe <- read.table('example_data.txt',header=T)

## Ensure that the strands are flipped appropriately
## (this does not need to be done if you have already done it)
df_updated <- strand.update(dataframe,
                        'exposure_effect_allele',
                        'exposure_non_effect_allele',
                        'exposure_beta',
                        'outcome_effect_allele',
                        'outcome_non_effect_allele',
                        'outcome_beta')

## Run the all.regressions function which will output a dataframe of the regressions
my.regression <- all.regressions(df_updated$exposure_beta,
                                 df_updated$outcome_beta,
                                 df_updated$exposure_se,
                                 df_updated$outcome_se)

## Run the make.plot to make the interactive plot (should appear in the viewer in Rstudio)
## You can input your own regressions here (using the same format) and the plot should still work
fig <- make.plot(df_updated$exposure_beta,
                 df_updated$outcome_beta,
                 df_updated$exposure_se,
                 df_updated$outcome_se,
                 my.regression,
                 'Title')

## This can be used to save the document as an interactive html file
htmlwidgets::saveWidget(as_widget(fig), 'MR-Tool_Interactive_Plot.html')
