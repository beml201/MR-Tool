source('../src/MRVis.R')

#####
dataframe <- read.table('example_data.txt',header=T)

## Merge Data and flip -ve effect alleles on the exposures and make sure the associations are the correct way around
## (exposure == TSH/FT4 SNPs, outcome == summary SNPs)
df_updated <- strand.update(dataframe,
                        'exposure_effect_allele',
                        'exposure_non_effect_allele',
                        'exposure_beta',
                        'outcome_effect_allele',
                        'outcome_non_effect_allele',
                        'outcome_beta')

## Adjust SEs as 1.96*se for 95% CI of normalised values??
#df_updated$exposure_se <- 1.96*df_updated$exposure_se
#df_updated$outcome_se <- 1.96*df_updated$outcome_se

my.regression <- all.regressions(df_updated$exposure_beta,
                                 df_updated$outcome_beta,
                                 df_updated$exposure_se,
                                 df_updated$outcome_se)

fig <- make.plot(df_updated$exposure_beta,
                 df_updated$outcome_beta,
                 df_updated$exposure_se,
                 df_updated$outcome_se,
                 my.regression,
                 'Title')


htmlwidgets::saveWidget(as_widget(fig), 'MR-Tool_Interactive_Plot.html')
