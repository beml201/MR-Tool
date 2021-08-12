library(plotly)

## These functions were made from the basis of: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4849733/
## Bowden et al.

##Recommended to be used on a merged dataframe of exposures and outcomes to make sure all the alleles are the correct way around
strand.update <- function(merged_df, expo_A1_col, expo_A2_col, expo_effect_col, out_A1_col, out_A2_col, out_effect_col){
  nmax <- nrow(merged_df)
  
  ## Flip the exposure associations so that they are all positive
  for(i in 1:nmax){
    if(merged_df[i,expo_effect_col] < 0){
      allele_switch <- merged_df[i,expo_A1_col]
      merged_df[i,expo_A1_col] <- merged_df[i,expo_A2_col]
      merged_df[i,expo_A2_col] <- allele_switch
      merged_df[i,expo_effect_col] <- -merged_df[i,expo_effect_col]
    }
  }
  
  ## Flip the outcome associations so that they are the same way around as the exposures
  for(i in 1:nmax){
    if(toupper(merged_df[i,expo_A1_col]) != toupper(merged_df[i,out_A1_col])){
      allele_switch <- merged_df[i,out_A1_col]
      merged_df[i,out_A1_col] <- merged_df[i,out_A2_col]
      merged_df[i,out_A2_col] <- allele_switch
      merged_df[i,out_effect_col] <- -merged_df[i,out_effect_col]
    }
  }
  return(merged_df)
}

weighted.median <- function(x, weights) {
  betas <- x[order(x)]
  wj <- weights[order(x)]
  wj <- wj/sum(wj) ## Standardized so the sum should == 1
  sj <- cumsum(wj)
  pj <- 100*(sj - wj/2)
  n <- max(which(pj< 50)) ## Index of median
  beta_wm <- betas[n] + (betas[n+1]-betas[n])*(50-pj[n])/(pj[n+1]-pj[n])
  return(beta_wm) 
}

boot.strap = function(x, y, x_se, y_se, weighting, n){
  med = NULL
  for(i in 1:n){
    beta_x = rnorm(length(x), mean=x, sd=x_se)
    beta_y = rnorm(length(y), mean=y, sd=y_se)
    beta_IV  = beta_y/beta_x
    med[i] = weighted.median(beta_IV, weighting)
    }
  return(sd(med)) 
  }

t.score <- function(beta,se){beta/se}
p.val <- function(t, df, dist = 't'){
  if(dist == 't'){2*pt(t,df)}else if(dist == 'norm'){2*pnorm(t, 0, 1)}
  }

IVW <- function(x, y, x_se, y_se, weighting = function(xw,yw) yw^-2){
  lm(y ~ 0 + x, weights =  weighting(x_se, y_se))
  }

Egger <- function(x, y, x_se, y_se, weighting = function(xw,yw) yw^-2){
  lm(y ~ x, weights = weighting(x_se, y_se))
  }

## Median done by calculating the beta using median.weighted and se using bootstrap
## use the t_score and p_val function to calculate them
median.regression <- function(x, y, x_se, y_se, weighting = function(xw,yw) (yw/xw)^-2){
  betaIV <- y/x
  betaIV_se <- weighting(x, y_se)
  beta_est <- weighted.median(betaIV, betaIV_se)
  se_est <- boot.strap(x, y, x_se, y_se, betaIV_se, 1000)
  t_est <- t.score(beta_est, se_est)
  p_est <- p.val(t_est, 1, dist = 'norm') ## df doesn't matter for normal dist
  list_out <- list('beta'=beta_est,'se'=se_est,'tval'=t_est,'p'=p_est)
  return(list_out)
  }

weights.penalised <- function(x, y, x_se, y_se, betaIVW){
  betaj <- y/x
  weights <- (y_se/x)^-2
  Qj <- weights*(betaj - betaIVW)^2
  qj <- pchisq(Qj, df = 1 , lower.tail = F)
  w_pen <- weights*pmin(1, 20*qj) ## pmin returns a vector of weights instead of just a single value
  return(w_pen)
  }

## Combine the regressions together
all.regressions <- function(x, y, x_se, y_se){
  ## Get all the regressions
  IVW_lm <- IVW(x, y, x_se, y_se)
  Egger_lm <- Egger(x, y, x_se, y_se)
  WM_lm <- median.regression(x, y, x_se, y_se)
  SM_lm <- median.regression(x, y, x_se, y_se, weighting = function(xw,yw) rep(1,length(xw)))
  PM_lm <- median.regression(x, y, x_se, y_se, weighting = function(xw,yw) weights.penalised(x, y, xw, yw, IVW_lm$coefficients[1]))
  
  ## Tidy the models that come from 'lm'
  IVW_summ <- summary(IVW_lm)$coefficients
  IVW_list <- list('beta'=IVW_summ[1],'se'=IVW_summ[2],'tval'=IVW_summ[3],'p'=IVW_summ[4])
  Egger_summ <- summary(Egger_lm)$coefficients
  Egger_list <- list('beta'=Egger_summ[2,1],'se'=Egger_summ[2,2],'tval'=Egger_summ[2,3],'p'=Egger_summ[2,4])
  
  ## Combine the models into a dataframe
  df_summary <- data.frame(rbind(IVW_list, Egger_list, WM_lm, SM_lm, PM_lm))
  df_summary$type <- c('IVW','MR-Egger','Weighted_Median','Simple_Median','Penalised_Weighted_Median')
  df_summary$intercept <- 0
  df_summary$intercept[which(df_summary$type == 'MR-Egger')] <- Egger_summ[1,1]
  df_summary <- df_summary[,c('type','beta','se','tval','p','intercept')]
  for(column in c('beta','se','tval','p')) df_summary[,column] <- as.numeric(df_summary[,column])
  return(df_summary)
}


make.plot <- function(x, y, x_se, y_se, regression_summary, title='Mendelian Randomisation Plot'){
  my.axes <- data.frame(x=c(min(x-x_se,0),max(x+x_se),rep(NA,length(x)-2)))
  my.colors <- c('red','blue','green','grey','purple')
  
  regression.names <- regression_summary$type
  regression.beta <- regression_summary$beta
  regression.intercept <- regression_summary$intercept
  
  fig1 <- plot_ly()
  fig1 <- fig1 %>% add_trace(x= x, y= y, type= 'scatter', name= 'Data Points',
                             error_x= list(array= 1.96*x_se, color= '#000000'),
                             error_y= list(array= 1.96*y_se, color= '#000000'))
  for(i in 1:nrow(regression_summary)){
    my.axes$y <- regression.intercept[i] + my.axes$x*regression.beta[i]
    fig1 <- fig1 %>% add_trace(data = my.axes,
                               x=~x,
                               y=~y, 
                               type = 'scatter', 
                               mode= 'lines',
                               name=regression.names[i],
                               legendgroup=toString(i),
                               color=my.colors[i])
  }
  fig1 <- fig1 %>% layout(
    xaxis= list(range= c(min(x-x_se,0),max(x+x_se))),
    yaxis= list(range= c(min(y-y_se),max(y+y_se)))
  )
  
  fig2 <- plot_ly() %>% add_trace()
  for(i in 1:nrow(regression_summary)){
    fig2 <- fig2 %>% add_trace(data= regression_summary[i,],
                               x= ~beta, 
                               y= ~type, 
                               type= 'scatter', 
                               error_x= list(array=~se*1.96),
                               name= regression.names[i],
                               legendgroup= toString(i),
                               showlegend= F,
                               color= my.colors[i])
  }
  
  my.reg.round <- regression_summary
  my.reg.round[,2:6] <- round(my.reg.round[,2:6],5)
  fig3 <- plot_ly(type="table",
                  header=list(values=names(my.reg.round),
                              fill = list(color = 'lightgrey')), 
                  cells=list(values=unname(my.reg.round)),
                  domain = list(x=c(0,1), y=c(0,0.4)))
  
  s1 <- subplot(fig1,fig2)
  fig <- subplot(s1,fig3,nrows=2) %>%
    layout(yaxis = list(domain=c(0.5,0.9)),
           yaxis3 = list(domain=c(0,0.5)))
  fig <- fig %>% layout(title = list(text=title,y=0.97))
  
}
