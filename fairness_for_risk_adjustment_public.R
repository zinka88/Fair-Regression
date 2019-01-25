# title: "Fairness for Risk Adjustment Public"
# author: Anna Zink
# created: 08/13/2018
# updated: 01/25/2019
# description: Analysis code for simulated analysis data 

## Load required packages
library(dplyr)
library(CVXR)

set.seed(1234)

# source functions to be used for analysis 
source('fairness_functions.R')

# read in sample data 
train<-read.csv('jasa_simulated_data.csv')

# create 5-folds in train data 
nfolds<-5
train$folds<-cut(seq(1,nrow(train)),breaks=nfolds,labels=FALSE)
folds<-train$folds

# create empty dataset to store data from each fold
n<-nrow(train)/nfolds
tmp_data<-data.frame(fold=integer(n), y=double(n), pred_ols=double(n), pred_1=double(n), 
                     pred_2_1=double(n), pred_2_2=double(n), pred_2_3=double(n), pred_2_4=double(n), 
                     pred_3_1=double(n), pred_3_2=double(n), pred_3_3=double(n), pred_3_4=double(n), 
                     pred_4_1=double(n), pred_4_2=double(n), pred_4_3=double(n), pred_4_4=double(n),
                     pred_5_1=double(n), pred_5_2=double(n), pred_5_3=double(n), pred_5_4=double(n))

# create empty final dataset
final_data<-tmp_data[0,]

# loop through each fold to run the analysis
for (i in 1:nfolds){

  ### PREPROCESS DATA ###
  index<-which(folds==i)
  test_i<-train[index,]
  train_i<-train[-index,]
  
  # define design matrix X and outcome y for train and test datasets
  y<-train_i$totpay
  X<-train_i[,!(names(train_i) %in% c('totpay','mh','folds'))]
  
  # indicators for whether observation has MHSUD (grp) or not (ref)
  grp<-train_i$mh
  ref<-!(train_i$mh)
  n_grp<-sum(grp)
  n_ref<-sum(ref)
  
  y_test<-test_i$totpay
  X_test<-test_i[,!(names(test_i) %in% c('totpay','mh','folds'))]
  
  # Scale data before optimizing 
  y_scale<-scale(y)
  X_scale<-scale(X)
  
  # calculate avg costs for those with/without MHSUD & scale 
  mhpay<-train_i[train_i$mh==1,'totpay']
  mhsud_cost<-mean(mhpay)
  refpay<-train_i[train_i$mh==0,'totpay']
  ref_cost<-mean(refpay)
  mhsud_cost_scale<-(mhsud_cost-mean(y))/sd(y)
  ref_cost_scale<-(ref_cost-mean(y))/sd(y)
  
  # scale test data based on the train dataset (scale used for predictions)
  X_test_scale<-sweep(sweep(X_test, 2L,attr(X_scale, 'scaled:center')), 2, attr(X_scale, 'scaled:scale'), "/")
  y_test_scale<-(y_test-mean(y))/sd(y)
  
  ### ESTIIMATION METHODS ###

  # OLS
  model_ols<-lm(y_scale~X_scale+0)
  beta_ols = as.matrix(coef(model_ols))
  
  # set up preliminaries for CVXR package 
  k<-length(beta_ols)
  beta<-Variable(k)
  loss<-sum((y_scale-X_scale %*% beta)^2)
  obj<-loss
  
  # Estimation Method 1: Constrained Regression
  prob<-Problem(Minimize(obj),list((t(grp) %*% (X_scale %*% beta))/n_grp==mhsud_cost_scale))
  result<-solve(prob)
  beta_1<-result$getValue(beta)
  print(paste("constrained reg:",result$status))
  
  # Estimation Method 2: Weighted Constrained Regression
  alphas = c(.2, .4, .6, .8)
  ols_avg<-(t(grp) %*% (X_scale %*% cbind(beta_ols)))/n_grp
  wcosts<-(1-alphas)*mhsud_cost_scale+alphas*as.vector(ols_avg)
  constraint<-list((t(grp) %*% (X_scale %*% beta))/n_grp == wcosts)
  constraint<-lapply(wcosts, function(x) list((t(grp) %*% (X_scale %*% beta))/n_grp == x))
  
  prob<-Problem(Minimize(obj),constraint[[1]])
  result<-solve(prob)
  beta_2_1<-result$getValue(beta)
  print(paste("weighted avg .2", result$status))
  
  prob<-Problem(Minimize(obj),constraint[[2]])
  result<-solve(prob)
  beta_2_2<-result$getValue(beta)
  print(paste("weighted avg .4",result$status))
  
  prob<-Problem(Minimize(obj),constraint[[3]])
  result<-solve(prob)
  beta_2_3<-result$getValue(beta)
  print(paste("weighted avg .6", result$status))
  
  prob<-Problem(Minimize(obj),constraint[[4]])
  result<-solve(prob)
  beta_2_4<-result$getValue(beta)
  print(paste("weighted avg .8", result$status))
  
  # Estimation Method 3: Covariance
  cstar = cov(grp, y_scale-predict(model_ols))
  share_grp<-n_grp/length(y)
  share_ref<-n_ref/length(y)
  
  prob<-Problem(Minimize(obj),list(share_ref*(sum(t(grp) %*% (y_scale - X_scale %*% beta))) - share_grp*(sum(t(ref) %*% (y_scale - X_scale %*% beta))) < .2*cstar))
  result<-solve(prob)
  beta_3_1<-result$getValue(beta)
  print(paste("cov .2", result$status))
  
  prob<-Problem(Minimize(obj),list(share_ref*(sum(t(grp) %*% (y_scale - X_scale %*% beta))) - share_grp*(sum(t(ref) %*% (y_scale - X_scale %*% beta))) < .4*cstar))
  result<-solve(prob)
  beta_3_2<-result$getValue(beta)
  print(paste("cov .4", result$status))
  
  prob<-Problem(Minimize(obj),list(share_ref*(sum(t(grp) %*% (y_scale - X_scale %*% beta))) - share_grp*(sum(t(ref) %*% (y_scale - X_scale %*% beta))) < .6*cstar))
  result<-solve(prob)
  beta_3_3<-result$getValue(beta)
  print(paste("cov .6", result$status))
  
  prob<-Problem(Minimize(obj),list(share_ref*(sum(t(grp) %*% (y_scale - X_scale %*% beta))) - share_grp*(sum(t(ref) %*% (y_scale - X_scale %*% beta))) < .8*cstar))
  result<-solve(prob)
  beta_3_4<-result$getValue(beta)
  print(paste("cov .8", result$status))
  
  # Estimation Method 4: Net Compensation
  prob<-Problem(Minimize(loss+5000*penalty(beta)))
  result<-solve(prob)
  beta_4_1<-result$getValue(beta)
  print(paste("net compensation 5000", result$status))
  
  prob<-Problem(Minimize(loss+10000*penalty(beta)))
  result<-solve(prob)
  beta_4_2<-result$getValue(beta)
  print(paste("net compensation 10000", result$status))
  
  prob<-Problem(Minimize(loss+20000*penalty(beta)))
  result<-solve(prob)
  beta_4_3<-result$getValue(beta)
  print(paste("net compensation 20000", result$status))

  prob<-Problem(Minimize(loss+30000*penalty(beta)))
  result<-solve(prob)
  beta_4_4<-result$getValue(beta)
  print(paste("net compensation 30000", result$status))
    
  # Estimation Method 5: Mean Residual Difference
  prob<-Problem(Minimize(loss+5000*penalty2(beta)))
  result<-solve(prob)
  beta_5_1<-result$getValue(beta)
  print(paste("mrd 5000", result$status))
  
  prob<-Problem(Minimize(loss+10000*penalty2(beta)))
  result<-solve(prob)
  beta_5_2<-result$getValue(beta)
  print(paste("mrd 10000", result$status))
  
  prob<-Problem(Minimize(loss+20000*penalty2(beta)))
  result<-solve(prob)
  beta_5_3<-result$getValue(beta)
  print(paste("mrd 20000", result$status))
  
  prob<-Problem(Minimize(loss+30000*penalty2(beta)))
  result<-solve(prob)
  beta_5_4<-result$getValue(beta)
  print(paste("mrd 30000", result$status))
   
  # create data for results from this fold
  tmp_data$fold<-i
  tmp_data$y<-y_test
  tmp_data$grp<-test_i$mh
  
  # get predictions for test_i dataset 
  betas<-list(beta_ols, beta_1, beta_2_1, beta_2_2, beta_2_3, beta_2_4, beta_3_1, beta_3_2, beta_3_3, beta_3_4, beta_4_1, beta_4_2, beta_4_3,beta_4_4,beta_5_1, beta_5_2, beta_5_3, beta_5_4)
  test_preds<-lapply(betas, get_preds)

  # save predictions in tmp_data
  tmp_data$pred_ols<-test_preds[1][[1]]
  tmp_data$pred_1<-test_preds[2][[1]]
  tmp_data$pred_2_1<-test_preds[3][[1]]
  tmp_data$pred_2_2<-test_preds[4][[1]]
  tmp_data$pred_2_3<-test_preds[5][[1]]
  tmp_data$pred_2_4<-test_preds[6][[1]]
  tmp_data$pred_3_1<-test_preds[7][[1]]
  tmp_data$pred_3_2<-test_preds[8][[1]]
  tmp_data$pred_3_3<-test_preds[9][[1]]
  tmp_data$pred_3_4<-test_preds[10][[1]]
  tmp_data$pred_4_1<-test_preds[11][[1]]
  tmp_data$pred_4_2<-test_preds[12][[1]]
  tmp_data$pred_4_3<-test_preds[13][[1]]
  tmp_data$pred_4_4<-test_preds[14][[1]]
  tmp_data$pred_5_1<-test_preds[15][[1]]
  tmp_data$pred_5_2<-test_preds[16][[1]]
  tmp_data$pred_5_3<-test_preds[17][[1]]
  tmp_data$pred_5_4<-test_preds[18][[1]]
  
  # save tmp_data in final dataset
  final_data<-rbind(final_data, tmp_data)
}

# Get evaluation measures for each estimation method 

# define data needed for metrics
test_y<-final_data$y
flag_mh<-final_data$grp
mhsud_cost<-mean(final_data[final_data$grp==1,'y'])
ref_cost<-mean(final_data[final_data$grp==0,'y'])

m1<-all_metrics(test_y, final_data$pred_ols, 'ols')
m2<-all_metrics(test_y, final_data$pred_1, 'average.constrained')
m3<-all_metrics(test_y, final_data$pred_2_1, 'weighted.avg.constrained alpha=.2')
m4<-all_metrics(test_y, final_data$pred_2_2, 'weighted.avg.constrained alpha=.4')
m5<-all_metrics(test_y, final_data$pred_2_3, 'weighted.avg.constrained alpha=.6')
m6<-all_metrics(test_y, final_data$pred_2_4, 'weighted.avg.constrained alpha=.8')
m7<-all_metrics(test_y, final_data$pred_3_1, 'covariance m=.2')
m8<-all_metrics(test_y, final_data$pred_3_2, 'covariance m=.4')
m9<-all_metrics(test_y, final_data$pred_3_3, 'covariance m=.6')
m10<-all_metrics(test_y, final_data$pred_3_4, 'covariance m=.8')
m11<-all_metrics(test_y, final_data$pred_4_1, 'net.compensation l=5000')
m12<-all_metrics(test_y, final_data$pred_4_2, 'net.compensation l=10000')
m13<-all_metrics(test_y, final_data$pred_4_3, 'net.compensation l=20000')
m14<-all_metrics(test_y, final_data$pred_4_4, 'net.compensation l=30000')
m15<-all_metrics(test_y, final_data$pred_5_1, 'mean.residual.difference l=5000')
m16<-all_metrics(test_y, final_data$pred_5_2, 'mean.residual.difference l=10000')
m17<-all_metrics(test_y, final_data$pred_5_3, 'mean.residual.difference l=20000')
m18<-all_metrics(test_y, final_data$pred_5_4, 'mean.residual.difference l=30000')

# create final dataset with measures for each estimation method - Table 1 in paper
metrics<-as.data.frame(rbind(m1, m2, m3, m4, m5, m6, m7, m8,  m9, m10, m11, m12, m13, m14, m15, m16, m17, m18))


