# title: Fairness functions
# author: Anna Zink
# created: 07/17/2018
# updated: 01/25/2019 
# description: Functions used to run fairness analysis 

# Net compensation penalty - penalize models with rev less than cost
penalty<-function(beta){
  avg_rev_mh<-((t(grp) %*% (X_scale %*% beta))/n_grp)
  fair<-(mhsud_cost_scale - avg_rev_mh)
}

# Mean residual difference penalty 
penalty2<-function(beta){
  avg_rev_mh<-((t(grp) %*% (X_scale %*% beta))/n_grp)
  avg_rev_ref<-((t(ref) %*% (X_scale %*% beta))/n_ref)
  fair<-( (mhsud_cost_scale-ref_cost_scale) - (avg_rev_mh - avg_rev_ref) )^2
}

# rescale y values 
rescale<-function(y_scale,pred){
  newpred<-pred*attr(y_scale,'scaled:scale')+attr(y_scale,'scaled:center')
}

# get predictions for test dataset and rescale them
get_preds<-function(beta){
  pred_scaled<-as.matrix(X_test_scale) %*% beta
  pred<-as.data.frame(rescale(y_scale, pred_scaled))
  return(pred$V1)
}

# r-squared
rsquared<-function(y,predy){
  SSR = sum((y-predy)^2)
  SST = sum((y-mean(y))^2)
  R2 = 1-SSR/SST
  return(R2)
}

# mse 
mse<-function(y,predy){
  SSR = sum((y-predy)^2)
  MSE = SSR/length(y)
  return(MSE)
}

# calculate average revenue for the group 
grp_rev<-function(predy){
  tmp<-as.data.frame(cbind(flag_mh, predy))
  names(tmp)<-c('mh','pred')
  grp_rev<-mean(tmp[tmp$mh==1,'pred'])
  ref_rev<-mean(tmp[tmp$mh==0,'pred'])
  rev_list<-list("grp"=grp_rev, "ref"=ref_rev)
  return(rev_list)
}

# calculated net compensation
overunder<-function(predy){ 
  rev<-grp_rev(predy)
  mhsud_rev<-rev$grp
  ref_rev<-rev$ref
  grp_ou<-mhsud_rev - mhsud_cost
  ref_ou<-ref_rev - ref_cost
  ou_list<-list("grp"=grp_ou, "ref"=ref_ou)
  return(ou_list)
}

# predicted ratio
predratio<-function(y,predy) {
  rev<-grp_rev(predy)
  mhsud_rev<-rev$grp
  ref_rev<-rev$ref
  grp_pr<-mhsud_rev/mhsud_cost
  ref_pr<-ref_rev/ref_cost
  pr_list<-list("grp"=grp_pr, "ref"=ref_pr)
  return(pr_list)
}

# corr btw grp and error
grpcorr<-function(y,predy){
  cval = cor(flag_mh, y-predy)
}

# cov btw grp and error
grpcov<-function(y,predy){
  cval_cov = cov(flag_mh, y-predy)
}

# call evaluation metrics and return df with metrics
all_metrics<-function(y,ypred,model){
  r2<-round(rsquared(y,ypred),3)
  mse<-round(mse(y,ypred),3)
  ou<-overunder(ypred)
  ou_grp<-round(ou$grp,3)
  ou_ref<-round(ou$ref,3)
  pr<-predratio(y,ypred)
  pr_grp<-round(pr$grp,3)
  pr_ref<-round(pr$ref,3)
  gc<-round(grpcorr(y,ypred),3)
  gcov<-round(grpcov(y,ypred),3)
  name<-model
  # create data frame to return (can combine later for print)
  df<-cbind(model, r2, mse, ou_grp, ou_ref, pr_grp, pr_ref, gc, gcov)
  return(df)
}

# load coefficients - upload coefficients for measures
load_coef<-function(file) {
  beta_tmp<-read.csv(file)
  beta<-beta_tmp$V1
  pred_scaled<-X2_scale %*% beta
  pred<-rescale(pred_scaled)
  return(pred) 
}


