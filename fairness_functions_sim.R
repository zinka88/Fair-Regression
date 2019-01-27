# title: Fairness functions sim
# date: 07/17/2018
# updated: 01/01/2019
# description: Functions used for fairness simulation

# solve method using CVXR solve function -- catches any failed solve attempts
my_solve<-function(x,j,simn){
  tryCatch(expr={result<-solve(x)},
           error = function(e){
             message("Solver failed")
             message(paste("Sim:",j,"Datsetsize:",simn))
             return(NA)}
  )
}

# Net compensation penalty - penalize models with rev less than cost
penalty<-function(beta,grp,X_scale,n_grp,grp_cost_scale){
  avg_rev_grp<-((t(grp) %*% (X_scale %*% beta))/n_grp)
  fair<-(grp_cost_scale - avg_rev_grp)
}

# Mean residual difference penalty 
penalty2<-function(beta,grp,X_scale,n_grp,ref,n_ref,grp_cost_scale,ref_cost_scale){
  avg_rev_grp<-((t(grp) %*% (X_scale %*% beta))/n_grp)
  avg_rev_ref<-((t(ref) %*% (X_scale %*% beta))/n_ref)
  fair<-( (grp_cost_scale-ref_cost_scale) - (avg_rev_grp - avg_rev_ref) )^2
}

# rescale y values 
rescale<-function(y_scale,pred){
  newpred<-pred*attr(y_scale,'scaled:scale')+attr(y_scale,'scaled:center')
}

# get predictions for test dataset and rescale them
get_preds<-function(beta,X_test_scale,y_scale){
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
grp_rev<-function(predy,grpvar,grp2var){
  tmp<-as.data.frame(cbind(grpvar, predy,grp2var))
  names(tmp)<-c('grp','pred','grp2')
  grp_rev<-mean(tmp[tmp$grp==1,'pred'])
  ref_rev<-mean(tmp[tmp$grp==0,'pred'])
  grp2_rev<-mean(tmp[tmp$grp2==1,'pred'])
  rev_list<-list("grp"=grp_rev, "ref"=ref_rev,"grp2"=grp2_rev)
  return(rev_list)
}

# calculated net compensation
overunder<-function(predy,grpvar,grp2var,grp_cost,ref_cost,grp2_cost){ 
  rev<-grp_rev(predy,grpvar,grp2var)
  grp_rev<-rev$grp
  ref_rev<-rev$ref
  grp2_rev<-rev$grp2
  grp_ou<-grp_rev - grp_cost
  ref_ou<-ref_rev - ref_cost
  grp2_ou<-grp2_rev - grp2_cost
  ou_list<-list("grp"=grp_ou, "ref"=ref_ou, "grp2"=grp2_ou)
  return(ou_list)
}

# predicted ratio
predratio<-function(y,predy,grpvar,grp2var,grp_cost,ref_cost,grp2_cost){
  rev<-grp_rev(predy,grpvar,grp2var)
  grp_rev<-rev$grp
  ref_rev<-rev$ref
  grp2_rev<-rev$grp2
  grp_pr<-grp_rev/grp_cost
  ref_pr<-ref_rev/ref_cost
  grp2_pr<-grp2_rev/grp2_cost
  pr_list<-list("grp"=grp_pr, "ref"=ref_pr, "grp2"=grp2_pr)
  return(pr_list)
}

# correlation btw grp and error
grpcorr<-function(y,predy, grpvar){
  cval = cor(grpvar, y-predy)
}

# covariance btw grp and error
grpcov<-function(y,predy,grpvar){
  cval_cov = cov(grpvar, y-predy)
}

# call evaluation metrics and return df with metrics
all_metrics<-function(y,ypred,grpvar,grp2var, model){
  temp<-data.frame(y,grpvar,grp2var)
  grp_cost<-mean(temp[temp$grpvar==1,'y'])
  ref_cost<-mean(temp[temp$grpvar==0,'y'])
  grp2_cost<-mean(temp[temp$grp2var==1,'y'])
  r2<-round(rsquared(y,ypred),3)
  mse<-round(mse(y,ypred),3)
  ou<-overunder(ypred,grpvar,grp2var,grp_cost,ref_cost,grp2_cost)
  ou_grp<-round(ou$grp,3)
  ou_ref<-round(ou$ref,3)
  ou_grp2<-round(ou$grp2,3)
  pr<-predratio(y,ypred,grpvar,grp2var,grp_cost,ref_cost,grp2_cost)
  pr_grp<-round(pr$grp,3)
  pr_ref<-round(pr$ref,3)
  pr_grp2<-round(pr$grp2,3)
  gc<-round(grpcorr(y,ypred,grpvar),3)
  gcov<-round(grpcov(y,ypred,grpvar),3)
  name<-model
  # create data frame to return (can combine later for print)
  df<-cbind(model, r2, mse, ou_grp, ou_ref, pr_grp, pr_ref, gc, gcov, ou_grp2, pr_grp2)
  return(df)
}