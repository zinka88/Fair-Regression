# title: Run simulation 
# date: 11/4/2018
# updated: 01/27/2019
# description: Function to simulation 5-fold fairness analysis

# Function: Sim 
# Parameters: 
#  scenario - simulation scenario being run (1,2,3)
#  df - data frame name
#  simn - sample size (1000, 10000)
#  nsims - number of simulations
#  yvar - outcome variable Y 
#  avar - sensitive group
#  avar2 - second sensitive group
#  includea - include avar in the prediction
#  m_cov - list of parameters for the covariance method
#  lambdas_nc - list of parameters for the net compensation penalty 
#  lambdas_nc_constraint - list of parameters for the net compensation constraint 
#  lambdas_mrd - list of parameters for the mean residual difference penalty 

sim<-function(scenario,df, simn, nsims, yvar, avar, avar2,includea=0,m_cov=c(.2),lambdas_nc=c(100,1000,5000), lambdas_nc_constraint=c(.2, .6, 1),lambdas_mrd=c(100,1000,5000)) {
  
  # create dataset to store results of simulation 
  simresults<-data.frame(sim=integer(0), ng=integer(0), min_resid=double(0), max_resid=double(0), median_resid=double(0), mean_resid=double(0), model=character(0), r2=double(0), mse=double(0), ou_grp=double(0), ou_ref=double(0), pr_grp=double(0), pr_ref=double(0), gc=double(0), gcov=double(0))
  
  # create dataset of any failures during the simulation
  fails<-data.frame(scenario=integer(0),sim=integer(0),size=integer(0),method=character(0),penalty=integer(0))
  
  # create dataset of time it takes to run each method
  times<-data.frame(method=character(0), time=double(0))

  for (j in 1:nsims) {
    
    print(paste("Sim: ", j))
    print(Sys.time())
    
    # sample from dataset
    data<-sample_n(df, simn)
    
    # add in folds 
    nfolds<-5
    data$folds<-cut(seq(1,nrow(data)),breaks=nfolds,labels=FALSE)
    folds<-data$folds
    
    # calculate results for each fold 
    for (i in 1:nfolds){

      ### step 0: split data into train/test based on the fold 
      index<-which(folds==i)
      test_i<-data[index,]
      train_i<-data[-index,]
      
      # define design matrix and y for train and test 
      y<-train_i[,yvar]
      if (includea==0) {
        X<-train_i[,!(names(train_i) %in% c(yvar, avar, avar2, 'folds'))]
      } else {X<-train_i[,!(names(train_i) %in% c(yvar, avar2, 'folds'))]}
      grp<-train_i[,avar]
      ref<-!train_i[,avar]
      n_grp<-sum(grp)
      n_ref<-sum(ref)
      
      y_test<-test_i[,yvar]
      if (includea==0) {
        X_test<-test_i[,!(names(test_i) %in% c(yvar, avar,avar2,'folds'))]
      }else {
        X_test<-test_i[,!(names(test_i) %in% c(yvar, avar2,'folds'))]
      }
      
      # Scale data before optimizing 
      y_scale<-scale(y)
      X_scale<-scale(X)
      
      # calculate avg grp costs for each fold
      grppay<-subset(train_i, eval(as.name(avar))==1)
      grp_cost<-mean(grppay[,yvar])
      refpay<-subset(train_i, eval(as.name(avar))==0)
      ref_cost<-mean(refpay[,yvar])
      
      grp_cost_scale<-(grp_cost-mean(y))/sd(y)
      ref_cost_scale<-(ref_cost-mean(y))/sd(y)
      
      # scale test data based on train scale for predictions
      # Compliments of stack exchange "Standardization/Normalization test data in R"
      X_test_scale<-sweep(sweep(X_test, 2L,attr(X_scale, 'scaled:center')), 2, attr(X_scale, 'scaled:scale'), "/")
      y_test_scale<-(y_test-mean(y))/sd(y)
      
      ### STEP 1: TRAIN ###
      ### run models on the train_i dataset & get predictions 
      
      # OLS
      model_ols<-lm(y_scale~X_scale+0)
      beta_ols = as.matrix(coef(model_ols))

      # set up preliminary for CVXR
      k<-length(beta_ols)
      beta<-Variable(k)
      loss<-sum((y_scale-X_scale %*% beta)^2)
      obj<-loss
      
      # Alt Method 1: Average constrained regression - estimated grp rev = actual grp cost
      timestart<-Sys.time()
      prob<-Problem(Minimize(obj),list((t(grp) %*% (X_scale %*% beta))/n_grp==grp_cost_scale))
      result<-solve(prob)
      beta_1<-result$getValue(beta)
      timeend<-Sys.time()
      times<-rbind(times, cbind(method='avg.constrained',time=timeend-timestart)) 
      test_betas<-cbind(beta_ols,beta_1)
      name<-c('ols','avg.constrained')

      print(paste("constrained reg:",result$status))
      
      # Alt Method 2: Weighted Average constrained regression
      alphas = c(.2, .4, .6, .8)
      ols_avg<-(t(grp) %*% (X_scale %*% cbind(beta_ols)))/n_grp
      wcosts<-(1-alphas)*grp_cost_scale+alphas*as.vector(ols_avg)
      constraint<-list((t(grp) %*% (X_scale %*% beta))/n_grp == wcosts)
      constraint<-lapply(wcosts, function(x) list((t(grp) %*% (X_scale %*% beta))/n_grp == x))
      
      c=1
      for (a in alphas){
        timestart<-Sys.time()
        prob<-Problem(Minimize(obj),constraint[[c]])
        result<-solve(prob)
        timeend<-Sys.time()
        times<-rbind(times, cbind(method=paste("weighted.avg",a),time=timeend-timestart))
        assign(paste("beta_2_",c, sep=""), result$getValue(beta))
        print(paste("weighted avg",a, result$status))
        test_betas<-cbind(test_betas, get(paste0("beta_2_",c)))
        name<-c(name, paste("weighted.avg",a))
        c=c+1 
      }
      
      # Alt Method 3: Covariance 
      cstar = cov(grp, y_scale-predict(model_ols))
      print(paste("covariance", cstar))
      share_grp<-n_grp/length(y)
      share_ref<-n_ref/length(y)
      
      c=1
      for (m in m_cov){
        timestart<-Sys.time()
        prob<-Problem(Minimize(obj),list(share_ref*(sum(t(grp) %*% (y_scale - X_scale %*% beta))) - share_grp*(sum(t(ref) %*% (y_scale - X_scale %*% beta))) < m*cstar))
        result<-my_solve(prob,j,simn)
        if (is.na(result[1])) {
          fails<-rbind(fails,cbind(scenario=scenario,sim=j,size=simn,method='cov',penalty=m))
        }
        else {
          timeend<-Sys.time()
          times<-rbind(times, cbind(method=paste("cov",m),time=timeend-timestart))
          assign(paste("beta_3_",c, sep=""), result$getValue(beta))
          print(paste("cov",m, result$status))
          test_betas<-cbind(test_betas, get(paste0("beta_3_",c)))
          name<-c(name, paste("cov",m))
        }
        c=c+1
      }
      
      # Alt Method 4: Net compensation penalty 
      c=1
      for (l in lambdas_nc){
        timestart<-Sys.time()
        prob<-Problem(Minimize(loss+l*penalty(beta,grp,X_scale,n_grp,grp_cost_scale)))
        result<-my_solve(prob,j,simn)
        if (is.na(result[1])) {
          fails<-rbind(fails,cbind(scenario=scenario,sim=j,size=simn,method='net comp',penalty=l))
        }
        else {
          timeend<-Sys.time()
          times<-rbind(times, cbind(method=paste("net.comp",l),time=timeend-timestart))
          assign(paste("beta_4_",c, sep=""), result$getValue(beta))
          print(paste("net comp",l, result$status))
          test_betas<-cbind(test_betas, get(paste0("beta_4_",c)))
          name<-c(name, paste("net.comp",l))
        }
        c=c+1
      # end of forloop
      }

      # Alt Method 4: Net compensation constraint
      c=1
      for (l in lambdas_nc_constraint){
        timestart<-Sys.time()
        prob<-Problem(Minimize(obj),list((sum(t(grp) %*% (y_scale - X_scale %*% beta)))/n_grp < l))
        result<-my_solve(prob,j,simn)
        if (is.na(result[1])) {
          fails<-rbind(fails,c(scenario=scenario,sim=j,size=simn,method='net comp constraint',penalty=l))
        }
        else {
          timeend<-Sys.time()
          times<-rbind(times, cbind(method=paste("net.comp.constraint",l),time=timeend-timestart))
          assign(paste("beta_4_",c, "b", sep=""), result$getValue(beta))
          print(paste("net comp constraint", l, result$status))
          test_betas<-cbind(test_betas, get(paste0("beta_4_",c, "b")))
          name<-c(name, paste("net.comp.constraint",l))
        } 
        c=c+1
      }
      
      # Alt Method 5: Mean residual difference 
      c=1
      # test 
      l=100
      for (l in lambdas_mrd){
        timestart<-Sys.time()
        prob<-Problem(Minimize(loss+l*penalty2(beta,grp,X_scale,n_grp,ref,n_ref,grp_cost_scale,ref_cost_scale)))
        result<-my_solve(prob,j,simn)
        if (is.na(result[1])) {
          fails<-rbind(fails,cbind(scenario=scenario,sim=j,size=simn,method='mrd',penalty=l))
        }
        else {
          timeend<-Sys.time()
          times<-rbind(times,cbind(method=paste("mrd",l),time=timeend-timestart))
          assign(paste("beta_5_",c, sep=""), result$getValue(beta))
          print(paste("mrd",l, result$status))
          test_betas<-cbind(test_betas, get(paste0("beta_5_",c)))
          name<-c(name, paste("mrd",l))
        } 
        c=c+1
      }        
      
      if (exists("alltimes")) { 
        alltimes<-rbind(alltimes,times)}
      else{alltimes<-times}
      
      # get predictions for test_i dataset
      test_preds<-apply(test_betas, 2, get_preds, X_test_scale=X_test_scale,y_scale=y_scale)
      test_preds<-as.data.frame(test_preds)
      names(test_preds)<-name
      
      # add in predictions for each row 
      tmp_data<-cbind(test_i,test_preds)
      tmp_data$y<-test_i[,yvar]
      tmp_data$grp<-test_i[,avar]
      tmp_data$grp2<-test_i[,avar2]
      
      # save data into full dataset
      if (exists("final_data")) { 
        final_data<-rbind(final_data,tmp_data)}
      else{final_data<-tmp_data}
      
      # end of fold 
    }

    for (n in name){
      pred<-final_data[,n]
      # calculate measures for each method 
      result<-all_metrics(final_data$y, pred, final_data$grp, final_data$grp2, n)
      # save results + information on the residual for the undercompensated group 
      grp<-final_data[(final_data$grp==1),]
      residual<-grp$y-grp[,n]
      simresults<-rbind(simresults, cbind(sim=j,ng=n_grp, min_resid=min(residual), max_resid=max(residual),median_resid=median(residual), mean_resid=mean(residual), result))
    }
    
    # end of iteration (of a sim)
  }
  return(simresults)
  # end of simulation function   
}
