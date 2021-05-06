ST.cluster.leroux.Wchange<-function(
  x,STy,offst,E,
  wpool,Q.w, detQ,
  tripletf.pool.T,
  neighbour_size.pool.T,
  tuning_para=1,
  beta.initial,
  theta.initial,
  phi.initial,
  tau2.initial,
  sigma2.initial,
  alpha.initial,
  w.initial,
  beta_mu=0,
  beta_v=10,
  beta.proposal.v=0.001,
  theta.proposal.v=0.001,
  phi.proposal.v=0.001,
  a=0.001,b=0.001, rho=0.99,tp,
  iterations,burnin,thinning=1){
  
  n<-nrow(STy)
  p<-ncol(x)
  max.class.num <- 80
  ### create matries to store simulated values
  beta_mc<-array(dim = c(iterations,p))
  theta_mc<-array(dim=c(iterations,tp))
  phi_mc<-array(dim=c(iterations,n*tp))
  tau2_mc<-array(dim=c(iterations,tp))
  sigma2_mc<-rep(NA,iterations)
  alpha_mc<-rep(NA,iterations)
  w_mc<-array(dim=c(iterations,tp))  
  quadvalue<-numeric(tp)
  prop.quadvalue<-numeric(tp)
  ### set initial values
  beta_mc[1,]<-beta.initial   
  theta_mc[1,]<-theta.initial
  phi_mc [1,]<-phi.initial 
  tau2_mc[1,]<-tau2.initial 
  sigma2_mc[1]<-sigma2.initial
  alpha_mc[1]<-alpha.initial
  w_mc[1,]<-w.initial 
  
  ### count the number of being accpeted for beta, random effects and rho
  cntbeta <- 0
  cnttheta1<-0
  cnttheta<-0
  cntphi <- 0
  cntw <- 0 
  

  for (t in 2:iterations){
 
    ###set the current values
    cur.beta<-beta_mc[t-1,]
    cur.theta<-theta_mc[t-1,]
    cur.phi<-phi_mc[t-1,] 
    ST.cur.phi<-matrix(cur.phi, nrow=n, ncol=tp, byrow=F) 
    cur.tau2<-tau2_mc[t-1,]
    cur.sigma2<-sigma2_mc[t-1]
    cur.alpha<-alpha_mc[t-1]
    
    cur.w.index<-w_mc[t-1,]
    cur.w<-list(t1=wpool[[1]][ ,,cur.w.index[1]],t2=wpool[[2]][ ,,cur.w.index[2]],t3=wpool[[3]][ ,,cur.w.index[3]],
         t4=wpool[[4]][ ,,cur.w.index[4]],t5=wpool[[5]][ ,,cur.w.index[5]],t6=wpool[[6]][ ,,cur.w.index[6]],
         t7=wpool[[7]][ ,,cur.w.index[7]])
    
    
    cur.Q.w<-list(t1=Q.w[[1]][[cur.w.index[1]]],t2=Q.w[[2]][[cur.w.index[2]]],t3=Q.w[[3]][[cur.w.index[3]]],
                  t4=Q.w[[4]][[cur.w.index[4]]],t5=Q.w[[5]][[cur.w.index[5]]],t6=Q.w[[6]][[cur.w.index[6]]],
                  t7=Q.w[[7]][[cur.w.index[7]]])
    

    
    cur.det.Q.w<- list(t1=detQ[[1]][cur.w.index[1]]*1,t2=detQ[[2]][cur.w.index[2]]*1,t3=detQ[[3]][cur.w.index[3]]*1,
                       t4=detQ[[4]][cur.w.index[4]]*1,t5=detQ[[5]][cur.w.index[5]]*1,t6=detQ[[6]][cur.w.index[6]]*1,
                       t7=detQ[[7]][cur.w.index[7]]*1)
    
   
    
    #update beta 
    #proposal_beta<-mvrnorm(1,mu=cur_beta, Sigma = diag(beta.proposal.v ,nrow =length(cur_beta)))
    
    prop.beta<- rnorm(1,mean =cur.beta,sd=sqrt(beta.proposal.v))
    betastore <- updatebeta(proposal_beta=prop.beta,cur_beta=cur.beta,
                            beta_mu=beta_mu,beta_v=beta_v,
                            n=n, tp=tp,x=x,offst=offst,y=STy, ST_cur_phi=ST.cur.phi, cur_theta=cur.theta)
    
    beta_mc[t,] <- betastore[[1]]
    cntbeta <- cntbeta + betastore[[2]]
    cur.beta<- beta_mc[t,]
    
    ###### update theta_1
    
    prop.theta1<-rnorm(1,mean =cur.theta[1],sd=sqrt(theta.proposal.v))
    theta1store <- updatetheta1(proposal_theta1=prop.theta1,
                                cur_theta1=cur.theta[1],
                                cur_theta1next=cur.theta[1+1], 
                                cur_alpha=cur.alpha,
                                cur_beta=cur.beta,
                                cur_sigma2=cur.sigma2, 
                                n=n, tp=tp,x=x,offst=offst,
                                y=STy, ST_cur_phi=ST.cur.phi)
    
    theta_mc[t,1] <- theta1store[[1]]
    cur.theta[1]<-theta_mc[t,1]
    cnttheta1 <- cnttheta1 + theta1store[[2]]
    ###### update theta_t, t=2,...tp
  
    prop.theta<-c(cur.theta[1], rnorm(tp-1,mean =cur.theta[2:tp],sd=sqrt(theta.proposal.v)))
    thetastore <- updatetheta(proposal_theta=prop.theta,cur_theta=cur.theta, 
                              cur_alpha=cur.alpha,cur_beta=cur.beta,
                              cur_sigma2=cur.sigma2, n=n, tp=tp,x=x,offst=offst,
                              y=STy, ST_cur_phi=ST.cur.phi)
    
    theta_mc[t,] <- thetastore[[1]]
    theta_mc[t,]<-theta_mc[t,]-mean(theta_mc[t,])
    cnttheta <- cnttheta + thetastore[[2]]
    cur.theta<- theta_mc[t,]
    
    ####### update phi_{it} individually
   phi.proposal.sd<-sqrt(phi.proposal.v)
    ST.prop.phi<- matrix(NA, ncol = tp,nrow = n)
    for (ii in 1:n){
      for (tt in 1:tp){
        ST.prop.phi[ii,tt]<-rnorm(1, mean=ST.cur.phi[ii,tt],sd=phi.proposal.sd)
      }
    }
    
  rcppphi<- matrix(NA,ncol = tp,nrow = n)  
  for (timepoint in 1:tp){
    tripletf<-tripletf.pool.T[[timepoint]][[cur.w.index[timepoint]]]
    neighbour_size<-neighbour_size.pool.T[[timepoint]][[cur.w.index[timepoint]]]
    
      updatephis<- updatephi2(proposal_phi=ST.prop.phi,
                            ST_cur_phi=ST.cur.phi,
                            cur_beta=cur.beta,
                            cur_tau2=cur.tau2,
                            cur_theta=cur.theta,
                            rho=rho,
                            n=n,tp=tp,x=x,offst=offst,y=STy,
                            tripletf = tripletf,
                            neighbour_size = neighbour_size, t=timepoint-1)
      rcppphi[,timepoint]<-updatephis[[1]]
      cntphi<-cntphi+updatephis[[2]]
  }
  
    
    transphi<- apply(rcppphi,2,function(colva) colva-mean(colva))
    phi_mc[t,]<-as.vector(transphi)
    cur.phi<- phi_mc[t,]
    ST.cur.phi<-matrix(cur.phi, nrow=n, ncol=tp, byrow=F)
   
    #######update alpha
    updatealphastore<- updatealpha(cur_theta=cur.theta, n=n,tp=tp)
    mean.value<-updatealphastore[[1]]
    var.value <- cur.sigma2/updatealphastore[[2]] 
    alpha_mc[t]<-rtruncnorm(1,a=0,b=0.9999999999,mean=mean.value,sd=sqrt(var.value))  
    cur.alpha<-alpha_mc[t]
    
    ### update sigma2  
    theta_quad_sum<-cur.theta[1]^2
    for (thetaco in 2:tp){
      theta_quad_sum<-theta_quad_sum+(cur.theta[thetaco]-cur.alpha*cur.theta[thetaco-1])^2
    } 
    sigma2_mc[t]<- rinvgamma(1, shape=a+tp/2.0, rate= b+theta_quad_sum/2)
    cur.sigma2<-sigma2_mc[t]
    
    for (timepoint in 1:tp){
      ## update w first setp within the same clustering methods
      poss.nums <- c((cur.w.index[timepoint]-tuning_para):(cur.w.index[timepoint]+tuning_para))
      if (cur.w.index[timepoint] %in% c(1:10)){
        max.poss.nums<-10
        min.poss.nums<-1
      }else if (cur.w.index[timepoint] %in% c(11:20)){
        max.poss.nums<-20
        min.poss.nums<-11
      }else if (cur.w.index[timepoint] %in% c(21:30)){
        max.poss.nums<-30
        min.poss.nums<-21
      } else if (cur.w.index[timepoint] %in% c(31:40)){
        max.poss.nums<-40
        min.poss.nums<-31
      }else if (cur.w.index[timepoint] %in% c(41:50)){
        max.poss.nums<-50
        min.poss.nums<-41
      }else if (cur.w.index[timepoint] %in% c(51:60)){
        max.poss.nums<-60
        min.poss.nums<-51
      }else if (cur.w.index[timepoint] %in% c(61:70)){
        max.poss.nums<-70
        min.poss.nums<-61
      }else if (cur.w.index[timepoint] %in% c(71:80)){
        max.poss.nums<-80
        min.poss.nums<-71
      }
      
      valid.nums <- poss.nums[poss.nums>0 & poss.nums<=max.class.num & poss.nums<=max.poss.nums & poss.nums>=min.poss.nums & poss.nums!=cur.w.index[timepoint]]
      choice <- sample(x=1:length(valid.nums), size=1)
      prop.w.index <- valid.nums[choice]
      
      prop.w<-wpool[[timepoint]][,,prop.w.index]
      prop.Q.w<-Q.w[[timepoint]][[prop.w.index]]
      prop.det.Q.w<-detQ[[timepoint]][prop.w.index]*1
      
     
      quadvalue[timepoint]<-quadraticform(cur.Q.w[[timepoint]], ST.cur.phi[,timepoint],n)
      full.w<- cur.det.Q.w[[timepoint]] - 0.5*quadvalue[timepoint]/cur.tau2[timepoint]
      

      prop.quadvalue[timepoint]<-quadraticform(prop.Q.w, ST.cur.phi[,timepoint],n)
      
      full.prop.w<- prop.det.Q.w - 0.5*prop.quadvalue[timepoint]/cur.tau2[timepoint]
   
      ###Calculating reverse probability for Metropolis-Hastings###
      
      back.poss.nums <- c((prop.w.index-tuning_para):(prop.w.index+tuning_para))
      
      if (prop.w.index %in% c(1:10)){
        back.max.poss.nums<-10
        back.min.poss.nums<-1
      }else if(prop.w.index %in% c(11:20)){
        back.max.poss.nums<-20
        back.min.poss.nums<-11
      }else if(prop.w.index %in% c(21:30)){
        back.max.poss.nums<-30
        back.min.poss.nums<-21
      } else if(prop.w.index %in% c(31:40)){
        back.max.poss.nums<-40
        back.min.poss.nums<-31
      }else if(prop.w.index %in% c(41:50)){
        back.max.poss.nums<-50
        back.min.poss.nums<-41
      }else if(prop.w.index %in% c(51:60)){
        back.max.poss.nums<-60
        back.min.poss.nums<-51
      }else if(prop.w.index %in% c(61:70)){
        back.max.poss.nums<-70
        back.min.poss.nums<-61
      }else if(prop.w.index %in% c(71:80)){
        back.max.poss.nums<-80
        back.min.poss.nums<-71
      }
      
      back.valid.nums <- back.poss.nums[back.poss.nums>0 & back.poss.nums<=max.class.num & back.poss.nums<=back.max.poss.nums & back.poss.nums>=back.min.poss.nums & back.poss.nums!=cur.w.index[timepoint]]
      W.to.W.star <- 1/length(valid.nums)
      W.star.to.W <- 1/length(back.valid.nums)
      ratio.w <- as.numeric(exp(full.prop.w - full.w +log(W.star.to.W) - log(W.to.W.star)))
      #####
      if(runif(1,0,1) < ratio.w){ 
        cur.w.index[timepoint]<-prop.w.index 
        cur.w[[timepoint]]<-prop.w
        w_mc[t,timepoint]<-cur.w.index[timepoint]
        cur.Q.w[[timepoint]]<-prop.Q.w    
        cur.det.Q.w[[timepoint]]<-prop.det.Q.w
        cntw <- cntw+1
        quadvalue[timepoint]<-prop.quadvalue[timepoint]
      }else{
        w_mc[t, timepoint]<-cur.w.index[timepoint]
      }
      
      
      ############## update W second stetp acorss different clusteirng methods
      
      poss.nums <- unique(c(seq(cur.w.index[timepoint],max.class.num,by=10),seq(cur.w.index[timepoint],1,by=-10))) 
      valid.nums <- poss.nums[poss.nums>0 & poss.nums<=max.class.num & poss.nums!=cur.w.index[timepoint]] 
      choice <- sample(x=1:length(valid.nums), size=1)
      prop.w.index<-valid.nums[choice]
      prop.w<-wpool[[timepoint]][,,prop.w.index]
      prop.Q.w<-Q.w[[timepoint]][[prop.w.index]]
      prop.det.Q.w<-detQ[[timepoint]][prop.w.index]*1
      
      
      quadvalue[timepoint]<-quadraticform(cur.Q.w[[timepoint]], ST.cur.phi[,timepoint],n)
      full.w<- cur.det.Q.w[[timepoint]] - 0.5*quadvalue[timepoint]/cur.tau2[timepoint]
      
     
      prop.quadvalue[timepoint]<-quadraticform(prop.Q.w, ST.cur.phi[,timepoint],n)
      full.prop.w<- prop.det.Q.w - 0.5*prop.quadvalue[timepoint]/cur.tau2[timepoint]
      
      ratio.w <- as.numeric(exp(full.prop.w - full.w))
      
      if(runif(1,0,1) < ratio.w){
        cur.w.index[timepoint]<-prop.w.index
        cur.w[[timepoint]] <- prop.w
        w_mc[t,timepoint]<-cur.w.index[timepoint]
        cur.Q.w[[timepoint]]<-prop.Q.w
        cur.det.Q.w[[timepoint]]<-prop.det.Q.w
        cntw <- cntw+1
        quadvalue[timepoint]<-prop.quadvalue[timepoint]
      }else{
        w_mc[t, timepoint]<-cur.w.index[timepoint]
      }
      
}
    
    ##update tau2
    for (timepoint in 1:tp){
      tau2_mc[t,timepoint]<-rinvgamma(1, shape=a+n/2.0, rate= b+quadvalue[timepoint]/2)
    }
    cur.tau2<-tau2_mc[t,]
    ### checking acceptance rate every 100 steps to see whether lie in 30-60 %, adjust themselves automatically
    if (floor(t/100)==t/100){   
      step.phi.acceptance.rate<-cntphi/(t*n*tp)
      if(step.phi.acceptance.rate<0.3){
        phi.proposal.v<-phi.proposal.v*(1/2)
      }else if
      (step.phi.acceptance.rate>0.6){
        phi.proposal.v<-phi.proposal.v*2
      } 
      step.beta.acceptance.rate<-cntbeta/t
      if(step.beta.acceptance.rate<0.3){
        beta.proposal.v<-beta.proposal.v*(1/2)
      }else if (step.beta.acceptance.rate>0.6){
        beta.proposal.v<-beta.proposal.v*2
      } 
      
      step.theta.acceptance.rate<-(cnttheta+cnttheta1)/(t*tp)
      if(step.theta.acceptance.rate<0.3){
        theta.proposal.v<-theta.proposal.v*(1/2)
      }else if (step.theta.acceptance.rate>0.6){
        theta.proposal.v<-theta.proposal.v*2
      } 
    } 
  } 

  ## store acceptance rate
  beta.acceptance.rate<-cntbeta/iterations
  theta.acceptance.rate <-(cnttheta+cnttheta1)/(iterations*tp)
  phi.acceptance.rate <-cntphi/(iterations*n*tp)
  w.acceptance.rate<-cntw/(iterations*2*tp)
  ## discard the burnin period
  if(burnin==0){
    
  }else{
    beta_mc<-as.data.frame(beta_mc[-c(1:burnin),])
    theta_mc<-theta_mc[-c(1:burnin),]   
    phi_mc<- phi_mc[-c(1:burnin),]   
    alpha_mc<-as.data.frame(alpha_mc[-c(1:burnin)])
    tau2_mc<-as.data.frame(tau2_mc[-c(1:burnin),])
    sigma2_mc<-as.data.frame(sigma2_mc[-c(1:burnin)])
    w_mc<-as.data.frame(w_mc[-c(1:burnin),])
  }
  
  ## thinning
  if (thinning==0){
    beta_mc<-as.matrix(beta_mc)
  }else if (thinning!=0){
    ## thinning
    beta_mc<-as.matrix(beta_mc[seq(1,nrow(beta_mc) ,thinning),])
    theta_mc<- theta_mc[seq(1,nrow(theta_mc),thinning),]
    phi_mc<- phi_mc[seq(1,nrow(phi_mc),thinning),]
    alpha_mc<-alpha_mc[seq(1,nrow(alpha_mc),thinning),]
    tau2_mc<-tau2_mc[seq(1,nrow(tau2_mc),thinning),]
    sigma2_mc<-sigma2_mc[seq(1,nrow(sigma2_mc),thinning),]
    w_mc<-w_mc[seq(1,nrow(w_mc),thinning),]
  }
  
  ## save the results
  #### return a list
  output<-list(beta.chain=as.data.frame(beta_mc),
               theta.chain=as.data.frame(theta_mc),
               phi.chain=as.data.frame(phi_mc),
               alpha.chain=as.data.frame(alpha_mc),
               tau2.chain=as.data.frame(tau2_mc), 
               sigma2.chain=as.data.frame(sigma2_mc), 
               w.chain=as.data.frame(w_mc),
               tuning_para=tuning_para,
               beta.acceptance.rate=beta.acceptance.rate, 
               theta.acceptance.rate=theta.acceptance.rate,
               phi.acceptance.rate=phi.acceptance.rate,
               w.acceptance.rate=w.acceptance.rate,
               w.initial=w.initial
  )
  return(output)
}


###mode function
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#Compute the DIC and p.d
common.modelfit <- function(samples.loglike, deviance.fitted)
{
  #### DIC
  mean.deviance <- -2 * sum(samples.loglike, na.rm=TRUE) /  nrow(samples.loglike)
  p.d <- mean.deviance - deviance.fitted
  DIC <- deviance.fitted + 2 * p.d
  
  
  #### Model fit criteria
  modelfit <- c(DIC, p.d)
  names(modelfit) <- c("DIC", "p.d")
  return(modelfit)  
}

