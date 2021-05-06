ST.cluster.leroux<-function(
  x,STy,offst,E,
  wpool,Q.w, detQ,
  tripletf.pool,
  neighbour_size.pool,
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
  max.class.num <- dim(wpool)[3]
  ### create matries to store simulated values
  beta_mc<-array(dim = c(iterations,p))
  theta_mc<-array(dim=c(iterations,tp))
  phi_mc<-array(dim=c(iterations,n*tp))
  tau2_mc<-array(dim=c(iterations,tp))
  sigma2_mc<-rep(NA,iterations)
  alpha_mc<-rep(NA,iterations)
  w_mc<-rep(NA,iterations)
  quadvalue<-numeric(tp)
  prop.quadvalue<-numeric(tp)
  ### set initial values
  beta_mc[1,]<-beta.initial   
  theta_mc[1,]<-theta.initial
  phi_mc [1,]<-phi.initial 
  tau2_mc[1,]<-tau2.initial 
  sigma2_mc[1]<-sigma2.initial
  alpha_mc[1]<-alpha.initial
  w_mc[1]<-w.initial 
  
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
    cur.w.index<-w_mc[t-1]
    cur.w<-wpool[,,cur.w.index]
    cur.Q.w<-Q.w[[cur.w.index]]
    cur.det.Q.w<-detQ[cur.w.index]*tp
    
    
    #update beta 
    
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
 
    tripletf<-tripletf.pool[[cur.w.index]]
    neighbour_size<-neighbour_size.pool[[cur.w.index]]
    
    phi.proposal.sd<-sqrt(phi.proposal.v)
    ST.prop.phi<- matrix(NA, ncol = tp,nrow = n)
    for (ii in 1:n){
      for (tt in 1:tp){
        ST.prop.phi[ii,tt]<-rnorm(1, mean=ST.cur.phi[ii,tt],sd=phi.proposal.sd)
      }
    }
    
    updatephis<- updatephi2_initial(proposal_phi=ST.prop.phi,
                            ST_cur_phi=ST.cur.phi,
                            cur_beta=cur.beta,
                            cur_tau2=cur.tau2,
                            cur_theta=cur.theta,
                            rho=rho,
                            n=n,tp=tp,x=x,offst=offst,y=STy,
                            tripletf = tripletf,
                            neighbour_size = neighbour_size)
    
    ##add constraints to make sure the identifiablilty
    transphi<- apply(updatephis[[1]],2,function(colva) colva-mean(colva))
    phi_mc[t,]<-as.vector(transphi)
    cur.phi<- phi_mc[t,]
    ST.cur.phi<-matrix(cur.phi, nrow=n, ncol=tp, byrow=F)
    cntphi<-cntphi+updatephis[[2]]
    
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
    
    
    ## update w first setp within the same clustering methods
    poss.nums <- c((cur.w.index-tuning_para):(cur.w.index+tuning_para))
    if (cur.w.index %in% c(1:10)){
      max.poss.nums<-10
      min.poss.nums<-1
    }else if (cur.w.index %in% c(11:20)){
      max.poss.nums<-20
      min.poss.nums<-11
    }else if (cur.w.index %in% c(21:30)){
      max.poss.nums<-30
      min.poss.nums<-21
    } else if (cur.w.index %in% c(31:40)){
      max.poss.nums<-40
      min.poss.nums<-31
    }else if (cur.w.index %in% c(41:50)){
      max.poss.nums<-50
      min.poss.nums<-41
    }else if (cur.w.index %in% c(51:60)){
      max.poss.nums<-60
      min.poss.nums<-51
    }else if (cur.w.index %in% c(61:70)){
      max.poss.nums<-70
      min.poss.nums<-61
    }else if (cur.w.index %in% c(71:80)){
      max.poss.nums<-80
      min.poss.nums<-71
    }
    
    valid.nums <- poss.nums[poss.nums>0 & poss.nums<=max.class.num & poss.nums<=max.poss.nums & poss.nums>=min.poss.nums & poss.nums!=cur.w.index]
    choice <- sample(x=1:length(valid.nums), size=1)
    prop.w.index <- valid.nums[choice]
    
    prop.w<-wpool[,,prop.w.index]
    prop.Q.w<-Q.w[[prop.w.index]]
    prop.det.Q.w<-detQ[prop.w.index]*tp
    
    for (itau2 in 1:tp){
      quadvalue[itau2]<-quadraticform(cur.Q.w, ST.cur.phi[,itau2],n)
    }
    
    full.w<- cur.det.Q.w - 0.5*sum(quadvalue/cur.tau2)
    for (itau2 in 1:tp){
      prop.quadvalue[itau2]<-quadraticform(prop.Q.w, ST.cur.phi[,itau2],n)
    }
    full.prop.w<- prop.det.Q.w - 0.5*sum(prop.quadvalue/cur.tau2)
    
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
    
    back.valid.nums <- back.poss.nums[back.poss.nums>0 & back.poss.nums<=max.class.num & back.poss.nums<=back.max.poss.nums & back.poss.nums>=back.min.poss.nums & back.poss.nums!=cur.w.index]
    W.to.W.star <- 1/length(valid.nums)
    W.star.to.W <- 1/length(back.valid.nums)
    ratio.w <- as.numeric(exp(full.prop.w - full.w +log(W.star.to.W) - log(W.to.W.star)))
    
    
    #####
    if(runif(1,0,1) < ratio.w){ 
      cur.w.index<-prop.w.index 
      cur.w <- prop.w
      w_mc[t]<-cur.w.index
      cur.Q.w<-prop.Q.w    
      cur.det.Q.w<-prop.det.Q.w
      cntw <- cntw+1
      quadvalue<-prop.quadvalue
    }else{
      w_mc[t]<-cur.w.index
    }
    ############## update W second stetp acorss different clusteirng methods
    
    poss.nums <- unique(c(seq(cur.w.index,max.class.num,by=10),seq(cur.w.index,1,by=-10))) 
    valid.nums <- poss.nums[poss.nums>0 & poss.nums<=max.class.num & poss.nums!=cur.w.index] 
    choice <- sample(x=1:length(valid.nums), size=1)
    prop.w.index<-valid.nums[choice]
    prop.w<-wpool[,,prop.w.index]
    prop.Q.w<-Q.w[[prop.w.index]]
    prop.det.Q.w<-detQ[prop.w.index]*tp
    
    for (itau2 in 1:tp){
      quadvalue[itau2]<-quadraticform(cur.Q.w, ST.cur.phi[,itau2],n)
    }
    full.w<- cur.det.Q.w - 0.5*sum(quadvalue/cur.tau2)
    
    for (itau2 in 1:tp){
      prop.quadvalue[itau2]<-quadraticform(prop.Q.w, ST.cur.phi[,itau2],n)
    }
    
    full.prop.w<- prop.det.Q.w - 0.5*sum(prop.quadvalue/cur.tau2)
    
    ratio.w <- as.numeric(exp(full.prop.w - full.w))
    
    if(runif(1,0,1) < ratio.w){
      cur.w.index<-prop.w.index
      cur.w <- prop.w
      w_mc[t]<-cur.w.index
      cur.Q.w<-prop.Q.w
      cur.det.Q.w<-prop.det.Q.w
      cntw <- cntw+1
      quadvalue<-prop.quadvalue
    }else{
      w_mc[t]<-cur.w.index
    }
    
    # update tau2
    for (itau2 in 1:tp){
      tau2_mc[t,itau2]<-rinvgamma(1, shape=a+n/2.0, rate= b+quadvalue[itau2]/2)
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
  w.acceptance.rate<-cntw/(iterations*2)
  ## discard the burnin period
  if(burnin==0){
    
  }else{
    beta_mc<-as.data.frame(beta_mc[-c(1:burnin),])
    theta_mc<-theta_mc[-c(1:burnin),]   
    phi_mc<- phi_mc[-c(1:burnin),]   
    alpha_mc<-as.data.frame(alpha_mc[-c(1:burnin)])
    tau2_mc<-as.data.frame(tau2_mc[-c(1:burnin),])
    sigma2_mc<-as.data.frame(sigma2_mc[-c(1:burnin)])
    w_mc<-as.data.frame(w_mc[-c(1:burnin)])
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


