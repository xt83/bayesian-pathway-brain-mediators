library(mvtnorm)
library(invgamma)
library(SuppDists)
library(tensorr)
library(Matrix)
library(truncnorm)

## Input ADNI data

# Brain structural connectivity
structural_connectivity<-read.csv("matrix.csv")
structural_connectivity<-structural_connectivity[,-1]
structural_connectivity<-as.matrix(structural_connectivity)

# Demographic data
data<-read.csv("demo.csv")

# Covariates
covariates<-data[,8:10]
covariates<-as.matrix(covariates)
# Including intercept
x<-cbind(1,x1)

# Time to AD onset
outcome<-data$outcome
log_outcome<-log(outcome)

# APOE status
APOE4<-data$APOE4

# Censoring Status
# c1=1: uncensored
# c1=0: censored
censor_status<-data1$event

## MCMC Update Process


# Update coefficient for covariates in AFT model
update_betax<-function(xtrain,ztrain,ytrain,sigma0,sigmax,omega,eta,a1h,a2h,betaz,alpha_k,beta_k,h,k1,k2,p,Atrain){
  n1<-length(xtrain[1,])
  n2<-length(xtrain[,1])
  y_hat<-rep(0,n2)
  Sigma<-solve((1/sigma0)*t(xtrain)%*%xtrain+(1/sigmax)*diag(n1))
  for(i in 1:n2){
    b1<-p*(i-1)+1
    b2<-b1+p-1
    A_1<-Atrain[b1:b2,]
    y_t=ytrain[i]
    for(j in 1:k2){
      y_t=y_t-as.vector(omega[j]*t(beta_k[j,])%*%A_1%*%beta_k[j,])
    }
    y_hat[i]=y_t-as.vector(betaz*ztrain[i])
  }
  u<-(1/sigma0)*Sigma%*%t(xtrain)%*%y_hat
  betax<-rmvnorm(1,mean = u,sigma = Sigma)
  return(betax)
}
# Update coefficient for exposure in AFT model
update_betaz<-function(xtrain,ztrain,ytrain,sigma0,sigmax,omega,eta,a1h,a2h,betax,alpha_k,beta_k,h,k1,k2,p,sigmaz,Atrain){
  n2<-length(ztrain)
  y_hat<-rep(0,n2)
  Sigma<-((sum(ztrain^2)/sigma0)+(1/sigmaz))^(-1)
  for(i in 1:n2){
    b1<-p*(i-1)+1
    b2<-b1+p-1
    A_1<-Atrain[b1:b2,1:p]
    y_t=ytrain[i]
    for(j in 1:k2){
      y_t=y_t-as.vector(omega[j]*t(beta_k[j,])%*%A_1%*%beta_k[j,])
    }
    y_hat[i]=y_t-sum(betax*xtrain[i,])
  }
  u<-sum(y_hat*ztrain)*Sigma/sigma0
  betaz<-rnorm(1,mean = u,sd=sqrt(Sigma))
  return(betaz)
}

# Update estimated structural connectivity matrix
A_matrix<-function(a1h,a2h,eta,alpha_k,xi,p,h,k1,zi){
  A_1<-matrix(0,nrow=p,ncol=p)
  for(j1 in 1:h){
    A_1=A_1+sum(a2h[j1,]*xi)*a1h[j1,]%*%t(a1h[j1,])
  }
  for(j2 in 1:k1){
    A_1=A_1+eta[j2]*zi*alpha_k[j2,]%*%t(alpha_k[j2,])
  }
  diag(A_1)<-0
  return(A_1)
}

# Update a1h and a2h for the coefficient tensor to adjust for the effects from the covariates.
update_a1h<-function(Atrain,a1h,a2h,xtrain,ztrain,eta,h,p,k1,alpha_k,sigmae,sigma_a){
  n2<-length(xtrain[,1])
  for(i in 1:h){
    for(j in 1:p){
      q1=0
      q2=0
      a0<-a1h[i,-j]
      for(q in 1:n2){
        q1=q1+as.vector(t(a0)%*%a0)*(as.vector(t(a2h[i,])%*%xtrain[q,])^2)
         b1<-p*(q-1)+1
         b2<-b1+p-1
        A_1<-Atrain[b1:b2,1:p]-A_matrix(xi=xtrain[q,],zi=ztrain[q],a1h,a2h,eta,alpha_k,p,h,k1)+sum(a2h[i,]*xtrain[q,])*a1h[i,]%*%t(a1h[i,])
        diag(A_1)<-0
        A_bar<-A_1[j,-j]
        q2<-q2+as.vector(t(a2h[i,])%*%xtrain[q,])*as.vector(t(A_bar)%*%a0)
      }
      sigma_1<-((q1/sigmae)+(1/sigma_a))^(-1)
      u<-sigma_1*q2/sigmae
      a1h[i,j]<-rnorm(1,u,sqrt(sigma_1))
    }
  }
 return(a1h) 
}
update_a2h<-function(Atrain,a1h,a2h,xtrain,ztrain,eta,h,p,k1,alpha_k,sigmae,sigma_a){
  n1<-length(xtrain[1,])
  n2<-length(xtrain[,1])
  for(i in 1:h){
    A_t<-a1h[i,]%*%t(a1h[i,])
    #a1: upper diagonal element of a1h%O%a1h
    a1<-t(A_t)[lower.tri(t(A_t))]
    q1=matrix(0,nrow=n1,ncol=n1)
    q2=rep(0,n1)
    for(q in 1:n2){
      q1=q1+as.vector(t(a1)%*%a1)*xtrain[q,]%*%t(xtrain[q,])
      b1<-p*(q-1)+1
      b2<-b1+p-1
      A_1<-Atrain[b1:b2,1:p]-A_matrix(xi=xtrain[q,],zi=ztrain[q],a1h,a2h,eta,alpha_k,p,h,k1)+sum(a2h[i,]*xtrain[q,])*a1h[i,]%*%t(a1h[i,])
      #a2: upper diagonal element of Ai_bar
      a2<-t(A_1)[lower.tri(t(A_1))]
      q2=q2+as.vector(t(a1)%*%a2)*xtrain[q,]
    }
    sigma_1<-solve((1/sigmae)*q1+(1/sigma_a)*diag(n1))
    u<-(1/sigmae)*sigma_1%*%q2
    a2h[i,]<-rmvnorm(1,mean = u,sigma = sigma_1)
  }
  return(a2h)
}

# Update within subgraph effect of each node contribute to outcome

update_beta_k<-function(xtrain,ztrain,Atrain,ytrain,sigma0,sigmax,omega,betax,alpha_k,h,k1,k2,p,sigmaz,v,beta_k,betaz,v0,v1){
  n2<-length(ztrain)
  inval<-1:k2
  for(i in 1:k2){
    for(j in 1:p){
      q1=0
      q2=0
      y_bar=0
      for(q in 1:n2){
        b1<-p*(q-1)+1
        b2<-b1+p-1
        A_1<-Atrain[b1:b2,1:p]
        q1=q1+(2*sum(omega[i]*beta_k[i,-j]*A_1[j,-j]))^2
        y_bar=ytrain[q]-sum(betax*xtrain[q,])-betaz*ztrain[q]-as.vector(omega[i]*t(beta_k[i,-j])%*%A_1[-j,-j]%*%beta_k[i,-j])
        for(s in inval[-i]){
          y_bar=y_bar-as.vector(omega[s]*t(beta_k[s,])%*%A_1%*%beta_k[s,])
        }
        q2=q2+2*sum(omega[i]*beta_k[i,-j]*A_1[j,-j])*y_bar
      }
      sigma_beta_kp<-((1-v[i,j])*v0+v[i,j]*v1)^2
      sigma_1<-((q1/sigma0)+(1/sigma_beta_kp))^(-1)
      u<-q2*sigma_1/sigma0
      beta_k[i,j]<-rnorm(1,mean = u,sd=sqrt(sigma_1))
    }
  }
  return(beta_k)
}


# Update impact between subgraphs contribute to outcome: omega
update_omega<-function(xtrain,ztrain,Atrain,ytrain,sigma0,sigmax,omega,betax,alpha_k,h,k2,p,sigmaz,v,beta_k,betaz,tau_omega){
   n2<-length(ztrain)
   w<-matrix(0,nrow=n2,ncol = k2)
   for(i in 1:k2){
     for(q in 1:n2){
       b1<-p*(q-1)+1
       b2<-b1+p-1
       A_1<-Atrain[b1:b2,1:p]
       w[q,i]<-as.vector(t(beta_k[i,])%*%A_1%*%beta_k[i,])
     }
   }
   y_bar=ytrain-xtrain%*%t(betax)-betaz*ztrain
   sigma_1=solve((1/sigma0)*(t(w)%*%w)+(1/tau_omega)*diag(k2))
   u<-(1/sigma0)*t(y_bar)%*%w%*%sigma_1
   omega<-rmvnorm(1,mean = u,sigma = sigma_1)
   return(omega)
}

# Update impact between subgraphs from exposure to mediator: eta

update_eta<-function(xtrain,ztrain,Atrain,ytrain,sigma0,sigmax,omega,betax,alpha_k,h,k1,p,sigmaz,v,beta_k,betaz,sigmae,a1h,a2h,tau_eta){
  Q<-matrix(0,nrow=k1,ncol=k1)
  U<-matrix(0,nrow=1,ncol=k1)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      Q1<-form_Q(i,j,k1,ztrain,alpha_k)
      A1<-form_A(Atrain,xtrain,a1h,a2h,h,ztrain,i,j,p)
      Q<-Q+(t(Q1)%*%Q1)
      U<-U+t(A1)%*%Q1
    }
  }
  sigma_1<-solve((1/sigmae)*Q+(1/tau_eta)*diag(k1))
  u<-(1/sigmae)*U%*%sigma_1
  eta<-rmvnorm(1,mean = u,sigma = sigma_1)
}
form_Q<-function(i,j,k1,ztrain,alpha_k){
  n2<-length(ztrain)
  Q1<-matrix(0,nrow=n2,ncol = k1)
  for(w in 1:n2){
    for(c in 1:k1){
      Q1[w,c]<-ztrain[w]*alpha_k[c,i]*alpha_k[c,j]
    }
  }
  return(Q1)
}
form_A<-function(Atrain,xtrain,a1h,a2h,h,ztrain,i,j,p){
  n2<-length(ztrain)
  A0<-rep(0,n2)
  for(q in 1:n2){
    b1<-p*(q-1)+1
    b2<-b1+p-1
    A_1<-Atrain[b1:b2,1:p]
    for(w in 1:h){
      A_1=A_1-sum(a2h[w,]*xtrain[q,])*a1h[w,]%*%t(a1h[w,])
    }
    A0[q]<-A_1[i,j]
  }
  return(A0)
}
# Update variance for eta and omega
update_tau_eta<-function(v_eta,eta,k1){
  tau<-rep(0,k1)
  for(i in 1:k1){
    t<-rinvGauss(1,v_eta/abs(eta[i]),(v_eta)^2)
    tau[i]<-1/t
  }
  return(tau)
}
update_tau_omega<-function(v_omega,omega,k2){
  tau<-rep(0,k2)
  for(i in 1:k2){
    t<-rinvGauss(1,v_omega/abs(omega[i]),(v_omega)^2)
    tau[i]<-1/t
  }
  return(tau)
}

# update the latent selection indicators introduced to identify nonzero elements within alpha_k, beta_j
update_v<-function(v0,v1,mu,epi,k2,p,beta_k,t){
  v<-matrix(0,nrow=k2,ncol=p)
  for(i in 1:k2){
    for(j in 1:p){
      p0<-(1/v0)*exp(-((beta_k[i,j]^2)/(2*v0^2)))
      p1<-(1/v1)*exp(-((beta_k[i,j]^2)/(2*v1^2))-mu)
      if((p0==0)&(p1==0)){
        v[i,j]=sample(c(0,1),1,replace=TRUE,prob = c(0.5,0.5))
      }
      else{
        v[i,j]=sample(c(0,1),1,replace=TRUE,prob = c(p0,p1))
      }
    }
  }
  return(v)
}

update_t<-function(v0,v2,mu,epi,k1,p,alpha_k,v){
  t<-matrix(0,nrow=k1,ncol=p)
  for(i in 1:k1){
    for(j in 1:p){
      p0<-(1/v0)*exp(-((alpha_k[i,j]^2)/(2*v0^2)))
      p1<-(1/v2)*exp(-((alpha_k[i,j]^2)/(2*v2^2))-mu)
      if((p0==0)&(p1==0)){
        t[i,j]=sample(c(0,1),1,replace=TRUE,prob = c(0.5,0.5))
      }
      else{
        t[i,j]=sample(c(0,1),1,replace=TRUE,prob = c(p0,p1))
      }  
    }
  }
  return(t)
}
# Update structural matrix
update_A<-function(a1h,a2h,eta,alpha_k,p,h,k1,xtrain,ztrain){
  n2<-length(ztrain)
  A0<-matrix(0,nrow=p*n2,ncol=p)
  for(i in 1:n2){
    b1<-p*(i-1)+1
    b2<-b1+p-1
    A1<-A_matrix(a1h,a2h,eta,alpha_k,xi=xtrain[i,],p,h,k1,zi=ztrain[i])
    A0[b1:b2,1:p]<-A1
  }
  return(A0)
}

# Update within subgraph effect of each node from exposure to mediator
update_alpha_k<-function(xtrain,ztrain,Atrain,ytrain,sigma0,sigmax,omega,betax,alpha_k,h,k1,p,sigmaz,t,beta_k,betaz,v0,v2,eta,sigmae,a1h,a2h){
  n2<-length(ztrain)
  for(i in 1:k1){
    for(j in 1:p){
      q1=0
      q2=0
      a0<-alpha_k[i,-j]
      for(q in 1:n2){
        q1=q1+(ztrain[q]^2)*(eta[i]^2)*sum(a0*a0)
        b1<-p*(q-1)+1
        b2<-b1+p-1
        A_1<-Atrain[b1:b2,1:p]
        for(j1 in 1:h){
        A_1=A_1-sum(a2h[j1,]*xtrain[q,])*a1h[j1,]%*%t(a1h[j1,])
        }
        if(j>1&j<p){
        a1a<-A_1[1:(j-1),j]
        a1b<-A_1[j,(j+1):p]
        a1<-c(a1a,a1b)
        }
        else if(j==1){
          a1<-A_1[1,2:p]
        }
        else if(j==p){
          a1<-A_1[1:(p-1),p]
        }
        q2=q2+ztrain[q]*eta[i]*sum(a0*a1)
      }
      sigma_alpha_kp<-((1-t[i,j])*v0+t[i,j]*v2)^2
      sigma_1<-((q1/sigmae)+(1/sigma_alpha_kp))^(-1)
      u<-q2*sigma_1/sigmae
      alpha_k[i,j]<-rnorm(1,mean = u,sd=sqrt(sigma_1))
    }
  }
  return(alpha_k)
}
# update error variance for matrix response regression
update_sigma0<-function(ytrain,y1){
  s=sum((ytrain-y1)^2)
  n3<-length(y1)
  sigma0=rinvgamma(1,(3+n3)/2,(s+1)/2)
}
# update error variance for AFT model
update_sigmae<-function(An,Atrain,p,ztrain){
  n2<-length(ztrain)
  n3<-n2*p*p-n2*p
  n3<-n3/2
  s=sum((An-Atrain)^2)/2
  sigmae=rinvgamma(1,(1+n3)/2,(s+3)/2)
}
# Update estimated log survival time
update_y<-function(An,betax,xtrain,ztrain,beta_k,omega,k2,p,betaz){
  n2<-length(ztrain)
  y1<-rep(0,n2)
  for(i in 1:n2){
    b1<-p*(i-1)+1
    b2<-b1+p-1
    A1<-An[b1:b2,1:p]
    w=0
    for(j in 1:k2){
      w=w+as.vector(omega[j]*t(beta_k[j,])%*%A1%*%beta_k[j,])
    }
    y1[i]<-w+betaz*ztrain[i]+sum(betax*xtrain[i,])
  }
  return(y1)
}
## MCMC function

mcmc<-function(xtrain,ztrain,Atrain,ytrain,censor,sigmax=0.4,sigma_a=0.5,v0=sqrt(0.1),v1=2,v2=2,numite=300,burn_in=200,k1=3,k2=3,h=2,p=83,epi=0.2,mu=0.2,nonparametric=T,ind=T,int=0){
  y_old<-ytrain
  A_new<-list()
  B_new<-list()
  v_new<-list()
  t_new<-list()
  a1h_new<-list()
  a2h_new<-list()
  betax_new<-list()
  betaz_new<-list()
  sigma0_new<-list()
  sigmae_new<-list()
  # Initialization Process
  censor_1=which(censor==0)
  eta_train<-matrix(0,nrow=numite,ncol=k1)
  lambda_train<-matrix(0,nrow=numite,ncol=k2)
  n1<-length(xtrain[1,])
  pre_train<-matrix(0,nrow=numite,ncol=length(ztrain))
  sigma0<-0.4
  sigmaz<-0.4
  sigmae<-0.4
  betaz<--1
  v_lambda<-0.4
  v_eta<-0.7
  tau_lambda<-rep(0.5,k2)
  tau_eta<-rep(0.5,k1)
  lambda<-rmvnorm(1,rep(0,k2),tau_lambda*diag(k2))
  eta<-rmvnorm(1,rep(0,k1),tau_eta*diag(k1))
  a1h<-matrix(0,nrow=h,ncol=p)
  a2h<-matrix(0,nrow=h,ncol=n1)
  for(i in 1:h){
    a1h[i,]<-rmvnorm(1,rep(0,p),sigma_a*diag(p))
    a2h[i,]<-rmvnorm(1,rep(0,n1),sigma_a*diag(n1))
  }
  beta_k<-matrix(0,nrow=k2,ncol=p)
  alpha_k<-matrix(0,nrow=k1,ncol=p)
  v<-matrix(0,nrow=k2,ncol=p)
  t<-matrix(0,nrow=k1,ncol=p)
  for(i in 1:k1){
    for(j in 1:p){
      t[i,j]<-sample(c(0,1),1,replace=TRUE,prob = c(0.5,0.5))
      alpha_k[i,j]<-rnorm(1,0,0.8)#sigma_alpha_kp)
    }
  }
  for(i in 1:k2){
    for(j in 1:p){
      v[i,j]<-sample(c(0,1),1,replace=TRUE,prob = c(0.5,0.5))
      beta_k[i,j]<-rnorm(1,0,0.2)#sigma_alpha_kp)
    }
  }
  # Update process via Gibbs Sampler
  for(i in 1:numite){
    betax<-update_betax(xtrain,ztrain,ytrain,sigma0,sigmax,omega,eta,a1h,a2h,betaz,alpha_k,beta_k,h,k1,k2,p,Atrain)
    a1h<-update_a1h(Atrain,a1h,a2h,xtrain,ztrain,eta,h,p,k1,alpha_k,sigmae,sigma_a)
    a2h<-update_a2h(Atrain,a1h,a2h,xtrain,ztrain,eta,h,p,k1,alpha_k,sigmae,sigma_a)
    betaz<-update_betaz(xtrain,ztrain,ytrain,sigma0,sigmax,omega,eta,a1h,a2h,betax,alpha_k,beta_k,h,k1,k2,p,sigmaz,Atrain)
    beta_k<-update_beta_k(xtrain,ztrain,Atrain,ytrain,sigma0,sigmax,omega,betax,alpha_k,h,k1,k2,p,sigmaz,v,beta_k,betaz,v0,v1)
    v<-update_v(v0,v1,mu,epi,k2,p,beta_k,t)
    beta_k[nv]<-0
    alpha_k<-update_alpha_k(xtrain,ztrain,Atrain,ytrain,sigma0,sigmax,omega,betax,alpha_k,h,k1,p,sigmaz,t,beta_k,betaz,v0,v2,eta,sigmae,a1h,a2h)
    t<-update_t(v0,v2,mu,epi,k1,p,alpha_k,v)
    alpha_k[nt]<-0
    omega<-update_omega(xtrain,ztrain,Atrain,ytrain,sigma0,sigmax,omega,betax,alpha_k,h,k2,p,sigmaz,v,beta_k,betaz,tau_omega)
    tau_omega<-update_tau_omega(v_omega,omega,k2)
    eta<-update_eta(xtrain,ztrain,Atrain,ytrain,sigma0,sigmax,omega,betax,alpha_k,h,k1,p,sigmaz,v,beta_k,betaz,sigmae,a1h,a2h,tau_eta)
    tau_eta<-update_tau_eta(v_eta,eta,k1)
    sigmaz<-rinvgamma(1,4/2,(1+betaz^2)/2)
    An<-update_A(a1h,a2h,eta,alpha_k,p,h,k1,xtrain,ztrain)
    sigmae<-update_sigmae(An,Atrain,p,ztrain)
    An<-Atrain
    y1<-update_y(An,betax,xtrain,ztrain,beta_k,omega,k2,p,betaz)
    sigma0<-update_sigma0(ytrain,y1)
    for(i5 in censor_1){
    ytrain[i5]<-rtruncnorm(1, a=y_old[i5], b=Inf, mean = y1[i5], sd =sqrt(sigma0))
    }
    # Save Bayesian inference result in each iteration
    pre_train[i,]<-y1
    eta_train[i,]<-eta
    omega_train[i,]<-omega
    A_new[[i]]<-alpha_k
    B_new[[i]]<-beta_k
    v_new[[i]]<-v
    t_new[[i]]<-t
    a1h_new[[i]]<-a1h
    a2h_new[[i]]<-a2h
    betax_new[[i]]<-betax
    betaz_new[[i]]<-betaz
    sigma0_new[[i]]<-sigma0
    sigmae_new[[i]]<-sigmae
  }
  return(list(pre_train,A_new,B_new,v_new,t_new,eta_train,omega_train,a1h_new,a2h_new,betax_new,betaz_new,sigma0_new,sigmae_new))
}
result<-mcmc(xtrain=x,ytrain=log_outcome,ztrain=APOE4,Atrain=structural_connectivity,censor=censor_status,numite=10000,k1=a,k2=b)

