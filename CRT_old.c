#include <R.h>
#include <math.h>
#include <Rmath.h>

void urnings_adaptive_CRT_paired_update(int*u,int*v,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,double*cumsum,int*x,int*y,int*k1,int*k2,double*P,int*queue,int*LL,int*LLsum,int*Score,int*n_options,int*Upd,int*n_scores,double*Q){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 double Correct=1.00;
 double Incorrect=1.00;
 int A=1;
 int j=0;
 int oldU=0;
 int oldV=0;
 int newU=0;
 int newV=0;
 int jj=0;
 int s=0;
 int success=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      /* Compute the normalising constant for the selection probabilities and a vector used for sampling from a multinomial distribution*/
      Mp=0;		
      for(int s=0;s<K[0];s++){
        Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
        cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
      }
      /* sample an item given the selection probabilities*/
      p=runif(0,Mp);
      j=0;
      for(int s=1;s<K[0];s++){
        if(p>cumsum[s]){
          j=j+1;
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
	    /* generate the observed accuracy*/
      A=4;        /*weight for the accuracy*/
	    L=1/(1+exp((delta[j+rep*K[0]]-theta[i+rep*N[0]])*A)); /* true probability correct*/
      p=runif(0,1);
      x[2]=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x[2]*A;                             
      v[j]=v[j]+(1-x[2])*A;	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<A;s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+A-v[j]-s);
        Incorrect=Incorrect*(k1[0]+A-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);
      /* generate the simulated response */
      p=runif(0,1);
      y[2]=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y[2]*A;
      v[j]=v[j]-(1-y[2])*A;       
       
	    /* if response is correct, also do updates based on response times*/
      if(x[2]==1){
	      for(int S=1;S>(-1);S--){/* loop over pseudo-responses with weights 2,1*/
          /* specify the correct weight*/
		      if(S==1){A=2;}
		      if(S==0){A=1;}
          /* generate the observed pseudo-response*/
      	  L=1/(1+exp((delta[j+rep*K[0]]-theta[i+rep*N[0]])*A));
      	  p=runif(0,1);
      	  x[S]=1*(L>p);
          /* add the balls based on the observed pseudo-response to the urns*/                      
		      u[i]=u[i]+x[S]*A;                             
      	  v[j]=v[j]+(1-x[S])*A;	
          /* compute the probability of X=1 given the urn configurations*/
      	  Correct=1.00;
      	  Incorrect=1.00;   
      	  for(int s=0;s<A;s++){              
         		Correct=Correct*(u[i]-s)*(k2[0]+A-v[j]-s);
         		Incorrect=Incorrect*(k1[0]+A-u[i]-s)*(v[j]-s);
          }
          L=Correct/(Correct+Incorrect);
          /* generate the simulated pseudo-response */
          p=runif(0,1);
          y[S]=1*(L>p);
          /* remove the balls based on the simulated pseudo-response from the urns: These would be the proposed values*/
          u[i]=u[i]-y[S]*A;
      	  v[j]=v[j]-(1-y[S])*A;       
        }
      }

      /* These are the proposed values*/
      newV=v[j];
	    newU=u[i];
      
      /* go back to the old values and then decide on whether the proposed values (if different from old) should be accepted*/
      u[i]=oldU;
      v[j]=oldV;
      if(newU!=oldU){
        /*Compute the normalising constant for the selection probability given the proposed values*/
		    Mp1=P[newU+newV*(k1[0]+1)];
      	for(int s=0;s<K[0];s++){
          if(s!=j){
          	Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
          }
      	}
        /* compute the MH acceptance probability*/
		    L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
        /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
        p=runif(0,1);
		    if(p<L){
			    u[i]=newU;
			    v[j]=newV;
		    }
	    }
	    /*set item j to the old value*/
       newV=v[j];
	    v[j]=oldV;
      /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	    if(newV!=oldV){
			  /*Find which of the possible update values (saved in vector Score) has been proposed*/
			  for(int t=0;t<n_scores[0];t++){
				  if((newV-oldV)==Score[t]){
					  s=t;
				  }
			  }
        /*for every value of the paired update count how many items are in the queue and can be updated*/
	      for(int t=0;t<n_scores[0];t++){
				  LLsum[t]=0;
				  for(int ii=0;ii<K[0];ii++){
					  LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					  LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				  }
			  }
        success=0;  
			  for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				  if(success==0){/*this is done only if the paired update has not been found yet*/
					  /*check whether for the particular possible update there are enough items in the queue*/
            Q[d]=1;
					  for(int t=0;t<n_scores[0];t++){
						  Q[d]=Q[d]*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					  }
					  if(Q[d]>0){/*If there are enough items*/
						  success=1;/*set succes to 1, such that we will not look further*/
						  for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							  if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								  for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									  /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                    for(int k=0;k<K[0];k++){
										  cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									  }
                    /*select an item randomly among those that are in the queue (where LL is not 0)*/
									  p=runif(0,LLsum[t]);
									  jj=0;
									  for(int k=1;k<K[0];k++){
										  if(p>cumsum[k]){
											  jj=jj+1;
										  }
									  }	
									  v[jj]=v[jj]+Score[t];/*update the selected item*/
									  queue[jj]=0;/*remove the selected item from the queue*/
									  LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									  LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								  }
							  }
						  }
						  v[j]=v[j]+Score[s];	/*update item j*/
					  }
				  }
			  }
			  if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				  queue[j]=Score[s];
			  }  
      }
        
         
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
    }
  }
  PutRNGstate();
}

void urnings_simpleX_HT(int*adaptive,int*paired,int*u,int*v,double*theta,double*delta,int*N,int*K,int*Rep,int*k1,int*k2,double*P,double*cumsum,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,double*MSE,double*mse_baseline,int*HT){
  double L=0;
  double p=0;
  double Mp=0;
  double Mp1=0;
  int oldU=0;
  int oldV=0;
  int newV=0;
  int newU=0;
  int j=0;
  int x=0;
  int y=0;
  int success=0;
  int Q=0;
  int s=0;
  int jj=0;
  double Correct=0;
  double Incorrect=0;
  double dif=0;
  GetRNGstate();
  for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    if(MSE[0]>mse_baseline[0]){/* only if MSE is still about the baseline, do the iteration*/
      HT[0]=rep;/*the current iteration is saved as potential hitting time, if after this iteration MSE goes below the baseline then this would be the recorded HT.*/
      for(int i=0;i<N[0];i++){/* loop over persons*/
    if(adaptive[0]==0){/* sample an item with equal probabilities*/
    p=runif(0,1.00*K[0]);
      j=0;
      for(int s=1;s<K[0];s++){
        if(p>s){
          j=j+1;
        }
      }
    }
    if(adaptive[0]==1){      /*adaptive item selection, kernels are pre-specified in matrix P*/
    Mp=0;		
      for(int s=0;s<K[0];s++){
        Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
        cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
      }
      p=runif(0,Mp);
      j=0;
      for(int s=1;s<K[0];s++){
        if(p>cumsum[s]){
          j=j+1;
        }
      }
    }
    /* save the current values of the person and the item */
    oldU=u[i];
    oldV=v[j];
    /* generate the observed accuracy*/
    L=1/(1+exp((delta[j]-theta[i]))); /* true probability correct, W is discrimination parameter, 1 is default value*/
    p=runif(0,1.00);
    x=1*(L>p);
    /* add the balls based on the observed response to the urns*/                      
    u[i]=u[i]+x;                             
    v[j]=v[j]+(1-x);	
    /* compute the probability of X=1 given the urn configurations*/
    Correct=u[i]*(k2[0]-v[j]);
    Incorrect=v[j]*(k1[0]-u[i]);   
    L=Correct/(Correct+Incorrect);      
    /* generate the simulated response */
    p=runif(0,1.00);
    y=1*(L>p);
    /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
    u[i]=u[i]-y;
    v[j]=v[j]-(1-y);
    
    if(adaptive[0]==1){/*Metropolis step for when item selection is adaptive*/
    newU=u[i];
      newV=v[j];
      u[i]=oldU;
      v[j]=oldV;
      if(newU!=oldU){
        /*Compute the normalising constant for the selection probability given the proposed values*/
        Mp1=P[newU+newV*(k1[0]+1)];
        for(int s=0;s<K[0];s++){
          if(s!=j){
            Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
          }
        }
        /* compute the MH acceptance probability*/
        L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
        /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
        p=runif(0,1);
        if(p<L){
          u[i]=newU;
          v[j]=newV;
        }
      }
    }
    if(paired[0]==1){/*paired update for when it is on*/
        /*set item j to the old value*/
        newV=v[j];
      v[j]=oldV;
      /*If the proposed value is different from the old one, check if a paired update can be performed*/ 
      if(newV!=oldV){
        /*Find which of the possible update values (saved in vector Score) has been proposed*/
        for(int t=0;t<n_scores[0];t++){
          if((newV-oldV)==Score[t]){
            s=t;
          }
        }
        /*for every possible value of an update, count how many items are in the queue and can be updated*/
        for(int t=0;t<n_scores[0];t++){
          LLsum[t]=0;
          for(int ii=0;ii<K[0];ii++){
            LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
            LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
          }
        }
        success=0;  
        for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update (in the rasch model case there is only one possible option for a paired update (-1 if item j gets an update +1, and +1 ifitem j gets an update of -1))*/
        if(success==0){/*this is done only if the paired update has not been found yet*/
        /*check whether for the particular possible update there are enough items in the queue*/
        Q=1;
          for(int t=0;t<n_scores[0];t++){
            Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
          }
          if(Q>0){/*If there are enough items*/
        success=1;/*set succes to 1, such that we will not look further*/
        for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
        if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
        for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
          /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
          for(int k=0;k<K[0];k++){
            cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
          }
          /*select an item randomly among those that are in the queue (where LL is not 0)*/
          p=runif(0,1.00*LLsum[t]);
          jj=0;
          for(int k=1;k<K[0];k++){
            if(p>cumsum[k]){
              jj=jj+1;
            }
          }	
          v[jj]=v[jj]+Score[t];/*update the selected item*/
          queue[jj]=0;/*remove the selected item from the queue*/
          LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
          LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
        }
        }
        }
        v[j]=v[j]+Score[s];	/*update item j*/
          }
        }
        }
        if(success==0){/*If no option for the paired update was found, put item j in the queue*/
          queue[j]=Score[s];
        }  
      }
    }
      }
      /*compute the average MSE across persons*/
      MSE[0]=0;
      for(int i=0;i<N[0];i++){
        dif=1.00*u[i]/k1[0]-1.00/(1+exp(-theta[i]));/*difference between the rating and the true value on the [0,1] scale*/
        MSE[0]=MSE[0]+dif*dif/N[0];
      }
    }
  }
  
  PutRNGstate();
}  

void urnings_RT(int*adaptive,int*paired,int*rule,int*u,int*v,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,double*cumsum,int*x,int*y,int*k1,int*k2,double*P,int*queue,int*LL,int*LLsum,int*Score,int*n_options,int*Upd,int*n_scores,double*Q,int*n_upd,int*N_obs1,int*Obs1,int*Exp1,int*N_obs2,int*Obs2,int*Exp2,int*N_obs3,int*Obs3,int*Exp3,int*N_obs4,int*Obs4,int*Exp4,int*N_obs5,int*Obs5,int*Exp5,int*N_obs6,int*Obs6,int*Exp6,int*N_obs7,int*Obs7,int*Exp7,int*N_obs8,int*Obs8,int*Exp8,int*N_obs9,int*Obs9,int*Exp9,int*W){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 double Correct=1.00;
 double Incorrect=1.00;
 int A=1;
 int j=0;
 int oldU=0;
 int oldV=0;
 int newU=0;
 int newV=0;
 int jj=0;
 int s=0;
 int success=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      /* Compute the normalising constant for the selection probabilities and a vector used for sampling from a multinomial distribution*/
      if(adaptive[rep]==1){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
          cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
        }
      }
      if(adaptive[rep]==0){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+1.00;
          cumsum[s+1]=cumsum[s]+1.00;
        }
      }
      /* sample an item given the selection probabilities*/
      p=runif(0,Mp);
      j=0;
      for(int s=1;s<K[0];s++){
        if(p>cumsum[s]){
          j=j+1;
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
	    /* generate the observed accuracy*/
             /*weight for the accuracy*/
	    L=1/(1+exp((delta[j+rep*K[0]]-theta[i+rep*N[0]])*W[0])); /* true probability correct*/
      p=runif(0,1);
      x[2]=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x[2]*W[0];                             
      v[j]=v[j]+(1-x[2])*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);
      /* generate the simulated response */
      p=runif(0,1);
      y[2]=1*(L>p);
      if(adaptive[rep]==0){
        N_obs1[u[i]+(k1[0]+5)*v[j]]=1+N_obs1[u[i]+(k1[0]+5)*v[j]];
        Obs1[u[i]+(k1[0]+5)*v[j]]=x[2]+Obs1[u[i]+(k1[0]+5)*v[j]];
        Exp1[u[i]+(k1[0]+5)*v[j]]=y[2]+Exp1[u[i]+(k1[0]+5)*v[j]];
      }
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y[2]*W[0];
      v[j]=v[j]-(1-y[2])*W[0];       

      if(rule[0]>0){ /*if RTs are used*/
	      /* if response is correct, or regarldess the response for Rule 1*/
        if((x[2]==1)|(rule[0]==1)){
	        for(int S=1;S>(-1);S--){/* loop over pseudo-responses with weights 2,1*/
            /* specify the correct weight*/
		        if(S==1){A=2;}
		        if(S==0){A=1;}
            /* generate the observed pseudo-response*/
      	    L=1/(1+exp((delta[j+rep*K[0]]-theta[i+rep*N[0]])*A));
      	    p=runif(0,1);
      	    x[S]=1*(L>p);
            /* add the balls based on the observed pseudo-response to the urns*/                      
		        u[i]=u[i]+x[S]*A;                             
      	    v[j]=v[j]+(1-x[S])*A;	
            /* compute the probability of X=1 given the urn configurations*/
      	    Correct=1.00;
      	    Incorrect=1.00;   
      	    for(int s=0;s<A;s++){              
         		  Correct=Correct*(u[i]-s)*(k2[0]+A-v[j]-s);
         		  Incorrect=Incorrect*(k1[0]+A-u[i]-s)*(v[j]-s);
            }
            L=Correct/(Correct+Incorrect);
            /* generate the simulated pseudo-response */
            p=runif(0,1);
            y[S]=1*(L>p);
            if(adaptive[rep]==0){
              if(S==1){
                if(x[2]==1){
                  N_obs2[u[i]+(k1[0]+2)*v[j]]=1+N_obs2[u[i]+(k1[0]+2)*v[j]];
                  Obs2[u[i]+(k1[0]+2)*v[j]]=x[1]+Obs2[u[i]+(k1[0]+2)*v[j]];
                  Exp2[u[i]+(k1[0]+2)*v[j]]=y[1]+Exp2[u[i]+(k1[0]+2)*v[j]];
                }
                if(x[2]==0){
                  N_obs3[u[i]+(k1[0]+2)*v[j]]=1+N_obs3[u[i]+(k1[0]+2)*v[j]];
                  Obs3[u[i]+(k1[0]+2)*v[j]]=x[1]+Obs3[u[i]+(k1[0]+2)*v[j]];
                  Exp3[u[i]+(k1[0]+2)*v[j]]=y[1]+Exp3[u[i]+(k1[0]+2)*v[j]];
                }
                  N_obs8[u[i]+(k1[0]+3)*v[j]]=1+N_obs8[u[i]+(k1[0]+3)*v[j]];
                  Obs8[u[i]+(k1[0]+3)*v[j]]=x[1]+Obs8[u[i]+(k1[0]+3)*v[j]];
                  Exp8[u[i]+(k1[0]+3)*v[j]]=y[1]+Exp8[u[i]+(k1[0]+3)*v[j]];
              }
              
              if(S==0){
                if((x[2]==1)&(x[1]==1)){
                  N_obs4[u[i]+(k1[0]+1)*v[j]]=1+N_obs4[u[i]+(k1[0]+1)*v[j]];
                  Obs4[u[i]+(k1[0]+1)*v[j]]=x[0]+Obs4[u[i]+(k1[0]+1)*v[j]];
                  Exp4[u[i]+(k1[0]+1)*v[j]]=y[0]+Exp4[u[i]+(k1[0]+1)*v[j]];
                }
                if((x[2]==1)&(x[1]==0)){
                  N_obs5[u[i]+(k1[0]+1)*v[j]]=1+N_obs5[u[i]+(k1[0]+1)*v[j]];
                  Obs5[u[i]+(k1[0]+1)*v[j]]=x[0]+Obs5[u[i]+(k1[0]+1)*v[j]];
                  Exp5[u[i]+(k1[0]+1)*v[j]]=y[0]+Exp5[u[i]+(k1[0]+1)*v[j]];  
                }
                if((x[2]==0)&(x[1]==1)){
                  N_obs6[u[i]+(k1[0]+1)*v[j]]=1+N_obs6[u[i]+(k1[0]+1)*v[j]];
                  Obs6[u[i]+(k1[0]+1)*v[j]]=x[0]+Obs6[u[i]+(k1[0]+1)*v[j]];
                  Exp6[u[i]+(k1[0]+1)*v[j]]=y[0]+Exp6[u[i]+(k1[0]+1)*v[j]];  
                }
                if((x[2]==0)&(x[1]==0)){
                  N_obs7[u[i]+(k1[0]+1)*v[j]]=1+N_obs7[u[i]+(k1[0]+1)*v[j]];
                  Obs7[u[i]+(k1[0]+1)*v[j]]=x[0]+Obs7[u[i]+(k1[0]+1)*v[j]];
                  Exp7[u[i]+(k1[0]+1)*v[j]]=y[0]+Exp7[u[i]+(k1[0]+1)*v[j]];
                }
                N_obs9[u[i]+(k1[0]+2)*v[j]]=1+N_obs9[u[i]+(k1[0]+2)*v[j]];
                  Obs9[u[i]+(k1[0]+2)*v[j]]=x[0]+Obs9[u[i]+(k1[0]+2)*v[j]];
                  Exp9[u[i]+(k1[0]+2)*v[j]]=y[0]+Exp9[u[i]+(k1[0]+2)*v[j]];
              }
            } 
            /* remove the balls based on the simulated pseudo-response from the urns: These would be the proposed values*/
            u[i]=u[i]-y[S]*A;
      	    v[j]=v[j]-(1-y[S])*A;       
          }
        }
      }
      /* These are the proposed values*/
      newV=v[j];
	    newU=u[i];
      
      if(adaptive[rep]==0){
        v[j]=newV;
        u[i]=newU;  
      }
      if(adaptive[rep]==1){
        /* go back to the old values and then decide on whether the proposed values (if different from old) should be accepted*/
        u[i]=oldU;
        v[j]=oldV;
        if(newU!=oldU){
          /*Compute the normalising constant for the selection probability given the proposed values*/
		      Mp1=P[newU+newV*(k1[0]+1)];
      	  for(int s=0;s<K[0];s++){
            if(s!=j){
          	  Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
            }
      	  }
          /* compute the MH acceptance probability*/
		      L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
          /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
          p=runif(0,1);
		      if(p<L){
			      u[i]=newU;
			      v[j]=newV;
          }
	      }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV=v[j];
	      v[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores[0];t++){
				    if((newV-oldV)==Score[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores[0];t++){
				    LLsum[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					    LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q[d]=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q[d]=Q[d]*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q[d]>0){/*If there are enough items*/
                n_upd[d]=n_upd[d]+1;
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue[j]=Score[s];
			    }  
        }
      }  
         
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
    }
  }
  PutRNGstate();
}  

void urnings_simpleX(int*adaptive,int*paired,int*u,int*v,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,int*k1,int*k2,double*P,double*cumsum,int*W,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 int oldU=0;
 int oldV=0;
 int newV=0;
 int newU=0;
 int j=0;
 int x=0;
 int y=0;
 int success=0;
 int Q=0;
 int s=0;
 int jj=0;
 double Correct=0;
 double Incorrect=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      if(adaptive[0]==0){/* sample an item with equal probabilities*/
        p=runif(0,1.00*K[0]);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>s){
            j=j+1;
          }
        }
      }
      if(adaptive[0]==1){      /*adaptive item selection, kernels are pre-specified in matrix P*/
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
          cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
        }
        p=runif(0,Mp);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>cumsum[s]){
            j=j+1;
          }
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
     /* generate the observed accuracy*/
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct, W is discrimination parameter, 1 is default value*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x*W[0];                             
      v[j]=v[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y*W[0];
      v[j]=v[j]-(1-y)*W[0];

      if(adaptive[0]==1){/*Metropolis step for when item selection is adaptive*/
        newU=u[i];
        newV=v[j];
        u[i]=oldU;
        v[j]=oldV;
        if(newU!=oldU){
          /*Compute the normalising constant for the selection probability given the proposed values*/
		      Mp1=P[newU+newV*(k1[0]+1)];
      	  for(int s=0;s<K[0];s++){
            if(s!=j){
          	  Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
            }
      	  }
          /* compute the MH acceptance probability*/
		      L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
          /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
          p=runif(0,1);
		      if(p<L){
			      u[i]=newU;
			      v[j]=newV;
          }
	      }
      }
      if(paired[0]==1){/*paired update for when it is on*/
	      /*set item j to the old value*/
        newV=v[j];
	      v[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired update can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores[0];t++){
				    if((newV-oldV)==Score[t]){
					    s=t;
				    }
			    }
          /*for every possible value of an update, count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores[0];t++){
				    LLsum[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					    LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update (in the rasch model case there is only one possible option for a paired update (-1 if item j gets an update +1, and +1 ifitem j gets an update of -1))*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue[j]=Score[s];
			    }  
        }
      }
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
    }
  }
  PutRNGstate();
}  

void urnings_combined_RT(int*adaptive,int*paired,int*u,int*v,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,int*k1,int*k2,double*P,double*cumsum,int*W,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 int oldU=0;
 int oldV=0;
 int newV=0;
 int newU=0;
 int j=0;
 int x=0;
 int y=0;
 int success=0;
 int Q=0;
 int s=0;
 int jj=0;
 double Correct=0;
 double Incorrect=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      /* sample an item with equal probabilities*/
      if(adaptive[0]==0){
        p=runif(0,1.00*K[0]);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>s){
            j=j+1;
          }
        }
      }
      /*adaptive item selection*/
      if(adaptive[0]==1){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
          cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
        }
        p=runif(0,Mp);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>cumsum[s]){
            j=j+1;
          }
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
      /* generate the observed accuracy*/
      W[0]=4;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x*W[0];                             
      v[j]=v[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y*W[0];
      v[j]=v[j]-(1-y)*W[0];
      
      W[0]=2;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x*W[0];                             
      v[j]=v[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y*W[0];
      v[j]=v[j]-(1-y)*W[0];

      W[0]=1;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x*W[0];                             
      v[j]=v[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y*W[0];
      v[j]=v[j]-(1-y)*W[0];
      
      if(adaptive[0]==1){
        newU=u[i];
        newV=v[j];
        u[i]=oldU;
        v[j]=oldV;
        if(newU!=oldU){
          /*Compute the normalising constant for the selection probability given the proposed values*/
		      Mp1=P[newU+newV*(k1[0]+1)];
      	  for(int s=0;s<K[0];s++){
            if(s!=j){
          	  Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
            }
      	  }
          /* compute the MH acceptance probability*/
		      L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
          /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
          p=runif(0,1);
		      if(p<L){
			      u[i]=newU;
			      v[j]=newV;
          }
	      }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV=v[j];
	      v[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores[0];t++){
				    if((newV-oldV)==Score[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores[0];t++){
				    LLsum[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					    LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue[j]=Score[s];
			    }  
        }
      }
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
    }
  }
  PutRNGstate();
}  

void urnings_simple(int*adaptive,int*paired,int*u,int*v,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,int*k1,int*k2,double*P,double*cumsum,int*W,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 int oldU=0;
 int oldV=0;
 int newV=0;
 int newU=0;
 int j=0;
 int x=0;
 int y=0;
 int success=0;
 int Q=0;
 int s=0;
 int jj=0;
 double Correct=0;
 double Incorrect=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      /* sample an item with equal probabilities*/
      if(adaptive[0]==0){
        p=runif(0,1.00*K[0]);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>s){
            j=j+1;
          }
        }
      }
      /*adaptive item selection*/
      if(adaptive[0]==1){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
          cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
        }
        p=runif(0,Mp);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>cumsum[s]){
            j=j+1;
          }
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
      /* generate the observed accuracy*/
      W[0]=4;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x*W[0];                             
      v[j]=v[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y*W[0];
      v[j]=v[j]-(1-y)*W[0];
      
      if(adaptive[0]==1){
        newU=u[i];
        newV=v[j];
        u[i]=oldU;
        v[j]=oldV;
        if(newU!=oldU){
          /*Compute the normalising constant for the selection probability given the proposed values*/
		      Mp1=P[newU+newV*(k1[0]+1)];
      	  for(int s=0;s<K[0];s++){
            if(s!=j){
          	  Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
            }
      	  }
          /* compute the MH acceptance probability*/
		      L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
          /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
          p=runif(0,1);
		      if(p<L){
			      u[i]=newU;
			      v[j]=newV;
          }
	      }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV=v[j];
	      v[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores[0];t++){
				    if((newV-oldV)==Score[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores[0];t++){
				    LLsum[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					    LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue[j]=Score[s];
			    }  
        }
      }
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
    }
  }
  PutRNGstate();
}  

void urnings_separate_RT(int*adaptive,int*paired,int*u,int*v,int*u2,int*v2,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,int*U2,int*V2,int*k1,int*k12,int*k2,int*k22,double*P,double*cumsum,int*W,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,int*Score2,int*n_scores2,int*n_options2,int*Upd2,int*queue2,int*LL2,int*LLsum2){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 int oldU=0;
 int oldV=0;
 int newV=0;
 int newU=0;
 int j=0;
 int x=0;
 int y=0;
 int success=0;
 int Q=0;
 int s=0;
 int jj=0;
 double Correct=0;
 double Incorrect=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      /* sample an item with equal probabilities*/
      if(adaptive[0]==0){
        p=runif(0,1.00*K[0]);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>s){
            j=j+1;
          }
        }
      }
      /*adaptive item selection*/
      if(adaptive[0]==1){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
          cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
        }
        p=runif(0,Mp);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>cumsum[s]){
            j=j+1;
          }
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
      /* generate the observed accuracy*/
      W[0]=4;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x*W[0];                             
      v[j]=v[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y*W[0];
      v[j]=v[j]-(1-y)*W[0];
      
      if(adaptive[0]==1){
        newU=u[i];
        newV=v[j];
        u[i]=oldU;
        v[j]=oldV;
        if(newU!=oldU){
          /*Compute the normalising constant for the selection probability given the proposed values*/
		      Mp1=P[newU+newV*(k1[0]+1)];
      	  for(int s=0;s<K[0];s++){
            if(s!=j){
          	  Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
            }
      	  }
          /* compute the MH acceptance probability*/
		      L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
          /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
          p=runif(0,1);
		      if(p<L){
			      u[i]=newU;
			      v[j]=newV;
          }
	      }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV=v[j];
	      v[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores[0];t++){
				    if((newV-oldV)==Score[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores[0];t++){
				    LLsum[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					    LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue[j]=Score[s];
			    }  
        }
      }
      
      oldU=u2[i];
      oldV=v2[j];
      W[0]=2;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u2[i]=u2[i]+x*W[0];                             
      v2[j]=v2[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u2[i]-s)*(k22[0]+W[0]-v2[j]-s);
        Incorrect=Incorrect*(k12[0]+W[0]-u2[i]-s)*(v2[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u2[i]=u2[i]-y*W[0];
      v2[j]=v2[j]-(1-y)*W[0];
      
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV=v2[j];
	      v2[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores2[0];t++){
				    if((newV-oldV)==Score2[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores2[0];t++){
				    LLsum2[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL2[ii+t*K[0]]=1*(queue2[ii]==Score2[t])*(v2[ii]+Score2[t]<k22[0]+1)*(v2[ii]+Score2[t]>(-1))*(ii!=j);
					    LLsum2[t]=LLsum2[t]+LL2[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options2[s];d<n_options2[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores2[0];t++){
						    Q=Q*(1-(LLsum2[t]<Upd2[t+d*n_scores2[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores2[0];t++){/*loop over possible update values*/
							    if(Upd2[t+d*n_scores2[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd2[t+d*n_scores2[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL2[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum2[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v2[jj]=v2[jj]+Score2[t];/*update the selected item*/
									    queue2[jj]=0;/*remove the selected item from the queue*/
									    LLsum2[t]=LLsum2[t]-LL2[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL2[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v2[j]=v2[j]+Score2[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue2[j]=Score2[s];
			    }  
        }
      }
    
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
      U2[i+rep*N[0]]=u2[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
      V2[j+rep*K[0]]=v2[j];
    }
  }
  PutRNGstate();
}  

void urnings_separate_RT2_real_data(int*nresp,int*P,int*I,int*K,int*X1,int*X2,int*Y1,int*Y2,int*n_obs_p,int*n_obs_i,int*u,int*v,int*u2,int*v2,int*Ubefore,int*Vbefore,int*U2before,int*V2before,int*Uafter,int*Vafter,int*Vafter2,int*U2after,int*V2after,int*V2after2,int*k1,int*k12,int*k2,int*k22,double*cumsum,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,int*Score2,int*n_scores2,int*n_options2,int*Upd2,int*queue2,int*LL2,int*LLsum2,int*N_obs_p,int*N_obs_i){
 double L=0;
 double p=0;
 int oldV=0;
 int newV=0;
 int oldV2=0;
 int newV2=0;
 int y=0;
 int success=0;
 int Q=0;
 int s=0;
 int jj=0;
 double Correct=0;
 double Incorrect=0;
 int i=0;
 int j=0;
 GetRNGstate();
 for(int r=0;r<nresp[0];r++){/*loop over response*/
  i=P[r]-1;
  j=I[r]-1;
  N_obs_p[r]=n_obs_p[i];
  N_obs_i[r]=n_obs_i[j];
  oldV=v[j];
  oldV2=v2[j];
  Ubefore[r]=u[i];
  Vbefore[r]=v[j];
  U2before[r]=u2[i];
  V2before[r]=v2[j]; 

  /*add balls based on observed data*/                  
	u[i]=u[i]+X1[r]*2;                             
  v[j]=v[j]+(1-X1[r])*2;	
  /* compute the probability of X=1 given the urn configurations*/
  Correct=1.00;
  Incorrect=1.00;   
  for(int s=0;s<2;s++){              
    Correct=Correct*(u[i]-s)*(k2[0]+2-v[j]-s);
    Incorrect=Incorrect*(k1[0]+2-u[i]-s)*(v[j]-s);
  }
  L=Correct/(Correct+Incorrect);      
  /* generate the simulated response */
  p=runif(0,1.00);
  y=1*(L>p);
  /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
  u[i]=u[i]-y*2;
  v[j]=v[j]-(1-y)*2;
  Y1[r]=y;

  /*add balls based on observed data*/
  u2[i]=u2[i]+X2[r];                             
  v2[j]=v2[j]+(1-X2[r]);	
  /* compute the probability of X=1 given the urn configurations*/
  Correct=1.00*u2[i]*(k22[0]+1-v2[j]);
  Incorrect=1.00*(k12[0]+1-u2[i])*v2[j];
  L=Correct/(Correct+Incorrect);      
  /* generate the simulated response */
  p=runif(0,1.00);
  y=1*(L>p);
  /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
  u2[i]=u2[i]-y;
  v2[j]=v2[j]-(1-y);
  Y2[r]=y;
    
  Uafter[r]=u[i];
  U2after[r]=u2[i];
  Vafter[r]=v[j];
  V2after[r]=v2[j];

  /*set item j to the old value*/
  newV=v[j];
	v[j]=oldV;
  /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	if(newV!=oldV){
		/*Find which of the possible update values (saved in vector Score) has been proposed*/
		for(int t=0;t<n_scores[0];t++){
			if((newV-oldV)==Score[t]){
				s=t;
			}
		}
    /*for every value of the paired update count how many items are in the queue and can be updated*/
	  for(int t=0;t<n_scores[0];t++){
			LLsum[t]=0;
			for(int ii=0;ii<K[0];ii++){
				LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
				LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
			}
		}
    success=0;  
		for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
			if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
			}
		}
		if(success==0){/*If no option for the paired update was found, put item j in the queue*/
			queue[j]=Score[s];
		}  
  }

Vafter2[r]=v[j];

newV2=v2[j];
v2[j]=oldV2;

if(newV2!=oldV2){
		/*Find which of the possible update values (saved in vector Score) has been proposed*/
		for(int t=0;t<n_scores2[0];t++){
			if((newV2-oldV2)==Score2[t]){
				s=t;
			}
		}
    /*for every value of the paired update count how many items are in the queue and can be updated*/
	  for(int t=0;t<n_scores2[0];t++){
			LLsum2[t]=0;
			for(int ii=0;ii<K[0];ii++){
				LL2[ii+t*K[0]]=1*(queue2[ii]==Score2[t])*(v2[ii]+Score2[t]<k22[0]+1)*(v2[ii]+Score2[t]>(-1))*(ii!=j);
				LLsum2[t]=LLsum2[t]+LL2[ii+t*K[0]];
			}
		}
    success=0;  
		for(int d=n_options2[s];d<n_options2[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
			if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores2[0];t++){
						    Q=Q*(1-(LLsum2[t]<Upd2[t+d*n_scores2[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores2[0];t++){/*loop over possible update values*/
							    if(Upd2[t+d*n_scores2[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd2[t+d*n_scores2[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL2[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum2[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v2[jj]=v2[jj]+Score2[t];/*update the selected item*/
									    queue2[jj]=0;/*remove the selected item from the queue*/
									    LLsum2[t]=LLsum2[t]-LL2[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL2[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v2[j]=v2[j]+Score2[s];	/*update item j*/
					    }
			}
		}
		if(success==0){/*If no option for the paired update was found, put item j in the queue*/
			queue2[j]=Score2[s];
		}  
  }
      
  V2after2[r]=v2[j];
  n_obs_p[i]=n_obs_p[i]+1;
  n_obs_i[j]=n_obs_i[j]+1;
}
  PutRNGstate();
  
  }  

void urnings_combined_RT_real_data(int*nresp,int*P,int*I,int*K,int*X1,int*X2,int*Y1,int*Y2,int*n_obs_p,int*n_obs_i,int*u,int*v,int*Ubefore,int*Vbefore,int*U_middle,int*V_middle,int*Uafter,int*Vafter,int*Vafter2,int*k1,int*k2,double*cumsum,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,int*N_obs_p,int*N_obs_i){
 double L=0;
 double p=0;
 int oldV=0;
 int newV=0;
 int y=0;
 int success=0;
 int Q=0;
 int s=0;
 int jj=0;
 double Correct=0;
 double Incorrect=0;
 int i=0;
 int j=0;
 GetRNGstate();
 for(int r=0;r<nresp[0];r++){/*loop over response*/
  i=P[r]-1;
  j=I[r]-1;
  N_obs_p[r]=n_obs_p[i];
  N_obs_i[r]=n_obs_i[j];
  oldV=v[j];
  Ubefore[r]=u[i];
  Vbefore[r]=v[j];
 
  /*add balls based on observed data*/                  
	u[i]=u[i]+X1[r]*2;                             
  v[j]=v[j]+(1-X1[r])*2;	
  /* compute the probability of X=1 given the urn configurations*/
  Correct=1.00;
  Incorrect=1.00;   
  for(int s=0;s<2;s++){              
    Correct=Correct*(u[i]-s)*(k2[0]+2-v[j]-s);
    Incorrect=Incorrect*(k1[0]+2-u[i]-s)*(v[j]-s);
  }
  L=Correct/(Correct+Incorrect);      
  /* generate the simulated response */
  p=runif(0,1.00);
  y=1*(L>p);
  /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
  u[i]=u[i]-y*2;
  v[j]=v[j]-(1-y)*2;
  Y1[r]=y;

  U_middle[r]=u[i];
  V_middle[r]=v[j];
  /*add balls based on observed data*/
  u[i]=u[i]+X2[r];                             
  v[j]=v[j]+(1-X2[r]);	
  /* compute the probability of X=1 given the urn configurations*/
  Correct=1.00*u[i]*(k2[0]+1-v[j]);
  Incorrect=1.00*(k1[0]+1-u[i])*v[j];
  L=Correct/(Correct+Incorrect);      
  /* generate the simulated response */
  p=runif(0,1.00);
  y=1*(L>p);
  /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
  u[i]=u[i]-y;
  v[j]=v[j]-(1-y);
  Y2[r]=y;
    
  Uafter[r]=u[i];
  Vafter[r]=v[j];
  
  /*set item j to the old value*/
  newV=v[j];
	v[j]=oldV;
  /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	if(newV!=oldV){
		/*Find which of the possible update values (saved in vector Score) has been proposed*/
		for(int t=0;t<n_scores[0];t++){
			if((newV-oldV)==Score[t]){
				s=t;
			}
		}
    /*for every value of the paired update count how many items are in the queue and can be updated*/
	  for(int t=0;t<n_scores[0];t++){
			LLsum[t]=0;
			for(int ii=0;ii<K[0];ii++){
				LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
				LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
			}
		}
    success=0;  
		for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
			if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
			}
		}
		if(success==0){/*If no option for the paired update was found, put item j in the queue*/
			queue[j]=Score[s];
		}  
  }

Vafter2[r]=v[j];


  n_obs_p[i]=n_obs_p[i]+1;
  n_obs_i[j]=n_obs_i[j]+1;
}
  PutRNGstate();
  
  }  

void urnings_separate_RT2(int*adaptive,int*paired,int*u,int*v,int*u2,int*v2,int*u3,int*v3,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,int*U2,int*V2,int*U3,int*V3,int*k1,int*k12,int*k13,int*k2,int*k22,int*k23,double*P,double*cumsum,int*W,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,int*Score2,int*n_scores2,int*n_options2,int*Upd2,int*queue2,int*LL2,int*LLsum2,int*Score3,int*n_scores3,int*n_options3,int*Upd3,int*queue3,int*LL3,int*LLsum3){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 int oldU=0;
 int oldV=0;
 int newV=0;
 int newU=0;
 int oldU2=0;
 int oldV2=0;
 int newU2=0;
 int newV2=0;
 int oldU3=0;
 int oldV3=0;
 int newU3=0;
 int newV3=0;
 int j=0;
 int x=0;
 int y=0;
 int success=0;
 int Q=0;
 int s=0;
 int jj=0;
 double Correct=0;
 double Incorrect=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      /* sample an item with equal probabilities*/
      if(adaptive[0]==0){
        p=runif(0,1.00*K[0]);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>s){
            j=j+1;
          }
        }
      }
      /*adaptive item selection*/
      if(adaptive[0]==1){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+P[u[i]+u2[i]+u3[i]+(v[s]+v2[s]+v3[s])*(k1[0]+k12[0]+k13[0]+1)];
          cumsum[s+1]=cumsum[s]+P[u[i]+u2[i]+u3[i]+(v[s]+v2[s]+v3[s])*(k1[0]+k12[0]+k13[0]+1)];
        }
        p=runif(0,Mp);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>cumsum[s]){
            j=j+1;
          }
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
      /* generate the observed accuracy*/
      W[0]=4;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x*W[0];                             
      v[j]=v[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y*W[0];
      v[j]=v[j]-(1-y)*W[0];
      
    oldU2=u2[i];
    oldV2=v2[j];  
 
      W[0]=2;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u2[i]=u2[i]+x*W[0];                             
      v2[j]=v2[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u2[i]-s)*(k22[0]+W[0]-v2[j]-s);
        Incorrect=Incorrect*(k12[0]+W[0]-u2[i]-s)*(v2[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u2[i]=u2[i]-y*W[0];
      v2[j]=v2[j]-(1-y)*W[0];
    
    oldU3=u3[i];
    oldV3=v3[j];

    W[0]=1;
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u3[i]=u3[i]+x*W[0];                             
      v3[j]=v3[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u3[i]-s)*(k23[0]+W[0]-v3[j]-s);
        Incorrect=Incorrect*(k13[0]+W[0]-u3[i]-s)*(v3[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u3[i]=u3[i]-y*W[0];
      v3[j]=v3[j]-(1-y)*W[0];
      
    

    if(adaptive[0]==1){
        newU=u[i];
        newV=v[j];
        newU2=u2[i];
        newV2=v2[j];
        newU3=u3[i];
        newV3=v3[j];
        u[i]=oldU;
        v[j]=oldV;
        u2[i]=oldU2;
        v2[j]=oldV2;
        u3[i]=oldU3;
        v3[j]=oldV3;
        if((newU+newU2+newU3)!=(oldU+oldU2+oldU3)){
          /*Compute the normalising constant for the selection probability given the proposed values*/
		      Mp1=P[newU+newU2+newU3+(newV+newV2+newV3)*(k1[0]+k12[0]+k13[0]+1)];
      	  for(int s=0;s<K[0];s++){
            if(s!=j){
          	  Mp1=Mp1+P[newU+newU2+newU3+(v[s]+v2[s]+v3[s])*(k1[0]+k12[0]+k13[0]+1)];
            }
      	  }
          /* compute the MH acceptance probability*/
		      L=P[newU+newU2+newU3+(newV+newV2+newV3)*(k1[0]+k12[0]+k13[0]+1)]/Mp1*Mp/P[oldU+oldU2+oldU3+(oldV+oldV2+oldV3)*(k1[0]+k12[0]+k13[0]+1)];
          /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
          p=runif(0,1.00);
		      if(p<L){
			      u[i]=newU;
			      v[j]=newV;
            u2[i]=newU2;
            v2[j]=newV2;
            u3[i]=newU3;
            v3[j]=newV3;
          }
	      }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV=v[j];
	      v[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores[0];t++){
				    if((newV-oldV)==Score[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores[0];t++){
				    LLsum[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					    LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue[j]=Score[s];
			    }  
        }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV2=v2[j];
	      v2[j]=oldV2;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV2!=oldV2){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores2[0];t++){
				    if((newV2-oldV2)==Score2[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores2[0];t++){
				    LLsum2[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL2[ii+t*K[0]]=1*(queue2[ii]==Score2[t])*(v2[ii]+Score2[t]<k22[0]+1)*(v2[ii]+Score2[t]>(-1))*(ii!=j);
					    LLsum2[t]=LLsum2[t]+LL2[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options2[s];d<n_options2[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores2[0];t++){
						    Q=Q*(1-(LLsum2[t]<Upd2[t+d*n_scores2[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores2[0];t++){/*loop over possible update values*/
							    if(Upd2[t+d*n_scores2[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd2[t+d*n_scores2[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL2[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum2[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v2[jj]=v2[jj]+Score2[t];/*update the selected item*/
									    queue2[jj]=0;/*remove the selected item from the queue*/
									    LLsum2[t]=LLsum2[t]-LL2[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL2[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v2[j]=v2[j]+Score2[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue2[j]=Score2[s];
			    }  
        }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV3=v3[j];
	      v3[j]=oldV3;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV3!=oldV3){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores3[0];t++){
				    if((newV3-oldV3)==Score3[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores3[0];t++){
				    LLsum3[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL3[ii+t*K[0]]=1*(queue3[ii]==Score3[t])*(v3[ii]+Score3[t]<k23[0]+1)*(v3[ii]+Score3[t]>(-1))*(ii!=j);
					    LLsum3[t]=LLsum3[t]+LL3[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options3[s];d<n_options3[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores3[0];t++){
						    Q=Q*(1-(LLsum3[t]<Upd3[t+d*n_scores3[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores3[0];t++){/*loop over possible update values*/
							    if(Upd3[t+d*n_scores3[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd3[t+d*n_scores3[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL3[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum3[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v3[jj]=v3[jj]+Score3[t];/*update the selected item*/
									    queue3[jj]=0;/*remove the selected item from the queue*/
									    LLsum3[t]=LLsum3[t]-LL3[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL3[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v3[j]=v3[j]+Score3[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue3[j]=Score3[s];
			    }  
        }
      }
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
      U2[i+rep*N[0]]=u2[i];
      U3[i+rep*N[0]]=u3[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
      V2[j+rep*K[0]]=v2[j];
      V3[j+rep*K[0]]=v3[j];
    }
  }
  PutRNGstate();
}  

void urnings_simple_RT(int*adaptive,int*paired,int*u,int*v,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,int*k1,int*k2,double*P,double*cumsum,int*W,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,int*rule){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 int oldU=0;
 int oldV=0;
 int newV=0;
 int newU=0;
 int j=0;
 int x=0;
 int y=0;
 int success=0;
 int Q=0;
 int s=0;
 int jj=0;
 double Correct=0;
 double Incorrect=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      /* sample an item with equal probabilities*/
      if(adaptive[0]==0){
        p=runif(0,1.00*K[0]);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>s){
            j=j+1;
          }
        }
      }
      /*adaptive item selection*/
      if(adaptive[0]==1){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
          cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
        }
        p=runif(0,Mp);
        j=0;
        for(int s=1;s<K[0];s++){
          if(p>cumsum[s]){
            j=j+1;
          }
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
     /* generate the observed accuracy*/
      L=1/(1+exp(W[0]*(delta[j+rep*K[0]]-theta[i+rep*N[0]]))); /* true probability correct*/
      p=runif(0,1.00);
      x=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x*W[0];                             
      v[j]=v[j]+(1-x)*W[0];	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<W[0];s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+W[0]-v[j]-s);
        Incorrect=Incorrect*(k1[0]+W[0]-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);      
      /* generate the simulated response */
      p=runif(0,1.00);
      y=1*(L>p);
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y*W[0];
      v[j]=v[j]-(1-y)*W[0];
      if((rule[0]==1)|(x==1)){
      for(int S=2;S>0;S--){/* loop over pseudo-responses with weights 2,1*/
            /* generate the observed pseudo-response*/
      	    L=1/(1+exp((delta[j+rep*K[0]]-theta[i+rep*N[0]])*S));
      	    p=runif(0,1);
      	    x=1*(L>p);
            /* add the balls based on the observed pseudo-response to the urns*/                      
		        u[i]=u[i]+x*S;                             
      	    v[j]=v[j]+(1-x)*S;	
            /* compute the probability of X=1 given the urn configurations*/
      	    Correct=1.00;
      	    Incorrect=1.00;   
      	    for(int s=0;s<S;s++){              
         		  Correct=Correct*(u[i]-s)*(k2[0]+S-v[j]-s);
         		  Incorrect=Incorrect*(k1[0]+S-u[i]-s)*(v[j]-s);
            }
            L=Correct/(Correct+Incorrect);
            /* generate the simulated pseudo-response */
            p=runif(0,1);
            y=1*(L>p);
            /* remove the balls based on the simulated pseudo-response from the urns: These would be the proposed values*/
            u[i]=u[i]-y*S;
      	    v[j]=v[j]-(1-y)*S;       
      }
      }
      if(adaptive[0]==1){
        newU=u[i];
        newV=v[j];
        u[i]=oldU;
        v[j]=oldV;
        if(newU!=oldU){
          /*Compute the normalising constant for the selection probability given the proposed values*/
		      Mp1=P[newU+newV*(k1[0]+1)];
      	  for(int s=0;s<K[0];s++){
            if(s!=j){
          	  Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
            }
      	  }
          /* compute the MH acceptance probability*/
		      L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
          /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
          p=runif(0,1);
		      if(p<L){
			      u[i]=newU;
			      v[j]=newV;
          }
	      }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV=v[j];
	      v[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores[0];t++){
				    if((newV-oldV)==Score[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores[0];t++){
				    LLsum[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					    LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q=Q*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q>0){/*If there are enough items*/
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,1.00*LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue[j]=Score[s];
			    }  
        }
      }
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
    }
  }
  PutRNGstate();
}  

void urnings_RT_misspecification1(int*adaptive,int*paired,int*rule,int*u,int*v,double*theta,double*delta,int*N,int*K,int*Rep,int*U,int*V,double*cumsum,int*x,int*y,int*k1,int*k2,double*P,int*queue,int*LL,int*LLsum,int*Score,int*n_options,int*Upd,int*n_scores,double*Q,int*n_upd,int*N_obs1,int*Obs1,int*Exp1,int*N_obs2,int*Obs2,int*Exp2,int*N_obs3,int*Obs3,int*Exp3,int*N_obs4,int*Obs4,int*Exp4,int*N_obs5,int*Obs5,int*Exp5,int*N_obs6,int*Obs6,int*Exp6,int*N_obs7,int*Obs7,int*Exp7,int*N_obs8,int*Obs8,int*Exp8,int*N_obs9,int*Obs9,int*Exp9){
 double L=0;
 double p=0;
 double Mp=0;
 double Mp1=0;
 double Correct=1.00;
 double Incorrect=1.00;
 int A=1;
 int j=0;
 int oldU=0;
 int oldV=0;
 int newU=0;
 int newV=0;
 int jj=0;
 int s=0;
 int success=0;
 GetRNGstate();
 for(int rep=0;rep<Rep[0];rep++){/* loop over iterations*/
    for(int i=0;i<N[0];i++){/* loop over persons*/
      /* Compute the normalising constant for the selection probabilities and a vector used for sampling from a multinomial distribution*/
      if(adaptive[rep]==1){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+P[u[i]+v[s]*(k1[0]+1)];
          cumsum[s+1]=cumsum[s]+P[u[i]+v[s]*(k1[0]+1)];
        }
      }
      if(adaptive[rep]==0){
        Mp=0;		
        for(int s=0;s<K[0];s++){
          Mp=Mp+1.00;
          cumsum[s+1]=cumsum[s]+1.00;
        }
      }
      /* sample an item given the selection probabilities*/
      p=runif(0,Mp);
      j=0;
      for(int s=1;s<K[0];s++){
        if(p>cumsum[s]){
          j=j+1;
        }
      }
      /* save the current values of the person and the item */
      oldU=u[i];
      oldV=v[j];
	    /* generate the observed accuracy*/
      A=4;        /*weight for the accuracy*/
	    L=1/(1+exp((delta[j+rep*K[0]]-theta[i+rep*N[0]])*A)); /* true probability correct*/
      p=runif(0,1);
      x[2]=1*(L>p);
      /* add the balls based on the observed response to the urns*/                      
	    u[i]=u[i]+x[2]*A;                             
      v[j]=v[j]+(1-x[2])*A;	
      /* compute the probability of X=1 given the urn configurations*/
      Correct=1.00;
      Incorrect=1.00;   
      for(int s=0;s<A;s++){              
        Correct=Correct*(u[i]-s)*(k2[0]+A-v[j]-s);
        Incorrect=Incorrect*(k1[0]+A-u[i]-s)*(v[j]-s);
      }
      L=Correct/(Correct+Incorrect);
      /* generate the simulated response */
      p=runif(0,1);
      y[2]=1*(L>p);
      if(adaptive[rep]==0){
        N_obs1[u[i]+(k1[0]+5)*v[j]]=1+N_obs1[u[i]+(k1[0]+5)*v[j]];
        Obs1[u[i]+(k1[0]+5)*v[j]]=x[2]+Obs1[u[i]+(k1[0]+5)*v[j]];
        Exp1[u[i]+(k1[0]+5)*v[j]]=y[2]+Exp1[u[i]+(k1[0]+5)*v[j]];
      }
      /* remove the balls based on the simulated response from the urns: These would be the proposed values*/
      u[i]=u[i]-y[2]*A;
      v[j]=v[j]-(1-y[2])*A;       

      if(rule[0]>0){ /*if RTs are used*/
	      /* if response is correct, or regarldess the response for Rule 1*/
        if((x[2]==1)|(rule[0]==1)){
	        for(int S=1;S>(-1);S--){/* loop over pseudo-responses with weights 2,1*/
            /* specify the correct weight*/
		        if(S==1){A=2;}
		        if(S==0){A=1;}
            /* generate the observed pseudo-response*/
      	    L=1/(1+exp((delta[j+rep*K[0]]-theta[i+rep*N[0]]*(x[2]==1))*A));
      	    p=runif(0,1);
      	    x[S]=1*(L>p);
            /* add the balls based on the observed pseudo-response to the urns*/                      
		        u[i]=u[i]+x[S]*A;                             
      	    v[j]=v[j]+(1-x[S])*A;	
            /* compute the probability of X=1 given the urn configurations*/
      	    Correct=1.00;
      	    Incorrect=1.00;   
      	    for(int s=0;s<A;s++){              
         		  Correct=Correct*(u[i]-s)*(k2[0]+A-v[j]-s);
         		  Incorrect=Incorrect*(k1[0]+A-u[i]-s)*(v[j]-s);
            }
            L=Correct/(Correct+Incorrect);
            /* generate the simulated pseudo-response */
            p=runif(0,1);
            y[S]=1*(L>p);
            if(adaptive[rep]==0){
              if(S==1){
                if(x[2]==1){
                  N_obs2[u[i]+(k1[0]+2)*v[j]]=1+N_obs2[u[i]+(k1[0]+2)*v[j]];
                  Obs2[u[i]+(k1[0]+2)*v[j]]=x[1]+Obs2[u[i]+(k1[0]+2)*v[j]];
                  Exp2[u[i]+(k1[0]+2)*v[j]]=y[1]+Exp2[u[i]+(k1[0]+2)*v[j]];
                }
                if(x[2]==0){
                  N_obs3[u[i]+(k1[0]+2)*v[j]]=1+N_obs3[u[i]+(k1[0]+2)*v[j]];
                  Obs3[u[i]+(k1[0]+2)*v[j]]=x[1]+Obs3[u[i]+(k1[0]+2)*v[j]];
                  Exp3[u[i]+(k1[0]+2)*v[j]]=y[1]+Exp3[u[i]+(k1[0]+2)*v[j]];
                }
                  N_obs8[u[i]+(k1[0]+3)*v[j]]=1+N_obs8[u[i]+(k1[0]+3)*v[j]];
                  Obs8[u[i]+(k1[0]+3)*v[j]]=x[1]+Obs8[u[i]+(k1[0]+3)*v[j]];
                  Exp8[u[i]+(k1[0]+3)*v[j]]=y[1]+Exp8[u[i]+(k1[0]+3)*v[j]];
              }
              
              if(S==0){
                if((x[2]==1)&(x[1]==1)){
                  N_obs4[u[i]+(k1[0]+1)*v[j]]=1+N_obs4[u[i]+(k1[0]+1)*v[j]];
                  Obs4[u[i]+(k1[0]+1)*v[j]]=x[0]+Obs4[u[i]+(k1[0]+1)*v[j]];
                  Exp4[u[i]+(k1[0]+1)*v[j]]=y[0]+Exp4[u[i]+(k1[0]+1)*v[j]];
                }
                if((x[2]==1)&(x[1]==0)){
                  N_obs5[u[i]+(k1[0]+1)*v[j]]=1+N_obs5[u[i]+(k1[0]+1)*v[j]];
                  Obs5[u[i]+(k1[0]+1)*v[j]]=x[0]+Obs5[u[i]+(k1[0]+1)*v[j]];
                  Exp5[u[i]+(k1[0]+1)*v[j]]=y[0]+Exp5[u[i]+(k1[0]+1)*v[j]];  
                }
                if((x[2]==0)&(x[1]==1)){
                  N_obs6[u[i]+(k1[0]+1)*v[j]]=1+N_obs6[u[i]+(k1[0]+1)*v[j]];
                  Obs6[u[i]+(k1[0]+1)*v[j]]=x[0]+Obs6[u[i]+(k1[0]+1)*v[j]];
                  Exp6[u[i]+(k1[0]+1)*v[j]]=y[0]+Exp6[u[i]+(k1[0]+1)*v[j]];  
                }
                if((x[2]==0)&(x[1]==0)){
                  N_obs7[u[i]+(k1[0]+1)*v[j]]=1+N_obs7[u[i]+(k1[0]+1)*v[j]];
                  Obs7[u[i]+(k1[0]+1)*v[j]]=x[0]+Obs7[u[i]+(k1[0]+1)*v[j]];
                  Exp7[u[i]+(k1[0]+1)*v[j]]=y[0]+Exp7[u[i]+(k1[0]+1)*v[j]];
                }
                N_obs9[u[i]+(k1[0]+2)*v[j]]=1+N_obs9[u[i]+(k1[0]+2)*v[j]];
                  Obs9[u[i]+(k1[0]+2)*v[j]]=x[0]+Obs9[u[i]+(k1[0]+2)*v[j]];
                  Exp9[u[i]+(k1[0]+2)*v[j]]=y[0]+Exp9[u[i]+(k1[0]+2)*v[j]];
              }
            } 
            /* remove the balls based on the simulated pseudo-response from the urns: These would be the proposed values*/
            u[i]=u[i]-y[S]*A;
      	    v[j]=v[j]-(1-y[S])*A;       
          }
        }
      }
      /* These are the proposed values*/
      newV=v[j];
	    newU=u[i];
      
      if(adaptive[rep]==0){
        v[j]=newV;
        u[i]=newU;  
      }
      if(adaptive[rep]==1){
        /* go back to the old values and then decide on whether the proposed values (if different from old) should be accepted*/
        u[i]=oldU;
        v[j]=oldV;
        if(newU!=oldU){
          /*Compute the normalising constant for the selection probability given the proposed values*/
		      Mp1=P[newU+newV*(k1[0]+1)];
      	  for(int s=0;s<K[0];s++){
            if(s!=j){
          	  Mp1=Mp1+P[newU+v[s]*(k1[0]+1)];
            }
      	  }
          /* compute the MH acceptance probability*/
		      L=P[newU+newV*(k1[0]+1)]/Mp1*Mp/P[oldU+oldV*(k1[0]+1)];
          /* Generate a random uniform (0,1) number to decide whether to accept the proposal*/
          p=runif(0,1);
		      if(p<L){
			      u[i]=newU;
			      v[j]=newV;
          }
	      }
      }
      if(paired[0]==1){
	      /*set item j to the old value*/
        newV=v[j];
	      v[j]=oldV;
        /*If the proposed value is different from the old one, check if a paired updated can be performed*/ 
	      if(newV!=oldV){
			    /*Find which of the possible update values (saved in vector Score) has been proposed*/
			    for(int t=0;t<n_scores[0];t++){
				    if((newV-oldV)==Score[t]){
					    s=t;
				    }
			    }
          /*for every value of the paired update count how many items are in the queue and can be updated*/
	        for(int t=0;t<n_scores[0];t++){
				    LLsum[t]=0;
				    for(int ii=0;ii<K[0];ii++){
					    LL[ii+t*K[0]]=1*(queue[ii]==Score[t])*(v[ii]+Score[t]<k2[0]+1)*(v[ii]+Score[t]>(-1))*(ii!=j);
					    LLsum[t]=LLsum[t]+LL[ii+t*K[0]];
				    }
			    }
          success=0;  
			    for(int d=n_options[s];d<n_options[s+1];d++){/*Loop over possible options of a paired update for the particular proposed update*/
				    if(success==0){/*this is done only if the paired update has not been found yet*/
					    /*check whether for the particular possible update there are enough items in the queue*/
              Q[d]=1;
					    for(int t=0;t<n_scores[0];t++){
						    Q[d]=Q[d]*(1-(LLsum[t]<Upd[t+d*n_scores[0]]));
					    }
					    if(Q[d]>0){/*If there are enough items*/
                n_upd[d]=n_upd[d]+1;
						    success=1;/*set succes to 1, such that we will not look further*/
						    for(int t=0;t<n_scores[0];t++){/*loop over possible update values*/
							    if(Upd[t+d*n_scores[0]]>0){/*if this update value is part of the paired update option*/
								    for(int h=1;h<(Upd[t+d*n_scores[0]]+1);h++){
									    /*update the values in the cumsum vector which is used for sampling from a multinomial distribution*/
                      for(int k=0;k<K[0];k++){
										    cumsum[k+1]=cumsum[k]+1.00*LL[k+t*K[0]];
									    }
                      /*select an item randomly among those that are in the queue (where LL is not 0)*/
									    p=runif(0,LLsum[t]);
									    jj=0;
									    for(int k=1;k<K[0];k++){
										    if(p>cumsum[k]){
											    jj=jj+1;
										    }
									    }	
									    v[jj]=v[jj]+Score[t];/*update the selected item*/
									    queue[jj]=0;/*remove the selected item from the queue*/
									    LLsum[t]=LLsum[t]-LL[jj+t*K[0]];/*update the normalising constant for the selection probabilities*/	
									    LL[jj+t*K[0]]=0;/*set the selection probability for the item to 0 such that it will not be selected again*/
								    }
							    }
						    }
						    v[j]=v[j]+Score[s];	/*update item j*/
					    }
				    }
			    }
			    if(success==0){/*If no option for the paired update was found, put item j in the queue*/
				    queue[j]=Score[s];
			    }  
        }
      }  
         
      /* save the values of the person urns*/
      U[i+rep*N[0]]=u[i];
    }
    /* loop over the items and save the values a the end of the iteration*/
    for(int j=0;j<K[0];j++){
      V[j+rep*K[0]]=v[j];
    }
  }
  PutRNGstate();
}  
