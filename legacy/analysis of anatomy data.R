# load the data
logs=read.csv("empirical_data/answers.csv")

# compute the accuracy of the response
logs$accuracy=1*(logs$item_asked==logs$item_answered)

# select data with a single type of items and from a single country
data=logs[(logs$type=='t2d')&(logs$ip_country=='CZ'),]

# remove persons with missing IDs 
data=data[!is.na(data$user),]
# save selected data
#save(data,file='anatomy_subset.RData')


# order the data chronologically based on variable time
library(tidyverse)
#load("data_analysis/anatomy_subset.RData")  
options(digits.secs = 6)
data$time = ymd_hms(data$time)
data_ordered <- data[order(data$time), ]

# save ordered data
#saveRDS(data_ordered, "data_ordered.rds")

#read the file with readRDS
#data<-readRDS("data_ordered.rds")

# remove persons with less than 10 responses
responses_per_person=table(data$user)
persons=as.numeric(names(responses_per_person[responses_per_person>9])) # these are persons with at least 10 responses

data1=data[is.element(data$user,persons),]
length(unique(data1$item_asked))
length(unique(data1$user))

# remove items with less than 100 responses
responses_per_item=table(data1$item_asked)
items=as.numeric(names(responses_per_item[responses_per_item>99])) # these are items with at least 100 responses

data2=data1[is.element(data1$item_asked,items),]

# for every item determine what is fast and what is slow [1 - grand median, 2 - item-level median, 3 - half of the time in which 90% responses were given]

# grand median
threshold1=median(data2$response_time)

threshold2<-threshold3<-numeric(length(items))
for(j in 1:length(items)){
	threshold2[j]=median(data2$response_time[data2$item_asked==items[j]]) # item-level median of RT
	threshold3[j]=quantile(data2$response_time[data2$item_asked==items[j]],0.9)/2 # half of the time in which 90% responses were given
}

# item-specific proportion of fast responses
p_fast1<-p_fast3<-numeric(length(items))
for(j in 1:length(items)){
	p_fast1[j]=mean(data2$response_time[data2$item_asked==items[j]]<threshold1)
	p_fast3[j]=mean(data2$response_time[data2$item_asked==items[j]]<threshold3[j])
}

# instead of the original item indicator make an item ID from 1 to the number of items
data2$item_id=0
for(j in 1:length(items)){
data2$item_id[which(data2$item_asked==items[j])]=j
}

# Create variable indicating whether each of the given responses is fast (3 different variables with different definitions of 'fast')
data2$FS1<-1*(data2$response_time<threshold1)
data2$FS2<-1*(data2$response_time<threshold2[data2$item_id])
data2$FS3<-1*(data2$response_time<threshold3[data2$item_id])

# instead of the original person indicator make a person ID from 1 to the number of persons
data2$person_id=0
for(j in 1:length(persons)){
data2$person_id[data2$user==persons[j]]=j
}

# create pseudo-response time variables for each response (separately for the three ways of defining 'fast')
data2$Y1<-1*((data2$accuracy==1)&(data2$FS1==1))+1*((data2$accuracy==0)&(data2$FS1==0))
data2$Y2<-1*((data2$accuracy==1)&(data2$FS2==1))+1*((data2$accuracy==0)&(data2$FS2==0))
data2$Y3<-1*((data2$accuracy==1)&(data2$FS3==1))+1*((data2$accuracy==0)&(data2$FS3==0))

# now the data are ready for analysis with urnings and is saved
#save(data2,file='anatomy_for_urnings.RData')

# load C code. It includes functions will run urnings for empirical data and saves urning values before and after the update for every response
# here we only use 1 pseudoRT variable in addition to accuracy
dyn.load("CRT_old.so")

# Analysis with separate urns for accuracy and pseudoRT.

# set urn sizes
# for accuracy
k1=20 # persons
k2=20 # items
# for pseudoRT
k12=10# persons
k22=10# items

# set input variable for the c code
nresp=nrow(data2) # number of responses
P=data2$person_id # vector of person IDs for all responses
I=data2$item_id # vector of item IDs for all response
N=length(persons) # number of persons
K=length(items) # number of items
X1=data2$accuracy # accuracy data
X2=data2$Y1 # pseudoRT data 
long=rep(0,nresp) # empty vector with as many elements as there are responses, it is used to save all sort of things within the c code

# set initial urnings at (0.5 x urnsize)
# for accuracy
u=rep(k1/2,N) # persons
v=rep(k2/2,K) # items
# for pseudoRT
u2=rep(k12/2,N) # persons 
v2=rep(k22/2,K) # items

# ancillary objects for paired update
queue<-queue2<-rep(0,K)
Upd=t(matrix(c(0,1,1,0),nrow=2))
n_options=c(0,1,2)
Score=c(-2,2)
n_scores=2
Upd2=t(matrix(c(0,1,1,0),nrow=2))
n_options2=c(0,1,2)
Score2=c(-1,1)
n_scores2=2
LL<-LL2<-rep(0,n_scores*K)
LLsum<-LLsum2<-rep(0,n_scores)

#void urnings_separate_RT_real_data
#(int*nresp,int*P,int*I,int*K,int*X1,int*X2,int*Y1,int*Y2,int*n_obs_p,int*n_obs_i,int*u,int*v,int*u2,int*v2,
#int*Ubefore,int*Vbefore,int*U2before,int*V2before,int*Uafter,int*Vafter,int*Vafter2,int*U2after,int*V2after,int*V2after2,
#int*k1,int*k12,int*k2,int*k22,double*cumsum,
#int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,
#int*Score2,int*n_scores2,int*n_options2,int*Upd2,int*queue2,int*LL2,int*LLsum2,int*N_obs_p,int*N_obs_i){

tmp<-.C("urnings_separate_RT2_real_data",
	as.integer(nresp),
	as.integer(P),
	as.integer(I),
	as.integer(K),	
	as.integer(X1),	
	as.integer(X2),
	as.integer(long), # to save simulated accuracy
	as.integer(long), # to save simulated pseudoRT	 
	as.integer(rep(0,N)), # to save the number of observations per person
	as.integer(rep(0,K)), # to save the number of observations per item	
	as.integer(u), # acc urning persons
	as.integer(v), # acc urning items
	as.integer(u2),# pseudoRT urning persons
	as.integer(v2),# psuedoRT urning items
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(k1),
	as.integer(k12),
	as.integer(k2),
	as.integer(k22),
	as.double(rep(0,K+1)),
	as.integer(Score),
	as.integer(n_scores),
	as.integer(n_options),
	as.integer(Upd),
	as.integer(queue),
	as.integer(LL),
	as.integer(LLsum),
	as.integer(Score2),
	as.integer(n_scores2),
	as.integer(n_options2),
	as.integer(Upd2),
	as.integer(queue2),
	as.integer(LL2),
	as.integer(LLsum2),	
	as.integer(long),
	as.integer(long))

n_obs_p<-tmp[[44]] # for every response how many responses the person gave so far
n_obs_i<-tmp[[45]] # for every response how many responses the item had so far

simulated1<-tmp[[7]] # simulated accuracy
simulated2<-tmp[[8]] # simulated pseudoRT

Ubefore<-tmp[[15]] # acc urning of the person before the update
Vbefore<-tmp[[16]] # acc urning of the item before the update


U2before<-tmp[[17]]# pseudoRT urning of the person before the update
V2before<-tmp[[18]]# pseudoRT urning of the item before the update

# state of the urning after putting the balls in but before removing them
U=Ubefore+X1*2 
V=Vbefore+(1-X1)*2
U2=U2before+X2
V2=V2before+(1-X2)

# evaluating fit to the pseudoRT data
# for every combination of the person and item ratings compute the number of observations with this combination, the proportion of observed pseudoRT=1 and the proportion of simulated pseudoRT=1. The responses with less than 20 responses of the person and the item so far are not included

# this needs to be corrected 0:11
h=c(0:10)
Observed2<-Expected2<-N_obs2<-matrix(ncol=11,nrow=11)
for(i in 1:11){
	for(j in 1:11){
		N_obs2[i,j]=sum((n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j]))	
		Observed2[i,j]=mean(X2[(n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j])])
		Expected2[i,j]=mean(simulated2[(n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j])])
	}
}


# evaluating fit to the accuracy data
# this needs to be corrected to (0:11)*2
h=c(0:10)*2
Observed<-Expected<-N_obs<-matrix(ncol=11,nrow=11)
for(i in 1:11){
	for(j in 1:11){
		N_obs[i,j]=sum((n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j]))	
		Observed[i,j]=mean(X1[(n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j])])
		Expected[i,j]=mean(simulated1[(n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j])])
	}
}

# plotting the results (should be on the diagonal line)
par(mfrow=c(1,2))
plot(c(Observed),c(Expected),cex=N_obs/5000)
abline(a=0,b=1)
plot(c(Observed2),c(Expected2),cex=N_obs/5000)
abline(a=0,b=1)

# Evaluating prediction accuracy

# compute combined urning value (sum of the acc urning and pseudoRT urning) divided by the combined urn size
U=(Ubefore+U2before)/30
V=(Vbefore+V2before)/30

twoPL=function(theta,delta,a){
  theta=log(theta/(1-theta))
  delta=log(delta/(1-delta))
  P=1/(1+exp(-a*(theta-delta)))
  P[delta==theta]=0.5
  return(P)
}

# Predicted response accuracy (probability correct given the urnings before the response, either only based on accuracy or the sum of two urnings)
Predicted2=twoPL(U,V,2)
Predicted=twoPL(Ubefore/20,Vbefore/20,2)

# MSE for accuracy,  responses with less than 20 responses of the person and the item so far are not included
# only accuracy urnings
mean((X1[(n_obs_p>19)&(n_obs_i>19)]-Predicted[(n_obs_p>19)&(n_obs_i>19)])^2)
#[1] 0.1958685
# acc urnings + pseudoRT urnings
mean((X1[(n_obs_p>19)&(n_obs_i>19)]-Predicted2[(n_obs_p>19)&(n_obs_i>19)])^2)
#[1] 0.1905444

# Below the same analysis is repeated but for the second way of defining pseudoRT

k1=20
k2=20
k12=10
k22=10

nresp=nrow(data2)
P=data2$person_id
I=data2$item_id
N=length(persons)
K=length(items)
X1=data2$accuracy
X2=data2$Y2
long=rep(0,nresp)

u=rep(k1/2,N)
v=rep(k2/2,K)
u2=rep(k12/2,N)
v2=rep(k22/2,K)

queue<-queue2<-rep(0,K)

Upd=t(matrix(c(0,1,1,0),nrow=2))
n_options=c(0,1,2)
Score=c(-2,2)
n_scores=2

Upd2=t(matrix(c(0,1,1,0),nrow=2))
n_options2=c(0,1,2)
Score2=c(-1,1)
n_scores2=2

LL<-LL2<-rep(0,n_scores*K)
LLsum<-LLsum2<-rep(0,n_scores)

#void urnings_separate_RT_real_data
#(int*nresp,int*P,int*I,int*K,int*X1,int*X2,int*Y1,int*Y2,int*n_obs_p,int*n_obs_i,int*u,int*v,int*u2,int*v2,
#int*Ubefore,int*Vbefore,int*U2before,int*V2before,int*Uafter,int*Vafter,int*Vafter2,int*U2after,int*V2after,int*V2after2,
#int*k1,int*k12,int*k2,int*k22,double*cumsum,
#int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,
#int*Score2,int*n_scores2,int*n_options2,int*Upd2,int*queue2,int*LL2,int*LLsum2,int*N_obs_p,int*N_obs_i){

tmp<-.C("urnings_separate_RT2_real_data",
as.integer(nresp),as.integer(P),as.integer(I),as.integer(K),as.integer(X1),as.integer(X2),as.integer(long),as.integer(long),as.integer(rep(0,N)),as.integer(rep(0,K)),as.integer(u),as.integer(v),as.integer(u2),as.integer(v2),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(k1),as.integer(k12),as.integer(k2),as.integer(k22),as.double(rep(0,K+1)),as.integer(Score),as.integer(n_scores),as.integer(n_options),as.integer(Upd),as.integer(queue),as.integer(LL),as.integer(LLsum),as.integer(Score2),as.integer(n_scores2),as.integer(n_options2),as.integer(Upd2),as.integer(queue2),as.integer(LL2),as.integer(LLsum2),as.integer(long),as.integer(long))

n_obs_p<-tmp[[44]]
n_obs_i<-tmp[[45]]

simulated1<-tmp[[7]]
simulated2<-tmp[[8]]

Ubefore<-tmp[[15]]
Vbefore<-tmp[[16]]

U2before<-tmp[[17]]
V2before<-tmp[[18]]

U=Ubefore+X1*2
V=Vbefore+(1-X1)*2

U2=U2before+X2
V2=V2before+(1-X2)

h=c(0:10)

Observed2<-Expected2<-N_obs2<-matrix(ncol=11,nrow=11)
for(i in 1:11){
	for(j in 1:11){
		N_obs2[i,j]=sum((n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j]))	
		Observed2[i,j]=mean(X2[(n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j])])
		Expected2[i,j]=mean(simulated2[(n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j])])
	}
}


h=c(0:10)*2

Observed<-Expected<-N_obs<-matrix(ncol=11,nrow=11)
for(i in 1:11){
	for(j in 1:11){
		N_obs[i,j]=sum((n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j]))	
		Observed[i,j]=mean(X1[(n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j])])
		Expected[i,j]=mean(simulated1[(n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j])])
	}
}

par(mfrow=c(1,2))
plot(c(Observed),c(Expected),cex=N_obs/5000)
abline(a=0,b=1)
plot(c(Observed2),c(Expected2),cex=N_obs/5000)
abline(a=0,b=1)

# prediction accuracy

U=(Ubefore+U2before)/30
V=(Vbefore+V2before)/30

twoPL=function(theta,delta,a){
  theta=log(theta/(1-theta))
  delta=log(delta/(1-delta))
  P=1/(1+exp(-a*(theta-delta)))
  P[delta==theta]=0.5
  return(P)
}

Predicted2=twoPL(U,V,2)
Predicted=twoPL(Ubefore/20,Vbefore/20,2)
mean((X1[(n_obs_p>19)&(n_obs_i>19)]-Predicted[(n_obs_p>19)&(n_obs_i>19)])^2)
#[1] 0.1957852
mean((X1[(n_obs_p>19)&(n_obs_i>19)]-Predicted2[(n_obs_p>19)&(n_obs_i>19)])^2)
#[1] 0.1902013


# And the same analysis again with the third way of defining pseudoRT
k1=20
k2=20
k12=10
k22=10

nresp=nrow(data2)
P=data2$person_id
I=data2$item_id
N=length(persons)
K=length(items)
X1=data2$accuracy
X2=data2$Y3
long=rep(0,nresp)

u=rep(k1/2,N)
v=rep(k2/2,K)
u2=rep(k12/2,N)
v2=rep(k22/2,K)

queue<-queue2<-rep(0,K)

Upd=t(matrix(c(0,1,1,0),nrow=2))
n_options=c(0,1,2)
Score=c(-2,2)
n_scores=2

Upd2=t(matrix(c(0,1,1,0),nrow=2))
n_options2=c(0,1,2)
Score2=c(-1,1)
n_scores2=2

LL<-LL2<-rep(0,n_scores*K)
LLsum<-LLsum2<-rep(0,n_scores)

#void urnings_separate_RT_real_data
#(int*nresp,int*P,int*I,int*K,int*X1,int*X2,int*Y1,int*Y2,int*n_obs_p,int*n_obs_i,int*u,int*v,int*u2,int*v2,
#int*Ubefore,int*Vbefore,int*U2before,int*V2before,int*Uafter,int*Vafter,int*Vafter2,int*U2after,int*V2after,int*V2after2,
#int*k1,int*k12,int*k2,int*k22,double*cumsum,
#int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,
#int*Score2,int*n_scores2,int*n_options2,int*Upd2,int*queue2,int*LL2,int*LLsum2,int*N_obs_p,int*N_obs_i){

tmp<-.C("urnings_separate_RT2_real_data",
as.integer(nresp),as.integer(P),as.integer(I),as.integer(K),as.integer(X1),as.integer(X2),as.integer(long),as.integer(long),as.integer(rep(0,N)),as.integer(rep(0,K)),as.integer(u),as.integer(v),as.integer(u2),as.integer(v2),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(long),as.integer(k1),as.integer(k12),as.integer(k2),as.integer(k22),as.double(rep(0,K+1)),as.integer(Score),as.integer(n_scores),as.integer(n_options),as.integer(Upd),as.integer(queue),as.integer(LL),as.integer(LLsum),as.integer(Score2),as.integer(n_scores2),as.integer(n_options2),as.integer(Upd2),as.integer(queue2),as.integer(LL2),as.integer(LLsum2),as.integer(long),as.integer(long))

n_obs_p<-tmp[[44]]
n_obs_i<-tmp[[45]]

simulated1<-tmp[[7]]
simulated2<-tmp[[8]]

Ubefore<-tmp[[15]]
Vbefore<-tmp[[16]]

U2before<-tmp[[17]]
V2before<-tmp[[18]]

U=Ubefore+X1*2
V=Vbefore+(1-X1)*2

U2=U2before+X2
V2=V2before+(1-X2)

h=c(0:10)

Observed2<-Expected2<-N_obs2<-matrix(ncol=11,nrow=11)
for(i in 1:11){
	for(j in 1:11){
		N_obs2[i,j]=sum((n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j]))	
		Observed2[i,j]=mean(X2[(n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j])])
		Expected2[i,j]=mean(simulated2[(n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j])])
	}
}



h=c(0:10)*2

Observed<-Expected<-N_obs<-matrix(ncol=11,nrow=11)
for(i in 1:11){
	for(j in 1:11){
		N_obs[i,j]=sum((n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j]))	
		Observed[i,j]=mean(X1[(n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j])])
		Expected[i,j]=mean(simulated1[(n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j])])
	}
}

par(mfrow=c(1,2))
plot(c(Observed),c(Expected),cex=N_obs/5000)
abline(a=0,b=1)
plot(c(Observed2),c(Expected2),cex=N_obs/5000)
abline(a=0,b=1)

# prediction accuracy

U=(Ubefore+U2before)/30
V=(Vbefore+V2before)/30

twoPL=function(theta,delta,a){
  theta=log(theta/(1-theta))
  delta=log(delta/(1-delta))
  P=1/(1+exp(-a*(theta-delta)))
  P[delta==theta]=0.5
  return(P)
}

Predicted2=twoPL(U,V,2)
Predicted=twoPL(Ubefore/20,Vbefore/20,2)
mean((X1[(n_obs_p>19)&(n_obs_i>19)]-Predicted[(n_obs_p>19)&(n_obs_i>19)])^2)
#[1] 0.1958778
mean((X1[(n_obs_p>19)&(n_obs_i>19)]-Predicted2[(n_obs_p>19)&(n_obs_i>19)])^2)
#[1] 0.1867542

# The third definition gives the biggest improvement of prediction accuracy. Further analysis is done for it.

# Looking at the development of the ratings of the most active player
# 5453 most active player
Uafter<-tmp[[19]]
U2after<-tmp[[22]]


Uafter[P==2363]->u_2363
U2after[P==2363]->u2_2363

plot(u_2363)

par(mfrow=c(1,1))
plot(u_2363/20,ylim=c(0,1),type='l')
lines((u_2363+u2_2363)/30,col=2)


# Analysis with a single urn and multiple updates (for Y3)

# urn sizes
k1=30
k2=30

nresp=nrow(data2)
P=data2$person_id
I=data2$item_id
N=length(persons)
K=length(items)
X1=data2$accuracy
X2=data2$Y3
long=rep(0,nresp)

u=rep(k1/2,N)
v=rep(k2/2,K)

queue<-rep(0,K)

Upd=list()

Upd[[1]]=list( 3,
   c(2,1),
   c(1,1,1))

Upd[[2]]=list( 2,
   c(1,1))
   
Upd[[3]]=list(1)

for(j in 4:6){
	Upd[[j]]=list()
	for(k in 1:length(Upd[[7-j]])){
		Upd[[j]][[k]]=-Upd[[7-j]][[k]]
	}
}

Updates=Upd

scores=c(-3,-2,-1,1,2,3)

Upd=list()
for(j in 1:6){
	Upd[[j]]=matrix(0,ncol=6,nrow=length(Updates[[j]]))
      for(k in 1:nrow(Upd[[j]])){
		for(i in 1:6){
			Upd[[j]][k,i]=sum(Updates[[j]][[k]]==scores[i])
		}	
	}
}


n_options<-c(0,cumsum(sapply(Upd,nrow)))

Upd_matrix=matrix(ncol=6,nrow=sum(sapply(Upd,nrow)))

for(i in 1:6){
	Upd_matrix[(n_options[i]+1):n_options[i+1],]=Upd[[i]]
}

Upd<-t(Upd_matrix)

# possible values of the item update
Score=c(-3:-1,1:3)
n_scores=6

LL<-rep(0,n_scores*K)
LLsum<-rep(0,n_scores)


#void urnings_combined_RT_real_data
#(int*nresp,int*P,int*I,int*K,int*X1,int*X2,int*Y1,int*Y2,int*n_obs_p,int*n_obs_i,int*u,int*v,int*Ubefore,int*Vbefore,int*U_middle,int*V_middle,int*Uafter,int*Vafter,int*Vafter2,int*k1,int*k2,double*cumsum,int*Score,int*n_scores,int*n_options,int*Upd,int*queue,int*LL,int*LLsum,int*N_obs_p,int*N_obs_i)


tmp<-.C("urnings_combined_RT_real_data",
	as.integer(nresp),
	as.integer(P),
	as.integer(I),
	as.integer(K),
	as.integer(X1),
	as.integer(X2),
	as.integer(long),
	as.integer(long),
	as.integer(rep(0,N)),
	as.integer(rep(0,K)),
	as.integer(u),
	as.integer(v),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(long),
	as.integer(k1),
	as.integer(k2),
	as.double(rep(0,K+1)),
	as.integer(Score),
	as.integer(n_scores),
	as.integer(n_options),
	as.integer(Upd),
	as.integer(queue),
	as.integer(LL),
	as.integer(LLsum),
	as.integer(long),
	as.integer(long))

n_obs_p<-tmp[[30]]
n_obs_i<-tmp[[31]]

simulated1<-tmp[[7]]
simulated2<-tmp[[8]]

# urnings before the update
Ubefore<-tmp[[13]]
Vbefore<-tmp[[14]]

# urnings after accuracy update but before the RT update
Umiddle<-tmp[[15]]
Vmiddle<-tmp[[16]]

# to evaluate fit to pseudoRT data we add X2 to the 'middle' urnings
U2=Umiddle+X2
V2=Vmiddle+(1-X2)

# This needs to be corrected. It should be 0:31
h=c(0:30)

Observed2<-Expected2<-N_obs2<-matrix(ncol=length(h),nrow=length(h))
for(i in 1:length(h)){
	for(j in 1:length(h)){
		N_obs2[i,j]=sum((n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j]))	
		Observed2[i,j]=mean(X2[(n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j])])
		Expected2[i,j]=mean(simulated2[(n_obs_p>19)&(n_obs_i>19)&(U2==h[i])&(V2==h[j])])
	}
}

# to evaluate fit to the accuracy data we add X1*2 to the urning before the update
U=Ubefore+X1*2
V=Vbefore+(1-X1)*2

# this needs to be corrected to 0:32

h=c(0:30)

Observed<-Expected<-N_obs<-matrix(ncol=length(h),nrow=length(h))
for(i in 1:length(h)){
	for(j in 1:length(h)){
		N_obs[i,j]=sum((n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j]))	
		Observed[i,j]=mean(X1[(n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j])])
		Expected[i,j]=mean(simulated1[(n_obs_p>19)&(n_obs_i>19)&(U==h[i])&(V==h[j])])
	}
}

par(mfrow=c(1,2))
plot(c(Observed),c(Expected),cex=N_obs/1000)
abline(a=0,b=1)
plot(c(Observed2),c(Expected2),cex=N_obs/1000)
abline(a=0,b=1)

Predicted3=twoPL(Ubefore/30,Vbefore/30,2)
mean((X1[(n_obs_p>19)&(n_obs_i>19)]-Predicted3[(n_obs_p>19)&(n_obs_i>19)])^2)
#[1] 0.1863147

sum(log(Predicted[(X1==1)&(n_obs_p>19)&(n_obs_i>19)]))+sum(log(1-Predicted[(X1==0)&(n_obs_p>19)&(n_obs_i>19)]))