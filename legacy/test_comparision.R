# R code for separate urnings for RA and RT where adaptivity is based on both ra and rt urnings

dyn.load("CRT.so")

# adaptivity is on
adaptive=1
# paired update is on
paired=1
# items are not starting cold, and the majority of people is also not cold (only the last 8 people onto which we zoom in are new)
cold=0

# number response per replication
nIt=500
# number of replications (the same system is repeated nrep times)
nrep=100

W=4
mu_theta=0
sd_theta=1/2 

N=1000 # sample size 

# sample true values of the persons parameters from multivariate normal
theta=qnorm(seq(1/(N+1),N/(N+1),by=1/(N+1)),mu_theta,sd_theta)
# theta is actually theta/4, so in the initial parametrisation theta~N(0,2)
theta=sample(theta,N) # suffling the order to avoid order effects
theta=c(theta,-1,-3/4,-1/2,-1/4,1/4,1/2,3/4,1)
N=length(theta) # updated sample size after 8 more people were added
new=c(501:N)

K=200 # number of items

# true item difficulties
delta=qnorm(seq(1/(K+1),K/(K+1),length=K),0,1/W)*2

#Directly specify the values for each iteration separately such that change can also be included (in the c code iteration-specific values are used to generate data)
Theta=matrix(ncol=nIt,nrow=N)
Delta=matrix(ncol=nIt,nrow=K)

for(j in 1:nIt){
	Theta[,j]=theta
	Delta[,j]=delta
}

# urn sizes for the persons and items for computing the matrix with selection probabilities
k1=70
k2=350

# normal kernel for item selection is formulated based on probabilities, not logits
M=0.7
SD=0.1
Prob2=matrix(1,nrow=k1+1,ncol=k2+1)
for(i in 0:(k1)){
	for(j in 0:(k2)){
		l=log((i+1)/(k1-i+1))-log((j+1)/(k2-j+1))
		Prob2[i+1,j+1]=dnorm(1/(1+exp(-W*l)),M,SD)
	}
}

# list to save the results of different methods
Res=list()

# conditions are actually not conditions, but methods, 0 - use only RA with weight 4, 2 - use RT but with separate urns, 1 - use RT with a single urn
for(condition in c(0,1,2)){

# to save the results of the last 8 people
UU<-UU2<-UU3<-array(dim=c(nrep,length(new),nIt))

# urn sizes for when separate urns are used
if(condition==2){
k1=40
k2=200
k12=k1/2
k22=k2/2
k13=k1/4
k23=k2/4
}

# urn sizes for when only RA is used, or a single RA-RT urn is used
if(condition<2){
k1=70
k2=350
k12=0
k22=0
k13=0
k23=0
}

# these are possible updates based on RA alone
if((condition==2)|(condition==0)){
	Upd=t(matrix(c(0,1,1,0),nrow=2))
	n_options=c(0,1,2)
	Score=c(-1,1)*W
	n_scores=2	
}

# these are possible updates when a single RA-RT urn is used
if(condition==1){
Upd=list()
Upd[[1]]=list(7,
   c(6,1),
   c(5,2),
   c(5,1,1),
   c(4,3),
   c(4,2,1),
   c(4,1,1,1),
   c(3,3,1),
   c(3,2,2),
   c(3,2,1,1),
   c(3,1,1,1,1),
   c(2,2,2,1),
   c(2,2,1,1,1),
   c(2,1,1,1,1,1),
   c(1,1,1,1,1,1,1))

Upd[[2]]=list(6,
   c(5,1),
   c(4,2), 
   c(4,1,1),
   c(3,3),
   c(3,2,1),
   c(3,1,1,1),
   c(2,2,2),
   c(2,2,1,1),
   c(2,1,1,1,1),
   c(1,1,1,1,1,1))

Upd[[3]]=list( 5,
   c(4,1),
   c(3,2),
   c(3,1,1),
   c(2,2,1),
   c(2,1,1,1),
   c(1,1,1,1,1))

Upd[[4]]=list(c(-4,7,1),
   c(-4,6,2),
   c(-4,5,3),
   c(4),
   c(3,1), 
   c(2,2),
   c(2,1,1),
   c(1,1,1,1))

Upd[[5]]=list( 3,
   c(2,1),
   c(1,1,1))

Upd[[6]]=list( 2,
   c(1,1))
   
Upd[[7]]=list(1)

for(j in 8:14){
	Upd[[j]]=list()
	for(k in 1:length(Upd[[15-j]])){
		Upd[[j]][[k]]=-Upd[[15-j]][[k]]
	}
}

Updates=Upd

scores=c(-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7)

Upd=list()
for(j in 1:14){
	Upd[[j]]=matrix(0,ncol=14,nrow=length(Updates[[j]]))
      for(k in 1:nrow(Upd[[j]])){
		for(i in 1:14){
			Upd[[j]][k,i]=sum(Updates[[j]][[k]]==scores[i])
		}	
	}
}


n_options<-c(0,cumsum(sapply(Upd,nrow)))

Upd_matrix=matrix(ncol=14,nrow=sum(sapply(Upd,nrow)))

for(i in 1:14){
	Upd_matrix[(n_options[i]+1):n_options[i+1],]=Upd[[i]]
}

Upd<-t(Upd_matrix)

# possible values of the item update
Score=c(-7:-1,1:7)
n_scores=14

}

# possible updates for the 2nd and 3rd urns when separate urns are used
Upd2=t(matrix(c(0,1,1,0),nrow=2))
n_options2=c(0,1,2)
Score2=c(-2,2)
n_scores2=2

Upd3=t(matrix(c(0,1,1,0),nrow=2))
n_options3=c(0,1,2)
Score3=c(-1,1)
n_scores3=2

for(rep in 1:nrep){
T1=Sys.time()

# starting values 
u=rep(k1/2,N)
v=rep(k2/2,K)

# the items and the not new persons are starting from invariant distributions
if(cold==0){
	v[1:(K/2)]=rbinom(K/2,k2,1/(1+exp(-delta[1:(K/2)])))
	v[K:(K/2+1)]=k2-v[1:(K/2)]
	u[-new]=rbinom(N-length(new),k1,1/(1+exp(-theta[-new])))
}

u2=rep(k12/2,N)
v2=rep(k22/2,K)

u3=rep(k13/2,N)
v3=rep(k23/2,K)

if(cold==0){
	if(condition==2){
		v2[1:(K/2)]=rbinom(K/2,k22,1/(1+exp(-delta[1:(K/2)])))
		v2[K:(K/2+1)]=k22-v2[1:(K/2)]
		u2[-new]=rbinom(N-length(new),k12,1/(1+exp(-theta[-new])))
		v3[1:(K/2)]=rbinom(K/2,k23,1/(1+exp(-delta[1:(K/2)])))
		v3[K:(K/2+1)]=k23-v3[1:(K/2)]
		u3[-new]=rbinom(N-length(new),k13,1/(1+exp(-theta[-new])))	
	}
}

# objects to save the sampled values after each iteration
U<-U2<-U3<-matrix(0,ncol=nIt,nrow=N)
V<-V2<-V3<-matrix(0,ncol=nIt,nrow=K)

# item queue. In case of separate urns, they are separate for the three urns
queue<-queue2<-queue3<-rep(0,K)

# objects that are needed for the paired update
LL<-LL2<-LL3<-rep(0,n_scores*K)
LLsum<-LLsum2<-LLsum3<-rep(0,n_scores)

if(condition==0){
tmp<-.C("urnings_simple",
        as.integer(adaptive), #indicator of the use of adaptive item selection
        as.integer(paired), #indicator for the inclusion of paired update
        as.integer(u), #student starting values
        as.integer(v), #item starting values
        as.double(Theta), #nplayers x niteration matrix of student true values
        as.double(Delta), #n_items x niteration matrix of item true values
        as.integer(N), #number of students
        as.integer(K), #number of items
        as.integer(nIt), #number of games
        as.integer(U),as.integer(V), #containers for the results?
        as.integer(k1),as.integer(k2), #urn sizes for students and items
        as.double(Prob2), #normal kernel matrix
        as.double(rep(0,K+1)), #no idea, but probably a helper for calculating the normalising constant
        as.integer(W), #weights
        as.integer(Score), #possible updates
        as.integer(n_scores), #number of possible updates other than 0
        as.integer(n_options), #number of possible updates as a vector? 
        as.integer(Upd), #helper for the paired update I guess
        as.integer(queue), #queue for the paired update
        as.integer(LL), #helpers for the paired update again
        as.integer(LLsum))#and again


U=matrix(tmp[[10]],ncol=nIt) #output is every input value for some 
#V=matrix(tmp[[11]],ncol=nIt)
}

if(condition==1){

tmp<-.C("urnings_combined_RT",
        as.integer(adaptive),
        as.integer(paired),
        as.integer(u),as.integer(v),
        as.double(Theta),as.double(Delta),
        as.integer(N),as.integer(K),as.integer(nIt),
        as.integer(U),as.integer(V),
        as.integer(k1),as.integer(k2),
        as.double(Prob2),as.double(rep(0,K+1)),as.integer(W),
        as.integer(Score),as.integer(n_scores),
        as.integer(n_options),as.integer(Upd),as.integer(queue),
        as.integer(LL),as.integer(LLsum))
U = urnings_combined_RT(student_starting = u,
                                   item_starting = v,
                                   Theta = Theta,
                                   Delta = Delta,
                                   n_students = N,
                                   n_items = K,
                                   n_games = nIt,
                                   student_urn_size = k1,
                                   item_urn_size = k2,
                                   weight = W,
                                   adaptive = adaptive,
                                   m_adapt = M,
                                   sd_adapt = SD,
                                   paired = paired,
                                   returns = "simple",
                                   OS = "MAC")


U=matrix(tmp[[10]],ncol=nIt)
#V=matrix(tmp[[11]],ncol=nIt)



}

if(condition==2){

U_s = urnings_separate_RT(student_starting = rbind(u,u2,u3),
                                 item_starting = rbind(v,v2,v3),
                                 Theta = Theta,
                                 Delta = Delta,
                                 n_students = N,
                                 n_items = K,
                                 n_games = nIt,
                                 student_urn_size = c(k1,k12,k13),
                                 item_urn_size = c(k2, k22, k23),
                                 weight = W,
                                 adaptive = adaptive,
                                 m_adapt = M,
                                 sd_adapt = SD,
                                 paired = paired,
                                 returns = "simple",
                                 OS = "MAC")

tmp<-.C("urnings_separate_RT2",
        as.integer(adaptive),
        as.integer(paired),
        as.integer(u),
        as.integer(v),
        as.integer(u2),
        as.integer(v2),
        as.integer(u3),
        as.integer(v3),
        as.double(Theta),
        as.double(Delta),
        as.integer(N),
        as.integer(K),
        as.integer(nIt),
        as.integer(U),
        as.integer(V),
        as.integer(U2),
        as.integer(V2),
        as.integer(U3),
        as.integer(V3),
        as.integer(k1),
        as.integer(k12),
        as.integer(k13),
        as.integer(k2),
        as.integer(k22),
        as.integer(k23),
        as.double(Prob2),
        as.double(rep(0,K+1)),
        as.integer(W),
        as.integer(Score),as.integer(n_scores),as.integer(n_options),as.integer(Upd),as.integer(queue),
as.integer(LL),as.integer(LLsum),
as.integer(Score2),as.integer(n_scores2),as.integer(n_options2),as.integer(Upd2),as.integer(queue2),
as.integer(LL2),as.integer(LLsum2),
as.integer(Score3),as.integer(n_scores3),as.integer(n_options3),as.integer(Upd3),as.integer(queue3),
as.integer(LL3),as.integer(LLsum3))

U=matrix(tmp[[14]],ncol=nIt)
#V=matrix(tmp[[13]],ncol=nIt)
U2=matrix(tmp[[16]],ncol=nIt)
#V2=matrix(tmp[[15]],ncol=nIt)
U3=matrix(tmp[[18]],ncol=nIt)
}

UU[rep,,]=U[new,]
#VV[rep,,]=V[new,]
UU2[rep,,]=U2[new,]
#VV2[rep,,]=V2[new,]
UU3[rep,,]=U3[new,]

T2=Sys.time()
print(c(rep,T2-T1))

}

Res[[condition+1]]=UU
if(condition==2){
	Res[[condition+1]]=list(UU,UU2,UU3)
}

}

# Below are computing things to compare the three methods based on bias and MSE


Means1=apply(Res[[1]],c(2,3),mean)/70
Means2=apply(Res[[2]],c(2,3),mean)/70
Means3=apply(Res[[3]][[1]]+Res[[3]][[2]]+Res[[3]][[3]],c(2,3),mean)/70
Means4=apply(Res[[3]][[1]],c(2,3),mean)/40
Means5=apply(Res[[3]][[2]],c(2,3),mean)/20
Means6=apply(Res[[3]][[3]],c(2,3),mean)/10

# plotting the average urnings across replications and responses for the 8 new people
j=2
plot(Means1[2,],type='l',ylim=c(0,1),ylab=expression(1/(1+exp(-theta/4))))
for(j in 1:8){
lines(Means1[j,])
lines(Means2[j,],col=2)
lines(Means3[j,],col=3)
abline(h=1/(1+exp(-theta[new][j])))
}

# here everything is transformed based to the original theta scale
TMeans1=log(Means1/(1-Means1))*4
TMeans2=log(Means2/(1-Means2))*4
TMeans3=log(Means3/(1-Means3))*4
plot(TMeans1[2,],type='l',ylim=c(-4,4),ylab=expression(theta),xlab='Number of responses')
for(j in 1:8){
lines(TMeans1[j,])
lines(TMeans2[j,],col=2)
lines(TMeans3[j,],col=3)
abline(h=theta[new][j]*4)
}
legend("topright", legend=c("Urnings", "Multiple trackers", "Multiple updates"),
       col=c("black", "green", "red"), lty=1:2, cex=0.8)


k1=300

par(mfrow=c(1,3))

mse_ra<-mse_combined<-mse_separate<-numeric(nIt)
for(j in 1:nrep){
mse_combined=mse_combined+apply(Res[[2]][j,,],2,FUN=function(X){mean((X/k1-1/(1+exp(-theta[new])))^2)})
mse_separate=mse_separate+apply((Res[[3]][[1]][j,,]+Res[[3]][[2]][j,,]),2,FUN=function(X){mean((X/k1-1/(1+exp(-theta[new])))^2)})
mse_ra=mse_ra+apply(Res[[1]][j,,],2,FUN=function(X){mean((X/k1-1/(1+exp(-theta[new])))^2)})
}
mse_combined=mse_combined/nrep
mse_separate=mse_separate/nrep
mse_ra=mse_ra/nrep

plot(mse_ra,type='l')
lines(mse_separate,col=2)
lines(mse_combined,col=3)

mse1=apply(Res[[3]][[1]][j,,],2,FUN=function(X){mean((X/60-1/(1+exp(-theta[new])))^2)})
mse2=apply(Res[[3]][[2]][j,,],2,FUN=function(X){mean((X/30-1/(1+exp(-theta[new])))^2)})

base=numeric(500)
for(j in 1:500){
u=rbinom(N,k1,1/(1+exp(-theta)))/k1
base[j]=mean((u-1/(1+exp(-theta)))[new]^2)
}
abline(h=mean(base),col=4)

mse_ra<-mse_combined<-mse_separate<-numeric(nIt)
for(j in 1:nrep){
mse_combined=mse_combined+apply(Res[[2]][j,,],2,FUN=function(X){mean((X/k1-1/(1+exp(-theta)))[-new]^2)})
mse_separate=mse_separate+apply((Res[[3]][[1]][j,,]+Res[[3]][[2]][j,,]),2,FUN=function(X){mean((X/k1-1/(1+exp(-theta)))[-new]^2)})
mse_ra=mse_ra+apply(Res[[1]][j,,],2,FUN=function(X){mean((X/k1-1/(1+exp(-theta)))[-new]^2)})
}
mse_combined=mse_combined/nrep
mse_separate=mse_separate/nrep
mse_ra=mse_ra/nrep

plot(mse_ra,type='l')
lines(mse_separate,col=2)
lines(mse_combined,col=3)

base=numeric(500)
for(j in 1:500){
u=rbinom(N,k1,1/(1+exp(-theta)))/k1
base[j]=mean((u-1/(1+exp(-theta)))[-new]^2)
}
abline(h=mean(base),col=4)

bias_ra<-bias_combined<-bias_separate<-numeric(nIt)
for(j in 1:nIt){
bias_combined[j]=mean(abs(colMeans(Res[[2]][,,j])/k1-1/(1+exp(-theta[new]))))
bias_separate[j]=mean(abs(colMeans(Res[[3]][[1]][,,j]+Res[[3]][[2]][,,j])/k1-1/(1+exp(-theta[new]))))
bias_ra[j]=mean(abs(colMeans(Res[[1]][,,j])/k1-1/(1+exp(-theta[new]))))
}

plot(bias_ra,type='l')
lines(bias_separate,col=2)
lines(bias_combined,col=3)

base=numeric(500)
for(j in 1:500){
u=matrix(rbinom(N*nrep,k1,rep(1/(1+exp(-theta)),nrep)),nrow=N)/k1
base[j]=mean(abs(rowMeans(u)-1/(1+exp(-theta)))[new])
}
abline(h=mean(base),col=4)

bias_ra<-bias_combined<-bias_separate<-numeric(nIt)
for(j in 1:nIt){
bias_combined[j]=mean(abs(colMeans(Res[[2]][,,j])/k1-1/(1+exp(-theta)))[-new])
bias_separate[j]=mean(abs(colMeans(Res[[3]][[1]][,,j]+Res[[3]][[2]][,,j])/k1-1/(1+exp(-theta)))[-new])
bias_ra[j]=mean(abs(colMeans(Res[[1]][,,j])/k1-1/(1+exp(-theta)))[-new])
}

plot(bias_ra,type='l')
lines(bias_separate,col=2)
lines(bias_combined,col=3)

base=numeric(500)
for(j in 1:500){
u=matrix(rbinom(N*nrep,k1,rep(1/(1+exp(-theta)),nrep)),nrow=N)/k1
base[j]=mean(abs(rowMeans(u)-1/(1+exp(-theta)))[-new])
}
abline(h=mean(base),col=4)



bias_ra<-bias_combined<-bias_separate<-numeric(nIt)
for(j in 1:nIt){
bias_combined[j]=mean(abs((colMeans(Res[[2]][,,j])/k1-1/(1+exp(-theta[new])))/(1/(1+exp(-theta[new])))))
bias_separate[j]=mean(abs((colMeans(Res[[3]][[1]][,,j]+Res[[3]][[2]][,,j])/k1-1/(1+exp(-theta[new])))/(1/(1+exp(-theta[new])))))
bias_ra[j]=mean(abs((colMeans(Res[[1]][,,j])/k1-1/(1+exp(-theta[new])))/(1/(1+exp(-theta[new])))))
}

plot(bias_ra,type='l')
lines(bias_separate,col=2)
lines(bias_combined,col=3)



cor_ra<-cor_combined<-cor_separate<-numeric(nIt)
for(j in 1:nrep){
cor_combined=cor_combined+apply(Res[[2]][j,,],2,FUN=function(X){cor(X,1/(1+exp(-theta[new])))})
cor_separate=cor_separate+apply((Res[[3]][[1]][j,,]+Res[[3]][[2]][j,,]),2,FUN=function(X){cor(X,1/(1+exp(-theta[new])))})
cor_ra=cor_ra+apply(Res[[1]][j,,],2,FUN=function(X){cor(X,1/(1+exp(-theta[new])))})
}
cor_combined=cor_combined/nrep
cor_separate=cor_separate/nrep
cor_ra=cor_ra/nrep


plot(cor_ra,type='l')
lines(cor_separate,col=2)
lines(cor_combined,col=3)
