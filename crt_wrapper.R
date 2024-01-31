urnings_simple = function(student_starting,
                          item_starting,
                          Theta,
                          Delta,
                          n_students,
                          n_items,
                          n_games,
                          student_urn_size,
                          item_urn_size,
                          weight = 1,
                          adaptive = 0,
                          m_adapt = 0.5,
                          sd_adapt = 0.1,
                          paired = 1,
                          returns = "simple",
                          OS = "MAC"){
################################################################################  
#Parameters
################################################################################
# FUNCTION ARGUMENTS
  # student_starting : (vector<int>) starting values for the student 
  #  item_starting : (vector<int>) starting values for the items 
  #  Theta : matrix[n_students, n_games]<double> matrix of true values
  #  Delta : matrix[n_items, n_games]<double> matrix of true values
  #  n_students: <int> number of students in the system
  #  n_items: <int> number of items in the system
  #  n_games: <int> number of items solved by each student
  #  student_urn_size: <int> the size of the student urns
  #  item_urn_size: <int> the size of the item urns
  #  weight = 1: <int> size or the stake of an update
  #  adaptive = 0: <int>[0,1] indicator 1 = adaptive, 0 = uniform item selection
  #  m_adapt = 0.5: <double> mu parameter of the Normal Kernel Method (probability)
  #  sd_adapt = 0.1: <double> sigma parameter of the Normal Kernel Method (probability)
  #  paired = 1: <int>[0,1] indicator 1 = paired update, 0 = no paired update
  #  returns = "simple": <char>["simple", "total"], "simple" = return student urnings,
  #                                                 "total" = return the whole list
  #  OS = "MAC" : <char>["MAC", "WINDOWS", "LINUX"]: selects the right file for the OS
  
#RETURN:
  #if simple: matrix[n_students, n_games]<integer>
  #if total: list
  
################################################################################  
#loading the right compiled version of the C file 
################################################################################
  switch (OS,
    "MAC" = dyn.load("CRT.so"),
    "WINDOWS" = dyn.load("CRT.dll"),
    "LINUX" = dyn.load("CRT.o")
  )
################################################################################  
#create the normal kernel values for the possible combinations
################################################################################
  Prob2=matrix(1,nrow=student_urn_size+1,ncol=item_urn_size+1)
  if(adaptive == 1){ 
    for(i in 0:(student_urn_size)){
      for(j in 0:(item_urn_size)){
        l=log((i+1)/(student_urn_size-i+1))-log((j+1)/(item_urn_size-j+1))
        Prob2[i+1,j+1]=dnorm(1/(1+exp(-weight*l)),m_adapt,sd_adapt)
      }
    }  
  }
################################################################################
#create helpers for the paired update
################################################################################
  Upd=t(matrix(c(0,1,1,0),nrow=2))
  n_options=c(0,1,2)
  Score=c(-1,1)*weight
  n_scores=2
  queue<rep(0,n_items)
  LL=rep(0,n_scores*K)
  LLsum=rep(0,n_scores)

################################################################################
#create helpers for the paired update
################################################################################
student_container=matrix(0,ncol=n_games,nrow=n_students)
item_container=matrix(0,ncol=n_games,nrow=n_items)
################################################################################
#call the C function to perform the simulation
################################################################################
  tmp<-.C("urnings_simple",
          as.integer(adaptive), #indicator of the use of adaptive item selection
          as.integer(paired), #indicator for the inclusion of paired update
          as.integer(student_starting), #student starting values
          as.integer(item_starting), #item starting values
          as.double(Theta), #nplayers x niteration matrix of student true values
          as.double(Delta), #n_items x niteration matrix of item true values
          as.integer(n_students), #number of students
          as.integer(n_items), #number of items
          as.integer(n_games), #number of games
          as.integer(student_container), #container for students
          as.integer(item_container), #containers for items
          as.integer(student_urn_size), # urn size for students
          as.integer(item_urn_size), #urn sizes for items
          as.double(Prob2), #normal kernel matrix
          as.double(rep(0,n_items+1)), #no idea, but probably a helper for calculating the normalising constant
          as.integer(weight), #weights
          as.integer(Score), #possible updates
          as.integer(n_scores), #number of possible updates other than 0
          as.integer(n_options), #number of possible updates as a vector? 
          as.integer(Upd), #helper for the paired update I guess
          as.integer(queue), #queue for the paired update
          as.integer(LL), #helpers for the paired update again
          as.integer(LLsum))#and again

################################################################################
#returning the results
################################################################################
  if(returns == "simple"){
    U=matrix(tmp[[10]],ncol=n_games)
    return(U)
  } else {
    return(tmp)
  }
}


urnings_combined_RT = function(student_starting,
                               item_starting,
                               Theta,
                               Delta,
                               n_students,
                               n_items,
                               n_games,
                               student_urn_size,
                               item_urn_size,
                               weight = 1,
                               adaptive = 0,
                               m_adapt = 0.5,
                               sd_adapt = 0.1,
                               paired = 1,
                               returns = "simple",
                               OS = "MAC"){
################################################################################  
#Parameters
################################################################################
  # FUNCTION ARGUMENTS
  # student_starting : (vector<int>) starting values for the student 
  #  item_starting : (vector<int>) starting values for the items 
  #  Theta : matrix[n_students, n_games]<double> matrix of true values
  #  Delta : matrix[n_items, n_games]<double> matrix of true values
  #  n_students: <int> number of students in the system
  #  n_items: <int> number of items in the system
  #  n_games: <int> number of items solved by each student
  #  student_urn_size: <int> the size of the student urns
  #  item_urn_size: <int> the size of the item urns
  #  weight = 1: <int> size or the stake of an update
  #  adaptive = 0: <int>[0,1] indicator 1 = adaptive, 0 = uniform item selection
  #  m_adapt = 0.5: <double> mu parameter of the Normal Kernel Method (probability)
  #  sd_adapt = 0.1: <double> sigma parameter of the Normal Kernel Method (probability)
  #  paired = 1: <int>[0,1] indicator 1 = paired update, 0 = no paired update
  #  returns = "simple": <char>["simple", "total"], "simple" = return student urnings,
  #                                                 "total" = return the whole list
  #  OS = "MAC" : <char>["MAC", "WINDOWS", "LINUX"]: selects the right file for the OS
  
  #RETURN:
  #if simple: matrix[n_students, n_games]<integer>
  #if total: list 
  
################################################################################  
#loading the right compiled version of the C file 
################################################################################
  switch (OS,
          "MAC" = dyn.load("CRT.so"),
          "WINDOWS" = dyn.load("CRT.dll"),
          "LINUX" = dyn.load("CRT.o")
  )
################################################################################  
#create the normal kernel values for the possible combinations
################################################################################
  Prob2=matrix(1,nrow=student_urn_size+1,ncol=item_urn_size+1)
  if(adaptive == 1){ 
    for(i in 0:(student_urn_size)){
      for(j in 0:(item_urn_size)){
        l=log((i+1)/(student_urn_size-i+1))-log((j+1)/(item_urn_size-j+1))
        Prob2[i+1,j+1]=dnorm(1/(1+exp(-weight*l)),m_adapt,sd_adapt)
      }
    }  
  }
  
################################################################################
#create helpers for the paired update 
  #TODO: Generalise this ASAP because it is hard to look at
################################################################################
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
    queue<rep(0,n_items)
    LL=rep(0,n_scores*n_items)
    LLsum=rep(0,n_scores)
    
################################################################################
#create containers for the results
################################################################################
    student_container=matrix(0,ncol=n_games,nrow=n_students)
    item_container=matrix(0,ncol=n_games,nrow=n_items)

################################################################################
#call the C function to perform the simulation
################################################################################  
    tmp<-.C("urnings_combined_RT",
            as.integer(adaptive),
            as.integer(paired),
            as.integer(student_starting),
            as.integer(item_starting),
            as.double(Theta),
            as.double(Delta),
            as.integer(n_students),
            as.integer(n_items),
            as.integer(n_games),
            as.integer(student_container),
            as.integer(item_container),
            as.integer(student_urn_size),
            as.integer(item_urn_size),
            as.double(Prob2),
            as.double(rep(0,K+1)),
            as.integer(weight),
            as.integer(Score),
            as.integer(n_scores),
            as.integer(n_options),
            as.integer(Upd),
            as.integer(queue),
            as.integer(LL),
            as.integer(LLsum))
################################################################################
#returning the results
################################################################################
    if(returns == "simple"){
      U=matrix(tmp[[10]],ncol=n_games)
      return(U)
    } else {
      return(tmp)
    }
  
}

urnings_separate_RT = function(student_starting,
                               item_starting,
                               Theta,
                               Delta,
                               n_students,
                               n_items,
                               n_games,
                               student_urn_size,
                               item_urn_size,
                               weight = 1,
                               adaptive = 0,
                               m_adapt = 0.5,
                               sd_adapt = 0.1,
                               paired = 1,
                               returns = "simple",
                               OS = "MAC"){
  ################################################################################  
  #Parameters
  ################################################################################
  # FUNCTION ARGUMENTS
  # student_starting : (matrix[3,n_students]<int>) starting values for the student 
  #  item_starting : (matrix[3, n_items]<int>) starting values for the items 
  #  Theta : matrix[n_students, n_games]<double> matrix of true values
  #  Delta : matrix[n_items, n_games]<double> matrix of true values
  #  n_students: <int> number of students in the system
  #  n_items: <int> number of items in the system
  #  n_games: <int> number of items solved by each student
  #  student_urn_size: vector[3]<int> the size of the student urns
  #  item_urn_size: vector[3]<int> the size of the item urns
  #  weight = 1: <int> size or the stake of an update
  #  adaptive = 0: <int>[0,1] indicator 1 = adaptive, 0 = uniform item selection
  #  m_adapt = 0.5: <double> mu parameter of the Normal Kernel Method (probability)
  #  sd_adapt = 0.1: <double> sigma parameter of the Normal Kernel Method (probability)
  #  paired = 1: <int>[0,1] indicator 1 = paired update, 0 = no paired update
  #  returns = "simple": <char>["simple", "total"], "simple" = return student urnings,
  #                                                 "total" = return the whole list
  #  OS = "MAC" : <char>["MAC", "WINDOWS", "LINUX"]: selects the right file for the OS
  
  #RETURN:
  #if simple: matrix[n_students, n_games]<integer>
  #if total: list 
  
  ################################################################################  
  #loading the right compiled version of the C file 
  ################################################################################
  switch(OS,
          "MAC" = dyn.load("CRT.so"),
          "WINDOWS" = dyn.load("CRT.dll"),
          "LINUX" = dyn.load("CRT.o")
  )
  ################################################################################  
  #create the normal kernel values for the possible combinations
  ################################################################################
  Prob2=matrix(1,nrow=student_urn_size[1]+1,ncol=item_urn_size[1]+1)
  if(adaptive == 1){ 
    for(i in 0:(student_urn_size[1])){
      for(j in 0:(item_urn_size[1])){
        l=log((i+1)/(student_urn_size[1]-i+1))-log((j+1)/(item_urn_size[1]-j+1))
        Prob2[i+1,j+1]=dnorm(1/(1+exp(-weight*l)),m_adapt,sd_adapt)
      }
    }  
  }
  
  ################################################################################
  #create helpers for the paired update 
  #TODO: Generalise this ASAP because it is hard to look at
  ################################################################################
  Upd=t(matrix(c(0,1,1,0),nrow=2))
  n_options=c(0,1,2)
  Score=c(-1,1)*W
  n_scores=2	
  
  # possible updates for the 2nd and 3rd urns when separate urns are used
  Upd2=t(matrix(c(0,1,1,0),nrow=2))
  n_options2=c(0,1,2)
  Score2=c(-2,2)
  n_scores2=2
  
  Upd3=t(matrix(c(0,1,1,0),nrow=2))
  n_options3=c(0,1,2)
  Score3=c(-1,1)
  n_scores3=2
  
  queue<-queue2<-queue3<-rep(0,n_items)
  LL<-LL2<-LL3<-rep(0,n_scores*n_items)
  LLsum<-LLsum2<-LLsum3<-rep(0,n_scores)
  ################################################################################
  #create containers for the results
  ################################################################################
  U<-U2<-U3<-matrix(0,ncol=n_games,nrow=n_students)
  V<-V2<-V3<-matrix(0,ncol=n_games,nrow=n_items)
  
  ################################################################################
  #call the C function to perform the simulation
  ################################################################################  
  tmp<-.C("urnings_separate_RT2",
          as.integer(adaptive),
          as.integer(paired),
          as.integer(student_starting[1,]),
          as.integer(item_starting[1,]),
          as.integer(student_starting[2,]),
          as.integer(item_starting[2,]),
          as.integer(student_starting[3,]),
          as.integer(item_starting[3,]),
          as.double(Theta),
          as.double(Delta),
          as.integer(n_students),
          as.integer(n_items),
          as.integer(n_games),
          as.integer(U),
          as.integer(V),
          as.integer(U2),
          as.integer(V2),
          as.integer(U3),
          as.integer(V3),
          as.integer(student_urn_size[1]),
          as.integer(student_urn_size[2]),
          as.integer(student_urn_size[3]),
          as.integer(item_urn_size[1]),
          as.integer(item_urn_size[2]),
          as.integer(item_urn_size[3]),
          as.double(Prob2),
          as.double(rep(0,K+1)),
          as.integer(weight),
          as.integer(Score),as.integer(n_scores),as.integer(n_options),as.integer(Upd),as.integer(queue),
          as.integer(LL),as.integer(LLsum),
          as.integer(Score2),as.integer(n_scores2),as.integer(n_options2),as.integer(Upd2),as.integer(queue2),
          as.integer(LL2),as.integer(LLsum2),
          as.integer(Score3),as.integer(n_scores3),as.integer(n_options3),as.integer(Upd3),as.integer(queue3),
          as.integer(LL3),as.integer(LLsum3))
  
  ################################################################################
  #returning the results
  ################################################################################
  if(returns == "simple"){
    U=matrix(tmp[[14]],ncol=nIt)
    U2=matrix(tmp[[16]],ncol=nIt)
    U3=matrix(tmp[[18]],ncol=nIt)
    res = list(U,U2,U3)
    return(res)
  } else {
    return(tmp)
  }
  
}

