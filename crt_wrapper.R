urnings_simple = function(student_starting,
                          item_starting,
                          Theta,
                          Delta,
                          n_students,
                          n_items,
                          n_games,
                          student_urn_size,
                          item_urn_size,
                          weight = 4,
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
  queue = rep(0,n_items)
  LL=rep(0,n_scores*n_items)
  LLsum=rep(0,n_scores)

################################################################################
#create helpers for the paired update
################################################################################
student_container=matrix(0,ncol=n_games,nrow=n_students)
item_container=matrix(0,ncol=n_games,nrow=n_items)
################################################################################
#call the C function to perform the simulation
################################################################################
  tmp=.C("urnings_simple",
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
                               weight = 4,
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
    
    
    n_options=c(0,cumsum(sapply(Upd,nrow)))
    
    Upd_matrix=matrix(ncol=14,nrow=sum(sapply(Upd,nrow)))
    
    for(i in 1:14){
      Upd_matrix[(n_options[i]+1):n_options[i+1],]=Upd[[i]]
    }
    
    Upd=t(Upd_matrix)
    
    # possible values of the item update
    Score=c(-7:-1,1:7)
    n_scores=14
    queue = rep(0,n_items)
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
    tmp=.C("urnings_combined_RT",
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
            as.double(rep(0,n_items+1)),
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
                               weight = 4,
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
  Prob2=matrix(1,nrow=sum(student_urn_size)+1,ncol=sum(item_urn_size)+1)
  if(adaptive == 1){ 
    for(i in 0:sum(student_urn_size)){
      for(j in 0:sum(item_urn_size)){
        l=log((i+1)/(sum(student_urn_size)-i+1))-log((j+1)/(sum(item_urn_size)-j+1))
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
  Score=c(-1,1)*weight
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
  
  queue=queue2=queue3=rep(0,n_items)
  LL=LL2=LL3=rep(0,n_scores*n_items)
  LLsum=LLsum2=LLsum3=rep(0,n_scores)
  ################################################################################
  #create containers for the results
  ################################################################################
  U=U2=U3=matrix(0,ncol=n_games,nrow=n_students)
  V=V2=V3=matrix(0,ncol=n_games,nrow=n_items)
  
  ################################################################################
  #call the C function to perform the simulation
  ################################################################################  
  tmp=.C("urnings_separate_RT2",
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
          as.double(rep(0,n_items+1)),
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
    U=matrix(tmp[[14]],ncol=n_games)
    U2=matrix(tmp[[16]],ncol=n_games)
    U3=matrix(tmp[[18]],ncol=n_games)
    res = list(U,U2,U3)
    return(res)
  } else {
    return(tmp)
  }
}

urnings_simple_misfit = function(misfit,
                                 ordering,
                          student_starting,
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
  # misfit: matrix[n_student, n_games]<int> true outcomes generated from different model
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
  queue = rep(0,n_items)
  LL=rep(0,n_scores*n_items)
  LLsum=rep(0,n_scores)
  
  ################################################################################
  #create helpers for the paired update
  ################################################################################
  student_container=matrix(0,ncol=n_games,nrow=n_students)
  item_container=matrix(0,ncol=n_games,nrow=n_items)
  ################################################################################
  #call the C function to perform the simulation
  ################################################################################
  tmp=.C("urnings_simple_misfit",
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
          as.integer(LLsum), #and again
          as.integer(misfit),
          as.integer(ordering))#
  
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

urnings_combined_RT_misfit = function(misfit,
                               student_starting,
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
  # misfit: array[n_student, n_games, 3] observed responses from other models for 3 dyadic expansion
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
  
  
  n_options=c(0,cumsum(sapply(Upd,nrow)))
  
  Upd_matrix=matrix(ncol=14,nrow=sum(sapply(Upd,nrow)))
  
  for(i in 1:14){
    Upd_matrix[(n_options[i]+1):n_options[i+1],]=Upd[[i]]
  }
  
  Upd=t(Upd_matrix)
  
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
  tmp=.C("urnings_combined_RT",
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
          as.integer(LLsum),
          as.integer(misfit[,,1]),
          as.integer(misfit[,,2]),
          as.integer(misfit[,,3]))
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


urnings_separate_RT_misfit = function(misfit,
                               student_starting,
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
  # misfit: array[n_student, n_games, 3] observed responses under different models for 3 dyadix exp
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
  
  queue=queue2=queue3=rep(0,n_items)
  LL=LL2=LL3=rep(0,n_scores*n_items)
  LLsum=LLsum2=LLsum3=rep(0,n_scores)
  ################################################################################
  #create containers for the results
  ################################################################################
  U=U2=U3=matrix(0,ncol=n_games,nrow=n_students)
  V=V2=V3=matrix(0,ncol=n_games,nrow=n_items)
  
  ################################################################################
  #call the C function to perform the simulation
  ################################################################################  
  tmp=.C("urnings_separate_RT2",
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
          as.integer(LL3),as.integer(LLsum3),
          as.integer(misfit[,,1]),
          as.integer(misfit[,,2]),
          as.integer(misfit[,,3]))
  
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


urnings_separate_RT_data = function(data,
                                    student_urn_size,
                                    item_urn_size,
                                    pRT_student_urn_size,
                                    pRT_item_urn_size,
                                    fast_indicator,
                                    OS = "MAC"){

  #-----------------------------------------------------------------------------
  # FUNCTION ARGUMENTS
  # data: data frame with the following columns:
    # - student_id: integer, unique identifier for the student (1:n_students)
    # - item_id: integer, unique identifier for the item (1:n_items)
    # - response: integer, 0 or 1, the response of the student (accuracy)
    # - pRT_1: integer 0 for slow 1 for fast response based on the grand median of all response times ("grand_median")
    # - pRT_2: integer 0 for slow 1 for fast response based on the median of the given  ("item_median")
    # - pRT_3: integer 0 for slow 1 for fast response based on the halving point of the 0-90th quantile of the response times of the given item  ("half_90_quantile")
  # student_urn_size: integer, urn size of the accuracy urn for the students
  # item_urn_size: integer, urn size of the accuracy urn for the items
  # pRT_student_urn_size: integer, urn size of the pseudo RT variable (half of the student urn size normally)
  # pRT_item_urn_size: integer, urn size of the pseudo RT variable (half of the item urn size normally)
  # fast_indicator: string, which definition of fast/slow responses is used see for the correct values above
  # OS: string, which operating system is used (for selecting the appropriate compiled C routine)
  
  # RETURN
  # list
    # see comments below 
  
  #-----------------------------------------------------------------------------
  # loading the compiled C routine
  switch(OS,
         "MAC" = dyn.load("CRT_old.so"),
         "WINDOWS" = dyn.load("CRT_old.dll"),
         "LINUX" = dyn.load("CRT_old.o")
  )
  
  #-----------------------------------------------------------------------------
  # setting up the arguments for the C wrapper
  
  #-----------------------------------------------------------------------------
  #disassembling the data to fit the C routine
  
  n_responses = nrow(data) 
  student_id = data$student_id 
  item_id=data$item_id 
  n_students=length(unique(student_id)) 
  n_items=length(unique(item_id)) 
  X1=data$response # accuracy data
  
  X2 = switch(fast_indicator, 
               "grand_median" = data$pRT_1,
               "item_median" = data$pRT_2,
               "half_90_quantile" = data$pRT_3)
  
  #-----------------------------------------------------------------------------
  # setting algorithm settings
  
  # set initial urnings at (0.5 x urnsize)
  # for accuracy
  student_starting = rep(student_urn_size/2,n_students) 
  item_starting =rep(item_urn_size/2,n_items) 
  # for pseudoRT
  pRT_student_starting = rep(pRT_student_urn_size/2,n_students)  
  pRT_item_starting=rep(pRT_item_urn_size/2,n_items) 
  
  #-----------------------------------------------------------------------------
  # auxiliary objects for paired update
  
  queue=queue2=rep(0,n_items) #queue where the items can rest :)
  Upd=t(matrix(c(0,1,1,0),nrow=2)) 
  n_options=c(0,1,2)
  Score=c(-2,2)
  n_scores=2
  Upd2=t(matrix(c(0,1,1,0),nrow=2))
  n_options2=c(0,1,2)
  Score2=c(-1,1)
  n_scores2=2
  LL=LL2=rep(0,n_scores*n_items)
  LLsum=LLsum2=rep(0,n_scores)
  
  #-----------------------------------------------------------------------------
  #containers for the results
  # empty vector with as many elements as there are responses
  #it is used to save all sort of things within the c code
  container_long = rep(0,n_responses) 
  
  tmp=.C("urnings_separate_RT2_real_data",
          as.integer(n_responses), #number of responses
          as.integer(student_id),
          as.integer(item_id),
          as.integer(n_items),	
          as.integer(X1),	
          as.integer(X2),
          as.integer(container_long), # to save simulated accuracy
          as.integer(container_long), # to save simulated pseudoRT	 
          as.integer(rep(0,n_students)), # to save the number of observations per person
          as.integer(rep(0,n_items)), # to save the number of observations per item	
          as.integer(student_starting), # acc urning persons
          as.integer(item_starting), # acc urning items
          as.integer(pRT_student_starting),# pseudoRT urning persons
          as.integer(pRT_item_starting),# psuedoRT urning items
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(student_urn_size),
          as.integer(pRT_student_urn_size),
          as.integer(item_urn_size),
          as.integer(pRT_item_urn_size),
          as.double(rep(0,n_items+1)),
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
          as.integer(container_long),
          as.integer(container_long))
  
  #-----------------------------------------------------------------------------
  # results
  n_obs_p=tmp[[44]] # for every response how many responses the person gave so far
  n_obs_i=tmp[[45]] # for every response how many responses the item had so far
  
  simulated1=tmp[[7]] # simulated accuracy
  simulated2=tmp[[8]] # simulated pseudoRT
  
  Ubefore=tmp[[15]] # acc urning of the person before the update
  Vbefore=tmp[[16]] # acc urning of the item before the update
  
  
  U2before=tmp[[17]]# pseudoRT urning of the person before the update
  V2before=tmp[[18]]# pseudoRT urning of the item before the update
  
  # state of the urning after putting the balls in but before removing them
  U=Ubefore+X1*2 
  V=Vbefore+(1-X1)*2
  U2=U2before+X2
  V2=V2before+(1-X2)
  
  #-----------------------------------------------------------------------------
  # return object
  res = list("student_urn_size" = student_urn_size,
             "pRT_student_urn_size" = pRT_student_urn_size,
             "X1" = X1,
             "X2" = X2,
             "item_urn_size" = item_urn_size,
             "pRT_item_urn_size" = pRT_item_urn_size,
             "student_u_after" = U,
             "item_u_after" = V,
             "student_pRT_u_after" = U2,
             "item_pRT_u_after" = V2,
             "student_u_before" = Ubefore,
             "item_u_before" = Vbefore,
             "student_pRT_u_before" = U2before,
             "item_pRT_u_before" = V2before,
             "simulated_accuarcy" = simulated1,
             "simulated_pRT" = simulated2,
             "n_observations_student" = n_obs_p,
             "n_observations_item" = n_obs_i)
  
  return(res)
}

urnings_combined_RT_data = function(data,
                                    student_urn_size,
                                    item_urn_size,
                                    fast_indicator,
                                    OS = "MAC"){
  
  #-----------------------------------------------------------------------------
  # FUNCTION ARGUMENTS
  # data: data frame with the following columns:
  # - student_id: integer, unique identifier for the student (1:n_students)
  # - item_id: integer, unique identifier for the item (1:n_items)
  # - response: integer, 0 or 1, the response of the student (accuracy)
  # - pRT_1: integer 0 for slow 1 for fast response based on the grand median of all response times ("grand_median")
  # - pRT_2: integer 0 for slow 1 for fast response based on the median of the given  ("item_median")
  # - pRT_3: integer 0 for slow 1 for fast response based on the halving point of the 0-90th quantile of the response times of the given item  ("half_90_quantile")
  # student_urn_size: integer, urn size of the accuracy urn for the students
  # item_urn_size: integer, urn size of the accuracy urn for the items
  # pRT_student_urn_size: integer, urn size of the pseudo RT variable (half of the student urn size normally)
  # pRT_item_urn_size: integer, urn size of the pseudo RT variable (half of the item urn size normally)
  # fast_indicator: string, which definition of fast/slow responses is used see for the correct values above
  # OS: string, which operating system is used (for selecting the appropriate compiled C routine)
  
  # RETURN
  # list
  # see comments below 
  
  #-----------------------------------------------------------------------------
  # loading the compiled C routine
  switch(OS,
         "MAC" = dyn.load("CRT_old.so"),
         "WINDOWS" = dyn.load("CRT_old.dll"),
         "LINUX" = dyn.load("CRT_old.o")
  )
  
  #-----------------------------------------------------------------------------
  # setting up the arguments for the C wrapper
  
  #-----------------------------------------------------------------------------
  #disassembling the data to fit the C routine
  
  n_responses = nrow(data) 
  student_id = data$student_id 
  item_id=data$item_id 
  n_students=length(unique(student_id)) 
  n_items=length(unique(item_id)) 
  X1=data$response # accuracy data
  
  X2 = switch(fast_indicator, 
              "grand_median" = data$pRT_1,
              "item_median" = data$pRT_2,
              "half_90_quantile" = data$pRT_3)
  
  #-----------------------------------------------------------------------------
  # setting algorithm settings
  
  # set initial urnings at (0.5 x urnsize)
  # for accuracy
  student_starting = rep(student_urn_size/2,n_students) 
  item_starting =rep(item_urn_size/2,n_items) 
  
  #-----------------------------------------------------------------------------
  # auxiliary objects for paired update
  
  queue<-rep(0,n_items)
  
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
  
  LL<-rep(0,n_scores*n_items)
  LLsum<-rep(0,n_scores)
  

  #-----------------------------------------------------------------------------
  #containers for the results
  # empty vector with as many elements as there are responses
  #it is used to save all sort of things within the c code
  container_long = rep(0,n_responses) 
  
  tmp<-.C("urnings_combined_RT_real_data",
          as.integer(n_responses),
          as.integer(student_id),
          as.integer(item_id),
          as.integer(n_items),
          as.integer(X1),
          as.integer(X2),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(rep(0,n_students)),
          as.integer(rep(0,n_items)),
          as.integer(student_starting),
          as.integer(item_starting),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(container_long),
          as.integer(student_urn_size),
          as.integer(item_urn_size),
          as.double(rep(0,n_items+1)),
          as.integer(Score),
          as.integer(n_scores),
          as.integer(n_options),
          as.integer(Upd),
          as.integer(queue),
          as.integer(LL),
          as.integer(LLsum),
          as.integer(container_long),
          as.integer(container_long))
  
  #-----------------------------------------------------------------------------
  # results
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
  
  #-----------------------------------------------------------------------------
  # return object
  res = list("student_urn_size" = student_urn_size,
             "X1" = X1,
             "X2" = X2,
             "item_urn_size" = item_urn_size,
             "student_pRT_u_after" = U2,
             "item_pRT_u_after" = V2,
             "student_u_before" = Ubefore,
             "item_u_before" = Vbefore,
             "student_u_middle" = Umiddle,
             "item_u_middle" = Vmiddle,
             "simulated_accuarcy" = simulated1,
             "simulated_pRT" = simulated2,
             "n_observations_student" = n_obs_p,
             "n_observations_item" = n_obs_i)
  
  return(res)
}

################################################################################
#Utility functions 
################################################################################
evaluate_fit_mu = function(res, 
                        threshold = 20){
  
  require(tidyverse)
  require(patchwork)
  #-----------------------------------------------------------------------------
  # pseudo RT
  
  #possible urning ratings
  h_student =c(0:res$pRT_student_urn_size)
  h_item = c(0:res$pRT_item_urn_size)
  all_u_student = length(h_student)
  all_u_item = length(h_item)
  
  #helpers to make the code human readable 
  student_obs_above_threshold = res$n_observations_student >= threshold
  item_obs_above_threshold = res$n_observations_item >= threshold
  
  #-----------------------------------------------------------------------------
  # fit plot for the pseudo RT variable
  
  Observed2=Expected2=N_obs2=matrix(ncol=all_u_item,nrow=all_u_student)
  for(i in 1:all_u_student){
    for(j in 1:all_u_item){
      #-------------------------------------------------------------------------
      # number of observation per urn size bins
      
      N_obs2[i,j]=sum(student_obs_above_threshold &
                      item_obs_above_threshold & 
                      (res$student_pRT_u_after==h_student[i]) & 
                      (res$item_pRT_u_after==h_item[j]))
      
      #-------------------------------------------------------------------------
      # observed frequency per urn size bins
      
      Observed2[i,j]=mean(res$X2[student_obs_above_threshold &
                                 item_obs_above_threshold & 
                                 (res$student_pRT_u_after == h_student[i])&
                                 (res$item_pRT_u_after == h_item[j])])
      
      #-------------------------------------------------------------------------
      # observed frequency per urn size bins
      
      Expected2[i,j]=mean(res$simulated_pRT[student_obs_above_threshold &
                                              item_obs_above_threshold & 
                                              (res$student_pRT_u_after == h_student[i])&
                                              (res$item_pRT_u_after == h_item[j])])
    }
  }
  
  #-----------------------------------------------------------------------------
  # fit plot for the accuracy data
  # this needs to be corrected to (0:11)*2
  
  h_student =c(0:res$pRT_student_urn_size) * 2
  h_item = c(0:res$pRT_item_urn_size) * 2
  
  Observed = Expected = N_obs = matrix(ncol=all_u_item,nrow=all_u_student)
  for(i in 1:all_u_student){
    for(j in 1:all_u_item){
      #-------------------------------------------------------------------------
      # number of observation per urn size bins
      
      N_obs[i,j]=sum(student_obs_above_threshold &
                        item_obs_above_threshold & 
                        (res$student_u_after==h_student[i]) & 
                        (res$item_u_after==h_item[j]))
      
      #-------------------------------------------------------------------------
      # observed frequency per urn size bins
      Observed[i,j]=mean(res$X1[student_obs_above_threshold &
                                   item_obs_above_threshold & 
                                   (res$student_u_after == h_student[i])&
                                   (res$item_u_after == h_item[j])])
      
      #-------------------------------------------------------------------------
      # observed frequency per urn size bins
      
      Expected[i,j]=mean(res$simulated_accuarcy[student_obs_above_threshold &
                                              item_obs_above_threshold & 
                                              (res$student_u_after == h_student[i])&
                                              (res$item_u_after == h_item[j])])
      
    }
  }
  
  #-----------------------------------------------------------------------------
  #plotting
  #plotting expected versus observed after flattening the matrices against each other as a scatter plot
  #the observation count sets the size of the dots
  
  plot_dat_acc = data.frame("observed" = c(Observed), "expected" = c(Expected), "n_obs" = c(N_obs) / 5000)
  plot_dat_pRT = data.frame("observed" = c(Observed2), "expected" = c(Expected2), "n_obs" = c(N_obs2) / 5000)
  
  p_acc = plot_dat_acc %>%
    ggplot(aes(x = observed, y = expected)) +
    geom_point(aes(size = n_obs), fill = NA, shape=1) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "red", linetype = "dotted") +
    jtools::theme_apa() +
    labs(x = "Observed",
         y = "Expected") + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") 
  
  p_pRT = plot_dat_pRT %>%
    ggplot(aes(x = observed, y = expected)) +
    geom_point(aes(size = n_obs), fill = NA, shape=1) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "red", linetype = "dotted") +
    jtools::theme_apa() +
    labs(x = "Observed",
         y = "Expected") + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") 
  
  return(p_acc + p_pRT)
}

evaluate_fit_ms = function(res, 
                           threshold = 20){
  
  require(tidyverse)
  require(patchwork)
  #-----------------------------------------------------------------------------
  # pseudo RT
  
  #possible urning ratings
  h_student =c(0:res$student_urn_size)
  h_item = c(0:res$item_urn_size)
  all_u_student = length(h_student)
  all_u_item = length(h_item)
  
  #helpers to make the code human readable 
  student_obs_above_threshold = res$n_observations_student >= threshold
  item_obs_above_threshold = res$n_observations_item >= threshold
  
  #-----------------------------------------------------------------------------
  # fit plot for the pseudo RT variable
  
  Observed2=Expected2=N_obs2=matrix(ncol=all_u_item,nrow=all_u_student)
  for(i in 1:all_u_student){
    for(j in 1:all_u_item){
      #-------------------------------------------------------------------------
      # number of observation per urn size bins
      
      N_obs2[i,j]=sum(student_obs_above_threshold &
                        item_obs_above_threshold & 
                        (res$student_pRT_u_after==h_student[i]) & 
                        (res$item_pRT_u_after==h_item[j]))
      
      #-------------------------------------------------------------------------
      # observed frequency per urn size bins
      
      Observed2[i,j]=mean(res$X2[student_obs_above_threshold &
                                   item_obs_above_threshold & 
                                   (res$student_pRT_u_after == h_student[i])&
                                   (res$item_pRT_u_after == h_item[j])])
      
      #-------------------------------------------------------------------------
      # expected frequency per urn size bins
      
      Expected2[i,j]=mean(res$simulated_pRT[student_obs_above_threshold &
                                              item_obs_above_threshold & 
                                              (res$student_pRT_u_after == h_student[i])&
                                              (res$item_pRT_u_after == h_item[j])])
    }
  }
  
  #-----------------------------------------------------------------------------
  # fit plot for the accuracy data
  U=res$student_u_before+res$X1*2
  V=res$item_u_before+(1-res$X1)*2
  
  Observed = Expected = N_obs = matrix(ncol=all_u_item,nrow=all_u_student)
  for(i in 1:all_u_student){
    for(j in 1:all_u_item){
      #-------------------------------------------------------------------------
      # number of observation per urn size bins
      
      N_obs[i,j]=sum(student_obs_above_threshold &
                       item_obs_above_threshold & 
                       (U==h_student[i]) & 
                       (V==h_item[j]))
      
      #-------------------------------------------------------------------------
      # observed frequency per urn size bins
      Observed[i,j]=mean(res$X1[student_obs_above_threshold &
                                  item_obs_above_threshold & 
                                  (U == h_student[i])&
                                  (V == h_item[j])])
      
      #-------------------------------------------------------------------------
      # expected frequency per urn size bins
      
      Expected[i,j]=mean(res$simulated_accuarcy[student_obs_above_threshold &
                                                  item_obs_above_threshold & 
                                                  (U == h_student[i])&
                                                  (V == h_item[j])])
      
    }
  }
  
  #-----------------------------------------------------------------------------
  #plotting
  #plotting expected versus observed after flattening the matrices against each other as a scatter plot
  #the observation count sets the size of the dots
  
  plot_dat_acc = data.frame("observed" = c(Observed), "expected" = c(Expected), "n_obs" = c(N_obs))
  plot_dat_pRT = data.frame("observed" = c(Observed2), "expected" = c(Expected2), "n_obs" = c(N_obs2))
  
  p_acc = plot_dat_acc %>%
    ggplot(aes(x = observed, y = expected)) +
    geom_point(aes(size = n_obs), fill = NA, shape=1) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "red", linetype = "dotted") +
    jtools::theme_apa() +
    labs(x = "Observed",
         y = "Expected") + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") 
  
  p_pRT = plot_dat_pRT %>%
    ggplot(aes(x = observed, y = expected)) +
    geom_point(aes(size = n_obs), fill = NA, shape=1) +
    geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), color = "red", linetype = "dotted") +
    jtools::theme_apa() +
    labs(x = "Observed",
         y = "Expected") + 
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none") 
  
  return(p_acc + p_pRT)
}


prediction_accuracy = function(res, 
                               is.multiple_urn = TRUE,
                               threshold = 20,
                               level = "RESPONSE",
                               data = NULL){
  
  #helpers to make the code human readable 
  student_obs_above_threshold = res$n_observations_student >= threshold
  item_obs_above_threshold = res$n_observations_item >= threshold
  source("util.R")
  if(is.multiple_urn){
    U=(res$student_u_before+res$student_pRT_u_before)/sum(res$student_urn_size, res$pRT_student_urn_size)
    V=(res$item_u_before+res$item_pRT_u_before)/sum(res$item_urn_size, res$pRT_item_urn_size)
  } else {
    U = res$student_u_before / res$student_urn_size
    V = res$item_u_before / res$item_urn_size
  }
  
  Predicted=twoPL(U,V,2)
  mse_pred = (res$X1[student_obs_above_threshold & 
                       item_obs_above_threshold] - 
                Predicted[student_obs_above_threshold & 
                            item_obs_above_threshold])^2
  #calculate mse per students 
  #mean over all mse values for each students
  if(level == "STUDENT"){
  dat = readRDS("empirical_data/emp_dat_analysis.rds")
  dat = dat$student_id
  mse_person = tapply(mse_pred, dat[student_obs_above_threshold & 
                                               item_obs_above_threshold], mean)
  }
  return(ifelse(level == "RESPONSE", mean(mse_pred), mean(mse_person)))
}

create_theta = function(mu_theta, sd_theta, n_games, n_students, new = 8, seed = 19186553, new_mean = log(0.7/0.3)){
  set.seed(seed)
  theta  = qnorm(seq(1/(n_students+1),n_students/(n_students+1),by=1/(n_students+1)),mu_theta,sd_theta)
  theta = sample(theta, n_students)
  theta = c(theta, (new_mean + seq(-2,2,length=new))/4) #8 new students with true values between -1 and 1, without 0 true value
  n_students = n_students + new
  Theta=matrix(0, ncol=n_games,nrow=n_students)
  for(j in 1:n_games){
    Theta[,j]=theta
  }
  return(Theta)
}

create_starting = function(param, n_agent, urn_size, items = FALSE, new = 8, cold_value = 0.5){
  if(items == TRUE){
    starting = rep(urn_size / 2, n_agent)
    true_probs = 1/(1+exp(-param[1:(n_agent/2)]))
    starting[1:(n_agent/2)] = rbinom(n_agent/2, urn_size, true_probs)
    starting[n_agent:(n_agent/2 + 1)] = urn_size - starting[1:(n_agent/2)]
  } else if(new > 0){
    starting = rep(round(cold_value*urn_size), n_agent + new)
    new_indicies = (n_agent+1):(n_agent+new)
    true_probs = 1/(1+exp(-param[-new_indicies]))
    starting[-new_indicies]=rbinom(n_agent,urn_size,true_probs)
  } else if(new == 0){
    true_probs = 1/(1+exp(-param))
    starting=rbinom(n_agent,urn_size,true_probs)
  }
  return(starting)
}


sim_observed_SRT = function(theta, delta, n_games, order){
  outcome = array(0,c(length(theta), n_games, 3))
  for(i in 1:length(theta)){
    curr_id = order[i,]
    curr_item = delta[curr_id]
    outcome[i,,1] = rbinom(n_games, size = 1, 1/(1+exp(4*(curr_item-theta[i]))))
    outcome[i,,2] = rbinom(n_games, size = 1, 1/(1+exp(2*(curr_item-theta[i]))))
    outcome[i,,3] = rbinom(n_games, size = 1, 1/(1+exp(1*(curr_item-theta[i]))))
  }
  return(outcome)
} 

sim_observed_CISRT = function(theta, delta, n_games, order){
  outcome = array(0,c(length(theta), n_games, 3))
  for(i in 1:length(theta)){
    curr_id = order[i, ]
    curr_item = delta[curr_id]
    outcome[i,,1] = acc = rbinom(n_games, size = 1, 1/(1+exp(4*(curr_item-theta[i]))))
    rev_acc = abs(acc-1)
    outcome[i,,2] = rbinom(n_games, size = 1, 1/(1+exp(2*(curr_item-theta[i])))) - rev_acc
    outcome[i,,3] = rbinom(n_games, size = 1, 1/(1+exp(1*(curr_item-theta[i])))) - rev_acc
  }
  return(outcome)
} 

sim_observed_LSRT = function(theta_ability, theta_speed, delta, n_games, order){
  n_students = length(theta_ability)
  outcome = array(0,c(n_students, n_games, 3))
  for(i in 1:n_students){
    curr_id = order[i, ]
    curr_item = delta[curr_id]
    prob = lapply(curr_item,
                  LSRT,
                  lower = 4,
                  upper = 7,
                  theta_1 = theta_ability[i],
                  theta_2 = theta_speed[i],
                  norm = 8)
    outcome[i,,1] = rbinom(n_games, size = 1, unlist(prob))
    prob = lapply(curr_item,
                   LSRT,
                   lower = 2,
                   upper = 3,
                   theta_1 = theta_ability[i],
                   theta_2 = theta_speed[i],
                   norm = 4)
    outcome[i,,2] = rbinom(n_games, size = 1, unlist(prob))
    prob = lapply(curr_item,
                  LSRT,
                  lower = 1,
                  upper = 1,
                  theta_1 = theta_ability[i],
                  theta_2 = theta_speed[i],
                  norm = 2)
    outcome[i,,2] = rbinom(n_games, size = 1, unlist(prob))
  }
  return(outcome)
}

LSRT = function(delta, lower,upper,theta_1, theta_2, norm = 8){
  c = list(0:7, c(1,0.5,0.25,0,0,0.25,0.5,1))
  norm_constant = numeric(norm)
  for(nc in 1:norm){
    norm_constant[nc] = exp(c[[1]][nc] * theta_1 + c[[2]][nc] * theta_2 - c[[1]][nc] * delta)
  }
  norm_constant = 1/sum(norm_constant)
  if(lower != upper){
    kernel = numeric(upper - lower + 1)
    counter = 1
    for(k in lower:upper){
      kernel[counter] = exp(c[[1]][k+1] * theta_1 + c[[2]][k+1] * theta_2 - c[[1]][k+1] * delta)
      counter = counter + 1
    }
    kernel = sum(kernel)
  } else {
    kernel = exp(c[[1]][lower+1] * theta_1 + c[[2]][lower+1] * theta_2 - c[[1]][lower+1] * delta)
  }

  return(norm_constant * kernel)  
}

