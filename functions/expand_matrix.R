expand_matrix <- function(mat, dist_et){
  m.length<-dim(mat)[1]
  m2.length <- nrow(dist_et)
  pop<- sum(dist_et$pop[m.length:m2.length])
  pop2 <- matrix(0,nrow = m2.length - m.length + 1,ncol = m2.length- m.length + 1 )
  for(i in m.length:m2.length){
    for(j in m.length:m2.length){
      pop2[i-m.length+1,j-m.length+1] <- dist_et$pop[i]*dist_et$pop[j]
    }
  }
  mat2 <- matrix(0,nrow = m2.length,ncol = m2.length)
  mat2[1:m.length,1:m.length] <- mat
  for(i in m.length:m2.length){
    mat2[1:(m.length-1),i] <- mat[1:(m.length-1),m.length]*dist_et$pop[i]/pop
    mat2[i,1:(m.length-1)] <- mat[m.length,1:(m.length-1)]*dist_et$pop[i]/pop
  }
  for(i in m.length:m2.length){
    for(j in m.length:m2.length){
      mat2[i,j] <- mat[m.length,m.length]*pop2[i-m.length+1,j-m.length+1]/sum(pop2)
    }
  }
  return(mat2)
}