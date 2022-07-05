library(rARPACK)
library(Matrix)
library(matlib)
calculate_R0 <- function(parameters,S = NULL, N = NULL){
  ## X_0 = (E, I, A, H), X_1 = (S,R,D)
  R0 <- with(parameters,{
    if(is.null(S)){ ## supposes that exists S_0 in parameters and compute R_0
      contact = diag(S_0)%*%(c%*%diag(1/S_0))
    } else{ ## then computes Rt based on the supplied S
      if(is.null(N)){## supposes that N == S_0 (might be wrong because of dead population)
        contact = diag(S)%*%(c%*%diag(1/S_0))
      } else{ ## computes the more correct Rt
        contact = diag(S)%*%(c%*%diag(1/N))
      }
    } 
    #### F matrix ####
    # E line
    F_EE <- beta*omega*contact
    F_EI <- beta*contact    
    F_EA <- beta*omega_a*contact
    F_EH <- beta*omega_s*contact
    F_E <- cbind(F_EE, F_EI, F_EA, F_EH)
    # other lines are zero
    F_mat <- rbind(F_E, matrix(0, nrow = 3*nrow(F_E),ncol = ncol(F_E)))
    #### T matrix #### writing in the lazy way because of readability
    # E line
    T_EE <- diag(1-gamma, nrow = 16, ncol = 16)
    T_EI <- T_EA <- T_EH <- matrix(0, ncol = 16, nrow = 16)
    T_E <- cbind(T_EE, T_EI, T_EA, T_EH)
    # I line
    T_IE <- diag((1-ihr)*gamma*(1-asymp))
    T_II <- diag(1-nu, nrow = 16, ncol = 16)
    T_IA <- T_IH <- matrix(0, nrow = 16, ncol = 16)
    T_I <- cbind(T_IE, T_II, T_IA, T_IH)
    # A line
    T_AE <- diag((1-ihr)*gamma*asymp)
    T_AA <- diag(1-nu, nrow = 16, ncol = 16)
    T_AI <- T_AH <- matrix(0, nrow = 16, ncol = 16)
    T_A <- cbind(T_AE, T_AI, T_AA, T_AH)
    # H line
    T_HE <- diag(ihr*gamma)
    T_HH <- diag(1-nus, nrow = 16, ncol = 16)
    T_HI <- T_HA <- matrix(0, nrow = 16, ncol = 16)
    T_H <- cbind(T_HE, T_HI, T_HA, T_HH)
    # final matrix
    T_mat <- rbind(T_E, T_I, T_A, T_H)
    # identity matrix
    Id_mat <- diag(1,nrow(T_mat),ncol(T_mat))
    # sparse because why not
    F_mat <- Matrix(F_mat, sparse = T)
    # automatically computes (I - T)^-1 and converts to sparse -- faster
    T_mat <- Matrix(inv(Id_mat - T_mat), sparse = T)
    # R0 is the maximum eigenvalue (that might be complex)
    R0 <- max(Re(eigs(F_mat%*%T_mat, 1, retvec = F)$values))
    # R0 <- eigs(F_mat%*%T_mat, 45)$values
    })
  return(R0)
}
