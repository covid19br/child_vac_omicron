#### T matrix #### writing in the lazy way because of readability

## Eu
# Eu.DxE cols
T_Eu.DEu.D <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_Eu.DEvA.D <- T_Eu.DEvP.D <- T_Eu.DEvC.D <- 
  T_Eu.DEwA.D <- T_Eu.DEwP.D <- T_Eu.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EuDxI cols
T_Eu.DIu.D <- T_Eu.DIvA.D <- T_Eu.DIvP.D <- T_Eu.DIvC.D <- 
  T_Eu.DIwA.D <- T_Eu.DIwP.D <- T_Eu.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EuDxA cols
T_Eu.DAu.D <- T_Eu.DAvA.D <- T_Eu.DAvP.D <- T_Eu.DAvC.D <- 
  T_Eu.DAwA.D <- T_Eu.DAwP.D <- T_Eu.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EuDxH cols
T_Eu.DHu.D <- T_Eu.DHvA.D <- T_Eu.DHvP.D <- T_Eu.DHvC.D <- 
  T_Eu.DHwA.D <- T_Eu.DHwP.D <- T_Eu.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_Eu.D <- cbind(T_Eu.DEu.D, T_Eu.DEvA.D, T_Eu.DEvP.D, T_Eu.DEvC.D, T_Eu.DEwA.D, T_Eu.DEwP.D, T_Eu.DEwC.D,
                T_Eu.DAu.D, T_Eu.DAvA.D, T_Eu.DAvP.D, T_Eu.DAvC.D, T_Eu.DAwA.D, T_Eu.DAwP.D, T_Eu.DAwC.D, 
                T_Eu.DIu.D, T_Eu.DIvA.D, T_Eu.DIvP.D, T_Eu.DIvC.D, T_Eu.DIwA.D, T_Eu.DIwP.D, T_Eu.DIwC.D, 
                T_Eu.DHu.D, T_Eu.DHvA.D, T_Eu.DHvP.D, T_Eu.DHvC.D, T_Eu.DHwA.D, T_Eu.DHwP.D, T_Eu.DHwC.D)
## EvA
# EvA.DxE cols
T_EvA.DEvA.D <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EvA.DEu.D <- T_EvA.DEvP.D <- T_EvA.DEvC.D <- 
  T_EvA.DEwA.D <- T_EvA.DEwP.D <- T_EvA.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvA.DxA cols
T_EvA.DAu.D <- T_EvA.DAvA.D <- T_EvA.DAvP.D <- T_EvA.DAvC.D <- 
  T_EvA.DAwA.D <- T_EvA.DAwP.D <- T_EvA.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvA.DxI cols
T_EvA.DIu.D <- T_EvA.DIvA.D <- T_EvA.DIvP.D <- T_EvA.DIvC.D <- 
  T_EvA.DIwA.D <- T_EvA.DIwP.D <- T_EvA.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvAxH cols
T_EvA.DHu.D <- T_EvA.DHvA.D <- T_EvA.DHvP.D <- T_EvA.DHvC.D <- 
  T_EvA.DHwA.D <- T_EvA.DHwP.D <- T_EvA.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EvA.D <- cbind(T_EvA.DEu.D, T_EvA.DEvA.D, T_EvA.DEvP.D, T_EvA.DEvC.D, T_EvA.DEwA.D, T_EvA.DEwP.D, T_EvA.DEwC.D,
                T_EvA.DAu.D, T_EvA.DAvA.D, T_EvA.DAvP.D, T_EvA.DAvC.D, T_EvA.DAwA.D, T_EvA.DAwP.D, T_EvA.DAwC.D, 
                T_EvA.DIu.D, T_EvA.DIvA.D, T_EvA.DIvP.D, T_EvA.DIvC.D, T_EvA.DIwA.D, T_EvA.DIwP.D, T_EvA.DIwC.D, 
                T_EvA.DHu.D, T_EvA.DHvA.D, T_EvA.DHvP.D, T_EvA.DHvC.D, T_EvA.DHwA.D, T_EvA.DHwP.D, T_EvA.DHwC.D)

## EvP
# EvP.DxE cols
T_EvP.DEvP.D <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EvP.DEu.D <- T_EvP.DEvA.D <- T_EvP.DEvC.D <- 
  T_EvP.DEwA.D <- T_EvP.DEwP.D <- T_EvP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvP.DxA cols
T_EvP.DAu.D <- T_EvP.DAvA.D <- T_EvP.DAvP.D <- T_EvP.DAvC.D <- 
  T_EvP.DAwA.D <- T_EvP.DAwP.D <- T_EvP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvP.DxI cols
T_EvP.DIu.D <- T_EvP.DIvA.D <- T_EvP.DIvP.D <- T_EvP.DIvC.D <- 
  T_EvP.DIwA.D <- T_EvP.DIwP.D <- T_EvP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvPxH cols
T_EvP.DHu.D <- T_EvP.DHvA.D <- T_EvP.DHvP.D <- T_EvP.DHvC.D <- 
  T_EvP.DHwA.D <- T_EvP.DHwP.D <- T_EvP.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EvP.D <- cbind(T_EvP.DEu.D, T_EvP.DEvA.D, T_EvP.DEvP.D, T_EvP.DEvC.D, T_EvP.DEwA.D, T_EvP.DEwP.D, T_EvP.DEwC.D,
                T_EvP.DAu.D, T_EvP.DAvA.D, T_EvP.DAvP.D, T_EvP.DAvC.D, T_EvP.DAwA.D, T_EvP.DAwP.D, T_EvP.DAwC.D, 
                T_EvP.DIu.D, T_EvP.DIvA.D, T_EvP.DIvP.D, T_EvP.DIvC.D, T_EvP.DIwA.D, T_EvP.DIwP.D, T_EvP.DIwC.D, 
                T_EvP.DHu.D, T_EvP.DHvA.D, T_EvP.DHvP.D, T_EvP.DHvC.D, T_EvP.DHwA.D, T_EvP.DHwP.D, T_EvP.DHwC.D)

## EvC
# EvC.DxE cols
T_EvC.DEvC.D <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EvC.DEu.D <- T_EvC.DEvA.D <- T_EvC.DEvP.D <- 
  T_EvC.DEwA.D <- T_EvC.DEwP.D <- T_EvC.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvC.DxA cols
T_EvC.DAu.D <- T_EvC.DAvA.D <- T_EvC.DAvP.D <- T_EvC.DAvC.D <- 
  T_EvC.DAwA.D <- T_EvC.DAwP.D <- T_EvC.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvC.DxI cols
T_EvC.DIu.D <- T_EvC.DIvA.D <- T_EvC.DIvP.D <- T_EvC.DIvC.D <- 
  T_EvC.DIwA.D <- T_EvC.DIwP.D <- T_EvC.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvCxH cols
T_EvC.DHu.D <- T_EvC.DHvA.D <- T_EvC.DHvP.D <- T_EvC.DHvC.D <- 
  T_EvC.DHwA.D <- T_EvC.DHwP.D <- T_EvC.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EvC.D <- cbind(T_EvC.DEu.D, T_EvC.DEvA.D, T_EvC.DEvP.D, T_EvC.DEvC.D, T_EvC.DEwA.D, T_EvC.DEwP.D, T_EvC.DEwC.D,
                T_EvC.DAu.D, T_EvC.DAvA.D, T_EvC.DAvP.D, T_EvC.DAvC.D, T_EvC.DAwA.D, T_EvC.DAwP.D, T_EvC.DAwC.D, 
                T_EvC.DIu.D, T_EvC.DIvA.D, T_EvC.DIvP.D, T_EvC.DIvC.D, T_EvC.DIwA.D, T_EvC.DIwP.D, T_EvC.DIwC.D, 
                T_EvC.DHu.D, T_EvC.DHvA.D, T_EvC.DHvP.D, T_EvC.DHvC.D, T_EvC.DHwA.D, T_EvC.DHwP.D, T_EvC.DHwC.D)

## EwA
# EwA.DxE cols
T_EwA.DEwA.D <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EwA.DEu.D <- T_EwA.DEvA.D <- T_EwA.DEvP.D <- 
  T_EwA.DEvC.D <- T_EwA.DEwP.D <- T_EwA.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwA.DxA cols
T_EwA.DAu.D <- T_EwA.DAvA.D <- T_EwA.DAvP.D <- T_EwA.DAvC.D <- 
  T_EwA.DAwA.D <- T_EwA.DAwP.D <- T_EwA.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwA.DxI cols
T_EwA.DIu.D <- T_EwA.DIvA.D <- T_EwA.DIvP.D <- T_EwA.DIvC.D <- 
  T_EwA.DIwA.D <- T_EwA.DIwP.D <- T_EwA.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwAxH cols
T_EwA.DHu.D <- T_EwA.DHvA.D <- T_EwA.DHvP.D <- T_EwA.DHvC.D <- 
  T_EwA.DHwA.D <- T_EwA.DHwP.D <- T_EwA.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EwA.D <- cbind(T_EwA.DEu.D, T_EwA.DEvA.D, T_EwA.DEvP.D, T_EwA.DEvC.D, T_EwA.DEwA.D, T_EwA.DEwP.D, T_EwA.DEwC.D,
                T_EwA.DAu.D, T_EwA.DAvA.D, T_EwA.DAvP.D, T_EwA.DAvC.D, T_EwA.DAwA.D, T_EwA.DAwP.D, T_EwA.DAwC.D, 
                T_EwA.DIu.D, T_EwA.DIvA.D, T_EwA.DIvP.D, T_EwA.DIvC.D, T_EwA.DIwA.D, T_EwA.DIwP.D, T_EwA.DIwC.D, 
                T_EwA.DHu.D, T_EwA.DHvA.D, T_EwA.DHvP.D, T_EwA.DHvC.D, T_EwA.DHwA.D, T_EwA.DHwP.D, T_EwA.DHwC.D)

## EwP
# EwP.DxE cols
T_EwP.DEwP.D <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EwP.DEu.D <- T_EwP.DEvA.D <- T_EwP.DEvP.D <- 
  T_EwP.DEvC.D <- T_EwP.DEwA.D <- T_EwP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwP.DxA cols
T_EwP.DAu.D <- T_EwP.DAvA.D <- T_EwP.DAvP.D <- T_EwP.DAvC.D <- 
  T_EwP.DAwA.D <- T_EwP.DAwP.D <- T_EwP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwP.DxI cols
T_EwP.DIu.D <- T_EwP.DIvA.D <- T_EwP.DIvP.D <- T_EwP.DIvC.D <- 
  T_EwP.DIwA.D <- T_EwP.DIwP.D <- T_EwP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwPxH cols
T_EwP.DHu.D <- T_EwP.DHvA.D <- T_EwP.DHvP.D <- T_EwP.DHvC.D <- 
  T_EwP.DHwA.D <- T_EwP.DHwP.D <- T_EwP.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EwP.D <- cbind(T_EwP.DEu.D, T_EwP.DEvA.D, T_EwP.DEvP.D, T_EwP.DEvC.D, T_EwP.DEwA.D, T_EwP.DEwP.D, T_EwP.DEwC.D,
                T_EwP.DAu.D, T_EwP.DAvA.D, T_EwP.DAvP.D, T_EwP.DAvC.D, T_EwP.DAwA.D, T_EwP.DAwP.D, T_EwP.DAwC.D, 
                T_EwP.DIu.D, T_EwP.DIvA.D, T_EwP.DIvP.D, T_EwP.DIvC.D, T_EwP.DIwA.D, T_EwP.DIwP.D, T_EwP.DIwC.D, 
                T_EwP.DHu.D, T_EwP.DHvA.D, T_EwP.DHvP.D, T_EwP.DHvC.D, T_EwP.DHwA.D, T_EwP.DHwP.D, T_EwP.DHwC.D)

## EwC
# EwC.DxE cols
T_EwC.DEwC.D <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EwC.DEu.D <- T_EwC.DEvA.D <- T_EwC.DEvP.D <- 
  T_EwC.DEvC.D <- T_EwC.DEwA.D <- T_EwC.DEwP.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwC.DxA cols
T_EwC.DAu.D <- T_EwC.DAvA.D <- T_EwC.DAvP.D <- T_EwC.DAvC.D <- 
  T_EwC.DAwA.D <- T_EwC.DAwP.D <- T_EwC.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwC.DxI cols
T_EwC.DIu.D <- T_EwC.DIvA.D <- T_EwC.DIvP.D <- T_EwC.DIvC.D <- 
  T_EwC.DIwA.D <- T_EwC.DIwP.D <- T_EwC.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwC.DxH cols
T_EwC.DHu.D <- T_EwC.DHvA.D <- T_EwC.DHvP.D <- T_EwC.DHvC.D <- 
  T_EwC.DHwA.D <- T_EwC.DHwP.D <- T_EwC.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EwC.D <- cbind(T_EwC.DEu.D, T_EwC.DEvA.D, T_EwC.DEvP.D, T_EwC.DEvC.D, T_EwC.DEwA.D, T_EwC.DEwP.D, T_EwC.DEwC.D,
                T_EwC.DAu.D, T_EwC.DAvA.D, T_EwC.DAvP.D, T_EwC.DAvC.D, T_EwC.DAwA.D, T_EwC.DAwP.D, T_EwC.DAwC.D, 
                T_EwC.DIu.D, T_EwC.DIvA.D, T_EwC.DIvP.D, T_EwC.DIvC.D, T_EwC.DIwA.D, T_EwC.DIwP.D, T_EwC.DIwC.D, 
                T_EwC.DHu.D, T_EwC.DHvA.D, T_EwC.DHvP.D, T_EwC.DHvC.D, T_EwC.DHwA.D, T_EwC.DHwP.D, T_EwC.DHwC.D)

#all E's together
T_E.D <- rbind(T_Eu.D, T_EvA.D, T_EvP.D, T_EvC.D, T_EwA.D, T_EwP.D, T_EwC.D)

# # A line
# T_AE <- diag((1-ihr)*gamma*asymp)
# T_AA <- diag(1-nu, nrow = 16, ncol = 16)
# T_AI <- T_AH <- matrix(0, nrow = 16, ncol = 16)
# T_A <- cbind(T_AE, T_AI, T_AA, T_AH)

## Au
# Au.DxE cols
T_Au.DEu.D <- diag((1-ihr)*gamma*asymp, ncol = age.bins, nrow = age.bins)
T_Au.DEvA.D <- T_Au.DEvP.D <- 
  T_Au.DEvC.D <- T_Au.DEwA.D <- T_Au.DEwP.D <- T_Au.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Au.DxA cols
T_Au.DAu.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_Au.DAvA.D <- T_Au.DAvP.D <- T_Au.DAvC.D <- 
  T_Au.DAwA.D <- T_Au.DAwP.D <- T_Au.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Au.DxI cols
T_Au.DIu.D <- T_Au.DIvA.D <- T_Au.DIvP.D <- T_Au.DIvC.D <- 
  T_Au.DIwA.D <- T_Au.DIwP.D <- T_Au.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Au.DxH cols
T_Au.DHu.D <- T_Au.DHvA.D <- T_Au.DHvP.D <- T_Au.DHvC.D <- 
  T_Au.DHwA.D <- T_Au.DHwP.D <- T_Au.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_Au.D <- cbind(T_Au.DEu.D, T_Au.DEvA.D, T_Au.DEvP.D, T_Au.DEvC.D, T_Au.DEwA.D, T_Au.DEwP.D, T_Au.DEwC.D,
                T_Au.DAu.D, T_Au.DAvA.D, T_Au.DAvP.D, T_Au.DAvC.D, T_Au.DAwA.D, T_Au.DAwP.D, T_Au.DAwC.D, 
                T_Au.DIu.D, T_Au.DIvA.D, T_Au.DIvP.D, T_Au.DIvC.D, T_Au.DIwA.D, T_Au.DIwP.D, T_Au.DIwC.D, 
                T_Au.DHu.D, T_Au.DHvA.D, T_Au.DHvP.D, T_Au.DHvC.D, T_Au.DHwA.D, T_Au.DHwP.D, T_Au.DHwC.D)
                
## AvA
# AvA.DxE cols
T_AvA.DEvA.D <- diag((1-ihr_vA.D)*gamma*asymp_vA.D, ncol = age.bins, nrow = age.bins)
T_AvA.DEu.D <- T_AvA.DEvP.D <- 
  T_AvA.DEvC.D <- T_AvA.DEwA.D <- T_AvA.DEwP.D <- T_AvA.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvA.DxA cols
T_AvA.DAvA.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AvA.DAu.D <- T_AvA.DAvP.D <- T_AvA.DAvC.D <- 
  T_AvA.DAwA.D <- T_AvA.DAwP.D <- T_AvA.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvA.DxI cols
T_AvA.DIu.D <- T_AvA.DIvA.D <- T_AvA.DIvP.D <- T_AvA.DIvC.D <- 
  T_AvA.DIwA.D <- T_AvA.DIwP.D <- T_AvA.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvA.DxH cols
T_AvA.DHu.D <- T_AvA.DHvA.D <- T_AvA.DHvP.D <- T_AvA.DHvC.D <- 
  T_AvA.DHwA.D <- T_AvA.DHwP.D <- T_AvA.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AvA.D <- cbind(T_AvA.DEu.D, T_AvA.DEvA.D, T_AvA.DEvP.D, T_AvA.DEvC.D, T_AvA.DEwA.D, T_AvA.DEwP.D, T_AvA.DEwC.D,
                T_AvA.DAu.D, T_AvA.DAvA.D, T_AvA.DAvP.D, T_AvA.DAvC.D, T_AvA.DAwA.D, T_AvA.DAwP.D, T_AvA.DAwC.D, 
                T_AvA.DIu.D, T_AvA.DIvA.D, T_AvA.DIvP.D, T_AvA.DIvC.D, T_AvA.DIwA.D, T_AvA.DIwP.D, T_AvA.DIwC.D, 
                T_AvA.DHu.D, T_AvA.DHvA.D, T_AvA.DHvP.D, T_AvA.DHvC.D, T_AvA.DHwA.D, T_AvA.DHwP.D, T_AvA.DHwC.D)

## AvP
# AvP.DxE cols
T_AvP.DEvP.D <- diag((1-ihr_vP.D)*gamma*asymp_vP.D, ncol = age.bins, nrow = age.bins)
T_AvP.DEu.D <- T_AvP.DEvA.D <- 
  T_AvP.DEvC.D <- T_AvP.DEwA.D <- T_AvP.DEwP.D <- T_AvP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvP.DxA cols
T_AvP.DAvP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AvP.DAu.D <- T_AvP.DAvA.D <- T_AvP.DAvC.D <- 
  T_AvP.DAwA.D <- T_AvP.DAwP.D <- T_AvP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvP.DxI cols
T_AvP.DIu.D <- T_AvP.DIvA.D <- T_AvP.DIvP.D <- T_AvP.DIvC.D <- 
  T_AvP.DIwA.D <- T_AvP.DIwP.D <- T_AvP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvP.DxH cols
T_AvP.DHu.D <- T_AvP.DHvA.D <- T_AvP.DHvP.D <- T_AvP.DHvC.D <- 
  T_AvP.DHwA.D <- T_AvP.DHwP.D <- T_AvP.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AvP.D <- cbind(T_AvP.DEu.D, T_AvP.DEvA.D, T_AvP.DEvP.D, T_AvP.DEvC.D, T_AvP.DEwA.D, T_AvP.DEwP.D, T_AvP.DEwC.D,
                T_AvP.DAu.D, T_AvP.DAvA.D, T_AvP.DAvP.D, T_AvP.DAvC.D, T_AvP.DAwA.D, T_AvP.DAwP.D, T_AvP.DAwC.D, 
                T_AvP.DIu.D, T_AvP.DIvA.D, T_AvP.DIvP.D, T_AvP.DIvC.D, T_AvP.DIwA.D, T_AvP.DIwP.D, T_AvP.DIwC.D, 
                T_AvP.DHu.D, T_AvP.DHvA.D, T_AvP.DHvP.D, T_AvP.DHvC.D, T_AvP.DHwA.D, T_AvP.DHwP.D, T_AvP.DHwC.D)
                
## AvC
# AvC.DxE cols
T_AvC.DEvC.D <- diag((1-ihr_vC.D)*gamma*asymp_vC.D, ncol = age.bins, nrow = age.bins)
T_AvC.DEu.D <- T_AvC.DEvA.D <- 
  T_AvC.DEvP.D <- T_AvC.DEwA.D <- T_AvC.DEwP.D <- T_AvC.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvC.DxA cols
T_AvC.DAvC.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AvC.DAu.D <- T_AvC.DAvA.D <- T_AvC.DAvP.D <- 
  T_AvC.DAwA.D <- T_AvC.DAwP.D <- T_AvC.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvC.DxI cols
T_AvC.DIu.D <- T_AvC.DIvA.D <- T_AvC.DIvP.D <- T_AvC.DIvC.D <- 
  T_AvC.DIwA.D <- T_AvC.DIwP.D <- T_AvC.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvC.DxH cols
T_AvC.DHu.D <- T_AvC.DHvA.D <- T_AvC.DHvP.D <- T_AvC.DHvC.D <- 
  T_AvC.DHwA.D <- T_AvC.DHwP.D <- T_AvC.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AvC.D <- cbind(T_AvC.DEu.D, T_AvC.DEvA.D, T_AvC.DEvP.D, T_AvC.DEvC.D, T_AvC.DEwA.D, T_AvC.DEwP.D, T_AvC.DEwC.D,
                T_AvC.DAu.D, T_AvC.DAvA.D, T_AvC.DAvP.D, T_AvC.DAvC.D, T_AvC.DAwA.D, T_AvC.DAwP.D, T_AvC.DAwC.D, 
                T_AvC.DIu.D, T_AvC.DIvA.D, T_AvC.DIvP.D, T_AvC.DIvC.D, T_AvC.DIwA.D, T_AvC.DIwP.D, T_AvC.DIwC.D, 
                T_AvC.DHu.D, T_AvC.DHvA.D, T_AvC.DHvP.D, T_AvC.DHvC.D, T_AvC.DHwA.D, T_AvC.DHwP.D, T_AvC.DHwC.D)

## AwA
# AwA.DxE cols
T_AwA.DEwA.D <- diag((1-ihr_wA.D)*gamma*asymp_wA.D, ncol = age.bins, nrow = age.bins)
T_AwA.DEu.D <- T_AwA.DEvP.D <- 
  T_AwA.DEvC.D <- T_AwA.DEvA.D <- T_AwA.DEwP.D <- T_AwA.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwA.DxA cols
T_AwA.DAwA.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AwA.DAu.D <- T_AwA.DAvP.D <- T_AwA.DAvC.D <- 
  T_AwA.DAvA.D <- T_AwA.DAwP.D <- T_AwA.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwA.DxI cols
T_AwA.DIu.D <- T_AwA.DIvA.D <- T_AwA.DIvP.D <- T_AwA.DIvC.D <- 
  T_AwA.DIwA.D <- T_AwA.DIwP.D <- T_AwA.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwA.DxH cols
T_AwA.DHu.D <- T_AwA.DHvA.D <- T_AwA.DHvP.D <- T_AwA.DHvC.D <- 
  T_AwA.DHwA.D <- T_AwA.DHwP.D <- T_AwA.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AwA.D <- cbind(T_AwA.DEu.D, T_AwA.DEvA.D, T_AwA.DEvP.D, T_AwA.DEvC.D, T_AwA.DEwA.D, T_AwA.DEwP.D, T_AwA.DEwC.D,
                T_AwA.DAu.D, T_AwA.DAvA.D, T_AwA.DAvP.D, T_AwA.DAvC.D, T_AwA.DAwA.D, T_AwA.DAwP.D, T_AwA.DAwC.D, 
                T_AwA.DIu.D, T_AwA.DIvA.D, T_AwA.DIvP.D, T_AwA.DIvC.D, T_AwA.DIwA.D, T_AwA.DIwP.D, T_AwA.DIwC.D, 
                T_AwA.DHu.D, T_AwA.DHvA.D, T_AwA.DHvP.D, T_AwA.DHvC.D, T_AwA.DHwA.D, T_AwA.DHwP.D, T_AwA.DHwC.D)

## AwP
# AwP.DxE cols
T_AwP.DEwP.D <- diag((1-ihr_wP.D)*gamma*asymp_wP.D, ncol = age.bins, nrow = age.bins)
T_AwP.DEu.D <- T_AwP.DEvA.D <- 
  T_AwP.DEvC.D <- T_AwP.DEwA.D <- T_AwP.DEvP.D <- T_AwP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwP.DxA cols
T_AwP.DAwP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AwP.DAu.D <- T_AwP.DAvA.D <- T_AwP.DAvC.D <- 
  T_AwP.DAwA.D <- T_AwP.DAvP.D <- T_AwP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwP.DxI cols
T_AwP.DIu.D <- T_AwP.DIvA.D <- T_AwP.DIvP.D <- T_AwP.DIvC.D <- 
  T_AwP.DIwA.D <- T_AwP.DIwP.D <- T_AwP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwP.DxH cols
T_AwP.DHu.D <- T_AwP.DHvA.D <- T_AwP.DHvP.D <- T_AwP.DHvC.D <- 
  T_AwP.DHwA.D <- T_AwP.DHwP.D <- T_AwP.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AwP.D <- cbind(T_AwP.DEu.D, T_AwP.DEvA.D, T_AwP.DEvP.D, T_AwP.DEvC.D, T_AwP.DEwA.D, T_AwP.DEwP.D, T_AwP.DEwC.D,
                T_AwP.DAu.D, T_AwP.DAvA.D, T_AwP.DAvP.D, T_AwP.DAvC.D, T_AwP.DAwA.D, T_AwP.DAwP.D, T_AwP.DAwC.D, 
                T_AwP.DIu.D, T_AwP.DIvA.D, T_AwP.DIvP.D, T_AwP.DIvC.D, T_AwP.DIwA.D, T_AwP.DIwP.D, T_AwP.DIwC.D, 
                T_AwP.DHu.D, T_AwP.DHvA.D, T_AwP.DHvP.D, T_AwP.DHvC.D, T_AwP.DHwA.D, T_AwP.DHwP.D, T_AwP.DHwC.D)
                
## AwC
# AwC.DxE cols
T_AwC.DEwC.D <- diag((1-ihr_wC.D)*gamma*asymp_wC.D, ncol = age.bins, nrow = age.bins)
T_AwC.DEu.D <- T_AwC.DEvA.D <- 
  T_AwC.DEvP.D <- T_AwC.DEwA.D <- T_AwC.DEwP.D <- T_AwC.DEvC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwC.DxA cols
T_AwC.DAwC.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AwC.DAu.D <- T_AwC.DAvA.D <- T_AwC.DAvP.D <- 
  T_AwC.DAwA.D <- T_AwC.DAwP.D <- T_AwC.DAvC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwC.DxI cols
T_AwC.DIu.D <- T_AwC.DIvA.D <- T_AwC.DIvP.D <- T_AwC.DIvC.D <- 
  T_AwC.DIwA.D <- T_AwC.DIwP.D <- T_AwC.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwC.DxH cols
T_AwC.DHu.D <- T_AwC.DHvA.D <- T_AwC.DHvP.D <- T_AwC.DHvC.D <- 
  T_AwC.DHwA.D <- T_AwC.DHwP.D <- T_AwC.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AwC.D <- cbind(T_AwC.DEu.D, T_AwC.DEvA.D, T_AwC.DEvP.D, T_AwC.DEvC.D, T_AwC.DEwA.D, T_AwC.DEwP.D, T_AwC.DEwC.D,
                T_AwC.DAu.D, T_AwC.DAvA.D, T_AwC.DAvP.D, T_AwC.DAvC.D, T_AwC.DAwA.D, T_AwC.DAwP.D, T_AwC.DAwC.D, 
                T_AwC.DIu.D, T_AwC.DIvA.D, T_AwC.DIvP.D, T_AwC.DIvC.D, T_AwC.DIwA.D, T_AwC.DIwP.D, T_AwC.DIwC.D, 
                T_AwC.DHu.D, T_AwC.DHvA.D, T_AwC.DHvP.D, T_AwC.DHvC.D, T_AwC.DHwA.D, T_AwC.DHwP.D, T_AwC.DHwC.D)

#all A's together
T_A.D <- rbind(T_Au.D, T_AvA.D, T_AvP.D, T_AvC.D, T_AwA.D, T_AwP.D, T_AwC.D)


# # I line
# T_IE <- diag((1-ihr)*gamma*(1-asymp))
# T_II <- diag(1-nu, nrow = 16, ncol = 16)
# T_IA <- T_IH <- matrix(0, nrow = 16, ncol = 16)
# T_I <- cbind(T_IE, T_II, T_IA, T_IH)

## Iu
# Iu.DxE cols
T_Iu.DEu.D <- diag((1-ihr)*gamma*(1-asymp), ncol = age.bins, nrow = age.bins)
T_Iu.DEvA.D <- T_Iu.DEvP.D <- 
  T_Iu.DEvC.D <- T_Iu.DEwA.D <- T_Iu.DEwP.D <- T_Iu.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Iu.DxI cols
T_Iu.DIu.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_Iu.DIvA.D <- T_Iu.DIvP.D <- T_Iu.DIvC.D <- 
  T_Iu.DIwA.D <- T_Iu.DIwP.D <- T_Iu.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Iu.DxA cols
T_Iu.DAu.D <- T_Iu.DAvA.D <- T_Iu.DAvP.D <- T_Iu.DAvC.D <- 
  T_Iu.DAwA.D <- T_Iu.DAwP.D <- T_Iu.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Iu.DxH cols
T_Iu.DHu.D <- T_Iu.DHvA.D <- T_Iu.DHvP.D <- T_Iu.DHvC.D <- 
  T_Iu.DHwA.D <- T_Iu.DHwP.D <- T_Iu.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_Iu.D <- cbind(T_Iu.DEu.D, T_Iu.DEvA.D, T_Iu.DEvP.D, T_Iu.DEvC.D, T_Iu.DEwA.D, T_Iu.DEwP.D, T_Iu.DEwC.D,
                T_Iu.DAu.D, T_Iu.DAvA.D, T_Iu.DAvP.D, T_Iu.DAvC.D, T_Iu.DAwA.D, T_Iu.DAwP.D, T_Iu.DAwC.D, 
                T_Iu.DIu.D, T_Iu.DIvA.D, T_Iu.DIvP.D, T_Iu.DIvC.D, T_Iu.DIwA.D, T_Iu.DIwP.D, T_Iu.DIwC.D, 
                T_Iu.DHu.D, T_Iu.DHvA.D, T_Iu.DHvP.D, T_Iu.DHvC.D, T_Iu.DHwA.D, T_Iu.DHwP.D, T_Iu.DHwC.D)
                
## IvA
# IvA.DxE cols
T_IvA.DEvA.D <- diag((1-ihr_vA.D)*gamma*(1-asymp_vA.D), ncol = age.bins, nrow = age.bins)
T_IvA.DEu.D <- T_IvA.DEvP.D <- 
  T_IvA.DEvC.D <- T_IvA.DEwA.D <- T_IvA.DEwP.D <- T_IvA.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvA.DxI cols
T_IvA.DIvA.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IvA.DIu.D <- T_IvA.DIvP.D <- T_IvA.DIvC.D <- 
  T_IvA.DIwA.D <- T_IvA.DIwP.D <- T_IvA.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvA.DxA cols
T_IvA.DAu.D <- T_IvA.DAvA.D <- T_IvA.DAvP.D <- T_IvA.DAvC.D <- 
  T_IvA.DAwA.D <- T_IvA.DAwP.D <- T_IvA.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvA.DxH cols
T_IvA.DHu.D <- T_IvA.DHvA.D <- T_IvA.DHvP.D <- T_IvA.DHvC.D <- 
  T_IvA.DHwA.D <- T_IvA.DHwP.D <- T_IvA.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IvA.D <- cbind(T_IvA.DEu.D, T_IvA.DEvA.D, T_IvA.DEvP.D, T_IvA.DEvC.D, T_IvA.DEwA.D, T_IvA.DEwP.D, T_IvA.DEwC.D,
                T_IvA.DAu.D, T_IvA.DAvA.D, T_IvA.DAvP.D, T_IvA.DAvC.D, T_IvA.DAwA.D, T_IvA.DAwP.D, T_IvA.DAwC.D, 
                T_IvA.DIu.D, T_IvA.DIvA.D, T_IvA.DIvP.D, T_IvA.DIvC.D, T_IvA.DIwA.D, T_IvA.DIwP.D, T_IvA.DIwC.D, 
                T_IvA.DHu.D, T_IvA.DHvA.D, T_IvA.DHvP.D, T_IvA.DHvC.D, T_IvA.DHwA.D, T_IvA.DHwP.D, T_IvA.DHwC.D)

## IvP
# IvP.DxE cols
T_IvP.DEvP.D <- diag((1-ihr_vP.D)*gamma*(1-asymp_vP.D), ncol = age.bins, nrow = age.bins)
T_IvP.DEu.D <- T_IvP.DEvA.D <- 
  T_IvP.DEvC.D <- T_IvP.DEwA.D <- T_IvP.DEwP.D <- T_IvP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvP.DxI cols
T_IvP.DIvP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IvP.DIu.D <- T_IvP.DIvA.D <- T_IvP.DIvC.D <- 
  T_IvP.DIwA.D <- T_IvP.DIwP.D <- T_IvP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvP.DxA cols
T_IvP.DAu.D <- T_IvP.DAvA.D <- T_IvP.DAvP.D <- T_IvP.DAvC.D <- 
  T_IvP.DAwA.D <- T_IvP.DAwP.D <- T_IvP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvP.DxH cols
T_IvP.DHu.D <- T_IvP.DHvA.D <- T_IvP.DHvP.D <- T_IvP.DHvC.D <- 
  T_IvP.DHwA.D <- T_IvP.DHwP.D <- T_IvP.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IvP.D <- cbind(T_IvP.DEu.D, T_IvP.DEvA.D, T_IvP.DEvP.D, T_IvP.DEvC.D, T_IvP.DEwA.D, T_IvP.DEwP.D, T_IvP.DEwC.D,
                T_IvP.DAu.D, T_IvP.DAvA.D, T_IvP.DAvP.D, T_IvP.DAvC.D, T_IvP.DAwA.D, T_IvP.DAwP.D, T_IvP.DAwC.D, 
                T_IvP.DIu.D, T_IvP.DIvA.D, T_IvP.DIvP.D, T_IvP.DIvC.D, T_IvP.DIwA.D, T_IvP.DIwP.D, T_IvP.DIwC.D, 
                T_IvP.DHu.D, T_IvP.DHvA.D, T_IvP.DHvP.D, T_IvP.DHvC.D, T_IvP.DHwA.D, T_IvP.DHwP.D, T_IvP.DHwC.D)
                
## IvC
# IvC.DxE cols
T_IvC.DEvC.D <- diag((1-ihr_vC.D)*gamma*(1-asymp_vC.D), ncol = age.bins, nrow = age.bins)
T_IvC.DEu.D <- T_IvC.DEvA.D <- 
  T_IvC.DEvP.D <- T_IvC.DEwA.D <- T_IvC.DEwP.D <- T_IvC.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvC.DxI cols
T_IvC.DIvC.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IvC.DIu.D <- T_IvC.DIvA.D <- T_IvC.DIvP.D <- 
  T_IvC.DIwA.D <- T_IvC.DIwP.D <- T_IvC.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvC.DxA cols
T_IvC.DAu.D <- T_IvC.DAvA.D <- T_IvC.DAvP.D <- T_IvC.DAvC.D <- 
  T_IvC.DAwA.D <- T_IvC.DAwP.D <- T_IvC.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvC.DxH cols
T_IvC.DHu.D <- T_IvC.DHvA.D <- T_IvC.DHvP.D <- T_IvC.DHvC.D <- 
  T_IvC.DHwA.D <- T_IvC.DHwP.D <- T_IvC.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IvC.D <- cbind(T_IvC.DEu.D, T_IvC.DEvA.D, T_IvC.DEvP.D, T_IvC.DEvC.D, T_IvC.DEwA.D, T_IvC.DEwP.D, T_IvC.DEwC.D,
                T_IvC.DAu.D, T_IvC.DAvA.D, T_IvC.DAvP.D, T_IvC.DAvC.D, T_IvC.DAwA.D, T_IvC.DAwP.D, T_IvC.DAwC.D, 
                T_IvC.DIu.D, T_IvC.DIvA.D, T_IvC.DIvP.D, T_IvC.DIvC.D, T_IvC.DIwA.D, T_IvC.DIwP.D, T_IvC.DIwC.D, 
                T_IvC.DHu.D, T_IvC.DHvA.D, T_IvC.DHvP.D, T_IvC.DHvC.D, T_IvC.DHwA.D, T_IvC.DHwP.D, T_IvC.DHwC.D)

## IwA
# IwA.DxE cols
T_IwA.DEwA.D <- diag((1-ihr_wA.D)*gamma*(1-asymp_wA.D), ncol = age.bins, nrow = age.bins)
T_IwA.DEu.D <- T_IwA.DEvP.D <- 
  T_IwA.DEvC.D <- T_IwA.DEvA.D <- T_IwA.DEwP.D <- T_IwA.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwA.DxI cols
T_IwA.DIwA.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IwA.DIu.D <- T_IwA.DIvP.D <- T_IwA.DIvC.D <- 
  T_IwA.DIvA.D <- T_IwA.DIwP.D <- T_IwA.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwA.DxA cols
T_IwA.DAu.D <- T_IwA.DAvA.D <- T_IwA.DAvP.D <- T_IwA.DAvC.D <- 
  T_IwA.DAwA.D <- T_IwA.DAwP.D <- T_IwA.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwA.DxH cols
T_IwA.DHu.D <- T_IwA.DHvA.D <- T_IwA.DHvP.D <- T_IwA.DHvC.D <- 
  T_IwA.DHwA.D <- T_IwA.DHwP.D <- T_IwA.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IwA.D <- cbind(T_IwA.DEu.D, T_IwA.DEvA.D, T_IwA.DEvP.D, T_IwA.DEvC.D, T_IwA.DEwA.D, T_IwA.DEwP.D, T_IwA.DEwC.D,
                T_IwA.DAu.D, T_IwA.DAvA.D, T_IwA.DAvP.D, T_IwA.DAvC.D, T_IwA.DAwA.D, T_IwA.DAwP.D, T_IwA.DAwC.D, 
                T_IwA.DIu.D, T_IwA.DIvA.D, T_IwA.DIvP.D, T_IwA.DIvC.D, T_IwA.DIwA.D, T_IwA.DIwP.D, T_IwA.DIwC.D, 
                T_IwA.DHu.D, T_IwA.DHvA.D, T_IwA.DHvP.D, T_IwA.DHvC.D, T_IwA.DHwA.D, T_IwA.DHwP.D, T_IwA.DHwC.D)

## IwP
# IwP.DxE cols
T_IwP.DEwP.D <- diag((1-ihr_wP.D)*gamma*(1-asymp_wP.D), ncol = age.bins, nrow = age.bins)
T_IwP.DEu.D <- T_IwP.DEvA.D <- 
  T_IwP.DEvC.D <- T_IwP.DEwA.D <- T_IwP.DEvP.D <- T_IwP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwP.DxI cols
T_IwP.DIwP.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IwP.DIu.D <- T_IwP.DIvA.D <- T_IwP.DIvC.D <- 
  T_IwP.DIwA.D <- T_IwP.DIvP.D <- T_IwP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwP.DxA cols
T_IwP.DAu.D <- T_IwP.DAvA.D <- T_IwP.DAvP.D <- T_IwP.DAvC.D <- 
  T_IwP.DAwA.D <- T_IwP.DAwP.D <- T_IwP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwP.DxH cols
T_IwP.DHu.D <- T_IwP.DHvA.D <- T_IwP.DHvP.D <- T_IwP.DHvC.D <- 
  T_IwP.DHwA.D <- T_IwP.DHwP.D <- T_IwP.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IwP.D <- cbind(T_IwP.DEu.D, T_IwP.DEvA.D, T_IwP.DEvP.D, T_IwP.DEvC.D, T_IwP.DEwA.D, T_IwP.DEwP.D, T_IwP.DEwC.D,
                T_IwP.DAu.D, T_IwP.DAvA.D, T_IwP.DAvP.D, T_IwP.DAvC.D, T_IwP.DAwA.D, T_IwP.DAwP.D, T_IwP.DAwC.D, 
                T_IwP.DIu.D, T_IwP.DIvA.D, T_IwP.DIvP.D, T_IwP.DIvC.D, T_IwP.DIwA.D, T_IwP.DIwP.D, T_IwP.DIwC.D, 
                T_IwP.DHu.D, T_IwP.DHvA.D, T_IwP.DHvP.D, T_IwP.DHvC.D, T_IwP.DHwA.D, T_IwP.DHwP.D, T_IwP.DHwC.D)
                
## IwC
# IwC.DxE cols
T_IwC.DEwC.D <- diag((1-ihr_wC.D)*gamma*(1-asymp_wC.D), ncol = age.bins, nrow = age.bins)
T_IwC.DEu.D <- T_IwC.DEvA.D <- 
  T_IwC.DEvP.D <- T_IwC.DEwA.D <- T_IwC.DEwP.D <- T_IwC.DEvC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwC.DxI cols
T_IwC.DIwC.D <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IwC.DIu.D <- T_IwC.DIvA.D <- T_IwC.DIvP.D <- 
  T_IwC.DIwA.D <- T_IwC.DIwP.D <- T_IwC.DIvC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwC.DxA cols
T_IwC.DAu.D <- T_IwC.DAvA.D <- T_IwC.DAvP.D <- T_IwC.DAvC.D <- 
  T_IwC.DAwA.D <- T_IwC.DAwP.D <- T_IwC.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwC.DxH cols
T_IwC.DHu.D <- T_IwC.DHvA.D <- T_IwC.DHvP.D <- T_IwC.DHvC.D <- 
  T_IwC.DHwA.D <- T_IwC.DHwP.D <- T_IwC.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IwC.D <- cbind(T_IwC.DEu.D, T_IwC.DEvA.D, T_IwC.DEvP.D, T_IwC.DEvC.D, T_IwC.DEwA.D, T_IwC.DEwP.D, T_IwC.DEwC.D,
                T_IwC.DAu.D, T_IwC.DAvA.D, T_IwC.DAvP.D, T_IwC.DAvC.D, T_IwC.DAwA.D, T_IwC.DAwP.D, T_IwC.DAwC.D, 
                T_IwC.DIu.D, T_IwC.DIvA.D, T_IwC.DIvP.D, T_IwC.DIvC.D, T_IwC.DIwA.D, T_IwC.DIwP.D, T_IwC.DIwC.D, 
                T_IwC.DHu.D, T_IwC.DHvA.D, T_IwC.DHvP.D, T_IwC.DHvC.D, T_IwC.DHwA.D, T_IwC.DHwP.D, T_IwC.DHwC.D)

#all I's together
T_I.D <- rbind(T_Iu.D, T_IvA.D, T_IvP.D, T_IvC.D, T_IwA.D, T_IwP.D, T_IwC.D)

# # H line
# T_HE <- diag(ihr*gamma)
# T_HH <- diag(1-nus, nrow = 16, ncol = 16)
# T_HI <- T_HA <- matrix(0, nrow = 16, ncol = 16)
# T_H <- cbind(T_HE, T_HI, T_HA, T_HH)

## Hu
# Hu.DxE cols
T_Hu.DEu.D <- diag(ihr*gamma, ncol = age.bins, nrow = age.bins)
T_Hu.DEvA.D <- T_Hu.DEvP.D <- 
  T_Hu.DEvC.D <- T_Hu.DEwA.D <- T_Hu.DEwP.D <- T_Hu.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Hu.DxH cols
T_Hu.DHu.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_Hu.DHvA.D <- T_Hu.DHvP.D <- T_Hu.DHvC.D <- 
  T_Hu.DHwA.D <- T_Hu.DHwP.D <- T_Hu.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Hu.DxI cols
T_Hu.DIu.D <- T_Hu.DIvA.D <- T_Hu.DIvP.D <- T_Hu.DIvC.D <- 
  T_Hu.DIwA.D <- T_Hu.DIwP.D <- T_Hu.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# Hu.DxA cols
T_Hu.DAu.D <- T_Hu.DAvA.D <- T_Hu.DAvP.D <- T_Hu.DAvC.D <- 
  T_Hu.DAwA.D <- T_Hu.DAwP.D <- T_Hu.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_Hu.D <- cbind(T_Hu.DEu.D, T_Hu.DEvA.D, T_Hu.DEvP.D, T_Hu.DEvC.D, T_Hu.DEwA.D, T_Hu.DEwP.D, T_Hu.DEwC.D,
                T_Hu.DAu.D, T_Hu.DAvA.D, T_Hu.DAvP.D, T_Hu.DAvC.D, T_Hu.DAwA.D, T_Hu.DAwP.D, T_Hu.DAwC.D, 
                T_Hu.DIu.D, T_Hu.DIvA.D, T_Hu.DIvP.D, T_Hu.DIvC.D, T_Hu.DIwA.D, T_Hu.DIwP.D, T_Hu.DIwC.D, 
                T_Hu.DHu.D, T_Hu.DHvA.D, T_Hu.DHvP.D, T_Hu.DHvC.D, T_Hu.DHwA.D, T_Hu.DHwP.D, T_Hu.DHwC.D)
                
## HvA
# HvA.DxE cols
T_HvA.DEvA.D <- diag(ihr_vA.D*gamma, ncol = age.bins, nrow = age.bins)
T_HvA.DEu.D <- T_HvA.DEvP.D <- 
  T_HvA.DEvC.D <- T_HvA.DEwA.D <- T_HvA.DEwP.D <- T_HvA.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvA.DxH cols
T_HvA.DHvA.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HvA.DHu.D <- T_HvA.DHvP.D <- T_HvA.DHvC.D <- 
  T_HvA.DHwA.D <- T_HvA.DHwP.D <- T_HvA.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvA.DxI cols
T_HvA.DIu.D <- T_HvA.DIvA.D <- T_HvA.DIvP.D <- T_HvA.DIvC.D <- 
  T_HvA.DIwA.D <- T_HvA.DIwP.D <- T_HvA.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvA.DxA cols
T_HvA.DAu.D <- T_HvA.DAvA.D <- T_HvA.DAvP.D <- T_HvA.DAvC.D <- 
  T_HvA.DAwA.D <- T_HvA.DAwP.D <- T_HvA.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HvA.D <- cbind(T_HvA.DEu.D, T_HvA.DEvA.D, T_HvA.DEvP.D, T_HvA.DEvC.D, T_HvA.DEwA.D, T_HvA.DEwP.D, T_HvA.DEwC.D,
                T_HvA.DAu.D, T_HvA.DAvA.D, T_HvA.DAvP.D, T_HvA.DAvC.D, T_HvA.DAwA.D, T_HvA.DAwP.D, T_HvA.DAwC.D, 
                T_HvA.DIu.D, T_HvA.DIvA.D, T_HvA.DIvP.D, T_HvA.DIvC.D, T_HvA.DIwA.D, T_HvA.DIwP.D, T_HvA.DIwC.D, 
                T_HvA.DHu.D, T_HvA.DHvA.D, T_HvA.DHvP.D, T_HvA.DHvC.D, T_HvA.DHwA.D, T_HvA.DHwP.D, T_HvA.DHwC.D)

## HvP
# HvP.DxE cols
T_HvP.DEvP.D <- diag(ihr_vP.D*gamma, ncol = age.bins, nrow = age.bins)
T_HvP.DEu.D <- T_HvP.DEvA.D <- 
  T_HvP.DEvC.D <- T_HvP.DEwA.D <- T_HvP.DEwP.D <- T_HvP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvP.DxH cols
T_HvP.DHvP.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HvP.DHu.D <- T_HvP.DHvA.D <- T_HvP.DHvC.D <- 
  T_HvP.DHwA.D <- T_HvP.DHwP.D <- T_HvP.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvP.DxI cols
T_HvP.DIu.D <- T_HvP.DIvA.D <- T_HvP.DIvP.D <- T_HvP.DIvC.D <- 
  T_HvP.DIwA.D <- T_HvP.DIwP.D <- T_HvP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvP.DxA cols
T_HvP.DAu.D <- T_HvP.DAvA.D <- T_HvP.DAvP.D <- T_HvP.DAvC.D <- 
  T_HvP.DAwA.D <- T_HvP.DAwP.D <- T_HvP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HvP.D <- cbind(T_HvP.DEu.D, T_HvP.DEvA.D, T_HvP.DEvP.D, T_HvP.DEvC.D, T_HvP.DEwA.D, T_HvP.DEwP.D, T_HvP.DEwC.D,
                T_HvP.DAu.D, T_HvP.DAvA.D, T_HvP.DAvP.D, T_HvP.DAvC.D, T_HvP.DAwA.D, T_HvP.DAwP.D, T_HvP.DAwC.D, 
                T_HvP.DIu.D, T_HvP.DIvA.D, T_HvP.DIvP.D, T_HvP.DIvC.D, T_HvP.DIwA.D, T_HvP.DIwP.D, T_HvP.DIwC.D, 
                T_HvP.DHu.D, T_HvP.DHvA.D, T_HvP.DHvP.D, T_HvP.DHvC.D, T_HvP.DHwA.D, T_HvP.DHwP.D, T_HvP.DHwC.D)
                
## HvC
# HvC.DxE cols
T_HvC.DEvC.D <- diag(ihr_vC.D*gamma, ncol = age.bins, nrow = age.bins)
T_HvC.DEu.D <- T_HvC.DEvA.D <- 
  T_HvC.DEvP.D <- T_HvC.DEwA.D <- T_HvC.DEwP.D <- T_HvC.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvC.DxH cols
T_HvC.DHvC.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HvC.DHu.D <- T_HvC.DHvA.D <- T_HvC.DHvP.D <- 
  T_HvC.DHwA.D <- T_HvC.DHwP.D <- T_HvC.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvC.DxI cols
T_HvC.DIu.D <- T_HvC.DIvA.D <- T_HvC.DIvP.D <- T_HvC.DIvC.D <- 
  T_HvC.DIwA.D <- T_HvC.DIwP.D <- T_HvC.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvC.DxA cols
T_HvC.DAu.D <- T_HvC.DAvA.D <- T_HvC.DAvP.D <- T_HvC.DAvC.D <- 
  T_HvC.DAwA.D <- T_HvC.DAwP.D <- T_HvC.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HvC.D <- cbind(T_HvC.DEu.D, T_HvC.DEvA.D, T_HvC.DEvP.D, T_HvC.DEvC.D, T_HvC.DEwA.D, T_HvC.DEwP.D, T_HvC.DEwC.D,
                T_HvC.DAu.D, T_HvC.DAvA.D, T_HvC.DAvP.D, T_HvC.DAvC.D, T_HvC.DAwA.D, T_HvC.DAwP.D, T_HvC.DAwC.D, 
                T_HvC.DIu.D, T_HvC.DIvA.D, T_HvC.DIvP.D, T_HvC.DIvC.D, T_HvC.DIwA.D, T_HvC.DIwP.D, T_HvC.DIwC.D, 
                T_HvC.DHu.D, T_HvC.DHvA.D, T_HvC.DHvP.D, T_HvC.DHvC.D, T_HvC.DHwA.D, T_HvC.DHwP.D, T_HvC.DHwC.D)

## HwA
# HwA.DxE cols
T_HwA.DEwA.D <- diag(ihr_wA.D*gamma, ncol = age.bins, nrow = age.bins)
T_HwA.DEu.D <- T_HwA.DEvP.D <- 
  T_HwA.DEvC.D <- T_HwA.DEvA.D <- T_HwA.DEwP.D <- T_HwA.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwA.DxH cols
T_HwA.DAwH.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HwA.DHu.D <- T_HwA.DHvP.D <- T_HwA.DHvC.D <- 
  T_HwA.DHvA.D <- T_HwA.DHwP.D <- T_HwA.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwA.DxI cols
T_HwA.DIu.D <- T_HwA.DIvA.D <- T_HwA.DIvP.D <- T_HwA.DIvC.D <- 
  T_HwA.DIwA.D <- T_HwA.DIwP.D <- T_HwA.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwA.DxA cols
T_HwA.DAu.D <- T_HwA.DAvA.D <- T_HwA.DAvP.D <- T_HwA.DAvC.D <- 
  T_HwA.DAwA.D <- T_HwA.DAwP.D <- T_HwA.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HwA.D <- cbind(T_HwA.DEu.D, T_HwA.DEvA.D, T_HwA.DEvP.D, T_HwA.DEvC.D, T_HwA.DEwA.D, T_HwA.DEwP.D, T_HwA.DEwC.D,
                T_HwA.DAu.D, T_HwA.DAvA.D, T_HwA.DAvP.D, T_HwA.DAvC.D, T_HwA.DAwA.D, T_HwA.DAwP.D, T_HwA.DAwC.D, 
                T_HwA.DIu.D, T_HwA.DIvA.D, T_HwA.DIvP.D, T_HwA.DIvC.D, T_HwA.DIwA.D, T_HwA.DIwP.D, T_HwA.DIwC.D, 
                T_HwA.DHu.D, T_HwA.DHvA.D, T_HwA.DHvP.D, T_HwA.DHvC.D, T_HwA.DHwA.D, T_HwA.DHwP.D, T_HwA.DHwC.D)

## HwP
# HwP.DxE cols
T_HwP.DEwP.D <- diag(ihr_wP.D*gamma, ncol = age.bins, nrow = age.bins)
T_HwP.DEu.D <- T_HwP.DEvA.D <- 
  T_HwP.DEvC.D <- T_HwP.DEwA.D <- T_HwP.DEvP.D <- T_HwP.DEwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwP.DxH cols
T_HwP.DHwP.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HwP.DHu.D <- T_HwP.DHvA.D <- T_HwP.DHvC.D <- 
  T_HwP.DHwA.D <- T_HwP.DHvP.D <- T_HwP.DHwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwP.DxI cols
T_HwP.DIu.D <- T_HwP.DIvA.D <- T_HwP.DIvP.D <- T_HwP.DIvC.D <- 
  T_HwP.DIwA.D <- T_HwP.DIwP.D <- T_HwP.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwP.DxA cols
T_HwP.DAu.D <- T_HwP.DAvA.D <- T_HwP.DAvP.D <- T_HwP.DAvC.D <- 
  T_HwP.DAwA.D <- T_HwP.DAwP.D <- T_HwP.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HwP.D <- cbind(T_HwP.DEu.D, T_HwP.DEvA.D, T_HwP.DEvP.D, T_HwP.DEvC.D, T_HwP.DEwA.D, T_HwP.DEwP.D, T_HwP.DEwC.D,
                T_HwP.DAu.D, T_HwP.DAvA.D, T_HwP.DAvP.D, T_HwP.DAvC.D, T_HwP.DAwA.D, T_HwP.DAwP.D, T_HwP.DAwC.D, 
                T_HwP.DIu.D, T_HwP.DIvA.D, T_HwP.DIvP.D, T_HwP.DIvC.D, T_HwP.DIwA.D, T_HwP.DIwP.D, T_HwP.DIwC.D, 
                T_HwP.DHu.D, T_HwP.DHvA.D, T_HwP.DHvP.D, T_HwP.DHvC.D, T_HwP.DHwA.D, T_HwP.DHwP.D, T_HwP.DHwC.D)
                
## HwC
# HwC.DxE cols
T_HwC.DEwC.D <- diag(ihr_wC.D*gamma, ncol = age.bins, nrow = age.bins)
T_HwC.DEu.D <- T_HwC.DEvA.D <- 
  T_HwC.DEvP.D <- T_HwC.DEwA.D <- T_HwC.DEwP.D <- T_HwC.DEvC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwC.DxH cols
T_HwC.DHwC.D <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HwC.DHu.D <- T_HwC.DHvA.D <- T_HwC.DHvP.D <- 
  T_HwC.DHwA.D <- T_HwC.DHwP.D <- T_HwC.DHvC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwC.DxI cols
T_HwC.DIu.D <- T_HwC.DIvA.D <- T_HwC.DIvP.D <- T_HwC.DIvC.D <- 
  T_HwC.DIwA.D <- T_HwC.DIwP.D <- T_HwC.DIwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwC.DxA cols
T_HwC.DAu.D <- T_HwC.DAvA.D <- T_HwC.DAvP.D <- T_HwC.DAvC.D <- 
  T_HwC.DAwA.D <- T_HwC.DAwP.D <- T_HwC.DAwC.D <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HwC.D <- cbind(T_HwC.DEu.D, T_HwC.DEvA.D, T_HwC.DEvP.D, T_HwC.DEvC.D, T_HwC.DEwA.D, T_HwC.DEwP.D, T_HwC.DEwC.D,
                T_HwC.DAu.D, T_HwC.DAvA.D, T_HwC.DAvP.D, T_HwC.DAvC.D, T_HwC.DAwA.D, T_HwC.DAwP.D, T_HwC.DAwC.D, 
                T_HwC.DIu.D, T_HwC.DIvA.D, T_HwC.DIvP.D, T_HwC.DIvC.D, T_HwC.DIwA.D, T_HwC.DIwP.D, T_HwC.DIwC.D, 
                T_HwC.DHu.D, T_HwC.DHvA.D, T_HwC.DHvP.D, T_HwC.DHvC.D, T_HwC.DHwA.D, T_HwC.DHwP.D, T_HwC.DHwC.D)
                
#all H's together
T_H.D <- rbind(T_Hu.D, T_HvA.D, T_HvP.D, T_HvC.D, T_HwA.D, T_HwP.D, T_HwC.D)


# final matrix
T_mat <- rbind(T_E.D, T_A.D, T_I.D, T_H.D)