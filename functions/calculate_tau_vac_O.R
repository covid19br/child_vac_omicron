#### T matrix #### writing in the lazy way because of readability
# # E line
# T_EE <- diag(1-gamma, nrow = 16, ncol = 16)
# T_EI <- T_EA <- T_EH <- matrix(0, ncol = 16, nrow = 16)
# T_E <- cbind(T_EE, T_EI, T_EA, T_EH)
## Eu
# Eu.OxE cols
T_Eu.OEu.O <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_Eu.OEvA.O <- T_Eu.OEvP.O <- T_Eu.OEvC.O <- 
  T_Eu.OEwA.O <- T_Eu.OEwP.O <- T_Eu.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EuDxI cols
T_Eu.OIu.O <- T_Eu.OIvA.O <- T_Eu.OIvP.O <- T_Eu.OIvC.O <- 
  T_Eu.OIwA.O <- T_Eu.OIwP.O <- T_Eu.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EuDxA cols
T_Eu.OAu.O <- T_Eu.OAvA.O <- T_Eu.OAvP.O <- T_Eu.OAvC.O <- 
  T_Eu.OAwA.O <- T_Eu.OAwP.O <- T_Eu.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EuDxH cols
T_Eu.OHu.O <- T_Eu.OHvA.O <- T_Eu.OHvP.O <- T_Eu.OHvC.O <- 
  T_Eu.OHwA.O <- T_Eu.OHwP.O <- T_Eu.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_Eu.O <- cbind(T_Eu.OEu.O, T_Eu.OEvA.O, T_Eu.OEvP.O, T_Eu.OEvC.O, T_Eu.OEwA.O, T_Eu.OEwP.O, T_Eu.OEwC.O,
                T_Eu.OAu.O, T_Eu.OAvA.O, T_Eu.OAvP.O, T_Eu.OAvC.O, T_Eu.OAwA.O, T_Eu.OAwP.O, T_Eu.OAwC.O, 
                T_Eu.OIu.O, T_Eu.OIvA.O, T_Eu.OIvP.O, T_Eu.OIvC.O, T_Eu.OIwA.O, T_Eu.OIwP.O, T_Eu.OIwC.O, 
                T_Eu.OHu.O, T_Eu.OHvA.O, T_Eu.OHvP.O, T_Eu.OHvC.O, T_Eu.OHwA.O, T_Eu.OHwP.O, T_Eu.OHwC.O)
## EvA
# EvA.OxE cols
T_EvA.OEvA.O <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EvA.OEu.O <- T_EvA.OEvP.O <- T_EvA.OEvC.O <- 
  T_EvA.OEwA.O <- T_EvA.OEwP.O <- T_EvA.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvA.OxA cols
T_EvA.OAu.O <- T_EvA.OAvA.O <- T_EvA.OAvP.O <- T_EvA.OAvC.O <- 
  T_EvA.OAwA.O <- T_EvA.OAwP.O <- T_EvA.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvA.OxI cols
T_EvA.OIu.O <- T_EvA.OIvA.O <- T_EvA.OIvP.O <- T_EvA.OIvC.O <- 
  T_EvA.OIwA.O <- T_EvA.OIwP.O <- T_EvA.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvAxH cols
T_EvA.OHu.O <- T_EvA.OHvA.O <- T_EvA.OHvP.O <- T_EvA.OHvC.O <- 
  T_EvA.OHwA.O <- T_EvA.OHwP.O <- T_EvA.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EvA.O <- cbind(T_EvA.OEu.O, T_EvA.OEvA.O, T_EvA.OEvP.O, T_EvA.OEvC.O, T_EvA.OEwA.O, T_EvA.OEwP.O, T_EvA.OEwC.O,
                T_EvA.OAu.O, T_EvA.OAvA.O, T_EvA.OAvP.O, T_EvA.OAvC.O, T_EvA.OAwA.O, T_EvA.OAwP.O, T_EvA.OAwC.O, 
                T_EvA.OIu.O, T_EvA.OIvA.O, T_EvA.OIvP.O, T_EvA.OIvC.O, T_EvA.OIwA.O, T_EvA.OIwP.O, T_EvA.OIwC.O, 
                T_EvA.OHu.O, T_EvA.OHvA.O, T_EvA.OHvP.O, T_EvA.OHvC.O, T_EvA.OHwA.O, T_EvA.OHwP.O, T_EvA.OHwC.O)

## EvP
# EvP.OxE cols
T_EvP.OEvP.O <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EvP.OEu.O <- T_EvP.OEvA.O <- T_EvP.OEvC.O <- 
  T_EvP.OEwA.O <- T_EvP.OEwP.O <- T_EvP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvP.OxA cols
T_EvP.OAu.O <- T_EvP.OAvA.O <- T_EvP.OAvP.O <- T_EvP.OAvC.O <- 
  T_EvP.OAwA.O <- T_EvP.OAwP.O <- T_EvP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvP.OxI cols
T_EvP.OIu.O <- T_EvP.OIvA.O <- T_EvP.OIvP.O <- T_EvP.OIvC.O <- 
  T_EvP.OIwA.O <- T_EvP.OIwP.O <- T_EvP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvPxH cols
T_EvP.OHu.O <- T_EvP.OHvA.O <- T_EvP.OHvP.O <- T_EvP.OHvC.O <- 
  T_EvP.OHwA.O <- T_EvP.OHwP.O <- T_EvP.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EvP.O <- cbind(T_EvP.OEu.O, T_EvP.OEvA.O, T_EvP.OEvP.O, T_EvP.OEvC.O, T_EvP.OEwA.O, T_EvP.OEwP.O, T_EvP.OEwC.O,
                T_EvP.OAu.O, T_EvP.OAvA.O, T_EvP.OAvP.O, T_EvP.OAvC.O, T_EvP.OAwA.O, T_EvP.OAwP.O, T_EvP.OAwC.O, 
                T_EvP.OIu.O, T_EvP.OIvA.O, T_EvP.OIvP.O, T_EvP.OIvC.O, T_EvP.OIwA.O, T_EvP.OIwP.O, T_EvP.OIwC.O, 
                T_EvP.OHu.O, T_EvP.OHvA.O, T_EvP.OHvP.O, T_EvP.OHvC.O, T_EvP.OHwA.O, T_EvP.OHwP.O, T_EvP.OHwC.O)

## EvC
# EvC.OxE cols
T_EvC.OEvC.O <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EvC.OEu.O <- T_EvC.OEvA.O <- T_EvC.OEvP.O <- 
  T_EvC.OEwA.O <- T_EvC.OEwP.O <- T_EvC.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvC.OxA cols
T_EvC.OAu.O <- T_EvC.OAvA.O <- T_EvC.OAvP.O <- T_EvC.OAvC.O <- 
  T_EvC.OAwA.O <- T_EvC.OAwP.O <- T_EvC.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvC.OxI cols
T_EvC.OIu.O <- T_EvC.OIvA.O <- T_EvC.OIvP.O <- T_EvC.OIvC.O <- 
  T_EvC.OIwA.O <- T_EvC.OIwP.O <- T_EvC.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EvCxH cols
T_EvC.OHu.O <- T_EvC.OHvA.O <- T_EvC.OHvP.O <- T_EvC.OHvC.O <- 
  T_EvC.OHwA.O <- T_EvC.OHwP.O <- T_EvC.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EvC.O <- cbind(T_EvC.OEu.O, T_EvC.OEvA.O, T_EvC.OEvP.O, T_EvC.OEvC.O, T_EvC.OEwA.O, T_EvC.OEwP.O, T_EvC.OEwC.O,
                T_EvC.OAu.O, T_EvC.OAvA.O, T_EvC.OAvP.O, T_EvC.OAvC.O, T_EvC.OAwA.O, T_EvC.OAwP.O, T_EvC.OAwC.O, 
                T_EvC.OIu.O, T_EvC.OIvA.O, T_EvC.OIvP.O, T_EvC.OIvC.O, T_EvC.OIwA.O, T_EvC.OIwP.O, T_EvC.OIwC.O, 
                T_EvC.OHu.O, T_EvC.OHvA.O, T_EvC.OHvP.O, T_EvC.OHvC.O, T_EvC.OHwA.O, T_EvC.OHwP.O, T_EvC.OHwC.O)

## EwA
# EwA.OxE cols
T_EwA.OEwA.O <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EwA.OEu.O <- T_EwA.OEvA.O <- T_EwA.OEvP.O <- 
  T_EwA.OEvC.O <- T_EwA.OEwP.O <- T_EwA.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwA.OxA cols
T_EwA.OAu.O <- T_EwA.OAvA.O <- T_EwA.OAvP.O <- T_EwA.OAvC.O <- 
  T_EwA.OAwA.O <- T_EwA.OAwP.O <- T_EwA.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwA.OxI cols
T_EwA.OIu.O <- T_EwA.OIvA.O <- T_EwA.OIvP.O <- T_EwA.OIvC.O <- 
  T_EwA.OIwA.O <- T_EwA.OIwP.O <- T_EwA.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwAxH cols
T_EwA.OHu.O <- T_EwA.OHvA.O <- T_EwA.OHvP.O <- T_EwA.OHvC.O <- 
  T_EwA.OHwA.O <- T_EwA.OHwP.O <- T_EwA.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EwA.O <- cbind(T_EwA.OEu.O, T_EwA.OEvA.O, T_EwA.OEvP.O, T_EwA.OEvC.O, T_EwA.OEwA.O, T_EwA.OEwP.O, T_EwA.OEwC.O,
                T_EwA.OAu.O, T_EwA.OAvA.O, T_EwA.OAvP.O, T_EwA.OAvC.O, T_EwA.OAwA.O, T_EwA.OAwP.O, T_EwA.OAwC.O, 
                T_EwA.OIu.O, T_EwA.OIvA.O, T_EwA.OIvP.O, T_EwA.OIvC.O, T_EwA.OIwA.O, T_EwA.OIwP.O, T_EwA.OIwC.O, 
                T_EwA.OHu.O, T_EwA.OHvA.O, T_EwA.OHvP.O, T_EwA.OHvC.O, T_EwA.OHwA.O, T_EwA.OHwP.O, T_EwA.OHwC.O)

## EwP
# EwP.OxE cols
T_EwP.OEwP.O <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EwP.OEu.O <- T_EwP.OEvA.O <- T_EwP.OEvP.O <- 
  T_EwP.OEvC.O <- T_EwP.OEwA.O <- T_EwP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwP.OxA cols
T_EwP.OAu.O <- T_EwP.OAvA.O <- T_EwP.OAvP.O <- T_EwP.OAvC.O <- 
  T_EwP.OAwA.O <- T_EwP.OAwP.O <- T_EwP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwP.OxI cols
T_EwP.OIu.O <- T_EwP.OIvA.O <- T_EwP.OIvP.O <- T_EwP.OIvC.O <- 
  T_EwP.OIwA.O <- T_EwP.OIwP.O <- T_EwP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwPxH cols
T_EwP.OHu.O <- T_EwP.OHvA.O <- T_EwP.OHvP.O <- T_EwP.OHvC.O <- 
  T_EwP.OHwA.O <- T_EwP.OHwP.O <- T_EwP.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EwP.O <- cbind(T_EwP.OEu.O, T_EwP.OEvA.O, T_EwP.OEvP.O, T_EwP.OEvC.O, T_EwP.OEwA.O, T_EwP.OEwP.O, T_EwP.OEwC.O,
                T_EwP.OAu.O, T_EwP.OAvA.O, T_EwP.OAvP.O, T_EwP.OAvC.O, T_EwP.OAwA.O, T_EwP.OAwP.O, T_EwP.OAwC.O, 
                T_EwP.OIu.O, T_EwP.OIvA.O, T_EwP.OIvP.O, T_EwP.OIvC.O, T_EwP.OIwA.O, T_EwP.OIwP.O, T_EwP.OIwC.O, 
                T_EwP.OHu.O, T_EwP.OHvA.O, T_EwP.OHvP.O, T_EwP.OHvC.O, T_EwP.OHwA.O, T_EwP.OHwP.O, T_EwP.OHwC.O)

## EwC
# EwC.OxE cols
T_EwC.OEwC.O <- diag(1-gamma, nrow = age.bins, ncol = age.bins)
T_EwC.OEu.O <- T_EwC.OEvA.O <- T_EwC.OEvP.O <- 
  T_EwC.OEvC.O <- T_EwC.OEwA.O <- T_EwC.OEwP.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwC.OxA cols
T_EwC.OAu.O <- T_EwC.OAvA.O <- T_EwC.OAvP.O <- T_EwC.OAvC.O <- 
  T_EwC.OAwA.O <- T_EwC.OAwP.O <- T_EwC.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwC.OxI cols
T_EwC.OIu.O <- T_EwC.OIvA.O <- T_EwC.OIvP.O <- T_EwC.OIvC.O <- 
  T_EwC.OIwA.O <- T_EwC.OIwP.O <- T_EwC.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# EwC.OxH cols
T_EwC.OHu.O <- T_EwC.OHvA.O <- T_EwC.OHvP.O <- T_EwC.OHvC.O <- 
  T_EwC.OHwA.O <- T_EwC.OHwP.O <- T_EwC.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_EwC.O <- cbind(T_EwC.OEu.O, T_EwC.OEvA.O, T_EwC.OEvP.O, T_EwC.OEvC.O, T_EwC.OEwA.O, T_EwC.OEwP.O, T_EwC.OEwC.O,
                T_EwC.OAu.O, T_EwC.OAvA.O, T_EwC.OAvP.O, T_EwC.OAvC.O, T_EwC.OAwA.O, T_EwC.OAwP.O, T_EwC.OAwC.O, 
                T_EwC.OIu.O, T_EwC.OIvA.O, T_EwC.OIvP.O, T_EwC.OIvC.O, T_EwC.OIwA.O, T_EwC.OIwP.O, T_EwC.OIwC.O, 
                T_EwC.OHu.O, T_EwC.OHvA.O, T_EwC.OHvP.O, T_EwC.OHvC.O, T_EwC.OHwA.O, T_EwC.OHwP.O, T_EwC.OHwC.O)

#all E's together
T_E.O <- rbind(T_Eu.O, T_EvA.O, T_EvP.O, T_EvC.O, T_EwA.O, T_EwP.O, T_EwC.O)

# # A line
# T_AE <- diag((1-ihr)*gamma*asymp)
# T_AA <- diag(1-nu, nrow = 16, ncol = 16)
# T_AI <- T_AH <- matrix(0, nrow = 16, ncol = 16)
# T_A <- cbind(T_AE, T_AI, T_AA, T_AH)

## Au
# Au.OxE cols
T_Au.OEu.O <- diag((1-ihr)*gamma*asymp, ncol = age.bins, nrow = age.bins)
T_Au.OEvA.O <- T_Au.OEvP.O <- 
  T_Au.OEvC.O <- T_Au.OEwA.O <- T_Au.OEwP.O <- T_Au.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Au.OxA cols
T_Au.OAu.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_Au.OAvA.O <- T_Au.OAvP.O <- T_Au.OAvC.O <- 
  T_Au.OAwA.O <- T_Au.OAwP.O <- T_Au.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Au.OxI cols
T_Au.OIu.O <- T_Au.OIvA.O <- T_Au.OIvP.O <- T_Au.OIvC.O <- 
  T_Au.OIwA.O <- T_Au.OIwP.O <- T_Au.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Au.OxH cols
T_Au.OHu.O <- T_Au.OHvA.O <- T_Au.OHvP.O <- T_Au.OHvC.O <- 
  T_Au.OHwA.O <- T_Au.OHwP.O <- T_Au.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_Au.O <- cbind(T_Au.OEu.O, T_Au.OEvA.O, T_Au.OEvP.O, T_Au.OEvC.O, T_Au.OEwA.O, T_Au.OEwP.O, T_Au.OEwC.O,
                T_Au.OAu.O, T_Au.OAvA.O, T_Au.OAvP.O, T_Au.OAvC.O, T_Au.OAwA.O, T_Au.OAwP.O, T_Au.OAwC.O, 
                T_Au.OIu.O, T_Au.OIvA.O, T_Au.OIvP.O, T_Au.OIvC.O, T_Au.OIwA.O, T_Au.OIwP.O, T_Au.OIwC.O, 
                T_Au.OHu.O, T_Au.OHvA.O, T_Au.OHvP.O, T_Au.OHvC.O, T_Au.OHwA.O, T_Au.OHwP.O, T_Au.OHwC.O)
                
## AvA
# AvA.OxE cols
T_AvA.OEvA.O <- diag((1-ihr_vA.O)*gamma*asymp_vA.O, ncol = age.bins, nrow = age.bins)
T_AvA.OEu.O <- T_AvA.OEvP.O <- 
  T_AvA.OEvC.O <- T_AvA.OEwA.O <- T_AvA.OEwP.O <- T_AvA.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvA.OxA cols
T_AvA.OAvA.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AvA.OAu.O <- T_AvA.OAvP.O <- T_AvA.OAvC.O <- 
  T_AvA.OAwA.O <- T_AvA.OAwP.O <- T_AvA.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvA.OxI cols
T_AvA.OIu.O <- T_AvA.OIvA.O <- T_AvA.OIvP.O <- T_AvA.OIvC.O <- 
  T_AvA.OIwA.O <- T_AvA.OIwP.O <- T_AvA.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvA.OxH cols
T_AvA.OHu.O <- T_AvA.OHvA.O <- T_AvA.OHvP.O <- T_AvA.OHvC.O <- 
  T_AvA.OHwA.O <- T_AvA.OHwP.O <- T_AvA.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AvA.O <- cbind(T_AvA.OEu.O, T_AvA.OEvA.O, T_AvA.OEvP.O, T_AvA.OEvC.O, T_AvA.OEwA.O, T_AvA.OEwP.O, T_AvA.OEwC.O,
                T_AvA.OAu.O, T_AvA.OAvA.O, T_AvA.OAvP.O, T_AvA.OAvC.O, T_AvA.OAwA.O, T_AvA.OAwP.O, T_AvA.OAwC.O, 
                T_AvA.OIu.O, T_AvA.OIvA.O, T_AvA.OIvP.O, T_AvA.OIvC.O, T_AvA.OIwA.O, T_AvA.OIwP.O, T_AvA.OIwC.O, 
                T_AvA.OHu.O, T_AvA.OHvA.O, T_AvA.OHvP.O, T_AvA.OHvC.O, T_AvA.OHwA.O, T_AvA.OHwP.O, T_AvA.OHwC.O)

## AvP
# AvP.OxE cols
T_AvP.OEvP.O <- diag((1-ihr_vP.O)*gamma*asymp_vP.O, ncol = age.bins, nrow = age.bins)
T_AvP.OEu.O <- T_AvP.OEvA.O <- 
  T_AvP.OEvC.O <- T_AvP.OEwA.O <- T_AvP.OEwP.O <- T_AvP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvP.OxA cols
T_AvP.OAvP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AvP.OAu.O <- T_AvP.OAvA.O <- T_AvP.OAvC.O <- 
  T_AvP.OAwA.O <- T_AvP.OAwP.O <- T_AvP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvP.OxI cols
T_AvP.OIu.O <- T_AvP.OIvA.O <- T_AvP.OIvP.O <- T_AvP.OIvC.O <- 
  T_AvP.OIwA.O <- T_AvP.OIwP.O <- T_AvP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvP.OxH cols
T_AvP.OHu.O <- T_AvP.OHvA.O <- T_AvP.OHvP.O <- T_AvP.OHvC.O <- 
  T_AvP.OHwA.O <- T_AvP.OHwP.O <- T_AvP.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AvP.O <- cbind(T_AvP.OEu.O, T_AvP.OEvA.O, T_AvP.OEvP.O, T_AvP.OEvC.O, T_AvP.OEwA.O, T_AvP.OEwP.O, T_AvP.OEwC.O,
                T_AvP.OAu.O, T_AvP.OAvA.O, T_AvP.OAvP.O, T_AvP.OAvC.O, T_AvP.OAwA.O, T_AvP.OAwP.O, T_AvP.OAwC.O, 
                T_AvP.OIu.O, T_AvP.OIvA.O, T_AvP.OIvP.O, T_AvP.OIvC.O, T_AvP.OIwA.O, T_AvP.OIwP.O, T_AvP.OIwC.O, 
                T_AvP.OHu.O, T_AvP.OHvA.O, T_AvP.OHvP.O, T_AvP.OHvC.O, T_AvP.OHwA.O, T_AvP.OHwP.O, T_AvP.OHwC.O)
                
## AvC
# AvC.OxE cols
T_AvC.OEvC.O <- diag((1-ihr_vC.O)*gamma*asymp_vC.O, ncol = age.bins, nrow = age.bins)
T_AvC.OEu.O <- T_AvC.OEvA.O <- 
  T_AvC.OEvP.O <- T_AvC.OEwA.O <- T_AvC.OEwP.O <- T_AvC.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvC.OxA cols
T_AvC.OAvC.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AvC.OAu.O <- T_AvC.OAvA.O <- T_AvC.OAvP.O <- 
  T_AvC.OAwA.O <- T_AvC.OAwP.O <- T_AvC.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvC.OxI cols
T_AvC.OIu.O <- T_AvC.OIvA.O <- T_AvC.OIvP.O <- T_AvC.OIvC.O <- 
  T_AvC.OIwA.O <- T_AvC.OIwP.O <- T_AvC.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AvC.OxH cols
T_AvC.OHu.O <- T_AvC.OHvA.O <- T_AvC.OHvP.O <- T_AvC.OHvC.O <- 
  T_AvC.OHwA.O <- T_AvC.OHwP.O <- T_AvC.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AvC.O <- cbind(T_AvC.OEu.O, T_AvC.OEvA.O, T_AvC.OEvP.O, T_AvC.OEvC.O, T_AvC.OEwA.O, T_AvC.OEwP.O, T_AvC.OEwC.O,
                T_AvC.OAu.O, T_AvC.OAvA.O, T_AvC.OAvP.O, T_AvC.OAvC.O, T_AvC.OAwA.O, T_AvC.OAwP.O, T_AvC.OAwC.O, 
                T_AvC.OIu.O, T_AvC.OIvA.O, T_AvC.OIvP.O, T_AvC.OIvC.O, T_AvC.OIwA.O, T_AvC.OIwP.O, T_AvC.OIwC.O, 
                T_AvC.OHu.O, T_AvC.OHvA.O, T_AvC.OHvP.O, T_AvC.OHvC.O, T_AvC.OHwA.O, T_AvC.OHwP.O, T_AvC.OHwC.O)

## AwA
# AwA.OxE cols
T_AwA.OEwA.O <- diag((1-ihr_wA.O)*gamma*asymp_wA.O, ncol = age.bins, nrow = age.bins)
T_AwA.OEu.O <- T_AwA.OEvP.O <- 
  T_AwA.OEvC.O <- T_AwA.OEvA.O <- T_AwA.OEwP.O <- T_AwA.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwA.OxA cols
T_AwA.OAwA.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AwA.OAu.O <- T_AwA.OAvP.O <- T_AwA.OAvC.O <- 
  T_AwA.OAvA.O <- T_AwA.OAwP.O <- T_AwA.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwA.OxI cols
T_AwA.OIu.O <- T_AwA.OIvA.O <- T_AwA.OIvP.O <- T_AwA.OIvC.O <- 
  T_AwA.OIwA.O <- T_AwA.OIwP.O <- T_AwA.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwA.OxH cols
T_AwA.OHu.O <- T_AwA.OHvA.O <- T_AwA.OHvP.O <- T_AwA.OHvC.O <- 
  T_AwA.OHwA.O <- T_AwA.OHwP.O <- T_AwA.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AwA.O <- cbind(T_AwA.OEu.O, T_AwA.OEvA.O, T_AwA.OEvP.O, T_AwA.OEvC.O, T_AwA.OEwA.O, T_AwA.OEwP.O, T_AwA.OEwC.O,
                T_AwA.OAu.O, T_AwA.OAvA.O, T_AwA.OAvP.O, T_AwA.OAvC.O, T_AwA.OAwA.O, T_AwA.OAwP.O, T_AwA.OAwC.O, 
                T_AwA.OIu.O, T_AwA.OIvA.O, T_AwA.OIvP.O, T_AwA.OIvC.O, T_AwA.OIwA.O, T_AwA.OIwP.O, T_AwA.OIwC.O, 
                T_AwA.OHu.O, T_AwA.OHvA.O, T_AwA.OHvP.O, T_AwA.OHvC.O, T_AwA.OHwA.O, T_AwA.OHwP.O, T_AwA.OHwC.O)

## AwP
# AwP.OxE cols
T_AwP.OEwP.O <- diag((1-ihr_wP.O)*gamma*asymp_wP.O, ncol = age.bins, nrow = age.bins)
T_AwP.OEu.O <- T_AwP.OEvA.O <- 
  T_AwP.OEvC.O <- T_AwP.OEwA.O <- T_AwP.OEvP.O <- T_AwP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwP.OxA cols
T_AwP.OAwP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AwP.OAu.O <- T_AwP.OAvA.O <- T_AwP.OAvC.O <- 
  T_AwP.OAwA.O <- T_AwP.OAvP.O <- T_AwP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwP.OxI cols
T_AwP.OIu.O <- T_AwP.OIvA.O <- T_AwP.OIvP.O <- T_AwP.OIvC.O <- 
  T_AwP.OIwA.O <- T_AwP.OIwP.O <- T_AwP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwP.OxH cols
T_AwP.OHu.O <- T_AwP.OHvA.O <- T_AwP.OHvP.O <- T_AwP.OHvC.O <- 
  T_AwP.OHwA.O <- T_AwP.OHwP.O <- T_AwP.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AwP.O <- cbind(T_AwP.OEu.O, T_AwP.OEvA.O, T_AwP.OEvP.O, T_AwP.OEvC.O, T_AwP.OEwA.O, T_AwP.OEwP.O, T_AwP.OEwC.O,
                T_AwP.OAu.O, T_AwP.OAvA.O, T_AwP.OAvP.O, T_AwP.OAvC.O, T_AwP.OAwA.O, T_AwP.OAwP.O, T_AwP.OAwC.O, 
                T_AwP.OIu.O, T_AwP.OIvA.O, T_AwP.OIvP.O, T_AwP.OIvC.O, T_AwP.OIwA.O, T_AwP.OIwP.O, T_AwP.OIwC.O, 
                T_AwP.OHu.O, T_AwP.OHvA.O, T_AwP.OHvP.O, T_AwP.OHvC.O, T_AwP.OHwA.O, T_AwP.OHwP.O, T_AwP.OHwC.O)
                
## AwC
# AwC.OxE cols
T_AwC.OEwC.O <- diag((1-ihr_wC.O)*gamma*asymp_wC.O, ncol = age.bins, nrow = age.bins)
T_AwC.OEu.O <- T_AwC.OEvA.O <- 
  T_AwC.OEvP.O <- T_AwC.OEwA.O <- T_AwC.OEwP.O <- T_AwC.OEvC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwC.OxA cols
T_AwC.OAwC.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_AwC.OAu.O <- T_AwC.OAvA.O <- T_AwC.OAvP.O <- 
  T_AwC.OAwA.O <- T_AwC.OAwP.O <- T_AwC.OAvC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwC.OxI cols
T_AwC.OIu.O <- T_AwC.OIvA.O <- T_AwC.OIvP.O <- T_AwC.OIvC.O <- 
  T_AwC.OIwA.O <- T_AwC.OIwP.O <- T_AwC.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# AwC.OxH cols
T_AwC.OHu.O <- T_AwC.OHvA.O <- T_AwC.OHvP.O <- T_AwC.OHvC.O <- 
  T_AwC.OHwA.O <- T_AwC.OHwP.O <- T_AwC.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_AwC.O <- cbind(T_AwC.OEu.O, T_AwC.OEvA.O, T_AwC.OEvP.O, T_AwC.OEvC.O, T_AwC.OEwA.O, T_AwC.OEwP.O, T_AwC.OEwC.O,
                T_AwC.OAu.O, T_AwC.OAvA.O, T_AwC.OAvP.O, T_AwC.OAvC.O, T_AwC.OAwA.O, T_AwC.OAwP.O, T_AwC.OAwC.O, 
                T_AwC.OIu.O, T_AwC.OIvA.O, T_AwC.OIvP.O, T_AwC.OIvC.O, T_AwC.OIwA.O, T_AwC.OIwP.O, T_AwC.OIwC.O, 
                T_AwC.OHu.O, T_AwC.OHvA.O, T_AwC.OHvP.O, T_AwC.OHvC.O, T_AwC.OHwA.O, T_AwC.OHwP.O, T_AwC.OHwC.O)

#all A's together
T_A.O <- rbind(T_Au.O, T_AvA.O, T_AvP.O, T_AvC.O, T_AwA.O, T_AwP.O, T_AwC.O)


# # I line
# T_IE <- diag((1-ihr)*gamma*(1-asymp))
# T_II <- diag(1-nu, nrow = 16, ncol = 16)
# T_IA <- T_IH <- matrix(0, nrow = 16, ncol = 16)
# T_I <- cbind(T_IE, T_II, T_IA, T_IH)

## Iu
# Iu.OxE cols
T_Iu.OEu.O <- diag((1-ihr)*gamma*(1-asymp), ncol = age.bins, nrow = age.bins)
T_Iu.OEvA.O <- T_Iu.OEvP.O <- 
  T_Iu.OEvC.O <- T_Iu.OEwA.O <- T_Iu.OEwP.O <- T_Iu.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Iu.OxI cols
T_Iu.OIu.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_Iu.OIvA.O <- T_Iu.OIvP.O <- T_Iu.OIvC.O <- 
  T_Iu.OIwA.O <- T_Iu.OIwP.O <- T_Iu.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Iu.OxA cols
T_Iu.OAu.O <- T_Iu.OAvA.O <- T_Iu.OAvP.O <- T_Iu.OAvC.O <- 
  T_Iu.OAwA.O <- T_Iu.OAwP.O <- T_Iu.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Iu.OxH cols
T_Iu.OHu.O <- T_Iu.OHvA.O <- T_Iu.OHvP.O <- T_Iu.OHvC.O <- 
  T_Iu.OHwA.O <- T_Iu.OHwP.O <- T_Iu.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_Iu.O <- cbind(T_Iu.OEu.O, T_Iu.OEvA.O, T_Iu.OEvP.O, T_Iu.OEvC.O, T_Iu.OEwA.O, T_Iu.OEwP.O, T_Iu.OEwC.O,
                T_Iu.OAu.O, T_Iu.OAvA.O, T_Iu.OAvP.O, T_Iu.OAvC.O, T_Iu.OAwA.O, T_Iu.OAwP.O, T_Iu.OAwC.O, 
                T_Iu.OIu.O, T_Iu.OIvA.O, T_Iu.OIvP.O, T_Iu.OIvC.O, T_Iu.OIwA.O, T_Iu.OIwP.O, T_Iu.OIwC.O, 
                T_Iu.OHu.O, T_Iu.OHvA.O, T_Iu.OHvP.O, T_Iu.OHvC.O, T_Iu.OHwA.O, T_Iu.OHwP.O, T_Iu.OHwC.O)
                
## IvA
# IvA.OxE cols
T_IvA.OEvA.O <- diag((1-ihr_vA.O)*gamma*(1-asymp_vA.O), ncol = age.bins, nrow = age.bins)
T_IvA.OEu.O <- T_IvA.OEvP.O <- 
  T_IvA.OEvC.O <- T_IvA.OEwA.O <- T_IvA.OEwP.O <- T_IvA.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvA.OxI cols
T_IvA.OIvA.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IvA.OIu.O <- T_IvA.OIvP.O <- T_IvA.OIvC.O <- 
  T_IvA.OIwA.O <- T_IvA.OIwP.O <- T_IvA.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvA.OxA cols
T_IvA.OAu.O <- T_IvA.OAvA.O <- T_IvA.OAvP.O <- T_IvA.OAvC.O <- 
  T_IvA.OAwA.O <- T_IvA.OAwP.O <- T_IvA.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvA.OxH cols
T_IvA.OHu.O <- T_IvA.OHvA.O <- T_IvA.OHvP.O <- T_IvA.OHvC.O <- 
  T_IvA.OHwA.O <- T_IvA.OHwP.O <- T_IvA.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IvA.O <- cbind(T_IvA.OEu.O, T_IvA.OEvA.O, T_IvA.OEvP.O, T_IvA.OEvC.O, T_IvA.OEwA.O, T_IvA.OEwP.O, T_IvA.OEwC.O,
                T_IvA.OAu.O, T_IvA.OAvA.O, T_IvA.OAvP.O, T_IvA.OAvC.O, T_IvA.OAwA.O, T_IvA.OAwP.O, T_IvA.OAwC.O, 
                T_IvA.OIu.O, T_IvA.OIvA.O, T_IvA.OIvP.O, T_IvA.OIvC.O, T_IvA.OIwA.O, T_IvA.OIwP.O, T_IvA.OIwC.O, 
                T_IvA.OHu.O, T_IvA.OHvA.O, T_IvA.OHvP.O, T_IvA.OHvC.O, T_IvA.OHwA.O, T_IvA.OHwP.O, T_IvA.OHwC.O)

## IvP
# IvP.OxE cols
T_IvP.OEvP.O <- diag((1-ihr_vP.O)*gamma*(1-asymp_vP.O), ncol = age.bins, nrow = age.bins)
T_IvP.OEu.O <- T_IvP.OEvA.O <- 
  T_IvP.OEvC.O <- T_IvP.OEwA.O <- T_IvP.OEwP.O <- T_IvP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvP.OxI cols
T_IvP.OIvP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IvP.OIu.O <- T_IvP.OIvA.O <- T_IvP.OIvC.O <- 
  T_IvP.OIwA.O <- T_IvP.OIwP.O <- T_IvP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvP.OxA cols
T_IvP.OAu.O <- T_IvP.OAvA.O <- T_IvP.OAvP.O <- T_IvP.OAvC.O <- 
  T_IvP.OAwA.O <- T_IvP.OAwP.O <- T_IvP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvP.OxH cols
T_IvP.OHu.O <- T_IvP.OHvA.O <- T_IvP.OHvP.O <- T_IvP.OHvC.O <- 
  T_IvP.OHwA.O <- T_IvP.OHwP.O <- T_IvP.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IvP.O <- cbind(T_IvP.OEu.O, T_IvP.OEvA.O, T_IvP.OEvP.O, T_IvP.OEvC.O, T_IvP.OEwA.O, T_IvP.OEwP.O, T_IvP.OEwC.O,
                T_IvP.OAu.O, T_IvP.OAvA.O, T_IvP.OAvP.O, T_IvP.OAvC.O, T_IvP.OAwA.O, T_IvP.OAwP.O, T_IvP.OAwC.O, 
                T_IvP.OIu.O, T_IvP.OIvA.O, T_IvP.OIvP.O, T_IvP.OIvC.O, T_IvP.OIwA.O, T_IvP.OIwP.O, T_IvP.OIwC.O, 
                T_IvP.OHu.O, T_IvP.OHvA.O, T_IvP.OHvP.O, T_IvP.OHvC.O, T_IvP.OHwA.O, T_IvP.OHwP.O, T_IvP.OHwC.O)
                
## IvC
# IvC.OxE cols
T_IvC.OEvC.O <- diag((1-ihr_vC.O)*gamma*(1-asymp_vC.O), ncol = age.bins, nrow = age.bins)
T_IvC.OEu.O <- T_IvC.OEvA.O <- 
  T_IvC.OEvP.O <- T_IvC.OEwA.O <- T_IvC.OEwP.O <- T_IvC.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvC.OxI cols
T_IvC.OIvC.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IvC.OIu.O <- T_IvC.OIvA.O <- T_IvC.OIvP.O <- 
  T_IvC.OIwA.O <- T_IvC.OIwP.O <- T_IvC.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvC.OxA cols
T_IvC.OAu.O <- T_IvC.OAvA.O <- T_IvC.OAvP.O <- T_IvC.OAvC.O <- 
  T_IvC.OAwA.O <- T_IvC.OAwP.O <- T_IvC.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IvC.OxH cols
T_IvC.OHu.O <- T_IvC.OHvA.O <- T_IvC.OHvP.O <- T_IvC.OHvC.O <- 
  T_IvC.OHwA.O <- T_IvC.OHwP.O <- T_IvC.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IvC.O <- cbind(T_IvC.OEu.O, T_IvC.OEvA.O, T_IvC.OEvP.O, T_IvC.OEvC.O, T_IvC.OEwA.O, T_IvC.OEwP.O, T_IvC.OEwC.O,
                T_IvC.OAu.O, T_IvC.OAvA.O, T_IvC.OAvP.O, T_IvC.OAvC.O, T_IvC.OAwA.O, T_IvC.OAwP.O, T_IvC.OAwC.O, 
                T_IvC.OIu.O, T_IvC.OIvA.O, T_IvC.OIvP.O, T_IvC.OIvC.O, T_IvC.OIwA.O, T_IvC.OIwP.O, T_IvC.OIwC.O, 
                T_IvC.OHu.O, T_IvC.OHvA.O, T_IvC.OHvP.O, T_IvC.OHvC.O, T_IvC.OHwA.O, T_IvC.OHwP.O, T_IvC.OHwC.O)

## IwA
# IwA.OxE cols
T_IwA.OEwA.O <- diag((1-ihr_wA.O)*gamma*(1-asymp_wA.O), ncol = age.bins, nrow = age.bins)
T_IwA.OEu.O <- T_IwA.OEvP.O <- 
  T_IwA.OEvC.O <- T_IwA.OEvA.O <- T_IwA.OEwP.O <- T_IwA.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwA.OxI cols
T_IwA.OIwA.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IwA.OIu.O <- T_IwA.OIvP.O <- T_IwA.OIvC.O <- 
  T_IwA.OIvA.O <- T_IwA.OIwP.O <- T_IwA.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwA.OxA cols
T_IwA.OAu.O <- T_IwA.OAvA.O <- T_IwA.OAvP.O <- T_IwA.OAvC.O <- 
  T_IwA.OAwA.O <- T_IwA.OAwP.O <- T_IwA.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwA.OxH cols
T_IwA.OHu.O <- T_IwA.OHvA.O <- T_IwA.OHvP.O <- T_IwA.OHvC.O <- 
  T_IwA.OHwA.O <- T_IwA.OHwP.O <- T_IwA.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IwA.O <- cbind(T_IwA.OEu.O, T_IwA.OEvA.O, T_IwA.OEvP.O, T_IwA.OEvC.O, T_IwA.OEwA.O, T_IwA.OEwP.O, T_IwA.OEwC.O,
                T_IwA.OAu.O, T_IwA.OAvA.O, T_IwA.OAvP.O, T_IwA.OAvC.O, T_IwA.OAwA.O, T_IwA.OAwP.O, T_IwA.OAwC.O, 
                T_IwA.OIu.O, T_IwA.OIvA.O, T_IwA.OIvP.O, T_IwA.OIvC.O, T_IwA.OIwA.O, T_IwA.OIwP.O, T_IwA.OIwC.O, 
                T_IwA.OHu.O, T_IwA.OHvA.O, T_IwA.OHvP.O, T_IwA.OHvC.O, T_IwA.OHwA.O, T_IwA.OHwP.O, T_IwA.OHwC.O)

## IwP
# IwP.OxE cols
T_IwP.OEwP.O <- diag((1-ihr_wP.O)*gamma*(1-asymp_wP.O), ncol = age.bins, nrow = age.bins)
T_IwP.OEu.O <- T_IwP.OEvA.O <- 
  T_IwP.OEvC.O <- T_IwP.OEwA.O <- T_IwP.OEvP.O <- T_IwP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwP.OxI cols
T_IwP.OIwP.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IwP.OIu.O <- T_IwP.OIvA.O <- T_IwP.OIvC.O <- 
  T_IwP.OIwA.O <- T_IwP.OIvP.O <- T_IwP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwP.OxA cols
T_IwP.OAu.O <- T_IwP.OAvA.O <- T_IwP.OAvP.O <- T_IwP.OAvC.O <- 
  T_IwP.OAwA.O <- T_IwP.OAwP.O <- T_IwP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwP.OxH cols
T_IwP.OHu.O <- T_IwP.OHvA.O <- T_IwP.OHvP.O <- T_IwP.OHvC.O <- 
  T_IwP.OHwA.O <- T_IwP.OHwP.O <- T_IwP.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IwP.O <- cbind(T_IwP.OEu.O, T_IwP.OEvA.O, T_IwP.OEvP.O, T_IwP.OEvC.O, T_IwP.OEwA.O, T_IwP.OEwP.O, T_IwP.OEwC.O,
                T_IwP.OAu.O, T_IwP.OAvA.O, T_IwP.OAvP.O, T_IwP.OAvC.O, T_IwP.OAwA.O, T_IwP.OAwP.O, T_IwP.OAwC.O, 
                T_IwP.OIu.O, T_IwP.OIvA.O, T_IwP.OIvP.O, T_IwP.OIvC.O, T_IwP.OIwA.O, T_IwP.OIwP.O, T_IwP.OIwC.O, 
                T_IwP.OHu.O, T_IwP.OHvA.O, T_IwP.OHvP.O, T_IwP.OHvC.O, T_IwP.OHwA.O, T_IwP.OHwP.O, T_IwP.OHwC.O)
                
## IwC
# IwC.OxE cols
T_IwC.OEwC.O <- diag((1-ihr_wC.O)*gamma*(1-asymp_wC.O), ncol = age.bins, nrow = age.bins)
T_IwC.OEu.O <- T_IwC.OEvA.O <- 
  T_IwC.OEvP.O <- T_IwC.OEwA.O <- T_IwC.OEwP.O <- T_IwC.OEvC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwC.OxI cols
T_IwC.OIwC.O <- diag(1-nu, ncol = age.bins, nrow = age.bins)
T_IwC.OIu.O <- T_IwC.OIvA.O <- T_IwC.OIvP.O <- 
  T_IwC.OIwA.O <- T_IwC.OIwP.O <- T_IwC.OIvC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwC.OxA cols
T_IwC.OAu.O <- T_IwC.OAvA.O <- T_IwC.OAvP.O <- T_IwC.OAvC.O <- 
  T_IwC.OAwA.O <- T_IwC.OAwP.O <- T_IwC.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# IwC.OxH cols
T_IwC.OHu.O <- T_IwC.OHvA.O <- T_IwC.OHvP.O <- T_IwC.OHvC.O <- 
  T_IwC.OHwA.O <- T_IwC.OHwP.O <- T_IwC.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_IwC.O <- cbind(T_IwC.OEu.O, T_IwC.OEvA.O, T_IwC.OEvP.O, T_IwC.OEvC.O, T_IwC.OEwA.O, T_IwC.OEwP.O, T_IwC.OEwC.O,
                T_IwC.OAu.O, T_IwC.OAvA.O, T_IwC.OAvP.O, T_IwC.OAvC.O, T_IwC.OAwA.O, T_IwC.OAwP.O, T_IwC.OAwC.O, 
                T_IwC.OIu.O, T_IwC.OIvA.O, T_IwC.OIvP.O, T_IwC.OIvC.O, T_IwC.OIwA.O, T_IwC.OIwP.O, T_IwC.OIwC.O, 
                T_IwC.OHu.O, T_IwC.OHvA.O, T_IwC.OHvP.O, T_IwC.OHvC.O, T_IwC.OHwA.O, T_IwC.OHwP.O, T_IwC.OHwC.O)

#all I's together
T_I.O <- rbind(T_Iu.O, T_IvA.O, T_IvP.O, T_IvC.O, T_IwA.O, T_IwP.O, T_IwC.O)

# # H line
# T_HE <- diag(ihr*gamma)
# T_HH <- diag(1-nus, nrow = 16, ncol = 16)
# T_HI <- T_HA <- matrix(0, nrow = 16, ncol = 16)
# T_H <- cbind(T_HE, T_HI, T_HA, T_HH)

## Hu
# Hu.OxE cols
T_Hu.OEu.O <- diag(ihr*gamma, ncol = age.bins, nrow = age.bins)
T_Hu.OEvA.O <- T_Hu.OEvP.O <- 
  T_Hu.OEvC.O <- T_Hu.OEwA.O <- T_Hu.OEwP.O <- T_Hu.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Hu.OxH cols
T_Hu.OHu.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_Hu.OHvA.O <- T_Hu.OHvP.O <- T_Hu.OHvC.O <- 
  T_Hu.OHwA.O <- T_Hu.OHwP.O <- T_Hu.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Hu.OxI cols
T_Hu.OIu.O <- T_Hu.OIvA.O <- T_Hu.OIvP.O <- T_Hu.OIvC.O <- 
  T_Hu.OIwA.O <- T_Hu.OIwP.O <- T_Hu.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# Hu.OxA cols
T_Hu.OAu.O <- T_Hu.OAvA.O <- T_Hu.OAvP.O <- T_Hu.OAvC.O <- 
  T_Hu.OAwA.O <- T_Hu.OAwP.O <- T_Hu.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_Hu.O <- cbind(T_Hu.OEu.O, T_Hu.OEvA.O, T_Hu.OEvP.O, T_Hu.OEvC.O, T_Hu.OEwA.O, T_Hu.OEwP.O, T_Hu.OEwC.O,
                T_Hu.OAu.O, T_Hu.OAvA.O, T_Hu.OAvP.O, T_Hu.OAvC.O, T_Hu.OAwA.O, T_Hu.OAwP.O, T_Hu.OAwC.O, 
                T_Hu.OIu.O, T_Hu.OIvA.O, T_Hu.OIvP.O, T_Hu.OIvC.O, T_Hu.OIwA.O, T_Hu.OIwP.O, T_Hu.OIwC.O, 
                T_Hu.OHu.O, T_Hu.OHvA.O, T_Hu.OHvP.O, T_Hu.OHvC.O, T_Hu.OHwA.O, T_Hu.OHwP.O, T_Hu.OHwC.O)
                
## HvA
# HvA.OxE cols
T_HvA.OEvA.O <- diag(ihr_vA.O*gamma, ncol = age.bins, nrow = age.bins)
T_HvA.OEu.O <- T_HvA.OEvP.O <- 
  T_HvA.OEvC.O <- T_HvA.OEwA.O <- T_HvA.OEwP.O <- T_HvA.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvA.OxH cols
T_HvA.OHvA.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HvA.OHu.O <- T_HvA.OHvP.O <- T_HvA.OHvC.O <- 
  T_HvA.OHwA.O <- T_HvA.OHwP.O <- T_HvA.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvA.OxI cols
T_HvA.OIu.O <- T_HvA.OIvA.O <- T_HvA.OIvP.O <- T_HvA.OIvC.O <- 
  T_HvA.OIwA.O <- T_HvA.OIwP.O <- T_HvA.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvA.OxA cols
T_HvA.OAu.O <- T_HvA.OAvA.O <- T_HvA.OAvP.O <- T_HvA.OAvC.O <- 
  T_HvA.OAwA.O <- T_HvA.OAwP.O <- T_HvA.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HvA.O <- cbind(T_HvA.OEu.O, T_HvA.OEvA.O, T_HvA.OEvP.O, T_HvA.OEvC.O, T_HvA.OEwA.O, T_HvA.OEwP.O, T_HvA.OEwC.O,
                T_HvA.OAu.O, T_HvA.OAvA.O, T_HvA.OAvP.O, T_HvA.OAvC.O, T_HvA.OAwA.O, T_HvA.OAwP.O, T_HvA.OAwC.O, 
                T_HvA.OIu.O, T_HvA.OIvA.O, T_HvA.OIvP.O, T_HvA.OIvC.O, T_HvA.OIwA.O, T_HvA.OIwP.O, T_HvA.OIwC.O, 
                T_HvA.OHu.O, T_HvA.OHvA.O, T_HvA.OHvP.O, T_HvA.OHvC.O, T_HvA.OHwA.O, T_HvA.OHwP.O, T_HvA.OHwC.O)

## HvP
# HvP.OxE cols
T_HvP.OEvP.O <- diag(ihr_vP.O*gamma, ncol = age.bins, nrow = age.bins)
T_HvP.OEu.O <- T_HvP.OEvA.O <- 
  T_HvP.OEvC.O <- T_HvP.OEwA.O <- T_HvP.OEwP.O <- T_HvP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvP.OxH cols
T_HvP.OHvP.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HvP.OHu.O <- T_HvP.OHvA.O <- T_HvP.OHvC.O <- 
  T_HvP.OHwA.O <- T_HvP.OHwP.O <- T_HvP.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvP.OxI cols
T_HvP.OIu.O <- T_HvP.OIvA.O <- T_HvP.OIvP.O <- T_HvP.OIvC.O <- 
  T_HvP.OIwA.O <- T_HvP.OIwP.O <- T_HvP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvP.OxA cols
T_HvP.OAu.O <- T_HvP.OAvA.O <- T_HvP.OAvP.O <- T_HvP.OAvC.O <- 
  T_HvP.OAwA.O <- T_HvP.OAwP.O <- T_HvP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HvP.O <- cbind(T_HvP.OEu.O, T_HvP.OEvA.O, T_HvP.OEvP.O, T_HvP.OEvC.O, T_HvP.OEwA.O, T_HvP.OEwP.O, T_HvP.OEwC.O,
                T_HvP.OAu.O, T_HvP.OAvA.O, T_HvP.OAvP.O, T_HvP.OAvC.O, T_HvP.OAwA.O, T_HvP.OAwP.O, T_HvP.OAwC.O, 
                T_HvP.OIu.O, T_HvP.OIvA.O, T_HvP.OIvP.O, T_HvP.OIvC.O, T_HvP.OIwA.O, T_HvP.OIwP.O, T_HvP.OIwC.O, 
                T_HvP.OHu.O, T_HvP.OHvA.O, T_HvP.OHvP.O, T_HvP.OHvC.O, T_HvP.OHwA.O, T_HvP.OHwP.O, T_HvP.OHwC.O)
                
## HvC
# HvC.OxE cols
T_HvC.OEvC.O <- diag(ihr_vC.O*gamma, ncol = age.bins, nrow = age.bins)
T_HvC.OEu.O <- T_HvC.OEvA.O <- 
  T_HvC.OEvP.O <- T_HvC.OEwA.O <- T_HvC.OEwP.O <- T_HvC.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvC.OxH cols
T_HvC.OHvC.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HvC.OHu.O <- T_HvC.OHvA.O <- T_HvC.OHvP.O <- 
  T_HvC.OHwA.O <- T_HvC.OHwP.O <- T_HvC.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvC.OxI cols
T_HvC.OIu.O <- T_HvC.OIvA.O <- T_HvC.OIvP.O <- T_HvC.OIvC.O <- 
  T_HvC.OIwA.O <- T_HvC.OIwP.O <- T_HvC.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HvC.OxA cols
T_HvC.OAu.O <- T_HvC.OAvA.O <- T_HvC.OAvP.O <- T_HvC.OAvC.O <- 
  T_HvC.OAwA.O <- T_HvC.OAwP.O <- T_HvC.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HvC.O <- cbind(T_HvC.OEu.O, T_HvC.OEvA.O, T_HvC.OEvP.O, T_HvC.OEvC.O, T_HvC.OEwA.O, T_HvC.OEwP.O, T_HvC.OEwC.O,
                T_HvC.OAu.O, T_HvC.OAvA.O, T_HvC.OAvP.O, T_HvC.OAvC.O, T_HvC.OAwA.O, T_HvC.OAwP.O, T_HvC.OAwC.O, 
                T_HvC.OIu.O, T_HvC.OIvA.O, T_HvC.OIvP.O, T_HvC.OIvC.O, T_HvC.OIwA.O, T_HvC.OIwP.O, T_HvC.OIwC.O, 
                T_HvC.OHu.O, T_HvC.OHvA.O, T_HvC.OHvP.O, T_HvC.OHvC.O, T_HvC.OHwA.O, T_HvC.OHwP.O, T_HvC.OHwC.O)

## HwA
# HwA.OxE cols
T_HwA.OEwA.O <- diag(ihr_wA.O*gamma, ncol = age.bins, nrow = age.bins)
T_HwA.OEu.O <- T_HwA.OEvP.O <- 
  T_HwA.OEvC.O <- T_HwA.OEvA.O <- T_HwA.OEwP.O <- T_HwA.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwA.OxH cols
T_HwA.OAwH.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HwA.OHu.O <- T_HwA.OHvP.O <- T_HwA.OHvC.O <- 
  T_HwA.OHvA.O <- T_HwA.OHwP.O <- T_HwA.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwA.OxI cols
T_HwA.OIu.O <- T_HwA.OIvA.O <- T_HwA.OIvP.O <- T_HwA.OIvC.O <- 
  T_HwA.OIwA.O <- T_HwA.OIwP.O <- T_HwA.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwA.OxA cols
T_HwA.OAu.O <- T_HwA.OAvA.O <- T_HwA.OAvP.O <- T_HwA.OAvC.O <- 
  T_HwA.OAwA.O <- T_HwA.OAwP.O <- T_HwA.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HwA.O <- cbind(T_HwA.OEu.O, T_HwA.OEvA.O, T_HwA.OEvP.O, T_HwA.OEvC.O, T_HwA.OEwA.O, T_HwA.OEwP.O, T_HwA.OEwC.O,
                T_HwA.OAu.O, T_HwA.OAvA.O, T_HwA.OAvP.O, T_HwA.OAvC.O, T_HwA.OAwA.O, T_HwA.OAwP.O, T_HwA.OAwC.O, 
                T_HwA.OIu.O, T_HwA.OIvA.O, T_HwA.OIvP.O, T_HwA.OIvC.O, T_HwA.OIwA.O, T_HwA.OIwP.O, T_HwA.OIwC.O, 
                T_HwA.OHu.O, T_HwA.OHvA.O, T_HwA.OHvP.O, T_HwA.OHvC.O, T_HwA.OHwA.O, T_HwA.OHwP.O, T_HwA.OHwC.O)

## HwP
# HwP.OxE cols
T_HwP.OEwP.O <- diag(ihr_wP.O*gamma, ncol = age.bins, nrow = age.bins)
T_HwP.OEu.O <- T_HwP.OEvA.O <- 
  T_HwP.OEvC.O <- T_HwP.OEwA.O <- T_HwP.OEvP.O <- T_HwP.OEwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwP.OxH cols
T_HwP.OHwP.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HwP.OHu.O <- T_HwP.OHvA.O <- T_HwP.OHvC.O <- 
  T_HwP.OHwA.O <- T_HwP.OHvP.O <- T_HwP.OHwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwP.OxI cols
T_HwP.OIu.O <- T_HwP.OIvA.O <- T_HwP.OIvP.O <- T_HwP.OIvC.O <- 
  T_HwP.OIwA.O <- T_HwP.OIwP.O <- T_HwP.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwP.OxA cols
T_HwP.OAu.O <- T_HwP.OAvA.O <- T_HwP.OAvP.O <- T_HwP.OAvC.O <- 
  T_HwP.OAwA.O <- T_HwP.OAwP.O <- T_HwP.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HwP.O <- cbind(T_HwP.OEu.O, T_HwP.OEvA.O, T_HwP.OEvP.O, T_HwP.OEvC.O, T_HwP.OEwA.O, T_HwP.OEwP.O, T_HwP.OEwC.O,
                T_HwP.OAu.O, T_HwP.OAvA.O, T_HwP.OAvP.O, T_HwP.OAvC.O, T_HwP.OAwA.O, T_HwP.OAwP.O, T_HwP.OAwC.O, 
                T_HwP.OIu.O, T_HwP.OIvA.O, T_HwP.OIvP.O, T_HwP.OIvC.O, T_HwP.OIwA.O, T_HwP.OIwP.O, T_HwP.OIwC.O, 
                T_HwP.OHu.O, T_HwP.OHvA.O, T_HwP.OHvP.O, T_HwP.OHvC.O, T_HwP.OHwA.O, T_HwP.OHwP.O, T_HwP.OHwC.O)
                
## HwC
# HwC.OxE cols
T_HwC.OEwC.O <- diag(ihr_wC.O*gamma, ncol = age.bins, nrow = age.bins)
T_HwC.OEu.O <- T_HwC.OEvA.O <- 
  T_HwC.OEvP.O <- T_HwC.OEwA.O <- T_HwC.OEwP.O <- T_HwC.OEvC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwC.OxH cols
T_HwC.OHwC.O <- diag(1-nus, ncol = age.bins, nrow = age.bins)
T_HwC.OHu.O <- T_HwC.OHvA.O <- T_HwC.OHvP.O <- 
  T_HwC.OHwA.O <- T_HwC.OHwP.O <- T_HwC.OHvC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwC.OxI cols
T_HwC.OIu.O <- T_HwC.OIvA.O <- T_HwC.OIvP.O <- T_HwC.OIvC.O <- 
  T_HwC.OIwA.O <- T_HwC.OIwP.O <- T_HwC.OIwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)
# HwC.OxA cols
T_HwC.OAu.O <- T_HwC.OAvA.O <- T_HwC.OAvP.O <- T_HwC.OAvC.O <- 
  T_HwC.OAwA.O <- T_HwC.OAwP.O <- T_HwC.OAwC.O <- matrix(0, ncol = age.bins, nrow = age.bins)

T_HwC.O <- cbind(T_HwC.OEu.O, T_HwC.OEvA.O, T_HwC.OEvP.O, T_HwC.OEvC.O, T_HwC.OEwA.O, T_HwC.OEwP.O, T_HwC.OEwC.O,
                T_HwC.OAu.O, T_HwC.OAvA.O, T_HwC.OAvP.O, T_HwC.OAvC.O, T_HwC.OAwA.O, T_HwC.OAwP.O, T_HwC.OAwC.O, 
                T_HwC.OIu.O, T_HwC.OIvA.O, T_HwC.OIvP.O, T_HwC.OIvC.O, T_HwC.OIwA.O, T_HwC.OIwP.O, T_HwC.OIwC.O, 
                T_HwC.OHu.O, T_HwC.OHvA.O, T_HwC.OHvP.O, T_HwC.OHvC.O, T_HwC.OHwA.O, T_HwC.OHwP.O, T_HwC.OHwC.O)
                
#all H's together
T_H.O <- rbind(T_Hu.O, T_HvA.O, T_HvP.O, T_HvC.O, T_HwA.O, T_HwP.O, T_HwC.O)


# final matrix
T_mat <- rbind(T_E.O, T_A.O, T_I.O, T_H.O)