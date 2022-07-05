update_R0 <- function(parameters,R0.D = NULL, R0.O = NULL){
  if(!is.null(R0.D)){
    rel.value <- R0.D/parameters$R0.D
    parameters$beta.D <- rel.value*parameters$beta.D
    parameters$beta_vA.D <- rel.value*parameters$beta_vA.D
    parameters$beta_wA.D <- rel.value*parameters$beta_wA.D
    parameters$beta_vP.D <- rel.value*parameters$beta_vP.D
    parameters$beta_wP.D <- rel.value*parameters$beta_wP.D
    parameters$beta_bP.D <- rel.value*parameters$beta_bP.D
    parameters$beta_vC.D <- rel.value*parameters$beta_vC.D
    parameters$beta_wC.D <- rel.value*parameters$beta_wC.D
    parameters$R0.D <- R0.D
  }
  if(!is.null(R0.O)){
    rel.value <- R0.O/parameters$R0.O
    parameters$beta.O <- rel.value*parameters$beta.O
    parameters$beta_vA.O <- rel.value*parameters$beta_vA.O
    parameters$beta_wA.O <- rel.value*parameters$beta_wA.O
    parameters$beta_vP.O <- rel.value*parameters$beta_vP.O
    parameters$beta_wP.O <- rel.value*parameters$beta_wP.O
    parameters$beta_bP.O <- rel.value*parameters$beta_bP.O
    parameters$beta_vC.O <- rel.value*parameters$beta_vC.O
    parameters$beta_wC.O <- rel.value*parameters$beta_wC.O
    parameters$R0.O <- R0.O
  }
  return(parameters)
}