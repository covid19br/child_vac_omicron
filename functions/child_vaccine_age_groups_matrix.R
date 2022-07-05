child_vaccine_age_groups_matrix <- function(mat, dist_et){
  mat2 <- mat
  ##these matrices describe how many daily contacts a person of line i has with a person of column j,
  ##thus, to change age proportions, if age bin X is reduced by 20%, the contacts with this age bin
  ##decreases 20%. Therefore, it reduces the contacts in column X by 20%. In the other hand, age bin Y
  ##receives 20% of the contacts with age bin X, thus increasing column Y by 20% of X.
  ##maintains 0-4 and 20-29, ....
  ##from 5-9 to 5-11, second agebin, 2/5 of third agebin goes to second
  mat2[,2] <- mat2[,2] + 0.40*mat2[,3]
  ##from 10-14 to 12-14, third agebin
  mat2[,3] <- mat2[,3] - 0.40*mat2[,3]
  ##from 12-14 to 12-17 third agebin
  mat2[,3] <- mat2[,3] + 0.60*mat2[,4]
  ##from 15-19 to 18-19 fourth agebin
  mat2[,4] <- mat2[,4] - 0.60*mat2[,4] 
  return(mat2)
}
