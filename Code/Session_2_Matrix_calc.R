#################################################################################
# Script for the course: Applied multivariate Statistics with R					#
# 					by Ralf B. Schäfer											#
# 					WS 2017/18													#
# The script provides matrix calculations for slides "mathematical basis"		#
#################################################################################

### first we create two matrices
mat1 <- matrix(c(4, 3, 4, 5, 3, 5, 3, 22, 2, 1), ncol = 2)
mat1
mat2 <- matrix(c(5, 3, 3, 5, 3, 4, 3, 45, 77, 3), ncol = 2)
mat2

# check number of columns and rows
ncol(mat2)
nrow(mat2)

# Addition of matrices
mat3 <- mat1 + mat2
mat3

# normal multiplication
mat4 <- mat1 * mat2
mat4

# matrix multiplication
mat5 <- t(mat1) %*% mat2
mat5
mat6 <- mat1 %*% t(mat2)
mat6
#  t(A) %*% B is not the same as A %*% t(B)

# crossprod gives same results
mat7 <- crossprod(mat1, mat2)
mat7
mat8 <- tcrossprod(mat1, mat2)
mat8

# computation of the inverse matrix
nrow(mat7) - qr(mat7)$rank
# inverse does exist because n = rank of matrix
# computation of inverse
mat7_inv <- solve(mat7)

mat7_inv %*% mat7
mat7 %*% mat7_inv
# both give the identity matrix

# identity matrix, matrix multiplication returns initial matrix
mat7
diag(x = 1, 2, 2)
mat7 %*% diag(x = 1, 2, 2)


# inverse does not exist because n > rank of matrix
nrow(mat8) - qr(mat8)$rank
# check:
solve(mat8)
# inverse does not exist!!!

# however, we can compute the generalized inverse
library(MASS)
ginv(mat8)

####################################################################
# use matrix analysis to derive betas for multiple regression      #
####################################################################
# create matrix with variables
reg_mat <- matrix(sample(seq(0, 9), 60, replace = TRUE), ncol = 4)
reg_mat

### add intercept term
inter <- rep(1, nrow(reg_mat))
full <- cbind(inter, reg_mat)

# create response variable
resp_var <- sample(seq(0, 9), 15, replace = TRUE)

# the equation for the betas for multiple regression is b = (t(X) X)^-1 t(X) Y
# we have learned how we can compute this:
betas <- solve(t(full) %*% full) %*% (t(full) %*% resp_var)
betas

# compare with multiple regression
coefficients(lm(resp_var ~ ., data = data.frame(full)))
# gives same results
