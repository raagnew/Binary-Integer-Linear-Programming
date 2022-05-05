bilp.max <- function(c,A,b,itermax){
# Simple Maximizing Binary Integer Linear Programming Code
# Uses collapsed dual approach
# Standard form: max c'x, subject to Ax <= b, x binary 0 or 1
# c n-vector, b m-vector, A m x n matrix
# Bob Agnew, raagnew1@gmail.com, raagnew.com
m <- length(b)
n <- length(c)
if (min(b) < 0){
stop("Initial Solution of Zeros Infeasible")
return}
# Collapsed relaxed LP dual minimization function
dual <- function(y){sum(b%*%y) + sum(pmax(0,c - y%*%A))}
# Optimization of collapsed dual with plain vanilla nonlinear solver
opt <- nlminb(rep(0,m),dual,control=list(iter.max=itermax),lower=rep(0,m),upper=Inf)
# Optimal dual solution is approximate, close but not exact
y <- opt$par
v <- c - y%*%A
ord <- order(v,decreasing=TRUE)
x <- rep(0,n) # Initial solution of zeros
# Sequential binary assignments to ensure constraints satisfied
for (j in 1:n){
x1 <- x
k <- ord[j]
x1[k] <- 1
ifelse(min(b - A%*%x1) >= 0,x[k] <- 1,break)}
result <- list(x,cbind(A%*%x,b),as.numeric(c%*%x),opt$objective)
names(result) <- c("solution","constraints","maximum","dual_bound")
return(result)}

bilp.min <- function(c,A,b,itermax){
# Simple Minimizing Binary Integer Linear Programming Code
# Uses collapsed dual approach
# Standard form: min c'x, subject to Ax >= b, x binary 0 or 1
# c n-vector, b m-vector, A m x n matrix
# Bob Agnew, raagnew1@gmail.com, raagnew.com
m <- length(b)
n <- length(c)
if (min(A%*%rep(1,n) - b) < 0){
stop("Initial Solution of Ones Infeasible")
return}
# Negative collapsed relaxed LP dual minimization function
dual <- function(y){-sum(b%*%y) + sum(pmax(0,y%*%A - c))}
# Optimization of negative collapsed dual with plain vanilla nonlinear solver
opt <- nlminb(rep(0,m),dual,control=list(iter.max=itermax),lower=rep(0,m),upper=Inf)
# Optimal dual solution is approximate, close but not exact
y <- opt$par
v <- c - y%*%A
ord <- order(v,decreasing=TRUE)
x <- rep(1,n) # Initial solution of ones
# Sequential binary assignments to ensure constraints satisfied
for (j in 1:n){
x1 <- x
k <- ord[j]
x1[k] <- 0
ifelse(min(A%*%x1 - b) >= 0,x[k] <- 0,break)}
result <- list(x,cbind(A%*%x,b),as.numeric(c%*%x),-opt$objective)
names(result) <- c("solution","constraints","minimum","dual_bound")
return(result)}

# Textbook problem
# max 8*x1 + 11*x2 + 6*x3 + 4*x4
# subject to:
# 5*x1 + 7*x2 + 3*x4 <= 14
# 8*x1 + 4*x3 + 4*x4 <= 12
# 2*x1 + 10*x2 + 6*x3 + 4*x4 <= 15
# x1, x2, x3, x4 in {0,1}
c <- c(8,11,6,4)
A <- c(5,7,0,3)
A <- c(A,8,0,4,4)
A <- c(A,2,10,6,4)
A <- matrix(A,nrow=3,ncol=4,byrow=TRUE)
b <- c(14,12,15)

# bilp.max solution
result1 <- bilp.max(c,A,b,1000)
# Binary solution - optimal
result1$solution
#constraints satisfied
result1$constraints
# Maximal value
result1$maximum
# Dual upper bound
result1$dual_bound

# lpSolve solution for comparison
library("lpSolve")
result2 <- lp("max",c,A,const.dir=rep("<=",3),b,all.bin=TRUE)
# Binary solution - optimal
result2$solution
# Constraints satisfied
cbind(A%*%result2$solution,b)
# Maximal value
result2$objval

# Another Textbook Problem
# min 10*x1 + 12*x2 + 12*x3 + 13*x4 + 11*x5
# + 9*x6 + 7*x7 + 8*x8 + 8*x9
# subject to:
# x1 + x2 + x3 + x7 >= 1
# x3 + x4 + x5 + x8 >= 1
# x1 + x4 + x6 + x8 >= 1
# x1 + x2 + x5 + x6 + x9 >= 1
# x2 + x3 + x4 + x5 + x6 + x9 >= 1
# x1, x2, x3, x4, x5, x6, x7, x8, x9 in {0,1}
c <- c(10,12,12,13,11,9,7,8,8)
A <- c(1,1,1,0,0,0,1,0,0)
A <- c(A,0,0,1,1,1,0,0,1,0)
A <- c(A,1,0,0,1,0,1,0,1,0)
A <- c(A,1,1,0,0,1,1,0,0,1)
A <- c(A,0,1,1,1,1,1,0,0,1)
A <- matrix(A,nrow=5,ncol=9,byrow=TRUE)
b <- rep(1,5)

# bilp.min solution
result3 <- bilp.min(c,A,b,1000)
# Binary solution - suboptimal
result3$solution
# Constraints satisfied
result3$constraints
# Minimal value - suboptimal
result3$minimum
# Dual lower bound 
result3$dual_bound

# lpSolve solution for comparison
library("lpSolve")
result4 <- lp("min",c,A,const.dir=rep(">=",5),b,all.bin=TRUE)
# Binary solution - optimal
result4$solution
# Constraints satisfied exactly
cbind(A%*%result4$solution,b)
# Minimal value - optimal
result4$objval

# Randomized inputs for large problems
set.seed(17)
c <- sample.int(100,50000,replace=TRUE,prob=NULL)
A <- sample.int(50,5000000,replace=TRUE,prob=NULL)
A <- matrix(A,nrow=100,ncol=50000,byrow=TRUE)
b <- sample.int(1250000,100,replace=TRUE,prob=NULL)

date() #Time

# bilp.max solution
result5 <- bilp.max(c,A,b,10000)
# First 100 elements of binary solution
result5$solution[1:100]
# First 100 constraint elements
cbind((A[1:100,]%*%result5$solution),b[1:100])
# Constraints satisfied by maximal solution
sum(A%*%result5$solution <= b)
# Maximal value - close to dual upper bound
result5$maximum
# Dual upper bound 
result5$dual_bound

date() #Time - efficient

# lpSolve solution for comparison - gagged on large problem, just hung for hours
# library("lpSolve")
# result6 <- lp("max",c,A,const.dir=rep("<=",100),b,all.bin=TRUE)
# result6$solution[1:100]
# result6$objval

# bilp.min solution
result7 <- bilp.min(c,A,b,10000)
# First 100 elements of binary solution
result7$solution[1:100]
# First 100 constraint elements
cbind((A[1:100,]%*%result7$solution),b[1:100])
# Constraints satisfied by minimal solution
sum(A%*%result7$solution >= b)
# Minimal value - close to dual lower bound
result7$minimum
# Dual lower bound 
result7$dual_bound

date() #Time - efficient


 














