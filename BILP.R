BILP.MAX.EXACT <- function(c,b,A){
# Exact Maximizing Binary Integer Linear Programming Code
# Uses conventional binary IP solver
# Limited to small problems
# Maximize total value subject to available resources
# Includes 0-1 multidimensional knapsack problems
# Standard form: max c'x, subject to Ax <= b, x binary 0-1 vector
# c positive n-vector, b positive m-vector, A nonnegative m x n matrix
# OUT <- BILP.MAX.EXACT(c,b,A) is a summary list of named results which can be queried
# E.G., The exact BILP binary solution vector is OUT$EXACT_BILP_SOLN
# Bob Agnew, raagnew1@gmail.com, raagnew.com
Time_In <- Sys.time()
library("lpSolve")
m <- length(b)
n <- length(c)
opt <- lp("max",c,A,const.dir=rep("<=",m),b,all.bin=TRUE)
df <- data.frame(as.vector(A%*%opt$solution),rep(" <= ",m),b)
names(df) <- c("Ax"," <= ","b")
Time_Out <- Sys.time()
del_time <- as.numeric(difftime(Time_Out,Time_In,units="secs")) # Imprecise for small deltas
result <- list(n,m,del_time,opt$solution,opt$objval,sum(as.vector(A%*%opt$solution) <= b),df)
names(result) <- c("NUMBER_VARIABLES","NUMBER_CONSTRAINTS","EXECUTION_SECONDS","EXACT_BILP_SOLN",
"EXACT_BILP_MAX","CONSTRAINTS_SATISFIED","CONSTRAINT_DETAILS")
return(result)}

BILP.MAX.LP <- function(c,b,A){
# Approximate Maximizing Binary Integer Linear Programming Code
# Relevant when exact solution is unavailable due to size limitation
# Employs Williamson primal-dual approximation
# BILP is relaxed to an LP
# Relaxed LP is solved explicitly, including dual variables
# BILP assignments are ordered primarily by shadow prices,
# secondarily by relaxed LP values for zero shadow prices
# Performs well on moderately large problems
# Maximize total value subject to available resources
# Includes 0-1 multidimensional knapsack problems
# Standard form: max c'x, subject to Ax <= b, x binary 0-1 vector
# c positive n-vector, b positive m-vector, A nonnegative m x n matrix
# Vector of zeros must be feasible
# OUT <- BILP.MAX.LP(c,b,A) is a summary list of named results which can be queried
# E.G., The approximate BILP binary solution vector is OUT$APPROX_BILP_SOLN
# Bob Agnew, raagnew1@gmail.com, raagnew.com
Time_In <- Sys.time()
library("lpSolve")
m <- length(b)
n <- length(c)
# Explicit solution of relaxed LP
b1 <- c(b,rep(1,n))
A1 <- rbind(A,diag(n))
opt <- lp("max",c,A1,const.dir=rep("<=",m+n),b1,compute.sens=1)
# Heuristic sort and binary assignment based on relaxed LP dual variables
# and relaxed LP solution for zero duals
z <- opt$duals[(m+1):(m+n)] + opt$duals[(m+n+1):(m+2*n)]
z1 <- opt$solution # Explicit relaxed LP primal solution
I <- (1:n)[z==0] # Indices of zero duals
J <- order(z) # Ordering of dual variables
K <- match(I,J) # Indices of I within J
J[K] <- I[order(z1[I])] # Reordering of zero duals based on primal LP solution
ord <- rev(J) # Reverse ordering of duals to descending
x <- rep(0,n) # Initial feasible solution of zeros
# Ordered sequential binary assignments to ensure feasibility
for (j in 1:n){
x1 <- x
k <- ord[j]
x1[k] <- 1
if (min(b - as.vector(A%*%x1)) >= 0){x[k] <- 1}}
obj <- as.numeric(c%*%x)
df <- data.frame(as.vector(A%*%x),rep(" <= ",m),b)
names(df) <- c("Ax"," <= ","b")
Time_Out <- Sys.time()
del_time <- as.numeric(difftime(Time_Out,Time_In,units="secs")) # Imprecise for small deltas
result <- list(n,m,del_time,x,obj,opt$objval,obj/opt$objval,sum(as.vector(A%*%x) <= b),df)
names(result) <- c("NUMBER_VARIABLES","NUMBER_CONSTRAINTS","EXECUTION_SECONDS","APPROX_BILP_SOLN",
"APPROX_BILP_MAX","RELAXED_LP_MAX","APPROXIMATION_RATIO","CONSTRAINTS_SATISFIED","CONSTRAINT_DETAILS")
return(result)}

BILP.MAX.APPROX <- function(c,b,A){
# Approximate Maximizing Binary Integer Linear Programming Code
# Relevant when even explicit relaxed LP solution is unavailable due to size limitation
# Employs Williamson primal-dual approximation
# BILP is relaxed to an LP
# Dual LP is solved in collapsed form and provides bound - avoids simplex
# BILP assignments are by ordered approximate shadow prices
# Performs well on many large problems
# Maximize total value subject to available resources
# Includes 0-1 multidimensional knapsack problems
# Standard form: max c'x, subject to Ax <= b, x binary 0-1 vector
# c n-vector, b positive m-vector, A nonnegative m x n matrix
# Vector of zeros must be feasible
# OUT <- BILP.MAX.APPROX(c,b,A) is a summary list of named results which can be queried
# E.G., The approximate BILP binary solution vector is OUT$APPROX_BILP_SOLN
# Bob Agnew, raagnew1@gmail.com, raagnew.com
Time_In <- Sys.time()
library("cmna")
m <- length(b)
n <- length(c)
# Collapsed relaxed LP dual minimization function
dual <- function(y){as.numeric(b%*%y) + sum(pmax(0,c - as.vector(y%*%A)))}
gr <- function(y){b - as.vector(A%*%as.vector(1*(c > as.vector(y%*%A))))}
# Golden section search for a good starting vector
f <- function(u){dual(rep(u,m))}
u <- goldsectmin(f,0,f(0)/sum(b),tol=1e-5,m=10000)
start <- rep(u,m) # Optimized start vector - generally better than zeros
# Optimization of collapsed dual with plain vanilla nonlinear solver
# nlminb nonlinear solver incorporates lower bounds
opt <- nlminb(start,dual,gradient=gr,control=list(iter.max=100000,eval.max=100000,
step.min=.01,step.max=.01,abs.tol=1e-20,x.tol=1e-20),lower=rep(0,m))
# Optimal dual solution is approximate, close but not exact
y <- opt$par # Formal constraint shadow prices
# Variable upper bound shadow prices, also called Reduced Costs
z <- c - as.vector(y%*%A) # Unrounded
# Heuristic sort and binary assignment based on approximate relaxed LP dual variables
ord <- order(z,decreasing=TRUE) # Descending order by dual-informed shadow prices
x <- rep(0,n) # Initial feasible solution of zeros
# Ordered sequential binary assignments to ensure feasibility
for (j in 1:n){
x1 <- x
k <- ord[j]
x1[k] <- 1
if (min(b - as.vector(A%*%x1)) >= 0){x[k] <- 1}}
obj <- as.numeric(c%*%x)
df <- data.frame(as.vector(A%*%x),rep(" <= ",m),b)
names(df) <- c("Ax"," <= ","b")
Time_Out <- Sys.time()
del_time <- as.numeric(difftime(Time_Out,Time_In,units="secs")) # Imprecise for small deltas
result <- list(n,m,del_time,opt$iterations,opt$evaluations,x,obj,opt$objective,
obj/opt$objective,sum(as.vector(A%*%x) <= b),df)
names(result) <- c("NUMBER_VARIABLES","NUMBER_CONSTRAINTS","EXECUTION_SECONDS","SOLVER_ITERATIONS",
"SOLVER_EVALUATIONS","APPROX_BILP_SOLN","APPROX_BILP_MAX","APPROX_DUAL_LP_MIN",
"APPROXIMATION_RATIO","CONSTRAINTS_SATISFIED","CONSTRAINT_DETAILS")
return(result)}

BILP.MIN.EXACT <- function(c,b,A){
# Exact Minimizing Binary Integer Linear Programming Code
# Uses conventional binary IP solver
# Oriented to small problems
# Standard form: min c'x, subject to Ax >= b, x binary 0-1 vector
# c positive n-vector, b positive m-vector, A m x n matrix
# Minimize cost while covering constraints
# Includes set covering problems
# OUT <- BILP.MIN.EXACT(c,b,A) is a summary list of named results which can be queried
# E.G., The exact BILP binary solution vector is OUT$EXACT_BILP_SOLN
# Bob Agnew, raagnew1@gmail.com, raagnew.com
Time_In <- Sys.time()
library("lpSolve")
m <- length(b)
n <- length(c)
opt <- lp("min",c,A,const.dir=rep(">=",m),b,all.bin=TRUE)
df <- data.frame(as.vector(A%*%opt$solution),rep(" >= ",m),b)
names(df) <- c("Ax"," >= ","b")
Time_Out <- Sys.time()
del_time <- as.numeric(difftime(Time_Out,Time_In,units="secs")) # Imprecise for small deltas
result <- list(n,m,del_time,opt$solution,opt$objval,sum(A%*%opt$solution >= b),df)
names(result) <- c("NUMBER_VARIABLES","NUMBER_CONSTRAINTS","EXECUTION_SECONDS","EXACT_BILP_SOLN",
"EXACT_BILP_MIN","CONSTRAINTS_SATISFIED","CONSTRAINT_DETAILS")
return(result)}

BILP.MIN.LP <- function(c,b,A){
# Simple Minimizing Binary Integer Linear Programming Code
# Standard form: min c'x, subject to Ax >= b, x binary 0 or 1
# c n-vector, b positive m-vector, A nonnegative m x n matrix
# Transformed to equivalent maximizing problem for solution
# Minimize cost while covering constraints
# Includes set covering problems, although roughly
# Vector of ones must be feasible
# OUT <- BILP.MIN.LP(c,b,A) is a summary list of named results which can be queried
# E.G., The approximate BILP binary solution vector is OUT$APPROX_BILP_SOLN
# Bob Agnew, raagnew1@gmail.com, raagnew.com
Time_In <- Sys.time()
m <- length(b)
n <- length(c)
ones <- rep(1,n) # Vector of ones, assumed feasible
obj_ones <- as.numeric(c%*%ones) # Objective function with all ones
out <- BILP.MAX.LP(c,as.vector(A%*%ones) - b,A)
x <- ones - out$APPROX_BILP_SOLN
obj <- obj_ones - out$APPROX_BILP_MAX
bnd <- obj_ones - out$RELAXED_LP_MAX
df <- data.frame(as.vector(A%*%x),rep(" >= ",m),b)
names(df) <- c("Ax"," >= ","b")
Time_Out <- Sys.time()
del_time <- as.numeric(difftime(Time_Out,Time_In,units="secs")) # Imprecise for small deltas
result <- list(n,m,del_time,x,obj,bnd,obj/bnd,sum(as.vector(A%*%x) >= b),df)
names(result) <- c("NUMBER_VARIABLES","NUMBER_CONSTRAINTS","EXECUTION_SECONDS","APPROX_BILP_SOLN",
"APPROX_BILP_MIN","RELAXED_LP_MIN","APPROXIMATION_RATIO","CONSTRAINTS_SATISFIED","CONSTRAINT_DETAILS")
return(result)}

BILP.MIN.APPROX <- function(c,b,A){
# Simple Minimizing Binary Integer Linear Programming Code
# Standard form: min c'x, subject to Ax >= b, x binary 0 or 1
# c n-vector, b postive m-vector, A nonnegative m x n matrix
# Transformed to equivalent maximizing problem for solution
# Minimize cost while covering constraints
# Includes set covering problems, although roughly
# Vector of ones must be feasible
# OUT <- BILP.MIN.APPROX(c,b,A) is a summary list of named results which can be queried
# E.G., The approximate BILP binary solution vector is OUT$APPROX_BILP_SOLN
# Bob Agnew, raagnew1@gmail.com, raagnew.com
Time_In <- Sys.time()
m <- length(b)
n <- length(c)
ones <- rep(1,n) # Vector of ones, assumed feasible
obj_ones <- as.numeric(c%*%ones) # Objective function with all ones
out <- BILP.MAX.APPROX(c,as.vector(A%*%ones) - b,A)
x <- ones - out$APPROX_BILP_SOLN
obj <- obj_ones - out$APPROX_BILP_MAX
bnd <- obj_ones - out$APPROX_DUAL_LP_MIN
df <- data.frame(as.vector(A%*%x),rep(" >= ",m),b)
names(df) <- c("Ax"," >= ","b")
Time_Out <- Sys.time()
del_time <- as.numeric(difftime(Time_Out,Time_In,units="secs")) # Imprecise for small deltas
result <- list(n,m,del_time,out$SOLVER_ITERATIONS,out$SOLVER_EVALUATIONS,x,obj,bnd,obj/bnd,sum(as.vector(A%*%x) >= b),df)
names(result) <- c("NUMBER_VARIABLES","NUMBER_CONSTRAINTS","EXECUTION_SECONDS","SOLVER_ITERATIONS",
"SOLVER_EVALUATIONS","APPROX_BILP_SOLN","APPROX_BILP_MIN","APPROX_DUAL_LP_MAX","APPROXIMATION_RATIO",
"CONSTRAINTS_SATISFIED","CONSTRAINT_DETAILS")
return(result)}

# Selected output indices
OUT1 <- c(1:3,5:8)
OUT2 <- c(1:5,7:10)

# Small Textbook Maximization Problem
# max 8*x1 + 11*x2 + 6*x3 + 4*x4
# subject to:
# 5*x1 + 7*x2 + 3*x4 <= 14
# 8*x1 + 4*x3 + 4*x4 <= 12
# 2*x1 + 10*x2 + 6*x3 + 4*x4 <= 15
# x1, x2, x3, x4 in {0,1}
c <- c(8,11,6,4)
b <- c(14,12,15)
A <- c(5,7,0,3)
A <- c(A,8,0,4,4)
A <- c(A,2,10,6,4)
A <- matrix(A,nrow=3,ncol=4,byrow=TRUE)

# Maximization summaries

BILP.MAX.EXACT(c,b,A) # Optimal

BILP.MAX.LP(c,b,A) # Optimal

BILP.MAX.APPROX(c,b,A) # Optimal

# Mirror-image minimization summaries

BILP.MIN.EXACT(c,b,A) # Optimal

BILP.MIN.LP(c,b,A) # Optimal

BILP.MIN.APPROX(c,b,A) # Optimal

# Multidimensional Knapsack Problems from J. E. Beasley OR-Library:
# people.brunel.ac.uk/~mastjjb/jeb/info.html

# Problem 7 in file mknap1.txt
v <- scan(file="c:/BILP/mknap1_7.txt",sep="")
c <- v[4:53]
b <- v[304:308]
A <- matrix(v[54:303],nrow=length(b),ncol=length(c),byrow=TRUE)

#Maximization summaries

BILP.MAX.EXACT(c,b,A) # Optimal

BILP.MAX.LP(c,b,A) # Close

BILP.MAX.APPROX(c,b,A) # Close

# Mirror-image minimization summaries

BILP.MIN.EXACT(c,b,A) # Optimal

BILP.MIN.LP(c,b,A) # Close

BILP.MIN.APPROX(c,b,A) # Close

# Problem 30 in file mknapcb9.txt
v <- scan(file="c:/BILP/mknapcb9_30.txt",sep="")
c <- v[4:503]
b <- v[15504:15533]
A <- matrix(v[504:15503],nrow=length(b),ncol=length(c),byrow=TRUE)

# Maximization summaries

# EXACT version didn't solve

BILP.MAX.LP(c,b,A)[OUT1] # Very close

BILP.MAX.APPROX(c,b,A)[OUT2] # Very close

# Chu & Beasley (1998) got an even better maximum of 300460 via a
# genetic algorithm but it was also much slower.

# Mirror-image minimization summaries

# EXACT version didn't solve

BILP.MIN.LP(c,b,A)[OUT1] # Very close

BILP.MIN.APPROX(c,b,A)[OUT2] # Very close

# Randomized Inputs for Large Randomized Problem (50000 variables, 300 constraints)
# Explicit LP versions wouldn't run.  Only approximations ran.
set.seed(17)
c <- sample.int(100,50000,replace=TRUE,prob=NULL)
b <- sample.int(1250000,300,replace=TRUE,prob=NULL)
A <- sample.int(50,15000000,replace=TRUE,prob=NULL)
# write(c(c,b,A),"c:/BILP/Large_Problem.txt") for Python version
A <- matrix(A,nrow=length(b),ncol=length(c),byrow=TRUE)

# Maximization summaries

# Exact version didn't run

# LP version didn't run

BILP.MAX.APPROX(c,b,A)[OUT2] # Extremely close

# Mirror-image minimization summaries

# Exact version didn't run

# LP version didn't run

BILP.MIN.APPROX(c,b,A)[OUT2] # Extremely close

# The following set covering problems show pitfalls of simple shadow price sorting.
# Sorting should work better with more differentiated values in c-vector.
# Bottom line, there are better approximate algorithms for these problems.
# Approximation ratios aren't particularly useful for these problems.

# Small Textbook Minimization Problem
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
b <- rep(1,5)
A <- c(1,1,1,0,0,0,1,0,0)
A <- c(A,0,0,1,1,1,0,0,1,0)
A <- c(A,1,0,0,1,0,1,0,1,0)
A <- c(A,1,1,0,0,1,1,0,0,1)
A <- c(A,0,1,1,1,1,1,0,0,1)
A <- matrix(A,nrow=5,ncol=9,byrow=TRUE)

# Minimization summaries

BILP.MIN.EXACT(c,b,A) # Optimal

BILP.MIN.LP(c,b,A) # Suboptimal

BILP.MIN.APPROX(c,b,A) # Suboptimal

# Mirror-image maximization summaries

BILP.MAX.EXACT(c,b,A) # Optimal

BILP.MAX.LP(c,b,A) # Optimal

BILP.MAX.APPROX(c,b,A) # Optimal

# Set covering problems from J. E. Beasley OR-Library:
# people.brunel.ac.uk/~mastjjb/jeb/info.html
# These problems are particularly challenging and
# approximation ratios are not useful.

# Problem in file scpe5.txt
v <- scan(file="c:/BILP/scpe5.txt",sep="")
m <- v[1]
n <- v[2]
c <- v[3:(n+2)]
b <- rep(1,m)
s <- n + 3
A <- NULL
for (i in 1:m){
x <- v[(s+1):(s+v[s])]
z <- y <- match(1:n,x)
z[is.na(y)] <- 0
z[!is.na(y)] <- 1
s <- s + v[s] + 1
A <- rbind(A,z)}

# Minimization summaries

BILP.MIN.EXACT(c,b,A) # Optimal

BILP.MIN.LP(c,b,A) # Suboptimal

BILP.MIN.APPROX(c,b,A) # Suboptimal

# Beasley (1987) achieved the exact minimum of
# 5 using some of the same ingredients, but again
# it took more time.

# Mirror-image maximization summaries

BILP.MAX.EXACT(c,b,A) # Optimal

BILP.MAX.LP(c,b,A) # Optimal

BILP.MAX.APPROX(c,b,A) # Optimal 

# Problem in file scpa5.txt
v <- scan(file="c:/BILP/scpa5.txt",sep="")
m <- v[1]
n <- v[2]
c <- v[3:(n+2)]
b <- rep(1,m)
s <- n + 3
A <- NULL
for (i in 1:m){
x <- v[(s+1):(s+v[s])]
z <- y <- match(1:n,x)
z[is.na(y)] <- 0
z[!is.na(y)] <- 1
s <- s + v[s] + 1
A <- rbind(A,z)}

# Minimization summaries

# Exact version didn't run

BILP.MIN.LP(c,b,A)[OUT1] # Suboptimal

BILP.MIN.APPROX(c,b,A)[OUT2] # Suboptimal

# Beasely (1987) indicates an exact minimum of 236 for this
# problem.  His algorithm achieved an upper bound of
# 240, better than our approximations.

# Mirror-image maximization summaries

# Exact version didn't run

BILP.MAX.LP(c,b,A)[OUT1] # Presumably suboptimal

BILP.MAX.APPROX(c,b,A)[OUT2] # Presumably suboptimal

# Problem in file scpnre5.txt
v <- scan(file="c:/BILP/scpnre5.txt",sep="")
m <- v[1]
n <- v[2]
c <- v[3:(n+2)]
b <- rep(1,m)
s <- n + 3
A <- NULL
for (i in 1:m){
x <- v[(s+1):(s+v[s])]
z <- y <- match(1:n,x)
z[is.na(y)] <- 0
z[!is.na(y)] <- 1
s <- s + v[s] + 1
A <- rbind(A,z)}

# Minimization summaries

# Exact version didn't run

BILP.MIN.LP(c,b,A)[OUT1] # Suboptimal

BILP.MIN.APPROX(c,b,A)[OUT2] # Suboptimal

# Beasley (1990) achieved an approximate minimum of 28,
# using a more complex algorithm, better than ours.

# Mirror-image maximization summaries

# Exact version didn't run 

BILP.MAX.LP(c,b,A)[OUT1] # Suboptimal

BILP.MAX.APPROX(c,b,A)[OUT2] # Presumably suboptimal

