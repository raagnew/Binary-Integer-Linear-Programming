import numpy as np
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from scipy.optimize import milp
from scipy.optimize import LinearConstraint
from scipy.optimize import Bounds
import time
import pandas as pd

def BILP_MAX_EXACT(c, b, A):
    # Exact Maximizing Binary Integer Linear Programming Code
    # Uses conventional binary IP solver
    # Limited to small problems
    # Maximize total value subject to available resources
    # Includes 0-1 multidimensional knapsack problems
    # Standard form: max c'x, subject to Ax <= b, x binary 0-1  vector
    # c positive n-vector, b positive m-vector, A nonnegative m x n matrix
    # OUT = BILP_MAX_EXACT(c,b,A) is a summary dictionary of named results which can be queried
    # E.G., The exact BILP binary solution vector is OUT["EXACT_BILP_SOLN"]
    # Bob Agnew, raagnew1@gmail.com, raagnew.com

    Time_In = time.time()

    m = len(b)
    n = len(c)
    
    c1 = np.negative(c)
    b_u = b
    
    integrality = np.ones_like(c)    
    bounds = Bounds(0,1)
    constraints = LinearConstraint(A,ub=b_u)    
    
    opt = milp(c=c1,integrality=integrality,bounds=bounds,
    constraints=constraints)

    solution = np.round(opt.x).astype(int)
    obj = -opt.fun
    Ax = np.dot(A, solution)
    constraints_satisfied = np.sum(Ax <= b)

    df = pd.DataFrame({
        "Ax": Ax,
        " <= ": [" <= "] * m,
        "b": b
    },index=list(range(1,(m+1))))
    
    Time_Out = time.time()
    del_time = Time_Out - Time_In

    result = {
        "NUMBER_VARIABLES": n,
        "NUMBER_CONSTRAINTS": m,
        "EXECUTION_SECONDS": del_time,
        "EXACT_BILP_SOLN": solution,
        "EXACT_BILP_MAX": obj,
        "CONSTRAINTS_SATISFIED": constraints_satisfied,
        "CONSTRAINT_DETAILS": df
    }

    return result

def BILP_MAX_APPROX(c, b, A):
    # Approximate Maximizing Binary Integer Linear Programming Code
    # Relevant for large problems where exact solution is unavailable
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
    # OUT = BILP_MAX_APPROX(c,b,A) returns a summary dictionary which can be queried
    # E.G., The approximate BILP binary solution vector is OUT['APPROX_BILP_SOLN'] 
    # Bob Agnew, raagnew1@gmail.com, raagnew.com

    Time_In = time.time()
    m = len(b)
    n = len(c)

    # Collapsed relaxed LP dual minimization function
    def dual(y):
        y = np.array(y)
        return np.dot(b, y) + np.sum(np.maximum(0, c - y @ A))

    def gr(y):
        y = np.array(y)
        return b - A @ (c > y @ A).astype(float)

    # Golden section search for a good starting vector
    def f(u):
        return dual([u] * m)

    res = minimize_scalar(f, bracket=(0, f(0)/np.sum(b)), method='golden',tol=1e-5,options={'maxiter':10000})
    u = res.x
    start = np.full(m, u)  # Optimized start vector
 
    bounds = [(0, None)] * m

    opt = minimize(dual, start, jac=gr, bounds=bounds,
    method='L-BFGS-B',options={'maxiter': 100000, 'ftol': 1e-20,   'gtol': 1e-20})

    y = opt.x  # Formal constraint shadow prices

    z = c - y @ A  # Unrounded variable upper bound shadow prices (Reduced Costs)

    ord_indices = np.argsort(-z)  # Descending order by dual-informed shadow prices

    x = np.zeros(n, dtype=int)  # Initial feasible solution of zeros

    for k in ord_indices:
        x1 = x.copy()
        x1[k] = 1
        if np.all(b - A @ x1 >= 0):
            x[k] = 1

    obj = np.dot(c, x)
    Ax = A @ x
    constraints_satisfied = np.sum(Ax <= b)

    df = pd.DataFrame({
        "Ax": Ax,
        " <= ": [" <= "] * m,
        "b": b
    },index=list(range(1,(m+1))))

    Time_Out = time.time()
    del_time = Time_Out - Time_In

    result = {
        "NUMBER_VARIABLES": n,
        "NUMBER_CONSTRAINTS": m,
        "EXECUTION_SECONDS": del_time,
        "SOLVER_ITERATIONS": opt.nit,
        "SOLVER_EVALUATIONS": opt.nfev,
        "APPROX_BILP_SOLN": x,
        "APPROX_BILP_MAX": obj,
        "APPROX_DUAL_LP_MIN": opt.fun,
        "APPROXIMATION_RATIO": obj / opt.fun,
        "CONSTRAINTS_SATISFIED": constraints_satisfied,
        "CONSTRAINT_DETAILS": df
    }

    return result

def BILP_MIN_EXACT(c, b, A):
    # Exact Minimizing Binary Integer Linear Programming Code
    # Uses conventional binary IP solver
    # Limited to small problems
    # Minimize total cost while covering stipulated constraints
    # Includes set covering problems
    # Standard form: min c'x, subject to Ax >= b, x binary 0-1  vector
    # c positive n-vector, b positive m-vector, A nonnegative m x n matrix
    # OUT = BILP_MIN_EXACT(c,b,A) is a summary dictionary of named results which can be queried
    # E.G., The exact BILP binary solution vector is OUT["EXACT_BILP_SOLN"]
    # Bob Agnew, raagnew1@gmail.com, raagnew.com

    Time_In = time.time()

    m = len(b)
    n = len(c)
    
    b_l = b
    
    integrality = np.ones_like(c)    
    bounds = Bounds(0,1)
    constraints = LinearConstraint(A,lb=b_l)    
    
    opt = milp(c,integrality=integrality,bounds=bounds,
    constraints=constraints)

    solution = np.round(opt.x).astype(int)
    obj = opt.fun
    Ax = A @ solution
    constraints_satisfied = np.sum(Ax >= b)

    df = pd.DataFrame({
        "Ax": Ax,
        " >= ": [" >= "] * m,
        "b": b
    },index=list(range(1,(m+1))))
    
    Time_Out = time.time()
    del_time = Time_Out - Time_In

    result = {
        "NUMBER_VARIABLES": n,
        "NUMBER_CONSTRAINTS": m,
        "EXECUTION_SECONDS": del_time,
        "EXACT_BILP_SOLN": solution,
        "EXACT_BILP_MIN": obj,
        "CONSTRAINTS_SATISFIED": constraints_satisfied,
        "CONSTRAINT_DETAILS": df
    }

    return result

def BILP_MIN_APPROX(c, b, A):
    # Simple Minimizing Binary Integer Linear Programming Code
    # Standard form: min c'x, subject to Ax >= b, x binary 0 or 1
    # c n-vector, b positive m-vector, A nonnegative m x n matrix
    # Transformed to equivalent maximizing problem for solution
    # Minimize cost while covering constraints
    # Includes set covering problems
    # Vector of ones must be feasible
    # OUT = BILP_MIN_APPROX(c,b,A) is a summary dictionary of named results which can be queried
    # E.G., The approximate BILP binary solution vector is OUT["APPROX_BILP_SOLN"]
    # Bob Agnew, raagnew1@gmail.com, raagnew.com

    Time_In = time.time()
    m = len(b)
    n = len(c)
    ones = [1] * n  # Vector of ones, assumed feasible
    
    obj_ones = np.dot(c,ones) # Objective function with all ones
    
    A_ones = np.dot(A,ones)
    out = BILP_MAX_APPROX(c, (A_ones - b), A)
    x = ones - out['APPROX_BILP_SOLN']
    obj = obj_ones - out['APPROX_BILP_MAX']
    bnd = obj_ones - out['APPROX_DUAL_LP_MIN']

    Ax = np.dot(A, x)
    df = pd.DataFrame({
        "Ax": Ax,
        " >= ": [" >= "] * m,
        "b": b
    },index=list(range(1,(m+1))))

    Time_Out = time.time()
    del_time = Time_Out - Time_In  # Imprecise for small deltas

    result = {
        "NUMBER_VARIABLES": n,
        "NUMBER_CONSTRAINTS": m,
        "EXECUTION_SECONDS": del_time,
        "SOLVER_ITERATIONS": out['SOLVER_ITERATIONS'],
        "SOLVER_EVALUATIONS": out['SOLVER_EVALUATIONS'],
        "APPROX_BILP_SOLN": x,
        "APPROX_BILP_MIN": obj,
        "APPROX_DUAL_LP_MAX": bnd,
        "APPROXIMATION_RATIO": obj / bnd,
        "CONSTRAINTS_SATISFIED": np.sum(Ax >= b),
      # "CONSTRAINTS_SATISFIED": sum(1 for val, bi in zip(Ax, b) if val >= bi),
        "CONSTRAINT_DETAILS": df
    }

    return result

print()
print()
print('Small Textbook Maximization Problem')
# Small Textbook Maximization Problem
# max 8*x1 + 11*x2 + 6*x3 + 4*x4
# subject to:
# 5*x1 + 7*x2 + 3*x4 <= 14
# 8*x1 + 4*x3 + 4*x4 <= 12
# 2*x1 + 10*x2 + 6*x3 + 4*x4 <= 15
# x1, x2, x3, x4 in {0,1}
c = [8, 11, 6, 4]
b = [14, 12, 15]
A = [
    [5, 7, 0, 3],
    [8, 0, 4, 4],
    [2, 10, 6, 4]
]
print()
print()
print('Maximization Solutions')

OUT = BILP_MAX_EXACT(c, b, A)

print()
print('Dictionary Summary for Exact Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

OUT = BILP_MAX_APPROX(c, b, A)

print()
print('Dictionary Summary for Approximate Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

print()
print('Mirror-Image Minimization Solutions')

OUT = BILP_MIN_EXACT(c, b, A)

print()
print('Dictionary Summary for Exact Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

OUT = BILP_MIN_APPROX(c, b, A)

print()
print('Dictionary Summary for Approximate Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

print()
print()
print('Multidimensional Knapsack Problems from J. E. Beasely OR-Library: people.brunel.ac.uk/~mastjjb/jeb/info.html')
print()
print()
print('Problem mknap1_7 in Beasley OR-Libary')
print()
with open('c:/BILP/mknap1_7.txt','r') as FILE:
    data_str = FILE.read()
    v = np.fromstring(data_str,sep=' ')
c = v[3:53]
b = v[303:309]
A = v[53:303].reshape(len(b),len(c))
print()
print('Maximization Solutions')

OUT = BILP_MAX_EXACT(c, b, A)

print()
print('Dictionary Summary for Exact Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

OUT = BILP_MAX_APPROX(c, b, A)

print()
print('Dictionary Summary for Approximate Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

print()
print('Mirror-Image Minimization Solutions')

OUT = BILP_MIN_EXACT(c, b, A)

print()
print('Dictionary Summary for Exact Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

OUT = BILP_MIN_APPROX(c, b, A)

print()
print('Dictionary Summary for Approximate Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()
 
print()
print('Problem mknapcb9_30 in Beasley OR-Libary')
print()
with open('c:/BILP/mknapcb9_30.txt','r') as FILE:
    data_str = FILE.read()
    v = np.fromstring(data_str,sep=' ')
c = v[3:503]
b = v[15503:15533]
A = v[503:15503].reshape(len(b),len(c))
print()
print('Maximization Solution')

OUT = BILP_MAX_APPROX(c, b, A)
del OUT['APPROX_BILP_SOLN']
del OUT['CONSTRAINT_DETAILS']

print()
print('Dictionary Summary for Approximate Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

print()
print('Mirror-Image Minimization Solution')

OUT = BILP_MIN_APPROX(c, b, A)
del OUT['APPROX_BILP_SOLN']
del OUT['CONSTRAINT_DETAILS']

print()
print('Dictionary Summary for Approximate Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

print()
print('Randomized Inputs for Large Randomized Problem (50000 variables, 300 constraints')
print('Using same randomized data as R')
print()
with open('c:/BILP/Large_Problem.txt','r') as FILE:
    data_str = FILE.read()
    v = np.fromstring(data_str,sep=' ')
c = v[0:50000]
b = v[50000:50300]
A = v[50300:15050300].reshape(len(b),len(c))

print('Maximization Solution')

OUT = BILP_MAX_APPROX(c, b, A)
del OUT['APPROX_BILP_SOLN']
del OUT['CONSTRAINT_DETAILS']

print()
print('Dictionary Summary for Approximate Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()

print()
print('Mirror-Image Minimization Solution')

OUT = BILP_MIN_APPROX(c, b, A)
del OUT['APPROX_BILP_SOLN']
del OUT['CONSTRAINT_DETAILS']

print()
print('Dictionary Summary for Approximate Solution')
print()
for key,value in OUT.items():
    print(f"{key}")
    print(f"{value}")
    print()



