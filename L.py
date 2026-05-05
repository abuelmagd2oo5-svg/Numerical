from sympy import *
import numpy as np


def bisection(a,b,expe_string,E):
    x =symbols('x')
    expr = lambdify(x,sympify(expe_string))
    f_A = expr(a)
    f_B = expr(b)
    A = 0
    B = 0
    if f_A*f_B > 0:
        print("f(a) and f(b) need  to have different signs")
        return
    if f_A < 0:
        A = a
        B = b
    else:
        A = b
        B = a
    while abs(B-A) > E:
        mid = (A+B)/2
        f_mid = expr(mid)
        print(f"{A}  {B} {mid} {f_mid} {abs(B-A)}")
        if  f_mid < 0:
            A = mid
        else:
            B = mid

# bisection(0,1,input("exp: "),10**-5)



def secant (x1,x2,expe_string,E):
    x =symbols('x')
    expr = lambdify(x,sympify(expe_string))
    Xi =x1
    Xi1 = x2
    n = 0
    while abs(Xi1 - Xi) > E:
         f_A = expr(Xi)
         print(f"{n} {Xi} {f_A} {abs(Xi1 - Xi)}")
         f_B = expr(Xi1)
         i1 =np.linalg.det(np.array([[Xi,f_A],[Xi1,f_B]]))/(f_B -f_A)
         Xi = Xi1
         Xi1 = i1
         n +=1

# x**3-2*(x**2)-5
# secant(1,4,input("exp: "),10**-5)

def newton_raphson(X,expe_string,E):
    x =symbols('x')
    expr = lambdify(x,sympify(expe_string))
    d1f = lambdify(x,diff(sympify(expe_string),x))
    d2f = lambdify(x,diff(sympify(expe_string),x,2))
    if abs(expr(X)*d2f(X)) <(d1f(X)**2):
        xi = X-(expr(X)/d1f(X))
        print(f"{xi} {abs(xi - X)}")
        while (abs(xi - X) > E):
            X = xi
            xi = X - (expr(X)/d1f(X))
            print(f"{xi} {abs(xi - X)}")
    elif input("f(x) may not converg will you continuo for n y,n: ") =='y':
        n = int(input("n = "))
        xi = X - (expr(X) / d1f(X))
        print(f"{xi} {abs(xi - X)}")
        for i in range(n):
            xi = xi - (expr(xi) / d1f(xi))
            print(f"{xi} {abs(xi - X)}")

# newton_raphson(1,input("exp: "),10*-3)
# # x**3+4*(x**2)-10


def midpoint_integration(expe_string, a, b, n=1):
    x = symbols('x')
    expr = lambdify(x, sympify(expe_string))
    h = (b - a) / n
    sum = 0
    for i in range(n):
        midpoint = a + (i + 0.5) * h
        sum += expr(midpoint)
    integ = h * sum
    print(f"integrating {a} to  {b} with  {n} intervals ,Result: {integ}\n")

midpoint_integration(input("expr: "), int(input("a: ")),int(input("b: ")), n=10)

#=========================================================================================================
def create_function(expr_str):
    x = symbols('x')
    expr = sympify(expr_str)
    return lambdify(x, expr, "math"), expr

def simple_iteration (g_expr_str, x0, tol):
    x = symbols('x')
    g_expr = sympify(g_expr_str)
    g = lambdify(x, g_expr, "math")
    g_derivative = lambdify(x, diff(g_expr, x), "math")

    # Convergence Check
    try:
        gp = abs(g_derivative(x0))
        print(f"|g'(x0)| = {gp}")

        if gp >= 1:
            print(" May NOT converge")
        else:
            print(" Likely converges")
    except:
        print(" Cannot evaluate convergence")

    i = 1
    while True:
        x1 = g(x0)
        print(f"Iter {i}: x={x1}, error={abs(x1-x0)}")

        if abs(x1 - x0) < tol:
            print("\n Root =", x1)
            return

        x0 = x1
        i += 1

        if i > 1000:
            print(" Not converging!")
            return

def jacobi():
    print("\nEnter coefficients for system:")

    print("Equation 1: a1*x1 + b1*x2 + c1*x3 = d1")
    a1, b1, c1, d1 = map(float, input("Enter a1 b1 c1 d1: ").split())

    print("Equation 2: a2*x1 + b2*x2 + c2*x3 = d2")
    a2, b2, c2, d2 = map(float, input("Enter a2 b2 c2 d2: ").split())

    print("Equation 3: a3*x1 + b3*x2 + c3*x3 = d3")
    a3, b3, c3, d3 = map(float, input("Enter a3 b3 c3 d3: ").split())

    # ===== Convergence Check =====
    if (abs(a1) > abs(b1) + abs(c1) and
        abs(b2) > abs(a2) + abs(c2) and
        abs(c3) > abs(a3) + abs(b3)):
        print(" System is diagonally dominant → likely converges")
    else:
        print(" Not diagonally dominant → may NOT converge")

    # Initial guesses
    x1 = x2 = x3 = 0.0
    tol = float(input("Enter tolerance: "))

    i = 1
    while True:
        x1_new = (d1 - b1*x2 - c1*x3) / a1
        x2_new = (d2 - a2*x1 - c2*x3) / b2
        x3_new = (d3 - a3*x1 - b3*x2) / c3

        print(f"Iter {i}: x1={x1_new}, x2={x2_new}, x3={x3_new}")

        # stopping condition
        if (abs(x1_new - x1) < tol and
            abs(x2_new - x2) < tol and
            abs(x3_new - x3) < tol):
            print("\n Solution:")
            print(f"x1 = {x1_new}, x2 = {x2_new}, x3 = {x3_new}")
            return

        x1, x2, x3 = x1_new, x2_new, x3_new
        i += 1

        if i > 1000:
            print(" Not converging!")
            return

def new_newton_raphson(f_expr_str, x0, tol):
    f_expr = sympify(f_expr_str)

    f = lambdify(x, f_expr, "math")
    f1 = lambdify(x, diff(f_expr, x), "math")
    f2 = lambdify(x, diff(f_expr, x, 2), "math")

    # Convergence Condition
    try:
        value = abs(f(x0) * f2(x0)) / (f1(x0)**2)
        print(f"|f(x)f''(x)/(f'(x))^2| = {value}")

        if value >= 1:
            print(" May NOT converge")
        else:
            print(" Likely converges")
    except:
        print(" Cannot evaluate condition")
    i = 1
    while True:
        if abs(f1(x0)) < 1e-12:
            print(" f'(x)=0 → Division error")
            return
        x1 = x0 - f(x0)/f1(x0)
        print(f"Iter {i}: x={x1}, error={abs(x1-x0)}")
        if abs(x1 - x0) < tol:
            print("\n Root =", x1)
            return
        x0 = x1
        i += 1
        if i > 1000:
            print(" Not converging!")
            return
