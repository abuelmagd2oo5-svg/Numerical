import math
from sympy import symbols, sympify, lambdify, diff

x = symbols('x')

# ================= CREATE FUNCTION =================
def create_function(expr_str):
    expr = sympify(expr_str)
    return lambdify(x, expr, "math"), expr

# ================= BISECTION =================
def bisection(f, a, b, tol):
    if f(a) * f(b) >= 0:
        print(" Invalid interval!")
        return

    i = 1
    while True:
        mid = (a + b) / 2
        print(f"Iter {i}: a={a}, b={b}, mid={mid}, f(mid)={f(mid)}, error={abs(b-a)}")

        if abs(b - a) < tol:
            print("\n Root =", mid)
            return

        if f(a) * f(mid) < 0:
            b = mid
        else:
            a = mid

        i += 1

# ================= SECANT =================
def secant(f, x0, x1, tol):
    i = 1
    while True:
        if abs(f(x1) - f(x0)) < 1e-12:
            print(" Division by zero!")
            return

        x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0))

        print(f"Iter {i}: x={x2}, f(x)={f(x2)}, error={abs(x2-x1)}")

        if abs(x2 - x1) < tol:
            print("\n Root =", x2)
            return

        x0, x1 = x1, x2
        i += 1

# ================= SIMPLE ITERATION =================
def fixed_point(g_expr_str, x0, tol):
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

# ================= NEWTON =================
def newton(f_expr_str, x0, tol):
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

#================== jacobi =============
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

# ================= MAIN (SMART MENU) =================
print("\n===== Numerical Methods Solver =====")

print("Choose method:")
print("1. Bisection")
print("2. Secant")
print("3. Simple Iteration")
print("4. Newton-Raphson")
print("5. Jacobi Method")

choice = int(input("Enter choice: "))

# ================= METHOD ROUTING =================

if choice == 1:
    tol = float(input("Enter tolerance: "))
    f_str = input("Enter f(x): ")
    f, _ = create_function(f_str)

    a, b = map(float, input("Enter a b: ").split())
    bisection(f, a, b, tol)

elif choice == 2:
    tol = float(input("Enter tolerance: "))
    f_str = input("Enter f(x): ")
    f, _ = create_function(f_str)

    x0, x1 = map(float, input("Enter x0 x1: ").split())
    secant(f, x0, x1, tol)

elif choice == 3:
    tol = float(input("Enter tolerance: "))
    g_str = input("Enter g(x): ")
    x0 = float(input("Enter initial guess x0: "))
    fixed_point(g_str, x0, tol)

elif choice == 4:
    tol = float(input("Enter tolerance: "))
    f_str = input("Enter f(x): ")
    x0 = float(input("Enter initial guess x0: "))
    newton(f_str, x0, tol)

elif choice == 5:
    tol = float(input("Enter tolerance: "))
    jacobi()

else:
    print(" Invalid choice!")