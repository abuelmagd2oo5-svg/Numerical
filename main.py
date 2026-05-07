import math
import re
from sympy import symbols, sympify, lambdify, diff, exp, E

x = symbols('x')


def safe_expr(expr_str):
    expr_str = expr_str.replace(" ", "")
    expr_str = re.sub(r"(\d)([a-zA-Z])", r"\1*\2", expr_str)  # Add implicit multiplication
    expr_str = expr_str.replace("^", "**")
    expr_str = expr_str.replace("e**", "exp")
    expr_str = expr_str.replace("e^", "exp")
    return expr_str


def create_function(expr_str):
    expr_str = safe_expr(expr_str)
    expr = sympify(expr_str, locals={"exp": exp, "E": E})
    return lambdify(x, expr, "math"), expr


def bisection(f, a, b, tol):
    if f(a) * f(b) >= 0:
        print("Invalid interval!")
        return
    i = 1
    while True:
        mid = (a + b) / 2
        print(f"Iter {i}: mid={mid}, f(mid)={f(mid)}, error={abs(b - a)}")
        if abs(b - a) < tol:
            print(f"\nRoot = {mid}")
            return
        if f(a) * f(mid) < 0:
            b = mid
        else:
            a = mid
        i += 1
        if i > 1000:
            print("Not converging!")
            return


def secant(f, x0, x1, tol):
    i = 1
    while True:
        f0 = f(x0)
        f1 = f(x1)
        if abs(f1 - f0) < 1e-14:
            print("Division by zero risk!")
            return
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        print(f"Iter {i}: x={x2}, error={abs(x2 - x1)}")
        if abs(x2 - x1) < tol:
            print(f"\nRoot = {x2}")
            return
        x0, x1 = x1, x2
        i += 1
        if i > 1000:
            print("Not converging!")
            return


def fixed_point(g_expr_str, x0, tol):
    g_expr = sympify(safe_expr(g_expr_str), locals={"exp": exp, "E": E})
    g = lambdify(x, g_expr, "math")
    g_prime = lambdify(x, diff(g_expr, x), "math")
    try:
        gp = abs(g_prime(x0))
        print(f"|g'(x0)| = {gp}")
        if gp >= 1:
            print("Warning: may NOT converge")
        else:
            print("Likely converges")
    except:
        pass
    i = 1
    while True:
        x1 = g(x0)
        print(f"Iter {i}: x={x1}, error={abs(x1 - x0)}")
        if abs(x1 - x0) < tol:
            print(f"\nRoot = {x1}")
            return
        x0 = x1
        i += 1
        if i > 1000:
            print("Not converging!")
            return


def newton(expr_str, x0, tol):
    expr = sympify(safe_expr(expr_str), locals={"exp": exp, "E": E})
    f = lambdify(x, expr, "math")
    f1 = lambdify(x, diff(expr, x), "math")
    f2 = lambdify(x, diff(expr, x, 2), "math")
    try:
        value = abs(f(x0) * f2(x0)) / (f1(x0) ** 2)
        print(f"|f f'' / (f')^2| = {value}")
        if value >= 1:
            print("Warning: may NOT converge")
    except:
        pass
    i = 1
    while True:
        if abs(f1(x0)) < 1e-12:
            print("Derivative zero!")
            return
        x1 = x0 - f(x0) / f1(x0)
        print(f"Iter {i}: x={x1}, error={abs(x1 - x0)}")
        if abs(x1 - x0) < tol:
            print(f"\nRoot = {x1}")
            return
        x0 = x1
        i += 1
        if i > 1000:
            print("Not converging!")
            return


def jacobi():
    """Generalized Jacobi method for n x n system"""
    n = int(input("Enter number of equations (n): "))

    A = []
    b = []

    print("\nEnter augmented matrix (coefficients and RHS):")
    for i in range(n):
        row = list(map(float, input(f"Row {i + 1}: ").split()))
        if len(row) != n + 1:
            print(f"Error: Expected {n + 1} values ({n} coefficients + RHS)")
            return
        A.append(row[:-1])
        b.append(row[-1])

    tol = float(input("\nEnter tolerance: "))
    x = list(map(float, input("Enter initial guess (space separated): ").split()))

    if len(x) != n:
        print(f"Error: Expected {n} initial values")
        return

    # Check for diagonal dominance
    diagonally_dominant = True
    for i in range(n):
        diag = abs(A[i][i])
        off_diag_sum = sum(abs(A[i][j]) for j in range(n) if j != i)
        if diag <= off_diag_sum:
            diagonally_dominant = False
            print(f"Warning: Row {i + 1} is not diagonally dominant")

    if not diagonally_dominant:
        proceed = input("Matrix not diagonally dominant. Continue? (y/n): ")
        if proceed.lower() != 'y':
            return

    # Check for zero diagonals
    for i in range(n):
        if abs(A[i][i]) < 1e-12:
            print(f"Error: A[{i + 1}][{i + 1}] = 0. Cannot solve with Jacobi")
            return

    iteration = 1
    while True:
        x_new = []

        # Calculate new values
        for i in range(n):
            s = 0
            for j in range(n):
                if j != i:
                    s += A[i][j] * x[j]
            x_new.append((b[i] - s) / A[i][i])

        # Print iteration results
        print(f"\nIter {iteration}:")
        for i in range(n):
            print(f"  x{i + 1} = {x_new[i]:.8f}", end=" ")
        print()

        # Check convergence
        max_error = max(abs(x_new[i] - x[i]) for i in range(n))
        print(f"  Max error = {max_error:.2e}")

        if max_error < tol:
            print("\n✓ Converged!")
            print("Solution:")
            for i in range(n):
                print(f"  x{i + 1} = {x_new[i]:.8f}")
            return

        x = x_new
        iteration += 1

        if iteration > 1000:
            print("\nNot converging within 1000 iterations!")
            print("Current approximation:", x)
            return


def gauss_seidel():
    """Generalized Gauss-Seidel method for n x n system"""
    n = int(input("Enter number of equations (n): "))

    A = []
    b = []

    print("\nEnter augmented matrix (coefficients and RHS):")
    for i in range(n):
        row = list(map(float, input(f"Row {i + 1}: ").split()))
        if len(row) != n + 1:
            print(f"Error: Expected {n + 1} values")
            return
        A.append(row[:-1])
        b.append(row[-1])

    tol = float(input("\nEnter tolerance: "))
    x = list(map(float, input("Enter initial guess (space separated): ").split()))

    if len(x) != n:
        print(f"Error: Expected {n} initial values")
        return

    # Check for zero diagonals
    for i in range(n):
        if abs(A[i][i]) < 1e-12:
            print(f"Error: Zero diagonal at row {i + 1}")
            return

    iteration = 1
    while True:
        x_old = x.copy()

        # Update using latest available values
        for i in range(n):
            s = 0
            for j in range(n):
                if j != i:
                    s += A[i][j] * x[j]
            x[i] = (b[i] - s) / A[i][i]

        # Print iteration results
        print(f"\nIter {iteration}:")
        for i in range(n):
            print(f"  x{i + 1} = {x[i]:.8f}", end=" ")
        print()

        # Check convergence
        max_error = max(abs(x[i] - x_old[i]) for i in range(n))
        print(f"  Max error = {max_error:.2e}")

        if max_error < tol:
            print("\n✓ Converged!")
            print("Solution:")
            for i in range(n):
                print(f"  x{i + 1} = {x[i]:.8f}")
            return

        iteration += 1

        if iteration > 1000:
            print("\nNot converging within 1000 iterations!")
            print("Current approximation:", x)
            return


def midpoint(expr_str, a, b, n):
    f = lambdify(x, sympify(safe_expr(expr_str)), "math")
    h = (b - a) / n
    total = 0
    for i in range(n):
        mid = a + (i + 0.5) * h
        total += f(mid)
    print(f"Midpoint Result = {h * total}")


def trapezoidal_integration(expr_str, a, b, n):
    expr_str = safe_expr(expr_str)
    expr = lambdify(x, sympify(expr_str), "math")
    h = (b - a) / n
    total = expr(a) + expr(b)

    for i in range(1, n):
        xi = a + i * h
        total += 2 * expr(xi)

    result = (h / 2) * total
    print(f"Trapezoidal Result = {result}")
    return result


def lagrange():
    n = int(input("Number of points: "))
    x_vals = []
    y_vals = []
    for i in range(n):
        xi, yi = map(float, input(f"Point {i + 1}: ").split())
        x_vals.append(xi)
        y_vals.append(yi)
    x_target = float(input("x to interpolate: "))
    result = 0.0
    for i in range(n):
        term = y_vals[i]
        for j in range(n):
            if i != j:
                if x_vals[i] == x_vals[j]:
                    print("Error: duplicate x values are not allowed!")
                    return
                term *= (x_target - x_vals[j]) / (x_vals[i] - x_vals[j])
        result += term
    print(f"Result = {result}")


# Main menu
while True:
    print("\n===== Numerical Methods =====")
    print("1. Bisection")
    print("2. Secant")
    print("3. Fixed Point")
    print("4. Newton")
    print("5. Jacobi (Generalized)")
    print("6. Gauss-Seidel (Generalized)")
    print("7. Midpoint")
    print("8. Trapezoidal")
    print("9. Lagrange")
    print("0. Exit")

    try:
        c = int(input("Choice: ").strip())
    except:
        print("Invalid input → enter number")
        continue

    try:
        if c == 0:
            print("Exiting...")
            break

        elif c == 1:
            f, _ = create_function(input("f(x): "))
            a, b = map(float, input("a b: ").split())
            tol = float(input("tol: "))
            bisection(f, a, b, tol)

        elif c == 2:
            f, _ = create_function(input("f(x): "))
            x0, x1 = map(float, input("x0 x1: ").split())
            tol = float(input("tol: "))
            secant(f, x0, x1, tol)

        elif c == 3:
            fixed_point(
                input("g(x): "),
                float(input("x0: ")),
                float(input("tol: "))
            )

        elif c == 4:
            newton(
                input("f(x): "),
                float(input("x0: ")),
                float(input("tol: "))
            )

        elif c == 5:
            jacobi()

        elif c == 6:
            gauss_seidel()

        elif c == 7:
            midpoint(
                input("f(x): "),
                float(input("a: ")),
                float(input("b: ")),
                int(input("n: "))
            )

        elif c == 8:
            trapezoidal_integration(
                input("f(x): "),
                float(input("a: ")),
                float(input("b: ")),
                int(input("n: "))
            )

        elif c == 9:
            lagrange()

        else:
            print("Invalid choice")

    except Exception as e:
        print(f"Error: {e}")