import numpy as np
from matplotlib import pyplot as plot
import scipy.integrate as integrate

def quadrature(f, a, b):
    return integrate.trapezoid([f(x) for x in np.linspace(a, b, 1000)], dx=(b-a)/1000)

def x_i(i, n):
    return 2.0 / n * i

def e_i(i, x, n):
    range_length = 2 / n
    if x_i(i - 1, n) <= x <= x_i(i, n):
        return (x - x_i(i - 1, n)) / range_length
    elif x_i(i, n) < x <= x_i(i + 1, n):
        return (x_i(i + 1, n) - x) / range_length
    return 0

def e_i_deriver(i, x, n):
    range_length = 2 / n
    if x_i(i - 1, n) <= x <= x_i(i, n):
        return 1 / range_length
    elif x_i(i, n) < x <= x_i(i + 1, n):
        return -1 / range_length
    return 0

def L_v(i, n):
    integral = lambda x: e_i(i, x, n) * (np.sin(x) - 2)
    return quadrature(integral, 0, 2) - 6 * e_i(i, 2, n)
        #∫v(sinx-2)dx - 6v(2)

def B_w_v(i, j, n):
    # integral1 = ∫v'w'dx
    integral1 = lambda x: e_i_deriver(i, x, n) * e_i_deriver(j, x, n)
    # integral2 = ∫vwdx
    integral2 = lambda x: e_i(i, x, n) * e_i(j, x, n)
    # B(w,v) = v(2)w(2) + ∫v'w'dx - ∫vwdx
    return e_i(i, 2, n) * e_i(j, 2, n) + quadrature(integral1, 0, 2) - quadrature(integral2, 0, 2)

def u(U_i, x, n):
    result = 0
    for i in range(1, len(U_i) + 1):
        result += U_i[i - 1] * e_i(i, x, n)

    return result - 2
    # u = w-2 # co wynika z u(0)=-2

def solve():
    n = int(input("podaj n: "))
    A = np.zeros((n, n))
    B = np.zeros(n)

    for i in range(1, n + 1):
        for j in range(1, n + 1):
            A[i - 1, j - 1] = B_w_v(j, i, n)

    for i in range(1, n + 1):
        B[i - 1] = L_v(i, n)

    U_i = np.linalg.solve(A, B)
    x = np.linspace(0, 2, 100)
    values = [u(U_i, xi, n) for xi in x]

    plot.plot(x, values, label="u(x)")
    plot.title("Wykres funkcji u(x)")
    plot.xlabel("x")
    plot.ylabel("u(x)")
    plot.grid()
    plot.legend()
    plot.show()

solve()
