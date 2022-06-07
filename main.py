import math
import operator as op
import numpy as np


def swap_rows(matrix, initial, final):
    if initial != final:
        tmp = np.copy(matrix)
        matrix[initial] = tmp[final]
        matrix[final] = tmp[initial]
    return matrix


def swap_cols(matrix, initial, final):
    if initial != final:
        tmp = np.copy(matrix)
        matrix[::, initial] = tmp[::, final]
        matrix[::, final] = tmp[::, initial]
    return matrix


# swap cols and rows in same time
def swap_cols_rows(matrix, init_col, final_col, init_row, final_row):
    matrix = swap_cols(matrix, init_col, final_col)
    matrix = swap_rows(matrix, init_row, final_row)
    return matrix


def pivoting(L, U, P, Q, x):
    abs_U = abs(U)
    if abs_U[x:, x:].sum() < pow(10, -14):
        return 0, L, U, P, Q

    i, j = np.where(abs_U[x:, x:] == abs_U[x:, x:].max())
    i[0] += x
    j[0] += x

    L = swap_cols_rows(L, j[0], x, i[0], x)
    U = swap_cols_rows(U, j[0], x, i[0], x)
    P = swap_rows(P, i[0], x)
    Q = swap_cols(Q, j[0], x)

    return 1, L, U, P, Q


def decomposition(matrix):
    U = np.copy(matrix)
    L = np.zeros((len(matrix), len(matrix)))
    P = np.eye(len(matrix))
    Q = np.eye(len(matrix))

    for i in range((len(matrix)) - 1):
        achived, L, U, P, Q = pivoting(L, U, P, Q, i)

        if achived == 0:
            break

        T = np.eye(len(matrix))
        for k in range(i + 1, len(matrix)):
            L[k, i] = U[k, i] / U[i, i]
            T[k, i] = (-1) * L[k, i]

        U = np.dot(T, U)

    L = L + np.eye(len(matrix))
    return L, U, P, Q


def solve_slae(A, b):
    L, U, P, Q = decomposition(A)
    x, y = np.zeros(len(b)), np.zeros(len(b))
    Pb = np.dot(P, b)
    # Ly=Pb
    for i in range(len(Pb)):
        y[i] = Pb[i]

        for j in range(i):
            y[i] -= y[j] * L[i, j]
    # Ux=y
    for i in range(len(Pb) - 1, -1, -1):
        x[i] = y[i]

        for j in range(i + 1, len(Pb)):
            x[i] -= x[j] * U[i, j]

        x[i] /= U[i, i]

    return np.dot(Q, x)


lim_a = 1.8
lim_b = 2.9
alpha = 0
beta = 4.0 / 7
exact = 57.48462064655285571820619434508191055583


def f(x):
    return 4 * math.cos(2.5*x) * math.exp(x*4.0/7) + 2.5 * math.sin(5.5*x) * math.exp(-3.0*x/5) + 4.3 * x

def moment_0(a: float, b: float):
    return (-7.0/3) * (math.pow((2.9 - b), 3.0/7) - math.pow((2.9 - a), 3.0/7))

def moment_1(a: float,b: float):
    return (-7.0/300) * (math.pow((2.9 - b), 3.0/7)*(30.0*b + 203.0) - math.pow((2.9 - a), 3.0/7)*(30.0*a + 203.0))

def moment_2(a: float,b: float):
    return (-7.0/25500) * (math.pow((2.9 - b), 3.0/7)*(1500.0*b*b + 6090.0*b + 41209.0) - math.pow((2.9 - a), 3.0/7)*(1500.0*a*a + 6090.0*a + 41209.0))

def moment_3(a: float,b: float):
    return (-7.0/2040000) * (math.pow((2.9-b), 3.0/7) * (85000.0*b*b*b + 304500.0*b*b + 1236270.0*b + 8365427.0) - math.pow((2.9-a), 3.07) * (85000.0*a*a*a + 304500.0*a*a + 1236270.0*a + 8365427.0))

def moment_4(a: float,b: float):
    return (-7.0/15810000) * (math.pow((2.9 - b), 3.0/7) * (5100000.0*b*b*b*b + 17255000.0*b*b*b + 61813500.0*b*b + 250962810.0*b + 1698181681) - math.pow((2.9 - a), 3.0/7) * (5100000.0*a*a*a*a + 17255000.0*a*a*a + 61813500.0*a*a + 250962810.0*a + 1698181681))

def moment_5(a: float,b: float):
    return (-7.0/12015600000) * (math.pow((2.9-b), 3.0/7) * (316200000.0*b*b*b*b*b + 1035300000.0*b*b*b*b + 3502765000.0*b*b*b + 12548140500.0*b*b + 50945450430.0*b + 344730881243.0) - math.pow((2.9-a), 3.0/7) * (316200000.0*a*a*a*a*a + 1035300000.0*a*a*a*a + 3502765000.0*a*a*a + 12548140500.0*a*a + 50945450430.0*a + 344730881243.0))

def cordano(coefs):
    a = coefs[2]
    b = coefs[1]
    c = coefs[0]

    p = b - (a ** 2) / 3
    q = c + (2 * (a ** 3)) / 27 - (a * b) / 3
    D = (q ** 2) / 4.0 + (p ** 3) / 27.0

    if D < 0:
        if q < 0:
            fi = math.atan(2.0 * math.sqrt(-D) / (-q))
        elif q > 0:
            fi = math.atan(2.0 * math.sqrt(-D) / (-q) + math.pi)
        else:
            fi = math.pi / 2.0

        x1 = 2.0 * math.sqrt(-p / 3.0) * math.cos(fi / 3.0) - a / 3.0
        x2 = 2.0 * math.sqrt(-p / 3.0) * math.cos(fi / 3.0 + 2.0 * math.pi / 3.0) - a / 3.0
        x3 = 2.0 * math.sqrt(-p / 3.0) * math.cos(fi / 3.0 + 4.0 * math.pi / 3.0) - a / 3.0
        return np.array([x1, x2, x3], float)

    elif D > 0:
        x1 = 0
        if (-q) / 2.0 + math.pow(D, 1.0 / 2.0) < 0:
            x1 += -math.pow((q) / 2.0 - math.pow(D, 1.0 / 2.0), 1.0 / 3.0)
        else:
            x1 += math.pow((-q) / 2.0 + math.pow(D, 1.0 / 2.0), 1.0 / 3.0)

        if (-q) / 2.0 - math.pow(D, 1.0 / 2.0) < 0:
            x1 += -math.pow(q / 2.0 + math.pow(D, 1.0 / 2.0), 1.0 / 3.0) - a / 3.0
        else:
            x1 += math.pow(-q / 2.0 - math.pow(D, 1.0 / 2.0), 1.0 / 3.0) - a / 3.0

        return np.array([x1], float)

    else:
        x1 = 2 * math.pow(-q / 2.0, 1.0 / 3.0) - a / 3.0
        x2 = -math.pow(-q / 2.0, 1.0 / 3.0) - a / 3.0

        return np.array([x1, x2], float)


def newton_cotes(a, b, parts):
    step = (b - a) / parts
    integral = 0.0
    for i in range(parts):
        b = a + step

        node_0 = a
        node_1 = a + (b - a) / 2
        node_2 = b
        nodes = np.array([node_0, node_1, node_2], float)

        moment0 = moment_0(a, b)
        moment1 = moment_1(a, b)
        moment2 = moment_2(a, b)
        moments = np.array([moment0, moment1, moment2], float)

        matrix = np.array([1, 1, 1], float)
        matrix = np.vstack((matrix, nodes))
        matrix = np.vstack((matrix, np.multiply(nodes, nodes)))

        solution = solve_slae(matrix, moments)
        for j in range(len(solution)):
            integral += solution[j] * f(nodes[j])

        a += step
    return integral


def gauss(a, b, step):
    lim1 = a
    integral = 0.0
    for i in range(1, math.ceil((lim_b - b) / step) + 2, 1):
        moment0 = moment_0(a, b)
        moment1 = moment_1(a, b)
        moment2 = moment_2(a, b)
        moment3 = moment_3(a, b)
        moment4 = moment_4(a, b)
        moment5 = moment_5(a, b)
        moments = np.array([moment0, moment1, moment2, moment3, moment4, moment5], float)

        node_0 = a
        node_1 = a + (b - a) / 2
        node_2 = b
        nodes = np.array([node_0, node_1, node_2], float)

        matrix = np.array([[moment0, moment1, moment2],
                           [moment1, moment2, moment3],
                           [moment2, moment3, moment4]], dtype=float)

        y = np.array([-moment3, -moment4, -moment5], dtype=float)

        solution = solve_slae(matrix, y)

        p = cordano(solution)

        p.sort()

        A = np.array([[1, 1, 1]], float)
        A = np.vstack((A, p))
        A = np.vstack((A, np.multiply(p, p)))

        solution = solve_slae(A, moments[:3])

        for j in range(len(solution)):
            integral += solution[j] * f(p[j])

        a = lim1 + (i) * step
        b = lim1 + (i + 1) * step

    return integral

np.set_printoptions(formatter={'float': '{: 0.4f}'.format})
print("--------Ньютон-Котс:--------\n")
method_mistake = 0.8
integral = newton_cotes(lim_a, lim_b, 1)
print(f"результат: {integral},\nточная погрешность: {abs(exact - integral)},\nметодическая погрешность : {method_mistake}\n")

a = lim_a
b = lim_b
e = 10
parts = 1
L = 2
degree = 4


print("--------Cоставная ИКФ c оптимальным шагом:--------\n")
step = (lim_b - lim_a) / 2
parts_opt = math.ceil((lim_b - lim_a) / (step * L * ((0.00001 / e) ** (1 / degree))))

a = lim_a
b = lim_b

integral_opt = newton_cotes(a, b, parts_opt)
e = abs((integral_opt - integral) / (math.pow(L, degree) - 1))

print(f"\nсоставная ИКФ с оптимальным шагом : {integral_opt},\nточная погрешность: {abs(exact - integral_opt)},\nошибка: {e}.\n ")


