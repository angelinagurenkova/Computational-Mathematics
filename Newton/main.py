import math
import numpy as np
import time
from LU import decomposition, solve_slae

Geps = 1e-3
solutionExists = 1


def F(x, mulBy=1):
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = x
    Fx = np.array([
        math.cos(x2 * x1) - math.exp(-3 * x3) + x4 * x5 ** 2 - x6 - math.sinh(2 * x8) * x9 + 2 * x10 + 2.000433974165385440,
        math.sin(x2 * x1) + x3 * x9 * x7 - math.exp(-x10 + x6) + 3 * x5 ** 2 - x6 * (x8 + 1) + 10.886272036407019994,
        x1 - x2 + x3 - x4 + x5 - x6 + x7 - x8 + x9 - x10 - 3.1361904761904761904,
        2 * math.cos(-x9 + x4) + x5 / (x3 + x1) - math.sin(x2 ** 2) + math.cos(x7 * x10) ** 2 - x8 - 0.1707472705022304757,
        math.sin(x5) + 2 * x8 * (x3 + x1) - math.exp(-x7 * (-x10 + x6)) + 2 * math.cos(x2) - 1.0 / (-x9 + x4) - 0.3685896273101277862,
        math.exp(x1 - x4 - x9) + x5 ** 2 / x8 + math.cos(3 * x10 * x2) / 2 - x6 * x3 + 2.0491086016771875115,
        x2 ** 3 * x7 - math.sin(x10 / x5 + x8) + (x1 - x6) * math.cos(x4) + x3 - 0.7380430076202798014,
        x5 * (x1 - 2 * x6) ** 2 - 2 * math.sin(-x9 + x3) + 0.15e1 * x4 - math.exp(x2 * x7 + x10) + 3.5668321989693809040,
        7 / x6 + math.exp(x5 + x4) - 2 * x2 * x8 * x10 * x7 + 3 * x9 - 3 * x1 - 8.4394734508383257499,
        x10 * x1 + x9 * x2 - x8 * x3 + math.sin(x4 + x5 + x6) * x7 - 0.78238095238095238096])

    return mulBy * Fx


def J(x):
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = x
    Jx = np.array([[-x2 * math.sin(x2 * x1), -x1 * math.sin(x2 * x1), 3 * math.exp(-3 * x3), x5 ** 2, 2 * x4 * x5,
                    -1, 0, -2 * math.cosh(2 * x8) * x9, -math.sinh(2 * x8), 2],
                   [x2 * math.cos(x2 * x1), x1 * math.cos(x2 * x1), x9 * x7, 0, 6 * x5,
                    -math.exp(-x10 + x6) - x8 - 1, x3 * x9, -x6, x3 * x7, math.exp(-x10 + x6)],
                   [1, -1, 1, -1, 1, -1, 1, -1, 1, -1],
                   [-x5 / (x3 + x1) ** 2, -2 * x2 * math.cos(x2 ** 2), -x5 / (x3 + x1) ** 2, -2 * math.sin(-x9 + x4),
                    1.0 / (x3 + x1), 0, -2 * math.cos(x7 * x10) * x10 * math.sin(x7 * x10), -1,
                    2 * math.sin(-x9 + x4), -2 * math.cos(x7 * x10) * x7 * math.sin(x7 * x10)],
                   [2 * x8, -2 * math.sin(x2), 2 * x8, 1.0 / (-x9 + x4) ** 2, math.cos(x5),
                    x7 * math.exp(-x7 * (-x10 + x6)), -(x10 - x6) * math.exp(-x7 * (-x10 + x6)), 2 * x3 + 2 * x1,
                    -1.0 / (-x9 + x4) ** 2, -x7 * math.exp(-x7 * (-x10 + x6))],
                   [math.exp(x1 - x4 - x9), -1.5 * x10 * math.sin(3 * x10 * x2), -x6, -math.exp(x1 - x4 - x9),
                    2 * x5 / x8, -x3, 0, -x5 ** 2 / x8 ** 2, -math.exp(x1 - x4 - x9), -1.5 * x2 * math.sin(3 * x10 * x2)],
                   [math.cos(x4), 3 * x2 ** 2 * x7, 1, -(x1 - x6) * math.sin(x4), x10 / x5 ** 2 * math.cos(x10 / x5 + x8),
                    -math.cos(x4), x2 ** 3, -math.cos(x10 / x5 + x8), 0, -1.0 / x5 * math.cos(x10 / x5 + x8)],
                   [2 * x5 * (x1 - 2 * x6), -x7 * math.exp(x2 * x7 + x10), -2 * math.cos(-x9 + x3), 1.5,
                    (x1 - 2 * x6) ** 2, -4 * x5 * (x1 - 2 * x6), -x2 * math.exp(x2 * x7 + x10), 0, 2 * math.cos(-x9 + x3),
                    -math.exp(x2 * x7 + x10)],
                   [-3, -2 * x8 * x10 * x7, 0, math.exp(x5 + x4), math.exp(x5 + x4),
                    -7.0 / x6 ** 2, -2 * x2 * x8 * x10, -2 * x2 * x10 * x7, 3, -2 * x2 * x8 * x7],
                   [x10, x9, -x8, math.cos(x4 + x5 + x6) * x7, math.cos(x4 + x5 + x6) * x7,
                    math.cos(x4 + x5 + x6) * x7, math.sin(x4 + x5 + x6), -x3, x2, x1]])

    return Jx


def norm_vec(x):
    return math.sqrt(sum([el ** 2 for el in x]))

def norm_inf_mat(A):
    return max([sum(row) for row in A])

x = [0.5, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5]
X = np.array(x)

def primitive_metod(x, maxIter=2000, eps=Geps):
    xc = x
    xp = np.zeros(len(x))
    operations = 0
    iterations = 0
    while (iterations < maxIter):
        L, U, P, ops = decomposition(J(xc))
        operations += ops
        dk, ops = solve_slae(L, U, P, F(xc, mulBy=-1).reshape(10, 1))
        operations += ops
        xc, xp = xc + dk, xc
        iterations += 1
        if (norm_vec(xc - xp) < eps):
            break

    print('Метод Ньютона:')
    print(f'X:  {[f"{x:.5f}" for x in xc]}')

    return (iterations, operations)

s_time = time.time()
its, ops = primitive_metod(X)
print(f"Итерации: {its}")
print(f"Операции: {ops}")
print(f"Время: {(time.time() - s_time):.1e}second\n")

def modified_metod(x, maxIter=2000, eps=Geps):
    xc = x
    xp = np.zeros(len(x))
    operations = 0
    L, U, P, ops = decomposition(J(xc))
    operations += ops
    iterations = 0
    while (iterations < maxIter):
        dk, ops = solve_slae(L, U, P, F(xc, mulBy=-1).reshape(10, 1))
        operations += ops
        xc, xp = xc + dk, xc
        iterations += 1
        if (norm_vec(xc - xp) < eps):
            break

    print('Модифицированный метод Ньютона:')
    print(f'X:  {[f"{x:.5f}" for x in xc]}')
    return (iterations, operations)

s_time = time.time()
its, ops = modified_metod(X)
print(f"Итерации: {its}")
print(f"Операции: {ops}")
print(f"Время: {(time.time() - s_time):.1e}second\n")

def combined_metod(x, k=7, maxIter=2000, eps=Geps):
    xc = x
    xp = np.zeros(len(x))
    operations = 0
    L, U, P = [[0] * len(x) for _ in range(3)]
    iterations = 0
    for _ in range(k):
        L, U, P, ops = decomposition(J(xc))
        operations += ops
        dk, ops = solve_slae(L, U, P, F(xc, mulBy=-1).reshape(10, 1))
        operations += ops
        xc, xp = xc + dk, xc
        iterations += 1
        if (norm_vec(dk) < eps):
            print('Combined newton:')
            print(f'X:  {[f"{x:.5f}" for x in xc]}')
            return (iterations, operations)

    while (iterations < maxIter):
        dk, ops = solve_slae(L, U, P, F(xc, mulBy=-1).reshape(10, 1))
        operations += ops
        xc, xp = xc + dk, xc
        iterations += 1
        if (norm_vec(xc - xp) < eps):
            break

    print('Комбинированный метод Ньютона:')
    print(f'X:  {[f"{x:.5f}" for x in xc]}')
    return (iterations, operations)

s_time = time.time()
its, ops = combined_metod(X)
print(f"Итерации: {its}")
print(f"Операции: {ops}")
print(f"Время: {(time.time() - s_time):.1e}second\n")

def auto_combined_metod(x, maxIter=2000, eps=Geps):
    xc = x
    xp = np.zeros(len(x))
    operations = 0
    L, U, P = [[0] * len(x) for _ in range(3)]
    iterations = 0
    while (iterations < maxIter):
        L, U, P, ops = decomposition(J(xc))
        operations += ops
        dk, ops = solve_slae(L, U, P, F(xc, mulBy=-1).reshape(10, 1))
        operations += ops
        xc, xp = xc + dk, xc
        iterations += 1
        if (norm_inf_mat(J(xc) - J(xp)) < eps):
            break

    while (iterations < maxIter):
        dk, ops = solve_slae(L, U, P, F(xc, mulBy=-1).reshape(10, 1))
        operations += ops
        xc, xp = xc + dk, xc
        iterations += 1
        if (norm_vec(xc - xp) < eps):
            break

    print('Автоматизированный комбинированный метод Ньютона:')
    print(f'X:  {[f"{x:.5f}" for x in xc]}')
    return (iterations, operations)

s_time = time.time()
its, ops = auto_combined_metod(X)
print(f"Итерации: {its}")
print(f"Операции: {ops}")
print(f"Время: {(time.time() - s_time):.1e}second\n")

def hybrid_metod(x, k=2, maxIter=2000, eps=Geps):
    xc = x
    xp = np.zeros(len(x))
    operations = 0
    L, U, P = [[0] * len(x) for _ in range(3)]
    iterations = 0
    while (iterations < maxIter):
        if (iterations % k == 0):  # primitive
            L, U, P, ops = decomposition(J(xc))
            operations += ops
        dk, ops = solve_slae(L, U, P, F(xc, mulBy=-1).reshape(10, 1))
        operations += ops
        xc, xp = xc + dk, xc
        iterations += 1
        if (norm_vec(xc - xp) < eps):
            break

    print('Гибридный метод Ньютона:')
    print(f'X:  {[f"{x:.5f}" for x in xc]}')
    return (iterations, operations)

s_time = time.time()
its, ops = hybrid_metod(X)
print(f"Итерации: {its}")
print(f"Операции: {ops}")
print(f"Время: {(time.time() - s_time):.1e}second\n")
