def E(n):
	return [
		[0] * i + [1] + [0] * (n - i - 1)
		for i in range(n)
	]


def decomposition(A):
	rows = len(A)
	oper = 0
	M = [[x for x in row]for row in A]

	P = E(rows)
	for i in range(rows - 1):  # N - 1 its
		elMaxRow = None
		for row in range(i + 1, rows):
			if elMaxRow is None or abs(M[row][i]) > abs(M[elMaxRow][i]):
				elMaxRow = row

		if (not i == elMaxRow):
			M[elMaxRow], M[i] = M[i], M[elMaxRow]
			P[elMaxRow], P[i] = P[i], P[elMaxRow]

		for j in range(i + 1, rows):
			M[j][i] = M[j][i] / M[i][i]
			oper += 1
			for k in range(i + 1, rows):
				M[j][k] = M[j][k] - M[j][i] * M[i][k]
				oper += 2

	L = [
		[(M[i][j] if i > j else 0) for j in range(rows)]
		for i in range(rows)
	]
	for i in range(rows):
		L[i][i] = 1

	U = [
		[(M[i][j] if i <= j else 0) for j in range(rows)]
		for i in range(rows)
	]

	return L, U, P, oper


def solve_slae(L, U, P, b):
	rows = len(L)
	oper = 0
	Pb = [
		b[P[rowIdx].index(1)]
		for rowIdx in range(rows)
	]

	y = [0] * rows
	for i in range(rows):
		y[i] = Pb[i][0] - sum([L[i][j] * y[j] for j in range(i)])
		oper += 2 * i

	x = [0] * rows
	for i in range(rows):
		idx = rows - i - 1
		x[idx] = (y[idx] - sum([x[j] * U[idx][j] for j in range(idx + 1, rows)])) / U[idx][idx]
		oper += 2 * i

	return x, oper