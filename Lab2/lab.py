import numpy as np
import matplotlib.pyplot as plt

def norm1(x):
    return np.linalg.norm(x, np.inf)

def norm2(x):
    return np.linalg.norm(x, 1)

def norm3(x):
    return np.linalg.norm(x, 3)

def calculate_LU(A):
    matrix_size = len(A)

    U = np.zeros((matrix_size, matrix_size))
    L = np.diag(np.ones(matrix_size))

    for i in range(0, matrix_size):
        for j in range(0, i):
            L[i, j] = (A[i, j] - sum(L[i, k] * U[k, j] for k in range(0, j))) / U[j, j]
        for j in range(i, matrix_size):
            U[i, j] = A[i, j] - sum(L[i, k] * U[k, j] for k in range(0, i))
    return L, U

def forward_substitution(L, f):
    matrix_size = len(L)
    y = np.zeros(matrix_size)

    for k in range(0, matrix_size):
        y[k] = (f[k] - sum(L[k, j] * y[j] for j in range(0, k))) / L[k, k]
    return y

def back_substitution(U, y):
    matrix_size = len(U)
    x = np.zeros(matrix_size)

    for k in range(matrix_size - 1, -1, -1):
        x[k] = ((y[k] - sum(U[k, j] * x[j] for j in range(k + 1, matrix_size)))) / U[k, k]
    return x

def LU_method(A, f):
    L, U = calculate_LU(A)
    y = forward_substitution(L, f)
    x = back_substitution(U, y)

    return x

def Gauss_method(A, f):
    matrix_size = len(A)
    x = np.zeros(matrix_size)

    for k in range(1, matrix_size):
        for j in range(k, matrix_size):
            m = A[j, k - 1] / A[k - 1, k - 1]

            for i in range(0, matrix_size):
                A[j, i] = A[j, i] - m * A[k - 1, i]

            f[j] -= m * f[k - 1]
    
    return back_substitution(A, f)

def Seidel_method(A, f, iter):
    D = np.diag(np.diag(A))
    L = np.tril(A) - D
    U = A - L - D

    R = -np.linalg.inv(L + D) @ U
    F = np.linalg.inv(L + D) @ f

    x = np.zeros(len(A))
    for i in range(0, iter):
        x = R @ x + F
    return x

def Jacobi_method(A, f, iter): #  not working
    D = np.diag(np.diag(A))
    L = np.tril(A, k = -1)
    U = A - L - D

    R = -np.linalg.inv(D) @ (L + U)
    F = np.linalg.inv(D) @ f

    x = np.zeros(len(A))
    for i in range(0, iter):
        x = R @ x + F
    return x

def relaxation_method(A, f, iter):
    D = np.diag(np.diag(A))
    L = np.tril(A, k = -1)
    U = A - L - D
    w = 1.1

    x = np.zeros(len(A))
    for i in range(0, iter):
        x = np.linalg.inv(D + w * L) @ (w * f - (w * U + (w - 1) * D) @ x)

    return x

def matrix_ctor():
    a = 10
    matrix_size = 100

    A = np.zeros((matrix_size, matrix_size))
    
    for i in range(0, matrix_size):
        for k in range(max(0, i - 4), i):
            A[i, k] = 1
        for j in range(i + 1, min(matrix_size, i + 5)):
            A[i, j] = 1
    
    A += np.diag(np.full(matrix_size, a, dtype = float))
    f = np.arange(1, matrix_size + 1, dtype = float)
    return A, f


def createResidualAxis(A, f, norm, method, iter_list):
    residual_axis = []
    
    for i in iter_list:
        x = method(A, f, i)
        residual_axis.append(norm(A @ x - f))
    return residual_axis

def createPlot(A, f, norm_list, method, iter_list, norm_name, method_name, file):
    iter_axis = iter_list
    residual_axis = []

    for norm in norm_list:
        residual_axis.append(createResidualAxis(A, f, norm, method, iter_list))
    residual_axis = np.abs(residual_axis)

    plt.xlabel('Iterations')
    plt.ylabel('Residual')
    plt.title(method_name)

    plt.semilogx()
    plt.semilogy()

    for i in residual_axis:
        plt.plot(iter_axis, i, marker='+', markersize=4)

    LegendList = norm_name
    plt.legend(LegendList, loc="upper right", fontsize="7")
    plt.savefig(file, dpi=800)
    plt.clf()


A, f = matrix_ctor()
epsilon = 10^-14

norm = [norm1, norm2, norm3]
norm_name = ["norm_1", "norm_2", "norm_3"]

x = LU_method(A.copy(), f.copy())

print("Calculate residual for LU_method:")
for n in range(0, 3):
    output = "residual of " + norm_name[n] + " equals to "
    residual = norm[n](A @ x - f)
    print(output, residual)

x = Gauss_method(A.copy(), f.copy())
print("Calculate residual for Gauss_method:")
for n in range(0, 3):
    output = "residual of " + norm_name[n] + " equals to "
    residual = norm[n](A @ x - f)
    print(output, residual)


iter_methods = [Seidel_method, Jacobi_method, relaxation_method]
methods_name = ["Seidel_method", "Jacobi_method", "relaxation_method"]

iter_list = np.arange(1, 100)

folder = 'graphs/'
for method in range(0, 3):
    createPlot(A.copy(), f.copy(), norm, iter_methods[method], iter_list, norm_name, methods_name[method], folder+methods_name[method]+".jpg")
