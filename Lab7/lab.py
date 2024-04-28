import sys
import os
sys.path.append(os.path.dirname(__file__) + '/..')
from Lab6.lab import expl_runge_kutta

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm

# y = [y, psi]
def f_y(x, y):
    return y[1]

def f_psi(x, y):
    return -x * np.cos(y[0]) * f_psi.p

def F(x_0: float, y_0: float, alpha: float, x_1: float, y_1: float, grid_size: int) -> float:
    F_list = [f_y, f_psi]
    return expl_runge_kutta(F_list, x_0, [y_0, alpha], x_1, grid_size)[-1][0] - y_1

def F_diff(x_0: float, y_0: float, alpha: float, x_1: float, y_1: float, delta: float, grid_size: int) -> float:
    return (F(x_0, y_0, alpha + delta, x_1, y_1, grid_size) - F(x_0, y_0, alpha, x_1, y_1, grid_size)) / delta

def shooting_method(x_0: float, y_0: float, alpha_0: float, x_1: float, y_1: float, grid_size: int):
    alpha = alpha_0
    delta = 1e-05
    
    result = F(x_0, y_0, alpha, x_1, y_1, grid_size)
    while (np.abs(result) > 1e-10):
        alpha -= result / F_diff(x_0, y_0, alpha, x_1, y_1, delta, grid_size)
        result = F(x_0, y_0, alpha, x_1, y_1, grid_size)

    return alpha

def get_points(x_0: float, y_0: float, alpha: float, x_1: float, y_1: float, grid_size: int):
    F_list = [f_y, f_psi]
    x_list = np.arange(x_0, x_1 + 1 / grid_size, 1 / grid_size)
    y_list = []
    for pair in expl_runge_kutta(F_list, x_0, [y_0, alpha], x_1, grid_size):
        y_list.append(pair[0])
    
    return x_list, y_list

def create_solution_plot(x_lists: list, y_lists: list, p_list: list, file: str):
    ax = plt.figure().add_subplot()

    for i in range(len(p_list)):
        ax.plot(x_lists[i], y_lists[i], label='p = ' + str(p_list[i]))
    ax.set_title(file)
    ax.legend(loc='lower right')
    ax.grid(alpha=0.1)
    
    plt.savefig('./graphs/' + file + '.jpg', dpi=800)
    plt.clf()

def quasi_method(p: float, grid_size: int):
    h = 1 / grid_size
    x_n = np.arange(0, h + 1, h)

    size = len(x_n)
    y = np.ones((1, size))
    y_0 = np.zeros(size)
    y = np.vstack([y, y_0])

    def solve_tridiagonal(a, b, c, d):
        c_ = np.zeros(size)
        d_ = np.zeros(size)

        c_[0] = c[0] / b[0]
        d_[0] = d[0] / b[0]

        for i in range(1, size, 1):
            c_[i] = c[i] / (b[i] - a[i] * c_[i - 1])
            d_[i] = (d[i] - a[i] * d_[i - 1]) / (b[i] - a[i] * c_[i - 1])


        x = np.zeros(size)
        x[-1] = d_[-1]
        for i in range(size - 2, -1, -1):
            x[i] = d_[i] - c_[i] * x[i + 1]
        
        return x

    iter = 1

    while norm(y[-1]-y[-2]) > 1e-10:
        a = [0] + [1] * (size - 1)
        c = (size - 1) * [1] + [0]
        b = -(2 + p * x_n * h**2 * np.sin(y[iter]))
        d = -p * x_n * h**2 * np.cos(y[iter]) - p * x_n * h**2 * y[iter] * np.sin(y[iter])

        x = solve_tridiagonal(a, b, c, d)
        y = np.vstack([y, x])

    return x_n, y[-1]


def tridiagonal_method(grid_size: int):
    def f(x):
        return np.cos(2 * np.pi * x)

    def P(x):
        return 10 + np.sin(2 * np.pi * x)
    
    h = 1 / grid_size
    x = np.arange(0, h + 1, h)
    size = len(x)
    y = np.zeros(size)

    a = np.zeros(size)
    b = np.zeros(size)
    c = np.zeros(size)
    ph = np.zeros(size)

    for i in range(size):
        a[i] = 1
        b[i] = 2 + P(x[i]) * h**2
        c[i] = 1
        ph[i] = f(x[i]) * h**2

    alpha = np.zeros(size)
    betta = np.zeros(size)
    gamma = np.zeros(size)

    alpha[1] = c[0] / b[0]
    betta[1] = -ph[0] / b[0]
    gamma[1] = a[0] / b[0]

    for i in range(1, size - 1):
        alpha[i + 1] = c[i] / (b[i] - alpha[i] * a[i])
        betta[i + 1] = (a[i] * betta[i] - ph[i]) / (b[i] - alpha[i] / a[i])
        gamma[i + 1] = (a[i] * gamma[i]) / (b[i] - alpha[i] * a[i])

    mu = np.zeros(size)
    nu = np.zeros(size)

    mu[size - 1] = -c[size - 1] / ((alpha[size - 1] + gamma[size - 1]) * a[size - 1] - b[size - 1])
    nu[size - 1] = (ph[size - 1] - a[size - 1] * betta[size - 1]) / ((alpha[size - 1] + gamma[size - 1]) * a[size - 1] - b[size - 1])

    for i in range(size - 1, 0, -1):
        mu[i - 1] = alpha[i] * mu[i] + gamma[i] * mu[size - 1]
        nu[i - 1] = betta[i] + alpha[i] * nu[i] + gamma[i] * nu[size - 1]

    y[0] = nu[0] / (1 - mu[0])
    y[size - 1] = mu[size - 1] * y[0] + nu[size - 1]

    for i in range(size - 1, 0, -1):
        y[i - 1] = alpha[i] * (mu[i] * y[0] + nu[i]) + betta[i] + gamma[i] * (mu[size - 1] * y[0] + nu[size - 1])

    return x, y

if __name__ == "__main__":
    grid_size = 1e+03
    x_lists = []
    y_lists = []

    x_lists_q = []
    y_lists_q = []

    p = [1, 3, 7, 25, 50, 100]
    for param in p:
        f_psi.p = param
        alpha = shooting_method(0, 0, 0, 1, 0, grid_size)
        x_list, y_list = get_points(0, 0, alpha, 1, 0, grid_size)

        x_lists.append(x_list)
        y_lists.append(y_list)

        x_list_q, y_list_q = quasi_method(param, grid_size)
        x_lists_q.append(x_list_q)
        y_lists_q.append(y_list_q)

    create_solution_plot(x_lists_q, y_lists_q, p, 'Newton')
    create_solution_plot(x_lists, y_lists, p, "Shooting")

    x_list_t, y_list_t = tridiagonal_method(grid_size)
    create_solution_plot([x_list_t], [y_list_t], ['Tridiagonal'], 'Tridiagonal')
