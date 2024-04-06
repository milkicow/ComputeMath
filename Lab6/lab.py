import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve
import scipy as sp

sigma = 10
r = 28
b = 8 / 3
grid_size = 100000

def f_x(y):
    return -sigma * (y[0] - y[1])

def f_y(y):
    return -y[0] * y[2] + r * y[0] - y[1]

def f_z(y):
    return y[0] * y[1] - b * y[2]

def calc_k_one_arg(f, arg_i: int, y: list, h: float) -> list: # classic 4th rang
    k = []
    k.append(f(y))
    k.append(f(y + [h / 2 * k[0] for i in y]))
    k.append(f(y + [h / 2 * k[1] for i in y]))
    k.append(f(y + [h * k[2] for i in y]))
    return k

def calc_k(F: list, y: list, h: float) -> list[list]:
    k = []
    for i in range(len(y)):
        k.append(calc_k_one_arg(F[i], i, y, h))
    return k

def runge_kutta_step(F: list, y: list, h: float):
    k = calc_k(F, y, h)
    y_next = []

    b = [1/6, 2/6, 2/6, 1/6]
    for i in range(len(y)):
        y_next.append(y[i] + h * sum(b_i * k_i for b_i, k_i in zip(b, k[i])))
    return y_next


def do_plot_expl(x0, y0, z0, x1, y1, z1, x2, y2, z2):
    ax = plt.figure().add_subplot(projection='3d')

    ax.plot(x0, y0, z0, label='runge_kutta')
    ax.plot(x1, y1, z1, label='adams')
    ax.plot(x2, y2, z2, label='backward_diff')
    ax.legend()
    plt.show()

def prepare_arrays_expl(points: list[list]):
    x = []
    y = []
    z = []

    for point in points:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])

    return x, y, z

def expl_runge_kutta(y_0: list, t_final: float):
    h = t_final / grid_size
    F = [f_x, f_y, f_z]
    
    y_list = [y_0]

    t = 0
    y_next = y_0
    while (t < t_final):
        y = y_next

        y_next = runge_kutta_step(F, y, h)
        t += h
        y_list.append(y_next)

    return y_list

def adams(y_0: list, t_final: float):
    h = t_final / grid_size
    F = [f_x, f_y, f_z]

    y_list = [y_0]

    y_prev = y_0
    y = y_0
    y_next = y_0
    t = 0
    while (t < t_final):
        y_prev = y
        y = y_next

        y_next = [0, 0, 0]
        for i in range(len(y)):
            y_next[i] = y[i] + h * (3 / 2 * F[i](y) - 1 / 2 * F[i](y_prev))
    
        y_list.append(y_next)
        t += h
    return y_list

def backward_diff(y_0: list, t_final: float):
    h = t_final / grid_size
    F = [f_x, f_y, f_z]

    y_list = [y_0]

    y_prev = y_0
    y = y_0
    y_next = y_0
    t = 0
    while (t < t_final):
        y_prev = y
        y = y_next
        
        y_next = [0, 0, 0]
        for i in range(len(y)):
            y_next[i] = y_prev[i] + 2 * h * F[i](y)

        y_list.append(y_next)
        t += 0.8
    return y_list

def explicit_methods():
    y_list = expl_runge_kutta([1, 1, 1], 50)
    y_list1 = adams([1, 1, 1], 50)
    y_list2 = backward_diff([1, 1, 1], 50)

    x0, y0, z0 = prepare_arrays_expl(y_list)
    x1, y1, z1 = prepare_arrays_expl(y_list1)
    x2, y2, z2 = prepare_arrays_expl(y_list2)
    do_plot_expl(x0, y0, z0, x1, y1, z1, x2, y2, z2)

alpha = 2.0
betta = 0.0015
gamma = 5.0
tetta0 = 3.0
fi0 = 0.0525
C = 5.0
k1 = 0.05
k2 = 0.35


def f_tetta(y): # y = [tetta, fi]
    return (alpha * y[0] * y[0]) / (y[0] + tetta0) - k1 * y[0] - gamma * y[0] * y[1]

def f_fi(y):
    return betta * y[0] * (1 - y[1] / C) * (1 + (y[1] / fi0) ^ 2) - k2 * y[1]

def F(y):
    return np.array([
        (alpha * y[0] * y[0]) / (y[0] + tetta0) - k1 * y[0] - gamma * y[0] * y[1],
         betta * y[0] * (1 - y[1] / C) * (1 + (y[1] / fi0) * (y[1] / fi0)) - k2 * y[1]
    ])

def impl_runge_kutta(f, y_0: list, t_final: float):
    h = t_final / grid_size

    y_list = [y_0]
    t_list = [0]

    t = 0
    y_next = y_0
    while (t < t_final):
        y = y_next
        k1_ = lambda k: k - f(y + 1/4 * h * k)
        k2_ = lambda k: k - f(y + 1/2 * h * k1 + 1/4 * h * k)

        k1 = fsolve(k1_, y)
        k2 = fsolve(k2_, y)

        y_next = y + h * (1/2 * k1 + 1/2 * k2)
        y_list.append(y_next)
        t_list.append(t)
        t += h
    return y_list, t_list


def rosenbrock(f, y_0: list, t_final: float):
    h = t_final / grid_size
    y_h = h

    y_list = [y_0]
    t_list = [0]

    t = 0
    y_next = y_0
    while (t < t_final):
        y = y_next

        J = (f(y + y_h) - f(y - y_h)) / (2 * y_h)
        E = np.eye(len(y_0))

        A = E - (1 + 1j) / 2 * h
        b = f(y)

        w = np.linalg.solve(A, b)
        y_next = y + h * w.real

        y_list.append(y_next)
        t_list.append(t)
        t += h
    return y_list, t_list

def nordsik(f, y_0: list, t_final: float):
    h = t_final / grid_size
    ord = 4

    y_list = [y_0]
    t_list = [0]

    e = np.array([0, 1, 0, 0, 0])
    l = np.array([251 / 720, 1, 11 / 12, 1 / 3, 1 / 24])
    P = sp.linalg.pascal(ord + 1, kind = 'upper')

    t = 0
    z_0 = np.zeros((ord + 1, len(y_0)))
    z_0[0] = y_0
    z_0[1] = h * f(y_0)

    y_next = y_0
    z_next = z_0

    while (t < t_final):
        y = y_next
        z = z_next

        Pz = P @ z

        # z_next_ = lambda zz: Pz + np.outer(l, (h * f(zz[0]) - np.dot(e, Pz)))
        # z_next = fsolve(z_next_, z)

        z_next = Pz + np.outer(l, (h * f(y) - np.dot(e, Pz)))
        y_next = z_next[0]

        y_list.append(y_next)
        t_list.append(t)

        t += h
    return y_list, t_list

def prepare_arrays_impl(points_lists: list[list]):
    tetta = []
    fi = []

    for points in points_lists:
        tetta.append(points[0])
        fi.append(points[1])

    return tetta, fi

def do_plot_impl(tetta0, fi0, tetta1, fi1, tetta2, fi2, t):
    ax = plt.figure().add_subplot()

    ax.plot(t, tetta0, label='tetta impl runge kutta')
    ax.plot(t, fi0, label='fi impl runge kutta')
    ax.plot(t, tetta1, label='tetta rosenbrock')
    ax.plot(t, fi1, label = 'fi rosenbrock')
    ax.plot(t, tetta2, label='tetta nordsik')
    ax.plot(t, fi2, label = 'fi nordsik')
    
    ax.legend()
    plt.show()

def implicit_methods():
    y_list, t_list = rosenbrock(F, np.array([tetta0, fi0]), 10)
    tetta1_, fi1_ = prepare_arrays_impl(y_list)

    y_list, t_list = impl_runge_kutta(F, np.array([tetta0, fi0]), 10)
    tetta0_, fi0_ = prepare_arrays_impl(y_list)

    y_list, t_list = nordsik(F, np.array([tetta0, fi0]), 10)
    tetta2_, fi2_ = prepare_arrays_impl(y_list)
    do_plot_impl(tetta0_, fi0_, tetta1_, fi1_, tetta2_, fi2_, t_list)


if __name__ == "__main__":
    explicit_methods()
    implicit_methods()