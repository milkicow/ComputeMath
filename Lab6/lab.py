import matplotlib.pyplot as plt
import numpy as np

sigma = 10
r = 28
b = 8 / 3
grid_size = 10000

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


def do_plot(x0, y0, z0, x1, y1, z1, x2, y2, z2):
    ax = plt.figure().add_subplot(projection='3d')

    ax.plot(x0, y0, z0, label='runge_kutta')
    ax.plot(x1, y1, z1, label='adams')
    ax.plot(x2, y2, z2, label='backward_diff')
    ax.legend()
    plt.show()

def prepare_arrays(points: list[list]):
    x = []
    y = []
    z = []

    for point in points:
        x.append(point[0])
        y.append(point[1])
        z.append(point[2])

    return x, y, z

def runge_kutta(y_0: list, t_final: float):
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


if __name__ == "__main__":
    y_list = runge_kutta([1, 1, 1], 50)
    y_list1 = adams([1, 1, 1], 50)
    y_list2 = backward_diff([1, 1, 1], 50)

    x0, y0, z0 = prepare_arrays(y_list)
    x1, y1, z1 = prepare_arrays(y_list1)
    x2, y2, z2 = prepare_arrays(y_list2)
    do_plot(x0, y0, z0, x1, y1, z1, x2, y2, z2)
