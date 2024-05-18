import matplotlib.pyplot as plt
import numpy as np

def do_plot(x, t, z):
    fig = plt.figure()
    ax = plt.axes(projection ='3d')
    x, t = np.meshgrid(x, t)
    ax.plot_surface(x, t, z)
    ax.set_title('Transport equastion')
    ax.set_xlabel('x Axis')
    ax.set_ylabel('t Axis')
    ax.set_zlabel('u Axis')
    plt.show()
    plt.close()


if __name__ == "__main__":
    with open("output.txt") as file:
        t_step, x_step = [float(x) for x in next(file).split()]
        numbers = np.array([[float(x) for x in line.split()] for line in file])


    print(t_step)
    print(x_step)
    print(numbers)

    x = np.arange(0, len(numbers[0]) * x_step, x_step)
    t = np.arange(0, len(numbers) * t_step, t_step)
    print(x)
    print(t)

    do_plot(x, t, numbers)