import numpy as np

x = np.array([0.25 * x for x in range(0, 9, 1)], dtype=float)
y = np.array([1.0, 0.989616, 0.958851, 0.908852,
      0.841471, 0.759188, 0.664997, 0.562278, 0.454649])

def trapezoid_method(x, y):
    return (x[1] - x[0]) * ((y[0] + y[-1]) / 2 + np.sum(y[1:-1]))

def Simpson_method(x, y):
    return 1/3 * (x[1] - x[0]) * (y[0] + y[-1] + 4 *  np.sum(y[1:-1:2]) + 2 * np.sum(y[2:-1:2]))

def Richardson_extrapolation(x, y, method, order):
    return method(x, y) + (method(x, y) - method(x[::2], y[::2])) / (2 ** order - 1)

def main():
    print(f'I by trapezoid method =  {trapezoid_method(x, y)}')
    print(f'I by Simpson method = {Simpson_method(x, y)}')
    print(f'I by trapezoid method with Richardson extrapolation = {Richardson_extrapolation(x, y, trapezoid_method, 2)}')

if __name__ == '__main__':
    main()