import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return 2 * np.power(x, 2) + 5 * x - 3

def func_der(x):
    return 4 * x + 5

def eq_func(x):
    return (-2 * np.power(x, 2) + 3) / 5

def eq_func_der(x):
    return -4 / 5 * x

def another_localize_solution(func, precision):
    step = np.power(10, -6, dtype = float)
    start_point = end_point = 0

    while func(start_point) * func(end_point) > -precision:
        start_point -= step
        end_point += step

    while func(start_point + step) * func(end_point) < -precision:
        start_point += step

    while func(start_point) * func(end_point - step) < -precision:
        end_point -= step

    return start_point, end_point

def localize_solution(func, precision, search_area, parts=10, find_solution=False):
    solution_range = [search_area[0], search_area[1]]

    while parts < np.power(2, 4):
        solution_range = np.linspace(solution_range[0], solution_range[1], parts)
        for i in range(0, len(solution_range) - 1):
            if func(solution_range[i]) * func(solution_range[i + 1]) < -precision:
                solution_range = solution_range[i:i+2]
                find_solution = True
                break
        solution_range = [solution_range[0], solution_range[-1]]
        parts *= 2

    if not find_solution:
        raise Exception("Can't find a solution!")
        exit()

    return solution_range


def check_convergence(eq_func_der, solution_range):
    step = np.power(10, -7, dtype = float)
    x_range = np.arange(solution_range[0], solution_range[1], step)
    for x in x_range:
        if eq_func_der(x) > 1:
            return False
    return True

def fixed_point_iteration(func, eq_func, eq_func_der, solution_range, precision):
    if not check_convergence(eq_func_der, solution_range):
        raise Exception("Method using eq_func_der don't converge")
    
    x = solution_range[0]
    # better take middle of solution range but because of localization 
    # it will be solution without iterations

    x_list = []
    while np.abs(func(x)) > precision:
        x = eq_func(x)
        x_list.append(x)
    return x, x_list

def Newton_iteraton(func, func_der, solution_range, precision):
    # the necessary condition for convergence is too large
    x = solution_range[0]
    # better take middle of solution range but because of localization 
    # it will be solution without iterations
    
    x_list = []
    while np.abs(func(x)) > precision:
        x = x - func(x) / func_der(x)
        x_list.append(x)
    return x, x_list

def create_equation_plot(x_fix, x_Newton, file):
    plt.xlabel('Iteration')
    plt.ylabel('x')
    plt.title('F(x) = 2x^2 +5x - 3')
    plt.xticks(np.arange(1, len(x_fix) + 1), fontsize = 6)

    plt.plot(np.arange(1, len(x_fix) + 1), x_fix, marker='+', markersize=4)
    plt.plot(np.arange(1, len(x_Newton) + 1), x_Newton, marker='*', markersize=4)

    LegendList = 'fixed_point_iteration', 'Newton_iteraton'
    plt.legend(LegendList, loc="upper right", fontsize="10")

    plt.savefig(file, dpi=800)
    plt.clf()

def create_system_plot(x_vals, y_vals, file):
    plt.xlabel('Iteration')
    plt.ylabel('x')
    plt.title('sin(x+1) - y = 1.2\n2x + cosy = 2')
    plt.xticks(np.arange(1, len(x_vals) + 1), fontsize = 6)

    plt.plot(np.arange(1, len(x_vals) + 1), x_vals, marker='+', markersize=4)
    plt.plot(np.arange(1, len(y_vals) + 1), y_vals, marker='*', markersize=4)

    LegendList = 'x', 'y'
    plt.legend(LegendList, loc="upper right", fontsize="10")

    plt.savefig(file, dpi=800)
    plt.clf()

def non_linear_equation(precision, solution_range):
    solution_fixed_point, x_fix = fixed_point_iteration(func, eq_func, eq_func_der, solution_range, precision)
    solution_Newton, x_Newton = Newton_iteraton(func, func_der, solution_range, precision)

    print('MPI solution    = ', solution_fixed_point)
    # print(x_fix)
    # print()
    print('Newton solution = ', solution_Newton)
    # print(x_Newton)

    create_equation_plot(x_fix, x_Newton, 'graphs/equation.jpg')

def reverse_Jacobian(x, y): #
    f1x = np.cos(x + 1)
    f1y = -1 
    f2x = 2
    f2y = -np.sin(y)

    J = np.array([[f1x, f1y], [f2x, f2y]])

    rev_J = np.linalg.inv(J)
    return rev_J

def sys_func(x, y):
    f1 = np.sin(x + 1) - y - 1.2
    f2 = 2 * x + np.cos(y) - 2
    return [f1, f2]

def non_linear_system(precision, solution_range):
    x = [1, 1]
    x_vals = []
    y_vals = []

    while np.abs(sys_func(x[0], x[1])[0] - sys_func(x[0], x[1])[1]) > precision:
        x = x - reverse_Jacobian(x[0], x[1]) @ sys_func(x[0], x[1])
        x_vals.append(x[0])
        y_vals.append(x[1])


    # print(x_vals)
    # print(y_vals)
    
    create_system_plot(x_vals, y_vals, 'graphs/system.jpg')
    print('System solution = ', x)

def set_precision_and_area():
    precision = np.power(10, -15, dtype = float) 
    search_area = (0, 3)
    return precision, search_area

def main():
    precision, search_area = set_precision_and_area()
    solution_range = localize_solution(func, precision, search_area)

    non_linear_equation(precision, solution_range)
    non_linear_system(precision, solution_range)

if __name__ == '__main__':
    main()