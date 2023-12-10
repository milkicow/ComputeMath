# ЧЕТВЕРТАЯ ЛАБА
#  сделать VI.9.32 + решить ее методом наименьших квадратов (целевую функцию для МНК выберете самостоятельно). 
# Для интереса, можете взять данные по годам для какой-нибудь другой страны из открытых источников по аналогии.

import numpy as np
import matplotlib.pyplot as plt

precision = 10**(-15)

years_USA = np.arange(1910, 2010, 10)
population_USA = np.array([92228496,
                        106021537,
                        123202624,
                        132164569,
                        151325798,
                        179323175,
                        203211926,
                        226545805, 
                        248709873,
                        281421906])

years_Russia = np.arange(1940, 2020, 10)
population_Russia = np.array([110098000, 
                              102067000, 
                              119045800, 
                              130079210, 
                              138126600, 
                              147665081, 
                              146890128,
                              142856536])

years_USSR = np.arange(1920, 2010, 10)
population_USSR = np.array([137727000,
                            157432000,
                            192598000,
                            179229000,
                            208808000,
                            241720000,
                            262436227,
                            294008571,
                            286877000,])

test_x = np.array([0.029, -0.08, -0.122, -0.185])
test_y = np.array([0.2, 0.25, 0.27, 0.3])

def diff_table(x_list: np.array, y_list: np.array):
    table_size = len(y_list)
    table = np.zeros([table_size, table_size])

    table[0] = y_list
    for i in range(1, table_size, 1):
        for j in range(0, table_size - i, 1):
            table[i][j] = (table[i - 1][j + 1] - table[i - 1][j]) / (x_list[j + i] - x_list[j])

    return table

def Newton_interpolation(x, x_list, y_list):
    b_list = np.zeros(len(y_list))
    diff = diff_table(x_list, y_list) 
    for i in range(0, len(y_list), 1):
        b_list[i] = diff[i][0]
    
    P = b_list[0]
    product = 1

    for i in range(1, len(b_list), 1):
        product *= (x - x_list[i - 1])
        P += b_list[i] * product

    return P

def prep_data(x_1_vals, y_1_vals, x_2_vals, y_2_vals, method):
    x_1 = np.arange(1920, 2011, 1)
    y_1 = np.array([])
    for x in x_1:
        y_1 = np.append(y_1, method(x, x_1_vals, y_1_vals))

    x_2 = np.arange(1910, 2011, 1)
    y_2 = np.array([])
    for x in x_2:
        y_2 = np.append(y_2, method(x, x_2_vals, y_2_vals))
    
    return x_1, y_1, x_2, y_2

def create_plot(x_1, y_1, x_2, y_2, method_name, file):
    plt.figure(figsize = (16, 9))
    plt.xlabel('Year')
    plt.xlim(1905, 2025)

    plt.ylabel('Population')
    plt.title('Population in Russia and USA' + ' using ' + method_name)

    plt.plot(x_1, y_1, color= "red")
    plt.plot(x_2, y_2, color= "blue")
    plt.grid(color = "black", linewidth = 0.45, linestyle = "dotted")
    plt.minorticks_on()
    plt.grid(which = "minor", color = "grey", linewidth = 0.25, linestyle = "dashed")

    plt.scatter(years_USSR, population_USSR, label= "+", color= "red", marker= "*", s=30) 
    plt.scatter(years_USA, population_USA, label= "*", color= "blue", marker= "*", s=30) 

    LegendList = 'Russia', 'USA'
    plt.legend(LegendList, loc="upper right", fontsize="10")

    plt.savefig(file, dpi=800)
    plt.clf()

def least_square_method(x, x_list, y_list):
    coef = get_coefficients(x_list, y_list)
    return coef[0] * x**2 + coef[1] * x + coef[2]

def least_square_matrix(x_list: np.array) -> np.array:
    return np.array([[sum_of_pows(x_list, 4), sum_of_pows(x_list, 3), sum_of_pows(x_list, 2)],
                     [sum_of_pows(x_list, 3), sum_of_pows(x_list, 2), sum_of_pows(x_list, 1)],
                     [sum_of_pows(x_list, 2), sum_of_pows(x_list, 1), sum_of_pows(x_list, 0)]])

def least_square_dependent(x_list: np.array, y_list: np.array):
    return np.array([np.sum(x_list**2 * y_list), np.sum(x_list * y_list), np.sum(y_list)])

def get_coefficients(x_list: np.array, y_list: np.array) -> np.array:
    A = least_square_matrix(x_list)
    b = least_square_dependent(x_list, y_list)
    return np.linalg.solve(A, b)

def sum_of_pows(x_list: np.array, pow: int):
    return np.sum(x_list ** pow)

def print_error(method, method_name):
    print('Russia population ' + method_name + f' error = {round(abs(291466400 - method(2010, years_USSR, population_USSR)) / (291466400) * 100, 3)} %')
    print('USA population ' + method_name + f' error    = {round(abs(308745538 - method(2010, years_USA, population_USA)) / (308745538) * 100, 3)} %')

def cubic_spline(x_list: np.array, y_list: np.array):
    a = y_list.copy()
    b = np.array(np.zeros(len(x_list) - 1))
    d = np.array(np.zeros(len(x_list) - 1))

    h = np.array([x_list[i + 1] - x_list[i] for i in range(0, len(x_list) - 1)])
    alpha = np.zeros(len(x_list) - 1)
    for i in range(1, len(x_list) - 1, 1):
        alpha[i] = 3 / h[i] * (a[i + 1] - a[i]) - 3 / h[i - 1] * (a[i] - a[i - 1])
    
    c = np.zeros(len(x_list))
    l = np.zeros(len(x_list))
    m = np.zeros(len(x_list))
    z = np.zeros(len(x_list))

    l[0] = 1
    m[0] = 0
    z[0] = 0

    for i in range(1, len(x_list) - 1):
        l[i] = 2 * (x_list[i + 1] - x_list[i - 1]) - h[i - 1] * m[i - 1]
        m[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[len(x_list) - 1] = 1
    z[len(x_list) - 1] = 0
    c[len(x_list) - 1] = 0

    for j in range(len(x_list) - 2, -1, -1):
        c[j] = z[j] - m[j] * c[j + 1]
        b[j] = (a[j + 1] - a[j]) / h[j] - (h[j] * (c[j + 1] + 2 * c[j])) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])
    
    ret_set = np.array([[0 for x in range(5)] for y in range(len(x_list) - 1)])
    for i in range(0, len(x_list) - 1):
        ret_set[i][0] = a[i]
        ret_set[i][1] = b[i]
        ret_set[i][2] = c[i]
        ret_set[i][3] = d[i]
        ret_set[i][4] = x_list[i]

    return ret_set


def choose_spline(x, separations):
    if (x - separations[0]) < precision:
        return 0
    
    for i in range(0, len(separations) - 1, 1):
        if ((separations[i+1] - x) > precision and (x - separations[i] > precision)) or separations[i] == x:
            return i
        
    if (x - separations[-1]) > precision or separations[-1] == x:
        return len(separations) - 1

def spline_method(x, x_list: np.array, y_list: np.array):
    set = cubic_spline(x_list, y_list)
    sep = x_list[:-1]

    i = choose_spline(x, sep)
    S = set[i][0] + set[i][1] * (x - set[i][4]) + set[i][2] * (x - set[i][4])**2 + set[i][3] * (x - set[i][4])**3
    return S

def main():

    print_error(Newton_interpolation, 'Newton interpolation')
    print_error(least_square_method, 'least square method')
    print_error(spline_method, 'spline method')

    prep_x_1, prep_y_1, prep_x_2, prep_y_2 = prep_data(years_USSR, population_USSR, years_USA, population_USA, Newton_interpolation)
    create_plot(prep_x_1, prep_y_1, prep_x_2, prep_y_2, 'Newton interpolation', 'graphs/Newton_interpolation.jpg')
    
    prep_x_1, prep_y_1, prep_x_2, prep_y_2 = prep_data(years_USSR, population_USSR, years_USA, population_USA, least_square_method)
    create_plot(prep_x_1, prep_y_1, prep_x_2, prep_y_2, 'least square method', 'graphs/least_square_method.jpg')

    prep_x_1, prep_y_1, prep_x_2, prep_y_2 = prep_data(years_USSR, population_USSR, years_USA, population_USA, spline_method)
    create_plot(prep_x_1, prep_y_1, prep_x_2, prep_y_2, 'spline method', 'graphs/spline_method.jpg')

if __name__ == '__main__':
    main()